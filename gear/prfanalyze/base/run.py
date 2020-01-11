#! /usr/bin/env python

from __future__ import print_function

import json, os, sys, csv, pimms

output_dir  = '/flywheel/v0/output'
input_dir   = '/flywheel/v0/input'
config_file = os.path.join(input_dir, 'config.json')
bids_dir    = os.path.join(input_dir, 'BIDS')
opts_file   = '/running/options.json'
verbose     = os.environ.get('VERBOSE', '0').strip() == '1'
force       = os.environ.get('FORCE', '0').strip() == '1'
solver_name = os.environ.get('PRF_SOLVER', None)
if solver_name is None:
    print("WARNING: The PRF_SOLVER environment variable is not set; using 'base'")
    solver_name = 'base'

    
# check for a separate config file
if len(sys.argv) > 1:
    config_file = sys.argv[1]

def die(*args):
    print(*args)
    sys.exit(1)
def note(*args):
    if verbose: print(*args)
    return None

if not os.path.isdir(bids_dir):
    die('no BIDS directory found!')

try:
    with open(config_file, 'r') as fl:
        conf = json.load(fl)
except Exception:
    die("Could not read config.json!")

if not pimms.is_map(conf):
    die("config.json must contain a single dictionary")
if 'subjectName' not in conf or not pimms.is_str(conf['subjectName']):
    die('config.json does not contain a valid "subjectName" entry')
if 'sessionName' not in conf or not pimms.is_str(conf['sessionName']):
    die('config.json does not contain a valid "sessionName" entry')

# we just have to find the relevant files then echo them for the calling script; in the case of the
# config file, we write out a new one in the /running directory
sub  = conf['subjectName']
ses  = conf['sessionName']
opts = conf.get('options', {})
with open(opts_file, 'w') as fl:
    json.dump(opts, fl)
note("Preparing solver \"%s\":" % (solver_name,))
note("  Subject: %s" % sub)
note("  Session: %s" % ses)
note("  Options: %s" % (opts,))

# find the relevant files in the BIDS dir; first, the BOLD image is easy to find:
func_dir = os.path.join(bids_dir, 'sub-' + sub, 'ses-' + ses, 'func')
bold_image = os.path.join(func_dir,
                          'sub-%s_ses-%s_task-prf_acq-normal_run-01_bold.nii.gz' % (sub, ses))
if not os.path.isfile(bold_image):
    die("BOLD image (%s) not found!" % bold_image)
# we get the stimulus filename from the events file:
events_file = os.path.join(func_dir, 'sub-%s_ses-%s_task-prf_events.tsv' % (sub, ses))
try:
    with open(events_file, 'r') as fl:
        rr = csv.reader(fl, delimiter='\t', quotechar='"')
        l0 = next(rr)
        if 'stim_file' not in l0:
            die('stim_file must be a column in the events file (%s)' % events_file)
        rows = [{k:v for (k,v) in zip(l0,r)} for r in rr]
except Exception:
    die("Could not load events file: %s" % events_file)
stim_file = set([r['stim_file'] for r in rows])
if len(stim_file) != 1:
    die("Multiple stimulus files found in events file (%s)" % events_File)
stim_file = os.path.join(bids_dir, 'stimuli', list(stim_file)[0])
if not os.path.isfile(stim_file):
    die("Stimulus file (%s) not found" % stim_file)
# we also need the stimulus json file, which is in the derivatives directory
stimjs_file = os.path.join(bids_dir, 'derivatives', 'prfsynth', 'sub-'+sub, 'ses-'+ses,
                           'sub-%s_ses-%s_task-prf_acq-normal_run-01_bold.json' % (sub, ses))
if not os.path.isfile(stimjs_file):
    die("Stimulus JSON file (%s) not found" % stimjs_file)

# Finally, we need to find the output directory (in the OUTPUT directory's BIDS directory)
# To figure out how we name the directory, we use the PRF_SOLVER environment variable
elif not solver_name.startswith('prfanalyze-'):
    solver_name = 'prfanalyze-' + solver_name
outbids_dir = os.path.join(output_dir, 'BIDS', 'derivatives', solver_name, 'sub-'+sub, 'ses-'+ses)
# if this directory already exists, it's safe to assume that we need a temporary directory (unless
# the force option is invoked):
if not force and os.path.isdir(outbids_dir):
    outbids_dir = None
    for k in range(1000):
        p = os.path.join(output_dir, 'BIDS', 'derivatives', solver_name + '_temp%03d' % k)
        if not os.path.isdir(p):
            # we've found it!
            outbids_dir = os.path.join(p, 'sub-'+sub, 'ses-'+ses)
            print("WARNING: Using temporary output directory: %s" % p)
            break
    if outbids_dir is None: die("Could not find a valid temporary directory!")
try:
    if not os.path.isdir(outbids_dir): os.makedirs(outbids_dir)
except Exception:
    die("Error creating output BIDS directory: %s" % outbids_dir)
note("Output BIDS directory: %s" % outbids_dir)

# Last thing before executing the script: we make a symlink from the output bids dir to /running
bids_link = '/running/output_bids'
try:
    if os.path.islink(bids_link): os.remove(bids_link)
    os.symlink(outbids_dir, bids_link)
except Exception:
    die("Could not create output link: %s" % bids_link)

# okay, we have the files; run the solver script!
try:
    pid = os.fork()
    if pid == 0:
        os.execl("/solve.sh", "/solve.sh",
                 opts_file, bold_image, stim_file, stimjs_file, outbids_dir)
    else:
        note("Beginning os.wait() for /solve.sh (child pid is %s)" % (pid,))
        os.wait()
except Exception:
    die("Failed to exec /solve.sh script!")

# If there are things to cleanup, specifically, if there is an estimates.mat file, we process that
estfl = '/running/output_bids/estimates.mat'
if os.path.isfile(estfl):
    note("Processing estimates.mat file...")
    from scipy.io import loadmat
    import nibabel as nib, numpy as np
    dat = loadmat(estfl)
    # decode the data...
    dat = dat['estimates'][0,0][4][0]
    (testdat, x0, y0, th, sigmin, sigmaj, pred) = [np.squeeze(u) for u in dat]
    # write out niftis
    nii_base = nib.load(bold_image)
    for (k,d) in zip(['testdata','modelpred'], [testdat, pred]):
        im = nib.Nifti1Image(np.reshape(d, (d.shape[0], 1, 1, d.shape[-1])),
                             nii_base.affine, nii_base.header)
        im.to_filename('/running/output_bids/' + k + '.nii.gz')
        note("  * " + k + ".nii.gz")
    for (k,d) in zip(['x0','y0','theta','sigmamajor','sigmaminor'], [x0,y0,th,sigmin,sigmaj]):
        im = nib.Nifti1Image(np.reshape(d, (-1, 1, 1, 1)), nii_base.affine, nii_base.header)
        im.to_filename('/running/output_bids/' + k + '.nii.gz')
        note("  * " + k + ".nii.gz")
else:
    note("No estimates.mat file found.")
    
# exit happily
sys.exit(0)
