#! /usr/bin/env python

import json, os, sys, csv, pimms
import nibabel as nib

output_dir  = '/flywheel/v0/output'
input_dir   = '/flywheel/v0/input'
config_file = os.path.join(input_dir, 'config.json')
bids_dir    = os.path.join(input_dir, 'BIDS')
# bids_link   = '/running/out'
# Fix so that it works in Singularity
bids_link   = os.path.join(output_dir,'out')
verbose     = os.environ.get('VERBOSE', '0').strip() == '1'
force       = os.environ.get('FORCE', '0').strip() == '1'
bids_fields = os.environ.get('FIELDS', 'task-prf_acq-normal')
solver_name = os.environ.get('PRF_SOLVER', None)
bids_fieldmap = [ss.split('-') for ss in bids_fields.split('_')]
bids_fields_noacq = '_'.join(['-'.join(ff) for ff in bids_fieldmap if ff[0] != 'acq'])
if bids_fields != '': bids_fields = '_' + bids_fields
if bids_fields_noacq != '': bids_fields_noacq = '_' + bids_fields_noacq
    
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
if 'isPRFSynthData' not in conf:
    note('Warning: "isPRFSynthData" not found in config JSON; assuming True.')
    conf['isPRFSynthData'] = True
    synthQ = True
else: synthQ = False
if 'solver' in conf:
    solver_name = conf['solver']

# now that we've read in the solver from the config we can process it.
if solver_name is None:
    print("WARNING: The PRF_SOLVER environment variable is not set; using 'base'")
    solver_name = 'base'
if not solver_name.startswith('prfanalyze-'):
    solver_name = 'prfanalyze-' + solver_name


# we just have to find the relevant files then echo them for the calling script; in the case of the
# config file, we write out a new one in the /running directory
sub  = conf['subjectName']
ses  = conf['sessionName']
opts = conf.get('options', {})
note("Preparing solver \"%s\":" % (solver_name,))
note("  Subject: %s" % sub)
note("  Session: %s" % ses)
note("  Options: %s" % (opts,))

# do all the processing for PRFSynthData
if conf['isPRFSynthData']:
    print("[base/run.py] Using synthetic data coming from prfsynthesize.")
    # find the relevant files in the BIDS dir; first, the BOLD image is easy to find:
    func_dir = os.path.join(bids_dir, 'sub-' + sub, 'ses-' + ses, 'func')
    # we get the stimulus filename from the events file:
    events_file = os.path.join(func_dir, 'sub-%s_ses-%s%s_events.tsv' % (sub, ses, bids_fields_noacq))
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
        die("Multiple stimulus files found in events file (%s)" % events_file)
    stim_file = os.path.join(bids_dir, 'stimuli', list(stim_file)[0])
    if not os.path.isfile(stim_file):
        die("Stimulus file (%s) not found" % stim_file)
    
    # Finally, we need to find the output directory (in the OUTPUT directory's BIDS directory)
    # To figure out how we name the directory, we use the PRF_SOLVER environment variable
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
    # # we make a symlink from the output bids dir to /running
    # try:
    #     if os.path.islink(bids_link): os.remove(bids_link)
    #     os.symlink(outbids_dir, bids_link)
    # except Exception:
    #     die("Could not create output link: %s" % bids_link)
    
    # Noahs solution didnt work, there should be a symlink for run.sh
    ## # we make a symlink from the output bids dir to /running
    try:
        if os.path.islink(bids_link): os.remove(bids_link)
        os.symlink(outbids_dir, bids_link)
    except Exception:
        die("Could not create output link: %s" % bids_link)
    # dump the options file in the output directory
    opts_file   = os.path.join(bids_link, 'options.json')
    with open(opts_file, 'w') as fl:
        json.dump(opts, fl)
    
    # We may have any number of runs, find them all:
    bold_prefix = 'sub-%s_ses-%s%s_run-' % (sub, ses, bids_fields)
    bold_suffix = '_bold.nii.gz'
    (pn,sn) = (len(bold_prefix), len(bold_suffix))
    processed = 0
    for flnm in os.listdir(func_dir):
        if not (flnm.startswith(bold_prefix) and flnm.endswith(bold_suffix)): continue
        runid = flnm[pn:-sn]
        bold_image = os.path.join(func_dir, flnm)
        # if the data are from prfsynth, we also need the stimulus json file, which is in the
        # derivatives directory
        stimjs_file = os.path.join(
            bids_dir, 'derivatives', 'prfsynth', 'sub-'+sub, 'ses-'+ses,
            'sub-%s_ses-%s%s_run-%s_bold.json' % (sub, ses, bids_fields, runid))
        if not os.path.isfile(stimjs_file):
            die("Stimulus JSON file (%s) not found" % stimjs_file)
        # okay, we have the files; run the solver script!
        try:
            pid = os.fork()
            if pid == 0:
                os.execl("/solve.sh", "/solve.sh",
                         opts_file, bold_image, stim_file, stimjs_file, outbids_dir)
            else:
                note("Beginning os.wait() for /solve.sh, run=%s (child pid is %s)" % (runid, pid))
                os.wait()
        except Exception:
            die("Failed to exec /solve.sh script!")
        
        nii_base = nib.load(bold_image)
        # If there are things to cleanup we do that; specifically, the estimates.json file:
        # estfl = os.path.join(bids_link, 'estimates.json')
        estfl = os.path.join(outbids_dir, 'estimates.json')
        if os.path.isfile(estfl):
            note("Processing estimates.json file...")
            import nibabel as nib, numpy as np
            with open(estfl, 'r') as fl:
                dat = json.load(fl)
            # decode the data...
            dat = {k: np.asarray([u[k] for u in dat]) for k in dat[0].keys()}
            for (k,v) in dat.items():
                if len(v.shape) == 2:            
                    im = nib.Nifti2Image(np.reshape(v, (v.shape[0], 1, 1, v.shape[-1])),
                                         nii_base.affine, nii_base.header)
                else:
                    im = nib.Nifti2Image(np.reshape(v, (-1, 1, 1, 1)),
                                         nii_base.affine, nii_base.header)
                print("Writing the estimates.json to nifti2 in outbids_dir: " + outbids_dir)
                # im.to_filename(os.path.join(bids_link, 'run-%s_%s.nii.gz' % (runid,k.lower())))
                im.to_filename(os.path.join(outbids_dir, 'run-%s_%s.nii.gz' % (runid,k.lower())))
            # os.rename(estfl, os.path.join(bids_link, 'run-%s_estimates.json' % (runid,)))
            os.rename(estfl, os.path.join(outbids_dir, 'run-%s_estimates.json' % (runid,)))
        else:
            note("No estimates.json file found.")
        # also rename results,.mat if it's there
        # resfli = os.path.join(bids_link, 'results.mat')
        # resflo = os.path.join(bids_link, 'run-%s_results.mat' % (runid,))
        resfli = os.path.join(outbids_dir, 'results.mat')
        resflo = os.path.join(outbids_dir, 'run-%s_results.mat' % (runid,))
        if os.path.isfile(resfli): os.rename(resfli, resflo)
        processed += 1

# if we have real data in prfprepare form
else:
    print("[base/run.py] Using real data, not coming from prfsynthesize.")
    
    # read additional bids fields for the filenames from the config.json
    if not 'tasks' in conf.keys(): die('Specify the tasks in the config file!')
    if not 'areas' in conf.keys(): die('Specify the areas in the config file!')
    tasks = conf['tasks'].split(']')[0].split('[')[-1].split(',')
    areas = conf['areas'].split(']')[0].split('[')[-1].split(',')
    
    if 'prfprepareAnalysis' in conf.keys():
        prfprep_analyis = conf['prfprepareAnalysis']
    else:
        print('You should consider adding an prfprepareAnalysisNumber in the config!')
        print('Automatically setting it to "01".')
        prfprep_analyis = '01'
    
    # define the prfprepare folder    
    prfprep_dir = os.path.join(bids_dir, 'derivatives', 'prfprepare', f'analysis-{prfprep_analyis}')
                               
    # find the relevant files in the BIDS dir; first, the BOLD image is easy to find:
    func_dir = os.path.join(prfprep_dir, f'sub-{sub}', f'ses-{ses}', 'func')
    
    # Finally, we need to find the output directory (in the OUTPUT directory's BIDS directory)
    # To figure out how we name the directory, we use the PRF_SOLVER environment variable
    p = os.path.join(output_dir, 'BIDS', 'derivatives', solver_name)
    outbids_dir = os.path.join(p, f'sub-{sub}', f'ses-{ses}')
    # if this directory already exists, it's safe to assume that we need a temporary directory (unless
    # the force option is invoked):
    if not force and os.path.isdir(outbids_dir):
        outbids_dir = None
        for k in range(1000):
            p = os.path.join(output_dir, 'BIDS', 'derivatives', solver_name + '_temp%03d' % k)
            if not os.path.isdir(p):
                # we've found it!
                outbids_dir = os.path.join(p, f'sub-{sub}', f'ses-{ses}')
                print("WARNING: Using temporary output directory: %s" % p)
                break
        if outbids_dir is None: die("Could not find a valid temporary directory!")
    try:
        if not os.path.isdir(outbids_dir): os.makedirs(outbids_dir)
    except Exception:
        die("Error creating output BIDS directory: %s" % outbids_dir)
    note("Output BIDS directory: %s" % outbids_dir)
    
    # dump the options file in the output directory
    opts_file = os.path.join(p, 'options.json')
    with open(opts_file, 'w') as fl:
        json.dump(opts, fl)
    
    # loop through specified areas and tasks
    for area in areas:
        for task in tasks:
            # We may have any number of runs, find them all:
            bold_prefix = f'sub-{sub}_ses-{ses}_task-{task}_run-'
            bold_suffix = f'_desc-{area}_bold.nii.gz'
            (pn,sn) = (len(bold_prefix), len(bold_suffix))
            processed = 0        
        
            for flnm in os.listdir(func_dir):
                if not (flnm.startswith(bold_prefix) and flnm.endswith(bold_suffix)): continue
                runid = flnm[pn:-sn]
                bold_image = os.path.join(func_dir, flnm)
                
                # we get the stimulus filename from the events file:
                events_file = os.path.join(func_dir, f'sub-{sub}_ses-{ses}_task-{task}_run-{runid}_events.tsv')
                # if ther is no run specific events file use the one withour run param
                if not os.path.exists(events_file):
                    events_file = os.path.join(func_dir, f'sub-{sub}_ses-{ses}_task-{task}_events.tsv')
            
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
                    die("Multiple stimulus files found in events file (%s)" % events_file)
                stims_dir = os.path.join(prfprep_dir, f'sub-{sub}', 'stimuli')
                stim_file = os.path.join(stims_dir, list(stim_file)[0])
                if not os.path.isfile(stim_file):
                    die("Stimulus file (%s) not found" % stim_file)
                
                
                if 'stimulus' not in conf:
                    die("In config.json, isPRFSynthData is False, but no stimulus settings were given.")
                stim = conf['stimulus']
                if not isinstance(stim, dict):
                    die('In config.json, stimulus data must be a dictionary')
                stim['isPRFSynthData'] = False
                    
                # make a temporary file
                import tempfile
                (fl, stimjs_file) = tempfile.mkstemp(suffix='.json', text=True)
                print("[base/run.py] This is the temp file with stim info: ")
                print(stimjs_file)
                print("[base/run.py] This is the content: ")
                print(stim)
                with open(stimjs_file, 'w') as json_data:
                       json.dump(stim, json_data)
                       
                # okay, we have the files; run the solver script!
                try:
                    pid = os.fork()
                    if pid == 0:
                        os.execl("/solve.sh", "/solve.sh",
                                 opts_file, bold_image, stim_file, stimjs_file, outbids_dir)
                    else:
                        note("Beginning os.wait() for /solve.sh, run=%s (child pid is %s)" % (runid, pid))
                        os.wait()
                except Exception:
                    die("Failed to exec /solve.sh script!")    
    
                nii_base = nib.load(bold_image)
                # If there are things to cleanup we do that; specifically, the estimates.json file:
                estfl = os.path.join(outbids_dir, 'estimates.json')
                if os.path.isfile(estfl):
                    note("Processing estimates.json file...")
                    import nibabel as nib, numpy as np
                    with open(estfl, 'r') as fl:
                        dat = json.load(fl)
                    # decode the data...
                    dat = {k: np.asarray([u[k] for u in dat]) for k in dat[0].keys()}
                    for (k,v) in dat.items():
                        if len(v.shape) == 2:            
                            im = nib.Nifti2Image(np.reshape(v, (v.shape[0], 1, 1, v.shape[-1])),
                                                 nii_base.affine, nii_base.header)
                        else:
                            im = nib.Nifti2Image(np.reshape(v, (-1, 1, 1, 1)),
                                                 nii_base.affine, nii_base.header)
                        print("Writing the estimates.json to nifti2 in outbids_dir: {outbids_dir}")
                        im.to_filename(os.path.join(outbids_dir, f'sub-{sub}_ses-{ses}_task-{task}_run-{runid}_desc-{area}_{k.lower()}.nii.gz'))
                    os.rename(estfl, os.path.join(outbids_dir, f'sub-{sub}_ses-{ses}_task-{task}_run-{runid}_desc-{area}_estimates.json'))
                else:
                    note("No estimates.json file found.")
                resfli = os.path.join(outbids_dir, 'results.mat')
                resflo = os.path.join(outbids_dir, f'sub-{sub}_ses-{ses}_task-{task}_run-{runid}_desc-{area}_results.mat')
                if os.path.isfile(resfli): os.rename(resfli, resflo)
                processed += 1
        
        
if processed == 0: die("No BOLD images found!")

# exit happily
sys.exit(0)