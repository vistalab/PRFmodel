import ctypes
import numpy as np
import popeye.og_hrf as og
import popeye.og as og_nohrf
import popeye.utilities as utils
import sharedmem, multiprocessing, json, os, sys, six, nibabel as nib, pimms
from popeye.visual_stimulus import VisualStimulus

# We should have a json file, a BOLD file, and a stimulus file
# First, read in the options
(opts_file, bold_file, stim_file, stimjs_file, outdir) = sys.argv[1:]

# import our files...
with open(opts_file, 'r') as fl:
    opts = json.load(fl)
bold_im = nib.load(bold_file)
stim_im = nib.load(stim_file)
with open(stimjs_file, 'r') as fl:
    stim_json = json.load(fl)

# seed random number generator so we get the same answers ...
np.random.seed(opts.get('seed', 2764932))
# if there is an HRF search or not...
fixed_hrf = opts.get('fixed_hrf', False)
if fixed_hrf is True: fixed_hrf = 0.0

# some other options we can extract
Ns = opts.get('grid_density', 3)
mps = opts.get('multiprocess', True)

# some post-processing
if mps == 'auto' or mps is True: mps = multiprocessing.cpu_count()
elif mps == 0: mps = 1

# Walk through each voxel/stim description
bold = np.reshape(np.asarray(bold_im.dataobj), (-1, bold_im.shape[-1]))
stim = np.squeeze(np.asarray(stim_im.dataobj))
if len(stim_json) != bold.shape[0]:
    raise ValueError('BOLD Image and Stimulus JSON do not have the same number of data points')
if fixed_hrf is not False: fields = ('theta', 'rho', 'sigma', 'beta', 'baseline')
else: fields = ('theta', 'rho', 'sigma', 'hrfdelay', 'beta', 'baseline')
def fit_voxel(tup):
    (ii, vx, js) = tup
    stdat = js['Stimulus']
    if pimms.is_list(stdat): stdat = stdat[0]
    height = stdat['fieldofviewVert']
    width = stdat['fieldofviewHorz']
    ### STIMULUS
    # First get a viewing distance and screen size
    dist = 100 # 100 cm is arbitrary
    stim_width = 2 * dist * np.tan(np.pi/180 * width/2)
    stimulus = VisualStimulus(stim, dist, stim_width, 1.0, float(js['TR']), ctypes.c_int16)
    if fixed_hrf is not False:
        model = og_nohrf.GaussianModel(stimulus, utils.double_gamma_hrf)
        model.hrf_delay = fixed_hrf
    else: model = og.GaussianModel(stimulus, utils.double_gamma_hrf)
    ### FIT
    ## define search grids
    # these define min and max of the edge of the initial brute-force search.
    x_grid = (-width/2,width/2)
    y_grid = (-height/2,height/2)
    s_grid = (1/stimulus.ppd + 0.25, 5.25)
    h_grid = (-1.0, 1.0)
    ## define search bounds
    # these define the boundaries of the final gradient-descent search.
    x_bound = (-width, width)
    y_bound = (-height, height)
    s_bound = (1/stimulus.ppd, 12.0) # smallest sigma is a pixel
    b_bound = (1e-8,None)
    u_bound = (None,None)
    h_bound = (-3.0,3.0)
    ## package the grids and bounds
    if fixed_hrf is not False:
        grids = (x_grid, y_grid, s_grid)
        bounds = (x_bound, y_bound, s_bound, b_bound, u_bound,)
    else:
        grids = (x_grid, y_grid, s_grid, h_grid)
        bounds = (x_bound, y_bound, s_bound, h_bound, b_bound, u_bound,)
    ## fit the response
    # auto_fit = True fits the model on assignment
    # verbose = 0 is silent
    # verbose = 1 is a single print
    # verbose = 2 is very verbose
    if fixed_hrf is not False:
        fit = og_nohrf.GaussianFit(model, vx, grids, bounds, Ns=Ns,
                                   voxel_index=(ii, 1, 1), auto_fit=True, verbose=2)
    else:
        fit = og.GaussianFit(model, vx, grids, bounds, Ns=Ns,
                             voxel_index=(ii, 1, 1), auto_fit=True, verbose=2)
    return (ii, vx) + tuple(fit.overloaded_estimate) + (fit.prediction,)
if mps == 1:
    voxs = [fit_voxel((ii, vx, js)) for (ii, vx, js) in zip(range(len(bold)), bold, stim_json)]
else:
    tups = list(zip(range(len(bold)), bold, stim_json))
    with sharedmem.Pool(np=mps) as pool:
        voxs = pool.map(fit_voxel, tups)
    voxs = list(sorted(voxs, key=lambda tup:tup[0]))
# Update the results to match the x0/y0, sigma style used by prfanalyze
all_fields = ('index','voxel') + fields + ('pred',)
res = {k:np.asarray([u[ii] for u in voxs]) for (ii,k) in enumerate(all_fields)}
rr = {}
rr['x0'] = np.cos(res['theta'])  * res['rho']
rr['y0'] = -np.sin(res['theta']) * res['rho']
rr['sigmamajor'] = res['sigma']
rr['sigmaminor'] = res['sigma']
rr['beta'] = res['beta']
rr['baseline'] = res['baseline']
if fixed_hrf is False: rr['hrfdelay'] = res['hrfdelay']
# Export the files
for (k,v) in six.iteritems(rr):
    im = nib.Nifti1Image(np.reshape(v, bold_im.shape[:-1]), bold_im.affine)
    im.to_filename(os.path.join(outdir, k + '.nii.gz'))
# Also export the prediction and testdata
im = nib.Nifti1Image(np.reshape(bold, bold_im.shape), bold_im.affine)
im.to_filename(os.path.join(outdir, 'testdata.nii.gz'))
im = nib.Nifti1Image(np.reshape(res['pred'], bold_im.shape), bold_im.affine)
im.to_filename(os.path.join(outdir, 'modelpred.nii.gz'))

# That's it!
print("Popeye finished succesfully.")
sys.exit(0)

    
    
