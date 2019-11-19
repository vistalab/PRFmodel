import ctypes
import numpy as np
import popeye.og_hrf as og
import popeye.utilities as utils
import json, os, sys, six, nibabel as nib, pimms
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

# some other options we can extract
Ns = opts.get('grid_density', 3)

# Walk through each voxel/stim description
bold = np.reshape(np.asarray(bold_im.dataobj), (-1, bold_im.shape[-1]))
stim = np.squeeze(np.asarray(stim_im.dataobj))
if len(stim_json) != bold.shape[0]:
    raise ValueError('BOLD Image and Stimulus JSON do not have the same number of data points')
fields = ('theta', 'rho', 'sigma', 'hrf_delay', 'beta', 'baseline')
res = {k:[] for k in fields}
for (ii, vx,js) in zip(range(len(bold)), bold, stim_json):
    stdat = js['Stimulus']
    if pimms.is_list(stdat): stdat = stdat[0]
    height = stdat['fieldofviewVert']
    width = stdat['fieldofviewHorz']
    ### STIMULUS
    # First get a viewing distance and screen size
    dist = 100 # 100 cm is arbitrary
    stim_width = 2 * dist * np.tan(np.pi/180 * width/2)
    stimulus = VisualStimulus(stim, dist, stim_width, 1.0, float(js['TR']), ctypes.c_int16)
    model = og.GaussianModel(stimulus, utils.double_gamma_hrf)
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
    grids = (x_grid, y_grid, s_grid, h_grid)
    bounds = (x_bound, y_bound, s_bound, h_bound, b_bound, u_bound,)
    ## fit the response
    # auto_fit = True fits the model on assignment
    # verbose = 0 is silent
    # verbose = 1 is a single print
    # verbose = 2 is very verbose
    fit = og.GaussianFit(model, vx, grids, bounds, Ns=Ns,
                         voxel_index=(ii, 1, 1), auto_fit=True, verbose=2)
    for (k,v) in zip(fields, fit.overloaded_estimate):
        res[k].append(v)

# Export the files
for (k,v) in six.iteritems(res):
    im = nib.Nifti1Image(np.reshape(v, bold_im.shape[:-1]), bold_im.affine)
    im.to_filename(os.path.join(outdir, k + '.nii.gz'))

# That's it!
print("Popeye finished succesfully.")
sys.exit(0)

    
    
