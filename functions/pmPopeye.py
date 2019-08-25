# %% markdown
# #Initialization
# %
import os, sys, six, h5py, pint, warnings, time, pandas, copy
os.getcwd()
import numpy       as np
import scipy       as sp
import scipy.stats as stats
import nibabel     as nib
import pyrsistent  as pyr

# Graphics libraries:
import ipyvolume         as ipv
import matplotlib        as mpl
import matplotlib.pyplot as plt
import matplotlib.tri    as tri
import seaborn           as sns

# Import neuropythy and pimms (related data-structures library)
import pimms
import neuropythy as ny
import popeye


%gui qt
%matplotlib inline

# Additional matplotlib preferences:
font_data = {'family':'sans-serif',
             'sans-serif':['Helvetica Neue', 'Helvetica', 'Arial'],
             'size': 10,
             'weight': 'light'}
mpl.rc('font',**font_data)
# we want relatively high-res images, especially when saving to disk.
mpl.rcParams['figure.dpi'] = 72*2
mpl.rcParams['savefig.dpi'] = 72*4


# %% markdown
# #Load Data
# %
# This is for the original 4x4 nifti
nii1 = ny.load('~/toolboxes/PRFmodel/data/examples/synthDataExample2_TR2.nii.gz', to='image')
data = np.reshape(nii1.dataobj, [nii1.shape[0]*nii1.shape[1], nii1.shape[-1]])

# This is for the 1D nifti
# With the 1D nifti I thing we will save lots of orientation problems. For example,
# Freesurfers MRIread read the data in different order than niftiRead, flipping x and y.
# With a 1D, we do a squeeze and that's it.
# nii1       = ny.load('~/toolboxes/PRFmodel/data/examples/synthDataExample3_1D_TR2.nii.gz', to='image')
# data       = np.squeeze(nii1.dataobj)
input_data = data[4,:]  # Noah said he used the 5th voxel for the example

# Stimulus related (we passed it as a nifti, to avoid sending matlab specific things)
nii2 = ny.load('~/toolboxes/PRFmodel/data/examples/Exp-103_binary-true_size-20x20.nii.gz', to='image')
stim = np.squeeze(nii2.dataobj)
# stim = np.roll(stim, -5, axis=2)
stimulus = popeye.visual_stimulus.VisualStimulus(stim.astype('int16'), 30, 10.57962, 0.5, 2, 'int16')


# %% markdown
# #Fit the Models
# ## Gaussian Model
# %
import popeye.og_hrf as og_hrf
import popeye.og as og

hrf   = popeye.utilities.double_gamma_hrf
model = og_hrf.GaussianModel(stimulus, hrf)
# model.hrf_delay = -0.25
# RSQarray = []
# for ii in np.array([-0.5,-0.25, 1,0.25, 0.5]):
model.hrf_delay = -0.25
fit = og_hrf.GaussianFit(model, input_data,
                     ((-10, 10), (-10, 10), (0.25, 6.25), (-3.0, 3.0)),
                     ((-10, 10), (-10, 10), (0.10, 12.0), (-6.0, 6.0)),
                     auto_fit=True, Ns=5, verbose=1)
    # RSQarray.append((ii,fit.RSQ))


# %
# find and show the parameters (x, y, sigma, n, beta, baseline)
sol = fit.overloaded_estimate
prediction = fit.prediction
for (k,v) in zip(['theta', 'rho', 'sigma', 'hrf_delay', 'beta', 'baseline'], tuple(sol)):
    print('%-10s'%k,':\t',v)

# %
t = np.linspace(0, 300, 150)
(m,b) = (1,0)
(m,b) = sp.stats.linregress(prediction, input_data)[:2]
for (y,sty,nm) in zip([input_data, prediction*m + b], ['k.', 'r-'], ['data','model']):
    plt.plot(t, y, sty, label=nm)
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('BOLD Response (AU)')
plt.title('Gaussian Model Fit');

# % markdown
# ##CSS Model
# __Notes:__
#
#   1. Popeye's CSS model does not fit the HRF delay parameter, so I used the delay found in the above fit.
#   2. This delay displayed above is 9.79; popeye adds 5 to it after the fitting (not sure why?) so I've used 4.79 here.

# %
import popeye.css as css

hrf = popeye.utilities.double_gamma_hrf
model = css.CompressiveSpatialSummationModel(stimulus, hrf)
model.hrf_delay = 4.79
fit = css.CompressiveSpatialSummationFit(
    model, input_data,
    ((-10, 10), (-10, 10), (0.25, 6.25), (0.01,  0.99)),
    ((-10, 10), (-10, 10), (0.10, 12.0), (0.001, 0.999)),
    auto_fit=True, Ns=5, verbose=1)


# %
# find and show the parameters (x, y, sigma, n, beta, baseline)
sol = fit.overloaded_estimate
prediction = fit.prediction
for (k,v) in zip(['theta', 'rho', 'sigma', 'n', 'beta', 'baseline'], tuple(sol)):
    print('%-10s'%k,':\t',v)

# %
t = np.linspace(0, 300, 150)
(m,b) = (1,0)
(m,b) = sp.stats.linregress(prediction, input_data)[:2]
for (y,sty,nm) in zip([input_data, prediction*m + b], ['k.', 'r-'], ['data','model']):
    plt.plot(t, y, sty, label=nm)
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('BOLD Response (AU)')
plt.title('CSS Model Fit');
