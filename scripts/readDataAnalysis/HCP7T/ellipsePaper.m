%% Explain


%% Download the 7T data, just 3 subjects, 
% Install neurophyty and other software
%{
    conda activate base
    conda install s3fs -c conda-forge
    conda install h5py
    pip install --user cpython
    pip install --user neuropythy
    export SUBJECTS_DIR=$PWD/fs
    export HCP_SUBJECTS_DIR=$PWD/hcp
    export NPYTHY_DATA_CACHE_ROOT=~/tsemp/npythy_cache
    export HCP_CREDENTIALS=~/.hcp-passwd
    export HCP_AUTO_DOWNLOAD=true
    jupyter console
%}

% In python
%{ 
import neuropythy as ny
import nibabel as nib
ny.hcp.download(536647)
ny.hcp.download(115017)
ny.hcp.download(164131)
%}

% Download time series from NYU. (1) vpn to NYU, and to local, (2) vpn to
% Stanford and upload to black

