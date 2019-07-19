function fname = pmForwardModelToNifti(DT, varargin)
% Writes BOLDnoise into a nifti file. The size of the nifti file will be related
% to the size of the table
% 
%  Inputs: table with one or several rows of parameters to calculate bold series
%          in fwd model
% 
%  Outputs: returns the same table with the last column of pm-s calculated
% 
%  See also: forwardModelToNifti
% 
%  GLU Vistalab 2019.05


% Read the inputs
varargin = mrvParamFormat(varargin);
p = inputParser;
p.addRequired ('DT', @istable);
p.addParameter('fname','BOLDnoise.nii.gz', @ischar);
p.addParameter('demean',false, @islogical);
p.parse(DT,varargin{:});

fname  = p.Results.fname;
demean = p.Results.demean;

%% Calculate
pm1    = DT.pm(1);
Ntime  = pm1.timePointsN;
Nvoxels = height(DT);

% We can do a flat nifti, but I think it id better to create just 1D.
% Otherwise I am not confident that reshape is returning the order I want
% dim1 = ceil(sqrt(Nvoxels));
% dim2 = ceil(Nvoxels/dim1);
dim1 = Nvoxels;
dim2 = 1;

% Create a cube of timeSeries
myCube = zeros(dim1,dim2,1,Ntime);
index = 0;
for ii=1:dim1
    % for jj=1:dim2
    jj = 1; 
        index = index + 1;
        if index <= Nvoxels
            tmp = DT.pm(index).BOLDnoise;
            if demean
                myCube(ii,jj,1,:) = (tmp-mean(tmp)) / mean(tmp);
            else
                myCube(ii,jj,1,:) = tmp;
            end
            assert(DT.pm(index).TR == pm1.TR, 'All BOLDnoise signals TRs should be equal')
        end
    % end
end
% Save it as a nifti
writeFileNifti(niftiCreate('data', myCube, 'tr', pm1.TR, ...
                            'fname',fname));

