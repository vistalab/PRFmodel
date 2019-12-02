function [fname,data] = pmForwardModelToNifti(DT, varargin)
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
p.addRequired ('DTDT'        , @istable);
p.addParameter('fname'     , 'BOLDnoise.nii.gz', @ischar);
p.addParameter('demean'    , false , @islogical);
p.addParameter('writenifti', true  , @islogical);
p.addParameter('maxrows'   , 32000 , @isnumeric);
p.parse(DT,varargin{:});

fname      = p.Results.fname;
demean     = p.Results.demean;
writenifti = p.Results.writenifti;
maxRows    = p.Results.maxrows;

% Check size
if height(DT) >= 32000^2
    error('File too large, maximum accepted is 1024e6')
end

%% Calculate
%% Check number of rows and if to big, make it 2D
% Nifti-s limit per dimension is about 32767, default set to 32000
    pm1     = DT.pm(1);
    Ntime   = pm1.timePointsN;
    Nvoxels = height(DT);
    
    % We can do a flat nifti, but I think it id better to create just 1D.
    % Otherwise I am not confident that reshape is returning the order I want
    % If the file is larger than 32000, we need to use other dims, but we
    % will fill first the first dimension and keep from there
    
    if Nvoxels < maxRows
        dim1 = Nvoxels;
        dim2 = 1;
        dim3 = 1;
    else
        dim1   = maxRows;
        dim2   = ceil(Nvoxels / maxRows);
        dim3   = 1;
    end
    
    % Create a cube of timeSeries
    data = zeros(dim1,dim2,1,Ntime);
    index = 0;
    for jj=1:dim2
        for ii=1:dim1
            index = index + 1;
            if index <= Nvoxels
                tmp = DT.pm(index).BOLDnoise;
                if demean
                    data(ii,jj,dim3,:) = (tmp-mean(tmp)) / mean(tmp);
                else
                    data(ii,jj,dim3,:) = tmp;
                end
                assert(DT.pm(index).TR == pm1.TR, 'All BOLDnoise signals TRs should be equal')
            end
        end
    end
    if writenifti
        % Save it as a nifti
        writeFileNifti(niftiCreate('data', data, 'tr', pm1.TR, ...
            'fname',fname));
    end
end 
