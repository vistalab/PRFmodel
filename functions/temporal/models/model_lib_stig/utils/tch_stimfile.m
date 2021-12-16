function [on, off, c, ims, ton, toff, tc, rd, cl] = tch_stimfile(stimfile)
% Compiles stimulus timing parameters for a given run of fMRI data.
%
% INPUT
%   stimfile: path to stimulus timing parameter file (*.txt)
%
% OUTPUTS
%   1) on: vector of image onset times (s)
%   2)off: vector of image offset times (s)
%   3) c: array of condition labels (one label per image)
%   4) ton: vector of trial onset times (s)
%   5) toff: vector of trial offset times (s)
%   6) tc: array of condition labels (one label per trial)
%   7) rd: run duration (s)
%   8) cl: array of all condition/trial types in experiment
%
% AS 2/2017

% check for file and open if found
if ~exist(stimfile, 'file')
    [on, off, c, ims, ton, toff, tc, rd, cl] = deal([]);
else
    [fid_par, message] = fopen(stimfile, 'r');
    % get experiment name and conditions
    ln = fgetl(fid_par);
    prts = strsplit(ln, ' ');
    exp = prts(1);
    cl = prts(3:end);
    cl = strrep(cl, ',', '');
    cl(strcmp('', cl)) = [];
    cl=sort(cl);
    % get run duration
    ln = fgetl(fid_par);
    prts = strsplit(ln, ' ');
    prts(strcmp('', prts)) = [];
    rd = str2num(prts{end});
    ln = fgetl(fid_par);
    ln = fgetl(fid_par);
    % code parameters of each stimulus
    itrials = []; c = {}; on = []; off = []; ims = {}; fcnt = 0;
    while ~feof(fid_par)
        ln = fgetl(fid_par);
        if isempty(ln) || isempty(findstr(sprintf(' '), ln)), break; end
        ln(ln==sprintf('\n')) = '';
        prts = strsplit(ln, ' ');
        fcnt = fcnt + 1;
        itrials(fcnt) = str2double(prts{1});
        c{fcnt} = prts{2};
        on(fcnt) = str2double(prts{3});
        off(fcnt) = on(fcnt) + str2double(prts{4});
        ims{fcnt} = strtok(prts{5}, '-');
    end
    fclose(fid_par);
    % code parameters of each trial
    ton = []; toff = []; tc = {}; bcnt = 0;
    for bb = 1:max(itrials)
        ii = find(itrials == bb, 1, 'first');
        zz = find(itrials == bb, 1, 'last');
        bcnt = bcnt + 1;
        ton(bcnt) = on(ii);
        toff(bcnt) = off(zz);
        tc{bcnt} = c{ii};
    end
    
    
end

end
