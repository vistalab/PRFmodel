function [ons, offs, dur, gap]=st_codestim(stiminput,fs)
% 
% if length(stiminput) > 1000 % this is rather an arbitary decision, if there is more than 1000 TRs ...                   
%      fs = 1000;             % the file is assumed to be in milliseconds         
% else
%      fs = 1;
% end

if notDefined('fs')
    fs = 1000;  
end
    
   
   f = find(diff([0,stiminput',0]~=0));
   ons = (f(1:2:end-1)) ./fs;  % Start indices
   offs = (f(2:2:end)-1) ./fs;  % Consecutive oneset™ counts
   dur = (f(2:2:end)-f(1:2:end-1)) ./fs;  % Consecutive onesâ€™ counts
   gap =  ons(2:end)- offs(1:end-1) ;


end
