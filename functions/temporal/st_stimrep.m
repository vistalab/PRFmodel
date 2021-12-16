%% returns how many times unique stimulus frame is repeated


function [reps, nrep,blankDur]= st_stimrep(stim)

flat_img = reshape(stim,size(stim,1)*size(stim,2),[]);
[uniqueImg]= unique(flat_img.','rows','stable').';
uniqueImg = reshape(uniqueImg,size(stim,1),size(stim,2),[]);
blankimg = uniqueImg(:,:,1);
% % search blanks
% blankimg= logical(zeros(size(uniqueImg(:,:,1))));
% counter = 1;
% for ii = 1:size(uniqueImg,3)
%        if isequal(blankimg, uniqueImg(:,:,ii))
%            blankidx(counter) = ii ;
%            counter = counter+1;
%        end
% end

uniqueImg=uniqueImg(:,:,2:end);

for uu =1:size(uniqueImg,3)
    idx=[];
    for ff = 1:size(stim,3)
        if  isequal(stim(:,:,ff),uniqueImg(:,:,uu))
            idx = [idx ff]; % uu is unique img, ff is index in stim
        end
    end
    reps{uu} = idx;
end


% clean up if repetition happens but in reality it is same stimulus,
% sweeping in different dirrection
% only accounts for doulbe direction for now
nreps =length(reps);
for rr = 1:nreps
    if max(diff(sort(reps{rr}))) ~= 1
        
        newseqidx = find(diff(sort(reps{rr}))~=1);
        warning('[st_stimconvert] Potential error in stimulus repitition')

        reps{end+1} = reps{rr}(newseqidx+1:end);
        reps{rr} = reps{rr}(1:newseqidx);


    end
    
    
end

nrep = length(reps{uu});


blankIDX=[];
for ff = 1:size(stim,3)
    if  isequal(stim(:,:,ff),blankimg(:,:,1))
    blankIDX(ff) = 1;
    end

end
blankIDX = find(blankIDX);
blankIDX = find(diff(sort(blankIDX))~=1);
blankIDX = blankIDX(2:end); % remove initial blank just incase
blankDur= diff(blankIDX);
if max(diff(blankDur)) == 0
    blankDur = mean(blankDur);
else
    error("cannot deal with variable blank durations")
end

warning(sprintf('[st_stimconvert] Does it have %d unqiue stims?',length(reps)))
warning(sprintf('[st_stimconvert] repeated %d times?',nrep))


end