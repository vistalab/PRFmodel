%% check if ROI is there and remove ROIs that does not have any voxels

function newtargetROIs = removeEmptyROIs(targetROIs)

    vw = initHiddenGray;
    roiFiles = getAllFiles('./3DAnatomy/ROIs/' ,'*_toon.mat',1);
    vw=loadROI(vw,roiFiles,[],[],1,0);

    for ii = 1:length(roiFiles)
        coords = vw.ROIs(ii).coords;
        [coordsIndex, coords] = roiIndices(vw, coords);
        if isempty(coordsIndex) || length(coordsIndex) < 2
            emptyROI{ii} = vw.ROIs(ii).name;
        else 
             emptyROI{ii} =[];
        end
        
        [~,name]=fileparts(roiFiles{ii});
        names{ii}=name;
    end

    emptyROI = emptyROI(~cellfun('isempty',emptyROI));

    Acommon = intersect(emptyROI,targetROIs);
    newtargetROIs = setxor(targetROIs,Acommon);

    % get actual files only
    newtargetROIs = intersect(names,newtargetROIs);


    
    disp(['Skipping ROIs: [' Acommon ']']) ;





end