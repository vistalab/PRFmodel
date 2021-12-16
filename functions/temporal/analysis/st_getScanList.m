function [trainSet, testSet]=st_getScanList(shuffled,wSplit)
load('ExpDesign.mat')

if notDefined('wSplit')
    wSplit = 1;
end

T = cell2table( ExpDesign ) ;
[C, ~, ic] = unique(T);


for ec = 1:height(C)
    if mod(ec,2) == 0
        scanList{ec} = find(ic == ec)';
    else
        scanList{ec} = find(ic == ec)';
    end
end

scanList=cell2mat(scanList');

seq = scanList(1:2:end,:);
shuf = scanList(2:2:end,:);
% 
% seq  = Shuffle(seq,2);
% shuf = Shuffle(shuf,2);
if wSplit == 1
    selectIdx =[1,2];
    selectIdx2 =3;
elseif wSplit == 2 
    selectIdx =[1,3];
    selectIdx2 =2;
elseif wSplit == 3
    selectIdx =[2,3];
    selectIdx2 =1;
elseif strcmp(wSplit,'avg')
    selectIdx = [1,2,3];
    selectIdx2 = [1,2,3];
end

seq_train = seq(:,selectIdx);
% seq_train = seq_train(:)';
seq_test = seq(:,selectIdx2);
% seq_test = seq_test(:)';

shuf_train = shuf(:,selectIdx);
% shuf_train = shuf_train(:)';
shuf_test = shuf(:,selectIdx2);
% shuf_test = shuf_test(:)';


if shuffled == 1 % shuffle
    trainSet = shuf_train';
    testSet = shuf_test';

else %seq
    trainSet = seq_train';
    testSet = seq_test';

end

