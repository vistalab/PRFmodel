
function str = toSetField(str, fieldName, fieldVal)
% function str = toSetField(str, fieldName, fieldVal)

for k = 1 : length(fieldName)
    str = setfield(str, fieldName{k}, fieldVal(k));
end
end