function [prm,fields]=getDefualtTemporalParams(temp_type)

c = Constants.getTemporalParams.temporalParams;
temp_type = strrep(temp_type,'_','-');

for i = 1:length(c)
    if strcmp(c{i}.type, temp_type)
        idx = i;
    end
end
fields = c{idx}.fields;
prm = c{idx}.prm;

end