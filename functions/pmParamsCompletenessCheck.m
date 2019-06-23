function checkedParams = pmParamsCompletenessCheck(passed,defaults)
% All params are defined in the defaults, and if no params are passed, the
% defaults params will be used. BUT, if the params that are passed are not
% complete, then we need to check if all the required values are there or not.
% If not, they will be filled with the values in defaults.
% If the field exist but it is empty, it will be edited with the default. 


    requiredFields  = fieldnames(defaults);
    checkedParams   = passed;
    for nf=1:length(requiredFields)
        if ~isfield(checkedParams, requiredFields{nf})
            checkedParams.(requiredFields{nf}) = defaults.(requiredFields{nf});
        end
        if isempty(checkedParams.(requiredFields{nf}))
            checkedParams.(requiredFields{nf}) = defaults.(requiredFields{nf});
        end
    end
end

