function output_cell = force_cell(input_var)
% Force input variable into cell array if not already.

if ~iscell(input_var)
    output_cell = {input_var};
else
    output_cell = input_var;
end


end

