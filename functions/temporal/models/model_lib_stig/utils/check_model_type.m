function optimize_flag = check_model_type(model_type)
% Function that checks if the model type argument passed as input to 
% tchModel object is a valid model implemented in the code. 
% INPUT
%   model_type: model name passed to tchModel (string)
% 
% OUTPUT
%   Error if input does not match a valid_model string, otherwise:
%   optimize_flag: 1 if model supports nonlinear optimization and 0 if not
% 
% AS 10/2017

valid_models = {'1ch-glm' '1ch-lin' '1ch-exp' '1ch-balloon' '2ch-lin-htd' ...
    '1ch-pow' '1ch-div' '1ch-dcts' '1ch-rect' '1ch-quad' '1ch-sig'...
    '2ch-lin-rect' '2ch-exp-rect' '2ch-lin-quad' '2ch-exp-quad' ...
    '2ch-pow-quad' '2ch-pow-rect' '2ch-lin-sig' '2ch-exp-sig'};

if sum(strcmp(model_type, valid_models)) == 0
    error('Invalid model');
end

optimize_flag = 1;
if contains(model_type, {'glm' 'balloon'})
    optimize_flag = 0;
end

end
