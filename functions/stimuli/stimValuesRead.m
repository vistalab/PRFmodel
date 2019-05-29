function stimValues = stimValuesRead(input)
% Check if it is a string, if it is reads the values and return values
% 
%  Inputs: path to a file or the values themselves
% 
%  Outputs: values
% 
% 

    if isstring(input) || ischar(input)
        s          = load(input);
        stimValues = s.stim;
    else
        stimValues = input;
    end


end

