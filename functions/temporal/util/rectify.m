function X_rect = rectify(X, polarity, thresh)
% Rectify data (X) in a given range determined by polarity and threshold.
% 
% INPUTS
%   1) X: numerical array
%   2) polarity: range of rectification
%       'positive' -- zeros numbers less than thresh (default)
%       'negative' -- zeros numbers greater than thresh
%       'abs' -- zeros numbers less than abs(thresh)
%   3) thresh: threshold for rectification
% 
% OUTPUT
%   X_rect: rectified version of input data X
% 
% AS 9/2017

if nargin < 2 || isempty(polarity); polarity = 'positive'; end
if nargin < 3 || isempty(thresh); thresh = 0; end
X_rect = X;
switch lower(polarity)
    case 'positive'
        X_rect(X < thresh) = 0;
    case 'negative'
        X_rect(X > thresh) = 0;
    case 'abs'
        X_rect(abs(X) < thresh) = 0;
    otherwise
        error('Polarity input must be positive, negative, or abs');
end

end
