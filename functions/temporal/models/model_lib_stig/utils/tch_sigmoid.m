function Y = tch_sigmoid(X, lambda_p, kappa_p, lambda_n, kappa_n)

if nargin < 4 || isempty(lambda_n); lambda_n = lambda_p; end
if nargin < 5 || isempty(kappa_n); kappa_n = kappa_p; end
X = rectify(X, 'abs', .001);
X_p = X; X_p(X < 0) = 0;
X_n = X; X_n(X_n > 0) = 0;

weibull_p = 1 - exp(-(X_p ./ lambda_p) .^ kappa_p);
weibull_n = 1 - exp(-(-X_n ./ lambda_n) .^ kappa_n);
Y = weibull_p + weibull_n;

end
