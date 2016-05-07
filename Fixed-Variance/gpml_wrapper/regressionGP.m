function [mu, var, lb, ub, amplitude, lengthscales] = regressionGP (X, y, Xtest, sigma2, amplitude, lengthscales)
% REGRESSIONGP Classification wrapper for the gpml library.
% [mu, lb, ub, amplitude, lengthscales] = REGRESSIONGP(X, y, Xtest, amplitude, lengthscales)
% [mu, lb, ub, amplitude, lengthscales] = REGRESSIONGP(X, y, Xtest)
%
%	X:             the training inputs  (NxD matrix: N instances, D dimensions)
%	y:             the training outputs (Nx1 matrix)
%	Xtest:         the test outputs     (MxD matrix: M instances, D dimensions)
%	amplitude:     the amplitude parameter of the RBF kernel
%	lengthscales:  the lengthscale parameters of the RBF kernel (1xD vector)
%	               If scalar, then one lengthscale will be used for all dimensions.
%
%	mu:  Mx1 vector of GP posterior mean
%	var: Mx1 vector of GP posterior variance
%	lb:  the estimated lower 95% quantile for the values in 'mu'
%	ub:  the estimated upper 95% quantile for the values in 'mu'
%	amplitude:     the amplitude parameter used (may be optimised)
%	lengthscales:  the lengthscale parameters used (may be optimised)
%
%	Note on hyperparameter optimisation:
%	If 'amplitude' or 'lengthscales' are not provided, 
%	then these will be set automatically (recommended)
%
% Dimitrios Milios (dmilios@inf.ed.ac.uk)

if nargin < 4
	error ('Required arguments: X, y, Xtest, sigma2');
end

N = size(X, 1);
M = size(Xtest, 1);
D = size(X, 2);

if size(y, 1) ~= N | ~isvector(y)
	error('y must be a Nx1 matrix');
end
if size(Xtest, 2) ~= D
	error('Xtest must have the same number of columns as X');
end


optimisation = false;
if nargin < 6
	lengthscales = ones(1, D);
	optimisation = true;
end
if nargin < 5
	amplitude = 1;
	optimisation = true;
end

if ~isscalar(amplitude) | amplitude <= 0
	error('amplitude has to be a positive scalar');
end
if isscalar(lengthscales)
	lengthscales = ones(1, D) .* lengthscales;
end
if ~isvector(lengthscales) | sum(lengthscales <= 0) > 0 | size(lengthscales, 1) ~= 1 | size(lengthscales, 2) ~= D
	error('lengthscales must be either a positive scalar, or a 1xD vector of positive values');
end



meanfunc = @meanConst; hyp.mean = 0;
covfunc = @covSEard;
hyp.cov = log([lengthscales amplitude]);
likfunc = @likGauss;
hyp.lik  = log(sigma2) / 2;
inferenceFunc = @infExact;


if optimisation
	hyp = minimize(hyp, @gp, -40, inferenceFunc, meanfunc, covfunc, likfunc, X, y);
end
amplitude = exp(hyp.cov(end));
lengthscales = exp(hyp.cov(1:D));


[ymu, ys2, fmu, fs2, lp] = gp(hyp, inferenceFunc, meanfunc, covfunc, likfunc, X, y, Xtest, ones(M,1));

mu = fmu;
var = fs2;
lb = mu - 2 * sqrt(var);
ub = mu + 2 * sqrt(var);

