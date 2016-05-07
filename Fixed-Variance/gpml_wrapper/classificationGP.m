function [p, lb, ub, amplitude, lengthscales] = classificationGP (X, y, Xtest, amplitude, lengthscales)
% CLASSIFICATIONGP Classification wrapper for the gpml library.
% [p, lb, ub, amplitude, lengthscales] = CLASSIFICATIONGP(X, y, Xtest, amplitude, lengthscales)
% [p, lb, ub, amplitude, lengthscales] = CLASSIFICATIONGP(X, y, Xtest)
%
%	X:             the training inputs  (NxD matrix: N instances, D dimensions)
%	y:             the training outputs (Nx1 matrix with entries: '+1' or '-1')
%	Xtest:         the test outputs     (MxD matrix: M instances, D dimensions)
%	amplitude:     the amplitude parameter of the RBF kernel
%	lengthscales:  the lengthscale parameters of the RBF kernel (1xD vector)
%	               If scalar, then one lengthscale will be used for all dimensions.
%
%	p:   Mx1 vector of class probabilities for the test points. 
%            An instance is classified:
%            as '+1' with probability 'p' and 
%            as '-1' with probability '1-p'.
%	lb:  the estimated lower 95% quantile for the values in 'p'
%	ub:  the estimated upper 95% quantile for the values in 'p'
%	amplitude:     the amplitude parameter used (may be optimised)
%	lengthscales:  the lengthscale parameters used (may be optimised)
%
%	Note on hyperparameter optimisation:
%	If 'amplitude' or 'lengthscales' are not provided, 
%	then these will be set automatically (recommended)
%
% Dimitrios Milios (dmilios@inf.ed.ac.uk)

if nargin < 3
	error ('Required arguments: X, y, Xtest');
end

N = size(X, 1);
M = size(Xtest, 1);
D = size(X, 2);

if size(y, 1) ~= N | ~isvector(y)
	error('y must be a Nx1 matrix with entries: +1 or -1 (N instances)');
end
for i = 1:N
	if y(i) ~= 1 & y(i) ~= -1
		error('y must be a Nx1 matrix with entries: +1 or -1 (N instances)');
	end
end
if size(Xtest, 2) ~= D
	error('Xtest must have the same number of columns as X');
end


optimisation = false;
if nargin < 5
	lengthscales = ones(1, D);
	optimisation = true;
end
if nargin < 4
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
likfunc = @likErf;
inferenceFunc = @infEP;


if optimisation
	hyp = minimize(hyp, @gp, -40, inferenceFunc, meanfunc, covfunc, likfunc, X, y);
end
amplitude = exp(hyp.cov(end));
lengthscales = exp(hyp.cov(1:D));


[ymu, ys2, fmu, fs2, lp] = gp(hyp, inferenceFunc, meanfunc, covfunc, likfunc, X, y, Xtest, ones(M,1));

p = normcdf(fmu ./ sqrt(1 + fs2)); % or equivalently:  exp(lp)
lb = normcdf((fmu - 2 * sqrt(fs2)) ./ sqrt(1 + fs2));
ub = normcdf((fmu + 2 * sqrt(fs2)) ./ sqrt(1 + fs2));

