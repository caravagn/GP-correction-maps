function [ amplitude, lengthscale, sigma2 ] = optimise_gpRBF( X, y, sigma2 )
% optimise_gpRBF Optimise the hyperparameters for a GP with isometric RBF kernel
% [ amplitude, lengthscale ] = optimise_gpRBF( X, y, sigma2 )
%
%	X:      the training set inputs 
% 		(i.e. NxD, where N is the number of instances and D the dimension)
%	y:      the training set outputs (i.e. Nx1)
%	sigma2: the variance parameter of the noise model (if not set, it will be optimised as well)
%
%	amplitude:   the optimised amplitude hyperparameter
%	lengthscale: the optimised amplitude lengthscale
%
% Dimitrios Milios (dmilios@inf.ed.ac.uk)

	optimiseNoise = false;
	if nargin < 3
		optimiseNoise = true;
	end

	trainingMean = mean(y);
	y = y - trainingMean;
	
	n = size(X, 1);
	nlog2pi = n * log(2 * pi) / 2;

	lower = [1e-7, 1e-7];
	upper = [1e7, 1e7];
	options = optimset('Display','off','Algorithm','active-set');

	initAmpl = max(y) - min(y);
	initLeng = mean(max(X) - min(X)) / 2;
	init = [initAmpl, initLeng];
	
	if optimiseNoise
		lower = [lower 1e-6];
		upper = [upper 1e3];
		init = [init initAmpl / 10];
	end
	
%  	init = log(init);
%  	lower = log(lower);
%  	upper = log(upper);
	[hbest,fbest,~,~,~,~,~] = fmincon(@negativeloglik, init, [],[],[],[],...
		lower,upper,[],options);
		
	for i = 1:10
		init = ([exprnd(initAmpl), exprnd(initLeng)]);
		if optimiseNoise
			init(3) = init(1) / 10;
		end
		
%  		init = log(init);
%  		lower = log(lower);
%  		upper = log(upper);
		[hval,fval,~,~,~,~,~] = fmincon(@negativeloglik, init, [],[],[],[],...
			lower,upper,[],options);
		if fval < fbest
			fbest = fval;
			hbest = hval;
		end
	end
	
%	hbest = exp(hbest);
	amplitude = hbest(1); 
	lengthscale = hbest(2); 
	if optimiseNoise
		sigma2 = hbest(3);
	end

	
	function [ minuslogl ] = negativeloglik( hyperparams )

%		hyperparams = exp(hyperparams);
	
		amplitude2 = hyperparams(1)^2;
		invLengthscale2 = 1 / hyperparams(2)^2;
		if optimiseNoise
			sigma2 = hyperparams(3);
		end
		
		if isscalar(sigma2)
			sigma2_D = eye(n) * sigma2;
		else
			sigma2_D = diag(sigma2);
		end
		Knn = covarianceRBF_vectorised(X, X, amplitude2, invLengthscale2);
		C = Knn + sigma2_D;
		[L, p] = chol(C,'lower');
		if p > 0
			disp('Non-P.D. covariance matrix during hyperparameter optimisation!')
			minuslogl = inf;
			return;
		end

		logdetC =  sum(log(diag(L))) * 2;
		invC_y = C\y;

		logl = - 0.5 * logdetC - 0.5 * y' * invC_y - nlog2pi;
		minuslogl = -logl;
	end
end
