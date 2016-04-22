function [ amplitude, lengthscale ] = optimise_gpRBF( X, y, sigma2 )
	trainingMean = mean(y);
	y = y - trainingMean;
	
	n = size(X, 1);
	nlog2pi = n * log(2 * pi) / 2;

	lower = [1e-5; 1e-5];
	upper = [1e5, 1e5];
	options = optimset('Display','off','Algorithm','active-set');

	initAmpl = var(y);
	initLeng = mean(max(X) - min(X)) / 2;
	init = [initAmpl, initLeng];
	[hbest,fbest,~,~,~,~,~] = fmincon(@negativeloglik, init, [],[],[],[],...
		lower,upper,[],options);
		
	for i = 1:5
		init = [exprnd(initAmpl), exprnd(initLeng)];
		[hval,fval,~,~,~,~,~] = fmincon(@negativeloglik, init, [],[],[],[],...
			lower,upper,[],options);
		if fval < fbest
			fbest = fval;
			hbest = hval;
		end
	end
		
	amplitude = hbest(1); 
	lengthscale = hbest(2); 


	
	function [ minuslogl ] = negativeloglik( hyperparams )
		
		amplitude2 = hyperparams(1)^2;
		invLengthscale2 = 1 / hyperparams(2)^2;
		
		sigma2_D = eye(n) * sigma2;
		Knn = covarianceRBF_vectorised(X, X, amplitude2, invLengthscale2);
		C = Knn + sigma2_D;
		U = chol(C);

		logdetC =  sum(log(diag(U))) * 2;
		invC_y = U\(U'\y);

		logl = - 0.5 * logdetC - 0.5 * y' * invC_y - nlog2pi;
		minuslogl = -logl;
	end
end
