function [U, invC_y] = train_gpRBF (X, y, amplitude, lengthscale, sigma2)
	trainingMean = mean(y);
	y = y - trainingMean;
	N = size(X, 1);     % training points
	
	amplitude2 = amplitude * amplitude;
	invLengthscale2 = 1 / (lengthscale * lengthscale);
	if isscalar(sigma2)
		sigma2_D = eye(N) * sigma2;
	else
		sigma2_D = diag(sigma2);
	end
	
	Knn = covarianceRBF_vectorised(X, X, amplitude2, invLengthscale2);
	C = Knn + sigma2_D;
	U = chol(C, 'upper');
	invC_y = U\(U'\y);
end
