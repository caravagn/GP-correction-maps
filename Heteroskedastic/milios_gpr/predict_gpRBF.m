function [gpMean, gpVar] = predict_gpRBF (U, invC_y, X, y, amplitude, lengthscale, Xtest)
	trainingMean = mean(y);
	y = y - trainingMean;
	
	amplitude2 = amplitude * amplitude;
	invLengthscale2 = 1 / (lengthscale * lengthscale);
	Kmm = covarianceRBF_vectorised(Xtest, Xtest, amplitude2, invLengthscale2);
	Kmn = covarianceRBF_vectorised(Xtest, X, amplitude2, invLengthscale2);
	Knm = Kmn';
	
	gpMean = Kmn * invC_y + trainingMean;
	if nargout > 1
		invC_Knm = U\(U'\Knm);
		gpVar = diag(Kmm) - sum(Kmn .* invC_Knm', 2);
	end
end
