function [gpMean, gpVar] = gpRBF(X, y, Xtest, amplitude, lengthscale, sigma2)
	[U, invC_y] = train_gpRBF (X, y, amplitude, lengthscale, sigma2);
	[gpMean, gpVar] = predict_gpRBF (U, invC_y, X, y, amplitude, lengthscale, Xtest);
end
