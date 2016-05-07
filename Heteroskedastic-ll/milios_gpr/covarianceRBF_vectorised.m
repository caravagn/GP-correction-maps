% this vectorised form of the covariance matrix calculation is much faster;
% comparisons are more fair this way
function [ cov_matrix ] = covarianceRBF_vectorised( X1, X2, amplitude2, invLengthscale2 )
	a = X1';
	b = X2';
	if (size(a,1) == 1)
		a = [a; zeros(1,size(a,2))];  
		b = [b; zeros(1,size(b,2))];  
	end 
	aa=sum(a.*a);
	bb=sum(b.*b);
	ab=a'*b;
	dist = sqrt(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab); 
	dist = real(dist);
	cov_matrix = amplitude2 * exp(-0.5 * dist.^2 * invLengthscale2);
end
