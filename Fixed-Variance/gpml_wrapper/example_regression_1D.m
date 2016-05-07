%% necessary to use the gpml library
addpath('gpml-matlab-v3.6-2015-07-07');
startup;



%% create data
n = 40;   % training points
m = 200;  % test points


% training inputs & outputs
X = linspace(-5, 5, n)';
y = sin(X) + randn(size(X)) * 0.1;

% test inputs
Xtest = linspace(-5, 5, m)';



%% GP regression
amplitude = 1;
lengthscale = 1;
sigma2 = 0.02;


% if amplitude & lengthscale are not set,
% then they will be internally optimised (recommended)
[mu, var, lb, ub] = regressionGP(X, y, Xtest, sigma2, amplitude, lengthscale);
% [mu, var, lb, ub, amplitude, lengthscale] = regressionGP(X, y, Xtest, sigma2);

[mu2, var2, lb2, ub2] = regressionGP(X, y, Xtest, sigma2);


%% Plotting
figure;
plot(X, y, 'o'); hold on;
plot(Xtest, mu, 'r-');
plot(Xtest, lb, 'c-');
plot(Xtest, ub, 'c-');


plot(Xtest, mu2, 'k-');
plot(Xtest, lb2, 'c-');
plot(Xtest, ub2, 'c-');

%hold off;
