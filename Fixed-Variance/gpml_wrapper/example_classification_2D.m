%% necessary to use the gpml library
addpath('gpml-matlab-v3.6-2015-07-07');
startup;



%% create data
n = 100;   % training points
m = 1600;  % test points

% training inputs & outputs
X = rand(n, 2) * 10;
y = (sin(sum(X, 2)) > 0) .* 2 - 1;

% test inputs
Xtest1 = linspace(1, 10, sqrt(m))';
Xtest2 = linspace(1, 10, sqrt(m))';
[Xtest1, Xtest2] = meshgrid(Xtest1, Xtest2);
Xtest = [reshape(Xtest1, m, 1), reshape(Xtest2, m, 1)];



%% GP probit regression (aka classification)
amplitude = 1;
lengthscale = 1;

% if amplitude & lengthscale are not set,
% then they will be internally optimised (recommended)
[p, lb, ub] = classificationGP(X, y, Xtest, amplitude, lengthscale);
% [p, lb, ub, amplitude, lengthscale] = classificationGP(X, y, Xtest);



%% Plotting
p_plot = reshape(p, sqrt(m), sqrt(m));
lb_plot = reshape(lb, sqrt(m), sqrt(m));
ub_plot = reshape(ub, sqrt(m), sqrt(m));

figure;
scatter3 (X(:, 1), X(:, 2), (y + 1) / 2, 'o');
hold on;
surf(Xtest1, Xtest2, p_plot);
title ('Class probabilities')
hold off;

figure;
scatter3 (X(:, 1), X(:, 2), (y + 1) / 2, 'o');
hold on;
surf(Xtest1, Xtest2, lb_plot);
title ('95% lower confidence interval for class probabilities')
hold off;

figure;
scatter3 (X(:, 1), X(:, 2), (y + 1) / 2, 'o');
hold on;
surf(Xtest1, Xtest2, ub_plot);
title ('95% upper confidence interval for class probabilities')
hold off;
