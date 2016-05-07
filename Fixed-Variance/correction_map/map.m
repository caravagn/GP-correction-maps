%% necessary to use the gpml library
addpath('../gpml_wrapper');
addpath('../gpml_wrapper/gpml-matlab-v3.6-2015-07-07');
addpath('../../Heteroskedastic/milios_gpr');
addpath('../../Heteroskedastic/plotting_extras');
startup;

%% create data
NUM_TRAINING_POINTS = 40;   % training points
NUM_TEST_POINTS = 200;  % test points

%% E(0) range [1, 100]
PAR_LEFT = 1;
PAR_RIGHT = 100;

% training inputs & outputs
TRAINING_POINTS = linspace(PAR_LEFT, PAR_RIGHT, NUM_TRAINING_POINTS)';
fprintf('TRAINING POINTS: ')

for m = 1:NUM_TRAINING_POINTS
	fprintf('%f ', TRAINING_POINTS(m))
end

fprintf('\nCREATING TRAINING SET\n')
TRAINING_POINTS_OBS = zeros(NUM_TRAINING_POINTS, 1);
TRAINING_POINTS_VALIDITY = zeros(NUM_TRAINING_POINTS, 1);
for m = 1:NUM_TRAINING_POINTS
	 TRAINING_POINTS_OBS(m,:) = model_metric(TRAINING_POINTS(m));
	 TRAINING_POINTS_VALIDITY(m,:) = MM_validity(TRAINING_POINTS(m));	

end

for m = 1:NUM_TRAINING_POINTS
	fprintf('%f ', TRAINING_POINTS_OBS(m,:))
end

% test inputs
TEST_POINTS = linspace(PAR_LEFT, PAR_RIGHT, NUM_TEST_POINTS)';

%% GP regression
sigma2 = 1;


[ amplitude, lengthscale ] = optimise_gpRBF( TRAINING_POINTS, TRAINING_POINTS_OBS, sigma2 )


% % if amplitude & lengthscale are not set,
% % then they will be internally optimised (recommended)
% [mu, var, lb, ub] = regressionGP(X, y, Xtest, sigma2, amplitude, lengthscale);
% % [mu, var, lb, ub, amplitude, lengthscale] = regressionGP(X, y, Xtest, sigma2);

[mu, variance, lb, ub] = regressionGP(TRAINING_POINTS, TRAINING_POINTS_OBS, TEST_POINTS, sigma2, amplitude, lengthscale);

 
% %% Plotting
%figure;
% g = subplot(2,1,1);

subplot(4, 1, [1 2 3]);

shadedErrorBar(TEST_POINTS, mu, [(ub-mu)'; (mu-lb)'], {'r'}, 1); hold on;
plot(TEST_POINTS, mu, 'r-', 'LineWidth', 2);
plot(TRAINING_POINTS, TRAINING_POINTS_OBS, 'o',  'MarkerFaceColor',[0 0 0], 'MarkerSize', 6); hold on;
set(gca,'xticklabel',[])


% plot(TEST_POINTS, mu, 'r-');
% plot(TEST_POINTS, lb, 'c-');
% plot(TEST_POINTS, ub, 'c-');
% title('Correction Map')
% xlabel('E0 - initial enzyme concentration');
% ylabel('Value');
ylim([-2 11])

% p = get(g,'position');
% p(4) = p(4)*1.30; 
% set(g, 'position', p);


% g = subplot(2,1,2) 
subplot(4, 1, 4);

plot(TRAINING_POINTS, TRAINING_POINTS_VALIDITY, 'o',  'MarkerFaceColor',[.4 .4 .4]+.3, 'MarkerSize', 4);
% title('QSSA validity (<< 1: valid)');
% xlabel('E0 - initial enzyme concentration');
ylim([0 1])
xlim([0 100])

% val_idx = find(abs(TRAINING_POINTS_VALIDITY - 0.4) < 0.1)
% val_idx = val_idx(1)
% line([TRAINING_POINTS(val_idx) TRAINING_POINTS(val_idx)], [0 TRAINING_POINTS_VALIDITY(val_idx)])
% rectangle('Position', [ 0 0 TRAINING_POINTS(val_idx) TRAINING_POINTS_VALIDITY(val_idx)])

% p = get(g,'position');
% p(4) = p(4)*.30; % Add 10 percent to height
% set(g, 'position', p);


hold off;
