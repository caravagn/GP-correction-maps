addpath('./milios_gpr');
addpath('./plotting_extras');

% If true, data is log-transformed
LOGTRANSFORM = true

% Regression statistics
NAME = 'averages.csv'

% Data generated in Python
TRAINING_POINTS = csvread(strcat('../Stochpy/results/ML_TRAINING_POINTS_', NAME))
TRAINING_POINTS_INPUT = csvread(strcat('../Stochpy/results/ML_TRAINING_POINTS_INPUT_', NAME))
TRAINING_POINTS_INPUT_VARIANCE = csvread(strcat('../Stochpy/results/ML_TRAINING_POINTS_VARIANCE_', NAME))


%% DEBUG
% REDUCTION = 25
% TRAINING_POINTS = TRAINING_POINTS(1:REDUCTION)
% TRAINING_POINTS_INPUT = TRAINING_POINTS_INPUT(1:REDUCTION)
% TRAINING_POINTS_INPUT_VARIANCE = TRAINING_POINTS_INPUT_VARIANCE(1:REDUCTION)


% Domains
PAR_LEFT = min(TRAINING_POINTS_INPUT)
PAR_RIGHT = max(TRAINING_POINTS_INPUT)

% Training size
NUM_TRAINING_POINTS = length(TRAINING_POINTS)

for m = 1:NUM_TRAINING_POINTS
 	fprintf('%f ', TRAINING_POINTS_INPUT(m))
 	fprintf('%f ', TRAINING_POINTS(m))
 	fprintf('%f \n', TRAINING_POINTS_INPUT_VARIANCE(m))
end


% Test size and points
NUM_TEST_POINTS = 200;  % test points
TEST_POINTS = linspace(PAR_LEFT, PAR_RIGHT, NUM_TEST_POINTS)';

% Options are
% VARIANCE = TRAINING_POINTS_INPUT_VARIANCE % Heteroschedastic regression
% VARIANCE = mean(TRAINING_POINTS_INPUT_VARIANCE) % Variance estimated across training points
% VARIANCE = 0.2 % Fixed variance


VARIANCE = 0.2

%%%%%%%%%%%%% GP regression

%% We estimate the hyperparameters from data
[ amplitude, lengthscale ] = optimise_gpRBF( TRAINING_POINTS_INPUT, TRAINING_POINTS, VARIANCE )
[ amplitude_h, lengthscale_h ] = optimise_gpRBF( TRAINING_POINTS_INPUT, TRAINING_POINTS, TRAINING_POINTS_INPUT_VARIANCE )

%% Regression
[gpMean, gpVar] = gpRBF(TRAINING_POINTS_INPUT, TRAINING_POINTS, TEST_POINTS, amplitude, lengthscale, VARIANCE);
[gpMean_h, gpVar_h] = gpRBF(TRAINING_POINTS_INPUT, TRAINING_POINTS, TEST_POINTS, amplitude_h, lengthscale_h, TRAINING_POINTS_INPUT_VARIANCE);


% Inverse map from log-transformation
if (LOGTRANSFORM)
	 gpMean = exp(gpMean);
	 gpMean_h = exp(gpMean_h);
	 TRAINING_POINTS = exp(TRAINING_POINTS);
end

% %% Plotting
figure;
%%subplot(2,1,1);
plot(TRAINING_POINTS_INPUT, TRAINING_POINTS, 'o'); hold on;
plot(TEST_POINTS, gpMean, 'r-');
plot(TEST_POINTS, gpMean_h, 'b-');

% Confidence intervals, inverse map from log-transformation
conf_int = sqrt(gpVar)*2
conf_int_h = sqrt(gpVar_h)*2

if (LOGTRANSFORM)
	plot(TEST_POINTS, exp(log(gpMean)+conf_int), 'c-');
	plot(TEST_POINTS, exp(log(gpMean)-conf_int), 'c-');
	plot(TEST_POINTS, exp(log(gpMean_h)+conf_int_h), 'g-');
	plot(TEST_POINTS, exp(log(gpMean_h)-conf_int_h), 'g-');

else
	plot(TEST_POINTS, gpMean+conf_int, 'c-');
	plot(TEST_POINTS, gpMean-conf_int, 'c-');
	plot(TEST_POINTS, gpMean_h+conf_int_h, 'g-');
	plot(TEST_POINTS, gpMean_h-conf_int_h, 'g-');
end



title('Correction Map for E[PR(t)]')
xlabel('Protein translation rate');
ylabel('Value');

hold off;

fig = figure;

ci_upper = [exp(log(gpMean)+conf_int)] - gpMean
ci_lower = gpMean - [exp(log(gpMean)-conf_int)]

ci_upper_h = [exp(log(gpMean_h)+conf_int_h)] - gpMean_h
ci_lower_h = gpMean_h - [exp(log(gpMean_h)-conf_int_h)]


h = shadedErrorBar(TEST_POINTS, gpMean_h, [ci_upper_h'; ci_lower_h'], {'r'}, 1); hold on;

h = plot(TEST_POINTS, gpMean_h, 'r-', 'LineWidth', 2);

h = plot(TEST_POINTS, gpMean, 'b-', 'LineWidth', 2);

if (LOGTRANSFORM)
	h = plot(TEST_POINTS, exp(log(gpMean)+conf_int), 'b-.');
	h = plot(TEST_POINTS, exp(log(gpMean)-conf_int), 'b-.');
else
	h = plot(TEST_POINTS, gpMean+conf_int, 'c-');
	h = plot(TEST_POINTS, gpMean-conf_int, 'c-');
	h = plot(TEST_POINTS, gpMean_h+conf_int_h, 'g-');
	ph = lot(TEST_POINTS, gpMean_h-conf_int_h, 'g-');
end

h = plot(TRAINING_POINTS_INPUT, TRAINING_POINTS, 'ko', 'MarkerFaceColor',[.49 1 .63], 'MarkerSize', 6); hold on;

hold off;
