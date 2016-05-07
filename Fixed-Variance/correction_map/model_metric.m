function r=model_metric(PARAM)

%% Max simulation time
t_final = 2;

%% property
t_ast = 1.5;

% parameters: 
% k1  = 2   1/(mole*second)
% k1r = 1   1/second
% k2  = 1.5 1/second
kinetics = [2 1 1.5];

%% Initial state -- E S ES P
y0 = [PARAM 60 0 0];

%% ODE solver parameters
vopt = odeset('InitialStep', 0.01, 'MaxStep', 0.01, 'RelTol', 0.01, 'AbsTol', 0.01);

% If I want to change a certain parameter (kinetics in this case)
fprintf('** PARAM E(0)=%s --> ', num2str(PARAM))

% Solve both models
fprintf('M || ')
[t_big, y_big] = ode45(@(t,y) ODE_MMenten(t,y, kinetics), [0 t_final], y0, vopt);

fprintf('m \n')
[t, y] = ode45(@(t,y) ODE_MMenten_QSSA(t,y, kinetics), [0 t_final], y0, vopt);

t_idx = find(abs(t-t_ast) < 0.01);
t_big_idx = find(abs(t_big-t_ast) < 0.01);

t_idx = t_idx(1);
t_big_idx = t_big_idx(1);

fprintf('\t Y = %s\n', num2str(y_big(t_big_idx, :)))
fprintf('\t y = %s\n', num2str(y(t_idx, :)))

r = norm(y_big(t_big_idx, [2 4]) - y(t_idx, [2 4]));
fprintf('\t norm S/P = %s\n', num2str(r))

 r = norm(y_big(t_big_idx, [4]) - y(t_idx, [4]));
 fprintf('\t norm P = %s\n', num2str(r))

%r = y_big(t_big_idx, 4) - y(t_idx, 4);




fprintf('\t r = %s\n', num2str(r))