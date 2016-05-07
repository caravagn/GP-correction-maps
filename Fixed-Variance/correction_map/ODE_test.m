function main(PARAM)

%% Max simulation time
t_final = 5;

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

t_idx = find(abs(t-t_ast) < 0.001);
t_big_idx = find(abs(t_big-t_ast) < 0.001);

fprintf('\t Y = %s\n', num2str(y_big(t_big_idx, :)))
fprintf('\t y = %s\n', num2str(y(t_idx, :)))

r = norm(y_big(t_big_idx, [2 4]) - y(t_idx, [2 4]));
fprintf('\t norm S/P = %s\n', num2str(r))

r = norm(y_big(t_big_idx, [4]) - y(t_idx, [4]));
fprintf('\t norm P = %s\n', num2str(r))

 % fprintf('EXACT MODEL')
 % disp([t_big y_big])

color_E =  [228 26 28]/256
color_S = [55 126 184]/256
color_ES = [77 175 74]/256
color_P = [152 78 163]/256

subplot(2,1,1)
plot(t_big, y_big(:, 1),  'LineWidth', 2, 'Color', color_E);
hold on;
plot(t_big, y_big(:, 2), 'b', 'LineWidth', 2, 'Color', color_S);
plot(t_big, y_big(:, 3), 'r', 'LineWidth', 2, 'Color', color_ES);
plot(t_big, y_big(:, 4), 'c', 'LineWidth', 2, 'Color', color_P);
ylim([0 70])

txt = text(t_big(t_big_idx), y_big(t_big_idx, 2), num2str( ceil(y_big(t_big_idx, 2))) );
txt.FontSize = 8;
txt = text(t_big(t_big_idx), y_big(t_big_idx, 4), num2str( ceil(y_big(t_big_idx, 4))) );
txt.FontSize = 8;

hold off;

ylabel('Copy Numbers');
legend('E', 'S', 'ES', 'P', 'Location', 'North');

% E_0 / (S_0+Km) << 1
S_0 = y0(2);
valid = PARAM / (S_0+ (kinetics(2) + kinetics(3))/kinetics(1));
fprintf('QSSA Validity (<1): %f\n', valid)

 % fprintf('QSSA MODEL')
 % disp([t y])

subplot(2,1,2) 
plot(t, y(:, 1), 'LineWidth', 2, 'Color', color_E);
hold on;
plot(t, y(:, 2),  'LineWidth', 2, 'Color', color_S);
plot(t, y(:, 3), 'LineWidth', 2, 'Color', color_ES);
plot(t, y(:, 4),  'LineWidth', 2, 'Color', color_P);
ylim([0 70])

txt = text(t(t_idx), y(t_idx, 2), num2str( ceil(y(t_idx, 2)) ) );
txt.FontSize = 8; 
txt = text(t(t_idx), y(t_idx, 4), num2str( ceil(y(t_idx, 4)) ) );
txt.FontSize = 8;

hold off;

title(sprintf('QSSA with E(0) = %.0f - Validity (<1): %f', PARAM, valid), 'FontName', 'Arial')
xlabel('time');
ylabel('Copy Numbers');
legend('E', 'S', 'ES', 'P', 'Location', 'North');

% %plot(t, y(:, 3), 'c');
% %plot(t, y(:, 4), 'r');
% hold off;
% xlabel('time (minutes)');
% ylabel('Population');
% %legend('E', 'S', 'ES', 'P', 'Location', 'North');
% legend('E_m', 'E_M', 'Location', 'North');