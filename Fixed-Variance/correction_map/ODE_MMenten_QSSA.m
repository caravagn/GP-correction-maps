%% Enzyme Reactions with Irreversible Henri-Michaelis-Menten Kinetics (m)
function ydot=ODE_MMenten(t,y,p)

k_1 = p(1);
k_1r = p(2);
k_2 = p(3);

% Total enzyme Etot = E_0 + ES_0
Etot = y(1) + y(3);

% Condition for QSSA
% E_0 / (S_0+Km) << 1
Km = (k_1r + k_2)/k_1;
Vmax = k_2 * y(1);

% QSSA
E = 0;
ES = 0;

% dS/dt = - Vmax S/(Km + S)
S  = - Vmax*y(2)/(Km + y(2));

% dP/dt = + dS/dt
P  = Vmax*y(2)/(Km + y(2));

ydot = [E; S; ES; P];
