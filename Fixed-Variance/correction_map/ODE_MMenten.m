%% Full Simple Model for Single Substrate Catalyzed Reactions (M)
function ydot=ODE_MMenten(t,y,p)

k_1 = p(1);
k_1r = p(2);
k_2 = p(3);

% dE/dt = k_1r ES - k1 E S + k2 ES
E  = k_1r * y(3)  -  k_1 * y(1) * y (2)  +  k_2 * y(3);

% dS/dt = k_1r ES - k_1 E S
S  = k_1r * y(3)  -  k_1 * y(1) * y (2);

% dES/dt = k_1 E S - k_2 ES - k_1r ES
ES = k_1 * y(1) * y (2)  -  k_2 * y(3)  -  k_1r * y(3);

% dP/dt = k_2 ES
P  = k_2 * y(3);

ydot = [E; S; ES; P];
