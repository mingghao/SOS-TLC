clear; clc;

syms x z theta dx dz dtheta F dF real    
syms x_d z_d r delta real                      
m = Parameter.m;                                   
g_acc = Parameter.g;                      
I = Parameter.I;

X = [x; z; theta; dx; dz; dtheta; F; dF];


f = [dx; 
     dz; 
     dtheta; 
     (F/m) * sin(theta); 
     (F/m) * cos(theta) - g_acc; 
     0;    
     dF;    
     0];

g1 = [0; 0; 0; 0; 0; 0; 0; 1];
g2 = [0; 0; 0; 0; 0; -1/I; 0; 0];
G = [g1, g2];

V = (x - x_d)^2 + (z - z_d)^2 -(r + delta)^2;

LfV = jacobian(V, X) * f;
Lf2V = jacobian(LfV, X) * f;
Lf3V = jacobian(Lf2V, X) * f;
Lf4V = jacobian(Lf3V, X) * f;

LgLf3V = jacobian(Lf3V, X) * G;

LfV = simplify(LfV);
Lf2V = simplify(Lf2V);
Lf3V = simplify(Lf3V);
Lf4V = simplify(Lf4V);
LgLf3V = simplify(LgLf3V);


matlabFunction(V, LfV, Lf2V, Lf3V, Lf4V, LgLf3V, ...
    'File', 'compute_TLC_derivatives', ...
    'Vars', {X, x_d, z_d, r, delta}, ...
    'Outputs', {'V', 'LfV', 'Lf2V', 'Lf3V', 'Lf4V', 'LgLf3V'});


%% create remainder function
                        
syms t u1 u2 real   

syms x0 z0 theta0 dx0 dz0 dtheta0 F0 dF0 real
X0 = [x0; z0; theta0; dx0; dz0; dtheta0; F0; dF0];



F_t = F0 + dF0*t + 0.5*u1*t^2;
dF_t = dF0 + u1*t;
dtheta_t = dtheta0 - (1/I)*u2*t;
theta_t = theta0 + dtheta0*t - 0.5*(1/I)*u2*t^2;

ax_t = (F_t / m) * (sqrt(3)/2);
az_t = (F_t / m) - g_acc;

dx_t = dx0 + int(ax_t, t, 0, t);
dz_t = dz0 + int(az_t, t, 0, t);
x_t  = x0 + int(dx_t, t, 0, t);
z_t  = z0 + int(dz_t, t, 0, t);

X_t = [x_t; z_t; pi/3; dx_t; dz_t; dtheta_t; F_t; dF_t];

Lf4V_t = subs(Lf4V, X, X_t);
LgLf3V_t = subs(LgLf3V, X, X_t);

% N_order = 4; 
% Lf4V_poly = taylor(Lf4V_t+LgLf3V_t*[u1;u2], t, 0, 'Order', N_order);

C_Lf4V = coeffs(Lf4V_t+LgLf3V_t*[u1;u2]-subs(LgLf3V, X, X0)*[u1;u2]-subs(Lf4V,X,X0), t, 'All');


matlabFunction(C_Lf4V, ...
    'File', 'compute_reaminder_poly_coeffs', ...
    'Vars', {X0, u1, u2, x_d, z_d}, ...
    'Outputs', {'out'});