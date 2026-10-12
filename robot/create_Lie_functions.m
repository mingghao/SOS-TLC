clear; clc;

syms x y vel acc theta dtheta real
syms x_o y_o r delta real

X = [x; y; vel; acc; theta; dtheta];

f = [vel*cos(theta);
     vel*sin(theta);
     acc;
     0;
     dtheta;
     0];

g1 = [0; 0; 0; 1; 0; 0];
g2 = [0; 0; 0; 0; 0; 1];
G = [g1, g2];

V = (x-x_o)^2 + (y-y_o)^2 - (r+delta)^2;

LfV = jacobian(V,X)*f;
Lf2V = jacobian(LfV,X)*f;
Lf3V = jacobian(Lf2V,X)*f;
LgLf2V = jacobian(Lf2V,X)*G;

LfV = simplify(LfV);
Lf2V = simplify(Lf2V);
Lf3V = simplify(Lf3V);
LgLf2V = simplify(LgLf2V);

matlabFunction(V,LfV,Lf2V,Lf3V,LgLf2V, ...
    'File','compute_TLC_derivatives', ...
    'Vars',{X,x_o,y_o,r,delta}, ...
    'Outputs',{'V','LfV','Lf2V','Lf3V','LgLf2V'});

%% The sampled-data remainder is bounded analytically.
syms jerk M t real  
syms x0 y0 vel0 acc0 theta0 dtheta0 real
X0 = [x0; y0; vel0; acc0; theta0; dtheta0];

f = [vel;
     vel;
     acc;
     0;
     dtheta;
     0];

LfV = jacobian(V,X)*f;
Lf2V = jacobian(LfV,X)*f;
Lf3V = jacobian(Lf2V,X)*f;
LgLf2V = jacobian(Lf2V,X)*G;

LfV = simplify(LfV);
Lf2V = simplify(Lf2V);
Lf3V = simplify(Lf3V);
LgLf2V = simplify(LgLf2V);


acc_t = acc0 + jerk*t;
vel_t = vel0 + acc0*t+ 0.5*jerk*t^2;
dtheta_t = dtheta0 + M*t;
theta_t = theta0 + dtheta0*t + 0.5*M*t^2;

x_t  = x0 + int(vel_t, t, 0, t);
y_t  = y0 + int(vel_t, t, 0, t);

X_t = [x_t; y_t; vel_t; acc_t; theta_t; dtheta_t];

Lf3V_t = subs(Lf3V, X, X_t);
LgLf2V_t = subs(LgLf2V, X, X_t);


C_Lf4V = coeffs(Lf3V_t+LgLf2V_t*[jerk;M]-subs(LgLf2V, X, X0)*[jerk;M]-subs(Lf3V,X,X0), t, 'All');


matlabFunction(C_Lf4V, ...
    'File', 'compute_reaminder_poly_coeffs', ...
    'Vars', {X0, jerk, M, x_o, y_o, r, delta}, ...
    'Outputs', {'out'});
