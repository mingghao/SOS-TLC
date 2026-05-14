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

