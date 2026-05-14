function dydt = quadcopter(t, y, u)
    % x = [x, z, theta, dx, dz, F,dF] u=[u_F,M]
    dydt = zeros(8, 1);
    dydt(1:3) = y(4:6);
    dydt(4) = sin(y(3))*y(7)/Parameter.m;
    dydt(5) = cos(y(3))*y(7)/Parameter.m-Parameter.g;
    dydt(6) = -u(2)/Parameter.I;
    dydt(7) = y(8);
    dydt(8) = u(1);
end

