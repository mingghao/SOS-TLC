function dydt = robot(t, x, u)
    % x = [x; y; vel; acc; theta; dtheta;],
    % u = [jerk; M]

    dydt = zeros(6, 1);
    dydt(1) = x(3)*cos(x(5));
    dydt(2) = x(3)*sin(x(5));
    dydt(3) = x(4);
    dydt(4) = u(1);
    dydt(5) = x(6);
    dydt(6) = u(2);
end
