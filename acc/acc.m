function dydt = acc(t, y, u)
    dydt = zeros(2, 1);
    F_r = Parameter.f_0.*sign(y(1)) + Parameter.f_1 .* y(1) + Parameter.f_2 .* y(1) .^ 2; 
    dydt(1) = -F_r ./ Parameter.m + u ./ Parameter.m; 
    dydt(2) = Parameter.v_0 - y(1);
end

