function solver = SDP_solver(Parameter)
        
        % Parameter Variables  
        p_v  = sdpvar(1, 1);
        p_v_2 = sdpvar(1, 1);
        p_z  = sdpvar(1, 1);
        p_Fr = sdpvar(1, 1);
        delta = sdpvar(1, 1);
        delta_2 = sdpvar(1, 1);

        % Decision Variables
        a = sdpvar(1, 1);
        slack = sdpvar(1, 1);
        Q0 = sdpvar(2, 2, 'symmetric');
        q3 = sdpvar(1, 1);

        Constraints = [];

        % control bounds
        Constraints = [Constraints, (-Parameter.g * Parameter.c_d <=a)&(a <= Parameter.g * Parameter.c_a)];

        % TLS
        v_bias = p_v - Parameter.v_d;
        A_tls = [2*v_bias*0.1, -1];
        b_tls = 2*v_bias*p_Fr*0.1/Parameter.m - v_bias^2;
        Constraints = [Constraints, A_tls * [a; slack] <= b_tls];

        % TLC
        b = p_z - Parameter.sigma;
        b_dot = Parameter.v_0 - p_v;
        remainder = (Parameter.f_1 * (-Parameter.g * Parameter.c_d - p_Fr/Parameter.m) * delta - Parameter.f_2 * p_v_2) / Parameter.m;

        c0 = b;
        c1 = b_dot;
        c2 = 0.5 * (p_Fr/Parameter.m - a + remainder);

        Constraints = [Constraints, Q0(1,1) == c0];
        Constraints = [Constraints, 2*Q0(1,2) + q3 * delta == c1];
        Constraints = [Constraints, Q0(2,2) - q3 == c2];
        Constraints = [Constraints, Q0 >= 0, q3 >= 0];

        Objective = (a - p_Fr/Parameter.m)^2 + Parameter.eps * slack^2;
        ops = sdpsettings('solver', 'mosek', 'verbose', 0);
        solver = optimizer(Constraints, Objective, ops, [p_v; p_v_2; p_z; p_Fr; delta;delta_2], a);
end

