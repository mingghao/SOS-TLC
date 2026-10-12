function solver = SDP_solver(Parameter)

    % Parameter Variables
    p_y = sdpvar(6, 1, 'full');
    Delta = sdpvar(1, 1, 'full');
    A_tls = sdpvar(2, 4, 'full');
    b_tls = sdpvar(2, 1, 'full');
    coeff_h = sdpvar(2, 4, 'full');
    coeff_hg = sdpvar(2, 2, 'full');
    remainders = sdpvar(2, 1, 'full');
    u_nom = sdpvar(2, 1, 'full');

    % Decision Variables: u(1) is acceleration; u(2) is angular rate
    u = sdpvar(2, 1, 'full');
    slack = sdpvar(2, 1, 'full');



    Constraints = [];
    Constraints = [Constraints, -Parameter.jerk_max <= u(1), u(1) <= Parameter.jerk_max];
    Constraints = [Constraints, -Parameter.M_max <= u(2), u(2) <= Parameter.M_max];

   
    Constraints = [Constraints, A_tls*[u; slack] <= b_tls];

   
    for obs_idx = 1:2
        Q1 = sdpvar(2, 2, 'symmetric');
        Q2 = sdpvar(2, 2, 'symmetric');
        coeff_c = [coeff_h(obs_idx,1);
                   coeff_h(obs_idx,2);
                   coeff_h(obs_idx,3)/2;
                   (coeff_h(obs_idx,4) + coeff_hg(obs_idx,1)*u(1) ...
                    + coeff_hg(obs_idx,2)*u(2) + remainders(obs_idx))/6];
        coeff_Q = [Delta*Q1(1,1);
                   2*Delta*Q1(1,2) - Q1(1,1) + Q2(1,1);
                   Delta*Q1(2,2) - 2*Q1(1,2) + 2*Q2(1,2);
                   -Q1(2,2) + Q2(2,2)];
        Constraints = [Constraints, coeff_c == coeff_Q];
        Constraints = [Constraints, Q1 >= 0, Q2 >= 0];
    end

    gamma1 = sdpvar(1,1);
    gamma2 = sdpvar(1,1);
    u_error = u - u_nom;
    Constraints = [Constraints, [gamma1, 0, u_error(1); ...
                                0, gamma1, u_error(2); ...
                                u_error(1), u_error(2), 1] >= 0];
    Constraints = [Constraints, [gamma2, 0, slack(1); ...
                                0, gamma2, slack(2); ...
                                slack(1), slack(2), 1] >= 0];
    weight_slack1 = 100;
    Objective = gamma1 + weight_slack1*gamma2;
    ops = sdpsettings('solver', 'mosek', 'verbose', 0);

    solver = optimizer(Constraints, Objective, ops, ...
        {p_y, Delta, A_tls, b_tls, coeff_h, coeff_hg, remainders, u_nom}, ...
        {u, slack});
end
