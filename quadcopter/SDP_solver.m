function solver = SDP_solver(Parameter)
    
    % Parameter Variables
    p_y = sdpvar(8, 1); 
    Delta = sdpvar(1, 1); 
    A_tls = sdpvar(2, 4);
    b_tls = sdpvar(2, 1);
    coeff_h = sdpvar(2,5);
    coeff_hg = sdpvar(2,2);
    remainders = sdpvar(2,1);
    % Decision Variables
    u = sdpvar(2, 1);
    slack = sdpvar(2, 1);
    
    % state
    theta = p_y(3); dtheta = p_y(6);
    F = p_y(7); dF = p_y(8);
    
    % input Constraint
    Constraints = [];
    Constraints = [Constraints, (-Parameter.M_max <=u(2))&(u(2) <= Parameter.M_max)];
    Constraints = [Constraints, (-Parameter.uF_max <=u(1))&(u(1) <= Parameter.uF_max)];

    % TLS
    Constraints = [Constraints, A_tls*[u;slack]<=b_tls];

    % theta Constraint
    h1 = Parameter.theta_max - theta;
    Lfh1 = -dtheta;
    LgLfh1 = u(2)/Parameter.I;
    Q1 = sdpvar(2, 2, 'symmetric');
    q1 = sdpvar(1, 1);

    coeff_c1 = [h1; Lfh1; LgLfh1/2];
    coeff_Q1 = [Q1(1,1); 2*Q1(1,2); Q1(2,2)];
    coeff_q1 = [0; Parameter.Delta_t*q1; -q1];
    Constraints = [Constraints, coeff_c1 == coeff_Q1 + coeff_q1];
    Constraints = [Constraints, Q1 >= 0, q1 >= 0];

    h2 = Parameter.theta_max + theta;
    Lfh2 = dtheta;
    LgLfh2 = -u(2)/Parameter.I;
    Q2 = sdpvar(2, 2, 'symmetric');
    q2 = sdpvar(1, 1);

    coeff_c2 = [h2; Lfh2; LgLfh2/2];
    coeff_Q2 = [Q2(1,1); 2*Q2(1,2); Q2(2,2)];
    coeff_q2 = [0; Parameter.Delta_t*q2; -q2];
    Constraints = [Constraints, coeff_c2 == coeff_Q2 + coeff_q2];
    Constraints = [Constraints, Q2 >= 0, q2 >= 0];


    % thrust Constraint
    h3 = Parameter.F_max - F;
    Lfh3 = -dF;
    LgLfh3 = -u(1);
    Q3 = sdpvar(2, 2, 'symmetric');
    q3 = sdpvar(1, 1);

    coeff_c3 = [h3; Lfh3; LgLfh3/2];
    coeff_Q3 = [Q3(1,1); 2*Q3(1,2); Q3(2,2)];
    coeff_q3 = [0; Parameter.Delta_t*q3; -q3];
    Constraints = [Constraints, coeff_c3 == coeff_Q3 + coeff_q3];
    Constraints = [Constraints, Q3 >= 0, q3 >= 0];

    h4 = F;
    Lfh4 = dF;
    LgLfh4 = u(1);
    Q4 = sdpvar(2, 2, 'symmetric');
    q4 = sdpvar(1, 1);

    coeff_c4 = [h4; Lfh4; LgLfh4/2];
    coeff_Q4 = [Q4(1,1); 2*Q4(1,2); Q4(2,2)];
    coeff_q4 = [0; Parameter.Delta_t*q4; -q4];
    Constraints = [Constraints, coeff_c4 == coeff_Q4 + coeff_q4];
    Constraints = [Constraints, Q4 >= 0, q4 >= 0];



    % safety constraints
    sos_slack = sdpvar(2,1);
    Q5 = sdpvar(3, 3, 'symmetric');
    q5 = sdpvar(2, 2, 'symmetric');
    coeff_c5 = [coeff_h(1,1)+sos_slack(1);
                coeff_h(1,2)* Delta;
                coeff_h(1,3)/2* Delta^2;
                coeff_h(1,4)/6* Delta^3;
                ((coeff_h(1,5)+coeff_hg(1,1)*u(1)+coeff_hg(1,2)*u(2))/24)* Delta^4];
    coeff_Q5 = [Q5(1,1);
                2*Q5(1,2);
                2*Q5(1,3) + Q5(2,2);
                2*Q5(2,3);
                Q5(3,3)];
    coeff_q5 = [0;
                q5(1,1);
                2 * q5(1,2) - q5(1,1);
                q5(2,2) - 2*q5(1,2);
                -q5(2,2)];
    Constraints = [Constraints, coeff_c5 == coeff_Q5 + coeff_q5];
    Constraints = [Constraints, Q5 >= 0, q5 >= 0];
    Q6 = sdpvar(3, 3, 'symmetric');
    q6 = sdpvar(2, 2, 'symmetric');
    coeff_c6 = [coeff_h(2,1)+sos_slack(2);
                coeff_h(2,2)* Delta;
                coeff_h(2,3)/2* Delta^2;
                coeff_h(2,4)/6* Delta^3;
                ((coeff_h(2,5)+coeff_hg(2,1)*u(1)+coeff_hg(2,2)*u(2))/24)* Delta^4];
    coeff_Q6 = [Q6(1,1);
                2*Q6(1,2);
                2*Q6(1,3) + Q6(2,2);
                2*Q6(2,3);
                Q6(3,3)];
    coeff_q6 = [0;
                q6(1,1);
                2 * q6(1,2) - q6(1,1);
                q6(2,2) - 2*q6(1,2);
                -q6(2,2)];
    Constraints = [Constraints, coeff_c6 == coeff_Q6 + coeff_q6];
    Constraints = [Constraints, Q6 >= 0, q6 >= 0];
    
    gamma1 = sdpvar(1,1);
    gamma2 = sdpvar(1,1);
    Constraints = [Constraints, [gamma1, 0, u(1); 0, gamma1, u(2); u(1), u(2), 1]>=0];
    Constraints = [Constraints, [gamma2, 0, slack(1); 0, gamma2, slack(2); slack(1), slack(2), 1]>=0];
    % Constraints = [Constraints, sos_slack(1)>=remainders(1)];
    % Constraints = [Constraints, sos_slack(2)>=remainders(2)];

    weight_slack1 = 100;
    weight_slack2 = 100;
    Objective = gamma1 + weight_slack1 * gamma2 + weight_slack2*(sos_slack(1) + sos_slack(2));
    ops = sdpsettings('solver', 'mosek', 'verbose', 0);
    
    solver = optimizer(Constraints, Objective, ops, {p_y, Delta, A_tls, b_tls, coeff_h, coeff_hg, remainders}, {u, slack,sos_slack});
end