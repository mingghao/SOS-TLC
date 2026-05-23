clear;
clc;
set(0,'DefaultTextInterpreter','latex');

% create SDP solver
sdp_solver =SDP_solver(Parameter);

%% main loop
ex=3; % ex=1, hocbf; 2, zoh-tlc; 3, sos_tlc

for sub_ex=1:3
    timestep = Parameter.dt;
    duration = 0:timestep:Parameter.T;
    timesteps = length(duration);
    u = zeros([2, timesteps - 1]);
    state_log = zeros([8, timesteps]);
    state_log(:, 1) = [Parameter.x_0; Parameter.z_0;Parameter.theta_0; 0; 0; 0; Parameter.g; 0];
    solve_time = zeros([1, timesteps - 1]);
    sos_slack_log = zeros(2,timesteps-1);
    remainders_log = zeros(2,timesteps-1);

    for i = 1:(timesteps - 1)
        if ex==1 
            tic;
            u(:,i) = compute_QP(state_log(:,i),@get_HOCBF,sub_ex);
            t=toc;
        elseif ex==2
            tic
            u(:,i) = compute_QP(state_log(:,i),@get_ZOHTLC,sub_ex);
            t=toc;
        else
            tic;
            [u(:,i),sos_slack_log(:,i),remainders_log(:,i)] = compute_SDP_SOSTLC(state_log(:,i),sdp_solver,sub_ex);
            t=toc;
        end
        solve_time(i) = t;
        [~, y] = ode45(@(t, y) quadcopter(t, y, u(:,i)),...
                    [0, timestep], state_log(:, i));
        y = y';
        state_log(:, i+1) = y(:, end);
    end  

    %% print info
    b1 = (state_log(1,:)-Parameter.x_o(1)).^2 + (state_log(2,:)-Parameter.z_o(1)).^2 - (Parameter.r(1)+Parameter.delta)^2;
    min_b1 = min(b1);
    b2 = (state_log(1,:)-Parameter.x_o(2)).^2 + (state_log(2,:)-Parameter.z_o(2)).^2 - (Parameter.r(2)+Parameter.delta)^2;
    min_b2 = min(b2);
    mean_time = mean(solve_time);
    sd_time = std(solve_time);

    disp([char(Parameter.exps(ex)), ': minimum b1:', num2str(min_b1),', minimum b2: ',num2str(min_b2),', mean time: ',num2str(mean_time),', standard deviation: ',num2str(sd_time)]);

    %% plot figures
    if sub_ex == 1
        figure('Position',[695, 494, 825, 578]);
        hold on;
        for plot_idx=1:length(Parameter.r)
            fill(Parameter.x_o(plot_idx)+Parameter.r(plot_idx)*cos(0:0.01:2*pi+0.01),Parameter.z_o(plot_idx)+Parameter.r(plot_idx)*sin(0:0.01:2*pi+0.01), 'b', 'FaceColor','#4DBEEE','EdgeColor', 'none');
            plot(Parameter.x_o(plot_idx)+(Parameter.r(plot_idx)+Parameter.delta)*cos(0:0.01:2*pi+0.01),Parameter.z_o(plot_idx)+(Parameter.r(plot_idx)+Parameter.delta)*sin(0:0.01:2*pi+0.01),'k--','LineWidth',1.5);
        end
        plot(state_log(1,:),state_log(2,:),'LineWidth',1.5,'Color',Parameter.color(ex),'LineStyle',Parameter.linestyle(sub_ex));
        axis equal;
        
        figure('Position',[1328, 647, 507, 350])
        subplot(2,1,1)
        plot([0 10], [Parameter.F_max,Parameter.F_max], 'k--',"LineWidth", 1.5);
        hold on;
        plot([0 10], [0,0], 'k--',"LineWidth", 1.5);
        plot(duration, state_log(7, :),"LineWidth", 1.5,'Color',Parameter.color(ex),'LineStyle',Parameter.linestyle(sub_ex));
        xlabel('$t/s$','Fontsize',20);
        ylabel('$F/N$','Fontsize',20)
        xlim([0, Parameter.T]);
        ylim([-1, 15]);
        grid on
        hold off;
        subplot(2,1,2)
        plot([0 10], [Parameter.theta_max,Parameter.theta_max], 'k--',"LineWidth", 1.5);
        hold on;
        plot([0 10], [-Parameter.theta_max,-Parameter.theta_max], 'k--',"LineWidth", 1.5);
        xlim([0, Parameter.T]);
        ylim([-Parameter.theta_max-0.2, Parameter.theta_max+0.2]);
        plot(duration, state_log(3, :),"LineWidth", 1.5,'Color',Parameter.color(ex),'LineStyle',Parameter.linestyle(sub_ex));
        xlabel('$t/s$','Fontsize',20);
        ylabel('$\theta/rad$','Fontsize',20)

        grid on
        hold off;



  
    else
        figure(1)
        hold on;
        plot(state_log(1,:),state_log(2,:),'LineWidth',1.5,'Color',Parameter.color(ex),'LineStyle',Parameter.linestyle(sub_ex));
        hold off;

        figure(2);
        subplot(2,1,1)
        hold on;
        plot(duration, state_log(7, :),"LineWidth", 1.5,'Color',Parameter.color(ex),'LineStyle',Parameter.linestyle(sub_ex));
        hold off;
        subplot(2,1,2)
        hold on;
        plot(duration, state_log(3, :),"LineWidth", 1.5,'Color',Parameter.color(ex),'LineStyle',Parameter.linestyle(sub_ex));
        hold off;
    end
end

figure(1);
hold on;
xlabel('$x/m$','Fontsize',20);
ylabel('$z/m$','Fontsize',20);
if ex == 1
    legend('','','','','p=2','p=4', 'p=6','Interpreter', 'latex','Fontsize',20);
else
    legend('','','','','$\Delta t =0.8$','$\Delta t =0.9$', '$\Delta t =1.0$','Interpreter', 'latex','Fontsize',20);
end

plot(0, 2, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 12, 'HandleVisibility', 'off');
text(0.2, 1.6, 'start', 'FontSize', 16, 'FontName', 'Times New Roman');
plot(10, 1, 'k^', 'MarkerFaceColor', 'k', 'MarkerSize', 14, 'HandleVisibility', 'off');
text(9.8, 1.6, 'target', 'FontSize', 16, 'FontName', 'Times New Roman', 'HorizontalAlignment', 'center');
xlim([-1,11]);
ylim([-1,5]);
legend('Box','off')
hold off;

figure(2);
hold on;
subplot(2,1,1)
% if ex == 1
%     legend('','','p=2','p=4', 'p=6','Interpreter', 'latex','Fontsize',20);
% else
%     legend('','','\Delta t =0.8','\Delta t =0.9', '\Delta t =1.0','Interpreter', 'latex','Fontsize',20);
% end
% legend('Box','off');
hold off;


function u = compute_QP(cur_state, safe_func,ex)
    [A_tls, b_tls] = get_TLS(cur_state);
    [A_tlc, b_tlc] = safe_func(cur_state,ex);
    H = diag([1,1,100,100]);
    F = [0; 0; 0; 0];
    A = [A_tls; A_tlc; 0, 1, 0, 0; 0, -1, 0, 0; 1, 0, 0, 0; -1, 0, 0, 0];
    b = [b_tls; b_tlc; Parameter.M_max; Parameter.M_max; Parameter.uF_max; Parameter.uF_max];
    [solution,~,~,~,~] = quadprog(H,F,A,b,[],[],[],[],[],Parameter.options);
    u = solution(1:2);
end

function [A_tls, b_tls] = get_TLS(state)
    ex = state(1) - Parameter.x_d;
    evx = state(4);
    ey = state(2) - Parameter.z_d;
    evy = state(5);
    
    Kp = 1.5; 
    Kd = 2.0; 
    
    acc_d = zeros(2,1);
    acc_d(1) = -Kp * ex - Kd * evx;
    acc_d(2) = -Kp * ey - Kd * evy;
    
    theta_d = atan(acc_d(1) / (acc_d(2) + Parameter.g));
    if theta_d > Parameter.theta_max
        theta_d = Parameter.theta_max;
    elseif theta_d < -Parameter.theta_max
        theta_d = -Parameter.theta_max;
    end
    
    F_d = sqrt(acc_d(1)^2 + (acc_d(2) + Parameter.g)^2) * Parameter.m;
    if F_d > Parameter.F_max
        F_d = Parameter.F_max;    
    end

    eF = state(7) - F_d;
    evF = state(8);
    etheta = state(3) - theta_d;
    evtheta = state(6);
    
    k_in = 0.2; 
    gamma = 20; 
    
    V1 = (eF + k_in*evF)^2;
    LfV1 = 2*(eF + k_in*evF)*evF;
    LgV1 = [2*(eF + k_in*evF)*k_in,  0];
    
    A1_tls = [LgV1, -1, 0];
    b1_tls = -gamma*V1 - LfV1;
    
    V2 = (etheta + k_in*evtheta)^2;
    LfV2 = 2*(etheta + k_in*evtheta)*evtheta;
    LgV2 = [0,  -2*(etheta + k_in*evtheta)*k_in/Parameter.I];
    
    A2_tls = [LgV2, 0, -1];
    b2_tls = -gamma*V2 - LfV2;
    
    A_tls = [A1_tls; A2_tls];
    b_tls = [b1_tls; b2_tls];
    
end


function [A_tlc, b_tlc] = get_ZOHTLC(state,ex)
    A1_tlc = [];
    b1_tlc = [];
    Delta = Parameter.tlc_t(ex);
    for obs_idx=1:length(Parameter.r)
        [V,LfV,Lf2V,Lf3V,Lf4V,LgLf3V] = compute_TLC_derivatives(state,Parameter.x_o(obs_idx),Parameter.z_o(obs_idx),Parameter.r(obs_idx),Parameter.delta);
        A1_tlc = [A1_tlc; -LgLf3V/factorial(4), 0, 0];
        b1_tlc = [b1_tlc; V/Delta^4 + LfV/Delta^3 + Lf2V/(factorial(2)*Delta^2) + Lf3V/(factorial(3)*Delta) + Lf4V/factorial(4)];
    end
    % theta constraints
    h1 = state(3)+Parameter.theta_max;
    dh1 = state(6);
    Lgh1 = [0,-1/Parameter.I];
    h2 = Parameter.theta_max-state(3);
    dh2 = -state(6);
    Lgh2 = [0,1/Parameter.I];
    A2_tlc = [-Parameter.Delta_t*Lgh1, 0, 0; -Parameter.Delta_t*Lgh2, 0, 0];
    b2_tlc = [h1+Parameter.Delta_t*dh1; h2+Parameter.Delta_t*dh2];


    % thrust constraints
    h1 = state(7);
    dh1 = state(8);
    Lgh1 = [1,0];
    h2 = Parameter.F_max-state(7);
    dh2 = -state(8);
    Lgh2 = [-1,0];
    A3_tlc = [-Parameter.Delta_t*Lgh1, 0, 0; -Parameter.Delta_t*Lgh2, 0, 0];
    b3_tlc = [h1+Parameter.Delta_t*dh1; h2+Parameter.Delta_t*dh2];


    A_tlc = [A1_tlc;A2_tlc;A3_tlc];
    b_tlc = [b1_tlc;b2_tlc;b3_tlc];


end
   
function [A_tlc, b_tlc] = get_HOCBF(state, ex)
    k1=Parameter.cbf_k(ex);k2=Parameter.cbf_k(ex);
    k3=Parameter.cbf_k(ex);k4=Parameter.cbf_k(ex);
    A1_tlc = [];
    b1_tlc = [];
    for obs_idx=1:length(Parameter.r)
        [V,LfV,Lf2V,Lf3V,Lf4V,LgLf3V] = compute_TLC_derivatives(state,Parameter.x_o(obs_idx),Parameter.z_o(obs_idx),Parameter.r(obs_idx),Parameter.delta);
        A1_tlc = [A1_tlc; -LgLf3V, 0, 0];
        b1_tlc = [b1_tlc; k1*k2*k3*k4*V + (k1*k2*k3+k1*k2*k4+k1*k3*k4+k2*k3*k4)*LfV + (k1*k2+k1*k3+k1*k4+k2*k3+k2*k4+k3*k4)*Lf2V + (k1+k2+k3+k4)*Lf3V + Lf4V];
    end
    % theta constraints
    k1 = 5;
    k2 =5;
    h1 = state(3)+Parameter.theta_max;
    dh1 = state(6);
    Lgh1 = [0,-1/Parameter.I];
    h2 = Parameter.theta_max-state(3);
    dh2 = -state(6);
    Lgh2 = [0,1/Parameter.I];
    A2_tlc = [-Lgh1, 0, 0; -Lgh2, 0, 0];
    b2_tlc = [k1*k2*h1+(k1+k2)*dh1; k1*k2*h2+(k1+k2)*dh2];


    % thrust constraints
    k1 = 5;
    k2 =5;
    h1 = state(7);
    dh1 = state(8);
    Lgh1 = [1,0];
    h2 = Parameter.F_max-state(7);
    dh2 = -state(8);
    Lgh2 = [-1,0];
    A3_tlc = [-Lgh1, 0, 0; -Lgh2, 0, 0];
    b3_tlc = [k1*k2*h1+(k1+k2)*dh1; k1*k2*h2+(k1+k2)*dh2];

    
    A_tlc = [A1_tlc;A2_tlc;A3_tlc];
    b_tlc = [b1_tlc;b2_tlc;b3_tlc];
end

function [A_tlc, b_tlc] = get_rTLC(state,ex)
    A1_tlc = [];
    b1_tlc = [];
    remainders = [0; 0];
    
    Delta = Parameter.tlc_t(ex);
    for obs_idx=1:length(Parameter.r)
        [V,LfV,Lf2V,Lf3V,Lf4V,LgLf3V] = compute_TLC_derivatives(state,Parameter.x_o(obs_idx),Parameter.z_o(obs_idx),Parameter.r(obs_idx),Parameter.delta);
        u_cases = [Parameter.uF_max, Parameter.M_max;
                   -Parameter.uF_max, Parameter.M_max;
                   Parameter.uF_max, -Parameter.M_max;
                   -Parameter.uF_max, -Parameter.M_max;];
        for cases = 1:4       
            remainder_coeffs = compute_reaminder_poly_coeffs(state, u_cases(cases,1), u_cases(cases,2), Parameter.x_o(obs_idx),Parameter.z_o(obs_idx));
            remainder_min = find_poly_min(remainder_coeffs, 0, Delta);
            remainders(obs_idx) = min(remainders(obs_idx), remainder_min*Delta^4/24);
        end


        A1_tlc = [A1_tlc; -LgLf3V/factorial(4), 0, 0];
        b1_tlc = [b1_tlc; V/Delta^4 + LfV/Delta^3 + Lf2V/(factorial(2)*Delta^2) + Lf3V/(factorial(3)*Delta) + Lf4V/factorial(4)+remainders(obs_idx)/Delta^4];
    end
    % theta constraints
    h1 = state(3)+Parameter.theta_max;
    dh1 = state(6);
    Lgh1 = [0,-1/Parameter.I];
    h2 = Parameter.theta_max-state(3);
    dh2 = -state(6);
    Lgh2 = [0,1/Parameter.I];
    A2_tlc = [-Parameter.Delta_t*Lgh1, 0, 0; -Parameter.Delta_t*Lgh2, 0, 0];
    b2_tlc = [h1+Parameter.Delta_t*dh1; h2+Parameter.Delta_t*dh2];


    % thrust constraints
    h1 = state(7);
    dh1 = state(8);
    Lgh1 = [1,0];
    h2 = Parameter.F_max-state(7);
    dh2 = -state(8);
    Lgh2 = [-1,0];
    A3_tlc = [-Parameter.Delta_t*Lgh1, 0, 0; -Parameter.Delta_t*Lgh2, 0, 0];
    b3_tlc = [h1+Parameter.Delta_t*dh1; h2+Parameter.Delta_t*dh2];

    
    A_tlc = [A1_tlc;A2_tlc;A3_tlc];
    b_tlc = [b1_tlc;b2_tlc;b3_tlc];
end


function [u,sos_slack,remainders] = compute_SDP_SOSTLC(state,solver,ex)
    
    [A_tls, b_tls] = get_TLS(state);
    Delta = Parameter.tlc_t(ex);
    remainders = [0; 0];
    for obs_idx=1:length(Parameter.r)
        u_cases = [Parameter.uF_max, Parameter.M_max;
                   -Parameter.uF_max, Parameter.M_max;
                   Parameter.uF_max, -Parameter.M_max;
                   -Parameter.uF_max, -Parameter.M_max];
        for cases = 1:4       
            remainder_coeffs = compute_reaminder_poly_coeffs(state, u_cases(cases,1), u_cases(cases,2), Parameter.x_o(obs_idx),Parameter.z_o(obs_idx));
            remainder_min = find_poly_min(remainder_coeffs, 0, Delta);
            remainders(obs_idx) = min(remainders(obs_idx), remainder_min*Delta^4/24);
        end
    end


    % theta constraints
    h1 = state(3)+Parameter.theta_max;
    dh1 = state(6);
    mu_bound = (-h1 - dh1*Parameter.Delta_t)/Parameter.Delta_t^2;
    p = [mu_bound, dh1, h1];
    [min_val1,~] = find_poly_min(p, 0, Parameter.Delta_t);

    h2 = Parameter.theta_max-state(3);
    dh2 = -state(6);
    mu_bound = (-h2 - dh2*Parameter.Delta_t)/Parameter.Delta_t^2;
    p = [mu_bound, dh2, h2];
    [min_val2,~] = find_poly_min(p, 0, Parameter.Delta_t);


    % thrust constraints
    h1 = state(7);
    dh1 = state(8);
    mu_bound = (-h1 - dh1*Parameter.Delta_t)/Parameter.Delta_t^2;
    p = [mu_bound, dh1, h1];
    [min_val3,~] = find_poly_min(p, 0, Parameter.Delta_t);

    h2 = Parameter.F_max-state(7);
    dh2 = -state(8);
    mu_bound = (-h2 - dh2*Parameter.Delta_t)/Parameter.Delta_t^2;
    p = [mu_bound, dh2, h2];
    [min_val4,~] = find_poly_min(p, 0, Parameter.Delta_t);

    % safe constraints
    [V1,LfV1,Lf2V1,Lf3V1,Lf4V1,LgLf3V1] = compute_TLC_derivatives(state,Parameter.x_o(1),Parameter.z_o(1),Parameter.r(1),Parameter.delta);
    mu_bound = (-V1 - remainders(1) - LfV1*Delta- Lf2V1*Delta^2/2 - Lf3V1*Delta^3/6)/Delta^4;
    p = [mu_bound, Lf3V1/6, Lf2V1/2, LfV1, V1+remainders(1)];
    [min_val5,~] = find_poly_min(p, 0, Delta);

    [V2,LfV2,Lf2V2,Lf3V2,Lf4V2,LgLf3V2] = compute_TLC_derivatives(state,Parameter.x_o(2),Parameter.z_o(2),Parameter.r(2),Parameter.delta);
    mu_bound = (-V2 - remainders(2) - LfV2*Delta- Lf2V2*Delta^2/2 - Lf3V2*Delta^3/6)/Delta^4;
    p = [mu_bound, Lf3V2/6, Lf2V2/2, LfV2, V2+remainders(2)];
    [min_val6,~] = find_poly_min(p, 0, Delta);
    min_val = min([min_val1,min_val2,min_val3,min_val4,min_val5, min_val6]);
    
    if min_val>=0
       u = compute_QP(state,@get_rTLC,ex);
       sos_slack = [0;0];
    else
        coeff_h = [V1,LfV1,Lf2V1,Lf3V1,Lf4V1;V2,LfV2,Lf2V2,Lf3V2,Lf4V2];
        coeff_hg = [LgLf3V1;LgLf3V2];
        solution = solver(state, Delta, A_tls, b_tls, coeff_h, coeff_hg, remainders);
        u = solution{1,1};
        sos_slack = solution{1,3};
    end
  
end


function [min_val, min_x] = find_poly_min(p, a, b)
    val_a = polyval(p, a);
    val_b = 0;
    
    dp = polyder(p); 
    roots_dp = roots(dp);
    
    real_roots = roots_dp(abs(imag(roots_dp)) < 1e-10);
    real_roots = real(real_roots); 
    
    valid_roots = real_roots(real_roots >= a & real_roots <= b);
    
    val_roots = polyval(p, valid_roots);
    
    all_x = [a; b; valid_roots];
    all_vals = [val_a; val_b; val_roots];
    
    [min_val, min_idx] = min(all_vals);
    min_x = all_x(min_idx);
end
