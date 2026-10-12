clear;
clc;
set(0,'DefaultTextInterpreter','latex');

% create SDP solver
sdp_solver = SDP_solver(Parameter);

%% main loop
ex = 3; % ex=1, HOCBF; 2, ZOH-TLC; 3, SOS-TLC

for sub_ex = 1:3
    timestep = Parameter.dt;
    duration = 0:timestep:Parameter.T;
    timesteps = length(duration);
    control_log = zeros(2,timesteps-1);
    state_log = zeros(6,timesteps);
    state_log(:,1) = [Parameter.x_0; Parameter.y_0; Parameter.vel_0; Parameter.acc_0; Parameter.theta_0; ...
                      Parameter.dtheta_0];

    for i = 1:(timesteps-1)
        t = duration(i);
        if ex == 1
            control_log(:,i) = compute_QP(state_log(:,i),@get_HOCBF,sub_ex,t);
        elseif ex == 2
            control_log(:,i) = compute_QP(state_log(:,i),@get_ZOHTLC,sub_ex,t);
        else
            control_log(:,i) = compute_SDP_SOSTLC(state_log(:,i),sdp_solver,sub_ex,t);
        end
        [~, trajectory] = ode45(@(time,state) robot(time,state,control_log(:,i)), ...
                                [0,timestep],state_log(:,i));
        state_log(:,i+1) = trajectory(end,:)';
    end

    %% print info
    for obs_idx = 1:length(Parameter.r)
        barrier = (state_log(1,:)-Parameter.x_o(obs_idx)).^2 + ...
                  (state_log(2,:)-Parameter.y_o(obs_idx)).^2 - ...
                  (Parameter.r(obs_idx)+Parameter.delta)^2;
        fprintf('%s, case %d: minimum b%d = %.6g\n', ...
                char(Parameter.exps(ex)),sub_ex,obs_idx,min(barrier));
    end

    %% plot figures
    if sub_ex == 1
        figure('Position',[695,494,825,578]);
        hold on;
        for obs_idx = 1:length(Parameter.r)
            circle_t = 0:0.01:2*pi+0.01;
            fill(Parameter.x_o(obs_idx)+Parameter.r(obs_idx)*cos(circle_t), ...
                 Parameter.y_o(obs_idx)+Parameter.r(obs_idx)*sin(circle_t), ...
                 'b','FaceColor','#4DBEEE','EdgeColor','none');
            plot(Parameter.x_o(obs_idx)+(Parameter.r(obs_idx)+Parameter.delta)*cos(circle_t), ...
                 Parameter.y_o(obs_idx)+(Parameter.r(obs_idx)+Parameter.delta)*sin(circle_t), ...
                 'k--','LineWidth',1.5);
        end
        figure('Position',[1328,647,507,350]);
        subplot(2,1,1);
        hold on;
        plot([0,Parameter.T],[Parameter.jerk_max,Parameter.jerk_max],'k--','LineWidth',1.5);
        plot([0,Parameter.T],-[Parameter.jerk_max,Parameter.jerk_max],'k--','LineWidth',1.5);
        ylabel('$u_F/(m/s^3)$','FontSize',20);
        grid on;
        subplot(2,1,2);
        hold on;
        plot([0,Parameter.T],[Parameter.M_max,Parameter.M_max],'k--','LineWidth',1.5);
        plot([0,Parameter.T],[-Parameter.M_max,-Parameter.M_max],'k--','LineWidth',1.5);
        ylabel('$M/(rad/s^2)$','FontSize',20);
        grid on;
    end

    figure(1);
    hold on;
    plot(state_log(1,:),state_log(2,:),'LineWidth',1.5, ...
         'Color',Parameter.color(ex),'LineStyle',Parameter.linestyle(sub_ex));
    axis equal;
    figure(2);
    subplot(2,1,1);
    hold on;
    plot(duration(1:end-1),control_log(1,:),'LineWidth',1.5, ...
         'Color',Parameter.color(ex),'LineStyle',Parameter.linestyle(sub_ex));
    subplot(2,1,2);
    hold on;
    plot(duration(1:end-1),control_log(2,:),'LineWidth',1.5, ...
         'Color',Parameter.color(ex),'LineStyle',Parameter.linestyle(sub_ex));
end

figure(1);
hold on;
ylim([0,6]);
xlim([0,11]);
plot(Parameter.x_0,Parameter.y_0,'ko','MarkerFaceColor','k', ...
     'MarkerSize',12,'HandleVisibility','off');
plot(Parameter.x_d,Parameter.y_d,'k^','MarkerFaceColor','k', ...
     'MarkerSize',14,'HandleVisibility','off');
xlabel('$x/m$','FontSize',20);
ylabel('$y/m$','FontSize',20);
if ex == 1
    legend('','','','','$p_1=p_2=p_3=1$','$p_1=p_2=p_3=2$','$p_1=p_2=p_3=3$', ...
        'Interpreter','latex','FontSize',20);
else
     legend('','','','','$\Delta t=0.4$','$\Delta t=0.5$','$\Delta t=0.6$', ...
        'Interpreter','latex','FontSize',20);
end
legend('Box','off');
hold off;

figure(2);
subplot(2,1,1);
xlabel('$t/s$','FontSize',20);
xlim([0,Parameter.T]);
ylim([-Parameter.jerk_max-1,Parameter.jerk_max+1])
subplot(2,1,2);
xlabel('$t/s$','FontSize',20);
xlim([0,Parameter.T]);
ylim([-Parameter.M_max-1,Parameter.M_max+1])


function u = compute_QP(cur_state,safe_func,ex,t)
    ref = straight_line_reference(t,[Parameter.x_0;Parameter.y_0], ...
                                      [Parameter.x_d;Parameter.y_d]);
    [A_tls,b_tls,u_nom] = get_TLS(cur_state,ref);
    [A_tlc,b_tlc] = safe_func(cur_state,ex);
    H = 2*diag([1,1,100,100]);
    F = [-2*u_nom;0;0];
    A = [A_tls;A_tlc; ...
         1,0,0,0;-1,0,0,0;0,1,0,0;0,-1,0,0];
    b = [b_tls;b_tlc;Parameter.jerk_max;Parameter.jerk_max; ...
         Parameter.M_max;Parameter.M_max];
    [solution,~,~] = quadprog(H,F,A,b,[],[],[],[],[],Parameter.options);
    u = solution(1:2);
end


function [A_tlc,b_tlc] = get_ZOHTLC(state,ex)
    A_tlc = [];
    b_tlc = [];
    Delta = Parameter.tlc_t(ex);
    for obs_idx = 1:length(Parameter.r)
        [V,LfV,Lf2V,Lf3V,LgLf2V] = compute_TLC_derivatives( ...
            state,Parameter.x_o(obs_idx),Parameter.y_o(obs_idx), ...
            Parameter.r(obs_idx),Parameter.delta);
        A_tlc = [A_tlc;-LgLf2V/factorial(3),0,0];
        b_tlc = [b_tlc;V/Delta^3+LfV/Delta^2+ ...
                  Lf2V/(2*Delta)+Lf3V/factorial(3)];
    end
end


function [A_tlc,b_tlc] = get_HOCBF(state,ex)
    k = Parameter.cbf_k(ex);
    A_tlc = [];
    b_tlc = [];
    for obs_idx = 1:length(Parameter.r)
        [V,LfV,Lf2V,Lf3V,LgLf2V] = compute_TLC_derivatives( ...
            state,Parameter.x_o(obs_idx),Parameter.y_o(obs_idx), ...
            Parameter.r(obs_idx),Parameter.delta);
        A_tlc = [A_tlc;-LgLf2V,0,0];
        b_tlc = [b_tlc;k^3*V+3*k^2*LfV+3*k*Lf2V+Lf3V];
    end
end


function [A_tlc,b_tlc] = get_rTLC(state,ex)
    A_tlc = [];
    b_tlc = [];
    Delta = Parameter.tlc_t(ex);
    remainders = remainder_bounds(state,Delta);
    for obs_idx = 1:length(Parameter.r)
        [V,LfV,Lf2V,Lf3V,LgLf2V] = compute_TLC_derivatives( ...
            state,Parameter.x_o(obs_idx),Parameter.y_o(obs_idx), ...
            Parameter.r(obs_idx),Parameter.delta);
        A_tlc = [A_tlc;-LgLf2V/factorial(3),0,0];
        b_tlc = [b_tlc;V/Delta^3+LfV/Delta^2+ ...
                  Lf2V/(2*Delta)+(Lf3V+remainders(obs_idx))/factorial(3)];
    end
    
end




function remainders = remainder_bounds(state,Delta)
    remainders = zeros(length(Parameter.r),1);
    u_cases = [Parameter.jerk_max,Parameter.M_max; ...
              -Parameter.jerk_max,Parameter.M_max; ...
               Parameter.jerk_max,-Parameter.M_max; ...
              -Parameter.jerk_max,-Parameter.M_max];
    for obs_idx = 1:length(Parameter.r)
        for case_idx = 1:size(u_cases,1)
            remainder_coeffs = compute_reaminder_poly_coeffs( ...
                state,u_cases(case_idx,1),u_cases(case_idx,2), ...
                Parameter.x_o(obs_idx),Parameter.y_o(obs_idx),Delta);
            remainder_min = find_poly_min(remainder_coeffs,0,Delta);
            remainders(obs_idx) = min(remainders(obs_idx),remainder_min);
        end
    end
end


function u = compute_SDP_SOSTLC(state,solver,ex,t)
    ref = straight_line_reference(t,[Parameter.x_0;Parameter.y_0], ...
                                      [Parameter.x_d;Parameter.y_d]);
    [A_tls,b_tls,u_nom] = get_TLS(state,ref);
    Delta = Parameter.tlc_t(ex);
    remainders = remainder_bounds(state,Delta);
    coeff_h = zeros(length(Parameter.r),4);
    coeff_hg = zeros(length(Parameter.r),2);
    for obs_idx = 1:length(Parameter.r)
        [V,LfV,Lf2V,Lf3V,LgLf2V] = compute_TLC_derivatives( ...
            state,Parameter.x_o(obs_idx),Parameter.y_o(obs_idx), ...
            Parameter.r(obs_idx),Parameter.delta);
        coeff_h(obs_idx,:) = [V,LfV,Lf2V,Lf3V];
        coeff_hg(obs_idx,:) = LgLf2V;
    end

    try
        candidate = compute_QP(state,@get_rTLC,ex,t);
        safe = true;
        for obs_idx = 1:length(Parameter.r)
            row = coeff_h(obs_idx,:);
            p = [(row(4)+coeff_hg(obs_idx,:)*candidate+remainders(obs_idx))/6, ...
                 row(3)/2,row(2),row(1)];
            safe = safe && find_poly_min(p,0,Delta) >= -1e-8;
        end
        if safe
            u = candidate;
            return;
        end
    catch ME
        if ~strcmp(ME.identifier,'compute_QP:Failed')
            rethrow(ME);
        end 
    end

    [solution,problem] = solver(state,Delta,A_tls,b_tls, ...
                                coeff_h,coeff_hg,remainders,u_nom);
    if problem ~= 0 || isempty(solution) || isempty(solution{1}) || ...
            any(~isfinite(solution{1}(:)))
        error('compute_SDP_SOSTLC:Failed', ...
              'SDP failed at t = %.4g s (problem=%d): %s', ...
              t,problem,yalmiperror(problem));
    end
    u = solution{1};
end


function [min_val,min_x] = find_poly_min(p,a,b)
    val_a = polyval(p,a);
    val_b = polyval(p,b);
    dp = polyder(p);
    roots_dp = roots(dp);
    real_roots = real(roots_dp(abs(imag(roots_dp)) < 1e-10));
    valid_roots = real_roots(real_roots >= a & real_roots <= b);
    all_x = [a;b;valid_roots(:)];
    all_vals = [val_a;val_b;polyval(p,valid_roots(:))];
    [min_val,min_idx] = min(all_vals);
    min_x = all_x(min_idx);
end
