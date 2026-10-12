clear;
clc;
set(0,'DefaultTextInterpreter','latex');

%% simulation time and state to be shown
sdp_solver =SDP_solver(Parameter);

%% rTLC
for ex=0:length(Parameter.delta)

    timestep = Parameter.dt;
    duration = 0:timestep:30;
    timesteps = length(duration);
    u = zeros([1, timesteps - 1]);
    state_log = zeros([2, timesteps]);
    state_log(:, 1) = [20; 90];
    solve_time = zeros([1, timesteps - 1]);
    %% paramters for QPs (not changed with time)
    H = 2 .* [1, 0; 0, Parameter.eps];

    %% main
    for i = 1:(timesteps - 1)
        if ex==0
            tic;
            u(i) = compute_QP(state_log(:,i),H,@get_HOCBF,ex);
            t=toc;
        elseif ex<=2
            tic
            u(i) = compute_QP(state_log(:,i),H,@get_ZOHTLC,ex);
            
            t=toc;
        elseif ex<=4
            tic;
            u(i) = compute_QP(state_log(:,i),H,@get_TLC,ex);
            t=toc;
        else
            tic;
            u(i) = compute_SDP_rTLC(state_log(:,i),H,sdp_solver,ex);
            t=toc;
        end
        solve_time(i) = t;
        [~, y] = ode45(@(t, y) acc(t, y, u(i)),...
                    [0, timestep], state_log(:, i));
        y = y';
        state_log(:, i+1) = y(:, end);
    end

    %% print info
    b = state_log(2,:)-Parameter.sigma;
    min_b = min(b);
    mean_time = mean(solve_time);
    sd_time = std(solve_time);
   
    disp([char(Parameter.exps(ex+1)), ': minmum b:', num2str(min_b),', mean time: ',num2str(mean_time),', standard deviation: ',num2str(sd_time)]);
    %% draw
    if ex==0
        figure('Position',[1000,231,1043,1007])
        ax1 = subplot(3,1,1);
        plot(ax1,[0 30], [Parameter.v_0,Parameter.v_0], 'k--',"LineWidth", 1.5);
        hold(ax1,'on');
        plot(ax1,[0 30], [Parameter.v_d,Parameter.v_d], 'k--',"LineWidth", 1.5);
        plot(ax1,duration, state_log(1, :),"LineWidth", 1.5,'Color',Parameter.color(ex+1));
        xlabel(ax1,'$t/s$','Fontsize',20);
        ylabel(ax1,'$v/(m/s)$','Fontsize',20)
        xlim(ax1,[0, 30]);
        ylim(ax1,[12, 25]);
        grid(ax1,"on");
  

        ax2 = subplot(3,1,2);
        plot(ax2,[0 30], [-0.7,-0.7], 'k--',"LineWidth", 1.5);
        hold(ax2,"on");
        plot(ax2,[0 30], [0.4,0.4], 'k--',"LineWidth", 1.5);
        plot(ax2,duration(1:end-1), u./Parameter.m./Parameter.g,"LineWidth", 1.5,'Color',Parameter.color(ex+1));
        xlabel(ax2,'$t/s$','Fontsize',20);
        ylabel(ax2,'$u/(Mg)$','Fontsize',20);
        xlim(ax2,[0, 30]);
        % ylim([-10, 100]);
        grid(ax2,"on");


        ax3 = subplot(3,1,3);
        plot(ax3, [0 30], [0.0, 0.0], 'k--',"LineWidth", 1.5);
        hold(ax3,"on");
       
        plot(ax3, duration, state_log(2,:)-Parameter.sigma,"LineWidth", 1.5,'Color',Parameter.color(ex+1));
        xlabel(ax3, '$t/s$','Fontsize',20);
        ylabel(ax3, '$h(x)$','Fontsize',20);
        xlim(ax3, [0, 30]);
        grid(ax3,"on");

        ax_in = axes('Position', [0.5, 0.19, 0.3, 0.12]); 
        plot(ax_in, [0 30], [0.0, 0.0], 'k--',"LineWidth", 1.5);
        hold(ax_in, 'on');
        plot(ax_in, duration, state_log(2,:)-Parameter.sigma,"LineWidth", 1.5,'Color',Parameter.color(ex+1));
        xlim(ax_in, [8, 11]); 
        ylim(ax_in, [-2, 5]); 
        grid(ax_in, 'on');
        box(ax_in, 'on');
    else

        plot(ax1, duration, state_log(1, :),"LineWidth", 1.5,'Color',Parameter.color(ex+1));
        plot(ax2, duration(1:end-1), u./Parameter.m./Parameter.g,"LineWidth", 1.5,'Color',Parameter.color(ex+1));
        plot(ax3, duration, state_log(2,:)-Parameter.sigma,"LineWidth", 1.5,'Color',Parameter.color(ex+1));
        plot(ax_in, duration, state_log(2,:)-Parameter.sigma,"LineWidth", 1.5,'Color',Parameter.color(ex+1));
        
    end
    figure(1)
    subplot(3,1,1)
    legend('','','$HOCBF$','$ZOH-TLC \quad \Delta t=1.5$','$ZOH-TLC \quad \Delta t=3$','$SD-SOS-TLC(QP) \quad \Delta t=1.5$','$SD-SOS-TLC(QP) \quad \Delta t=3$','$SD-SOS-TLC(SDP) \quad \Delta t=1.5$','$SD-SOS-TLC(SDP) \quad \Delta t=3$','Interpreter', 'latex','Fontsize',13)
    legend('Box','off')
    annotation('arrow', [0.38, 0.48], [0.17, 0.22], ...
            'Color', 'k', ...             
            'LineWidth', 1.5, ...        
            'HeadStyle', 'vback2', ...    
            'HeadLength', 10, ...        
            'HeadWidth', 8); 
end




%% functions
function u = compute_QP(cur_state, H, safe_func, ex)
    cur_vel = cur_state(1);
    cur_distance = cur_state(2);
    F_r = get_Fr(cur_vel);
    [A_tls, b_tls] = get_TLS(cur_vel, F_r);
    [A_tlc, b_tlc] = safe_func(F_r, cur_distance, cur_vel,ex);
    F = -2 .* [F_r/Parameter.m; 0];
    A = [A_tls; A_tlc; 1, 0;-1, 0];
    b = [b_tls; b_tlc; Parameter.g * Parameter.c_a; Parameter.g * Parameter.c_d];
    [solution,~,~,~,~] = quadprog(H,F,A,b,[],[],[],[],[],Parameter.options);
    u = solution(1)*Parameter.m;
end

function F_r = get_Fr(cur_vel)
    F_r = Parameter.f_0.*sign(cur_vel) + Parameter.f_1 .* cur_vel + ...
            Parameter.f_2 .* cur_vel .^ 2;
end

function [A_tls, b_tls] = get_TLS(cur_vel, F_r)
    v_bias = cur_vel - Parameter.v_d;
    A_tls = [2*v_bias*0.1, -1];
    b_tls = 2*v_bias*F_r*0.1/Parameter.m - v_bias^2;
end


function [A_tlc, b_tlc] = get_TLC(F_r, cur_distance, cur_vel,ex)  
    b = cur_distance - Parameter.sigma;
    b_dot = Parameter.v_0 - cur_vel;

    remainder = (Parameter.f_1*(-Parameter.g * Parameter.c_d-F_r/Parameter.m)*Parameter.delta(ex)-Parameter.f_2*cur_vel^2)/Parameter.m;
    mu_bound = (-b - b_dot*Parameter.delta(ex))/Parameter.delta(ex)^2-remainder;
    p = [mu_bound+remainder, b_dot, b];
    [min_val,~] = find_poly_min(p, 0, Parameter.delta(ex));
    if min_val>=0
        A_tlc = [1, 0];
        b_tlc = F_r/Parameter.m + 2*b/Parameter.delta(ex)^2 + 2*b_dot/Parameter.delta(ex) + remainder;
    else
        A_tlc = [b, 0];
        b_tlc = -b_dot^2./2 + b*(remainder + F_r / Parameter.m);
    end

end


function [A_tlc, b_tlc] = get_ZOHTLC(F_r, cur_distance, cur_vel,ex) 
    b = cur_distance - Parameter.sigma;
    b_dot = Parameter.v_0 - cur_vel;
    A_tlc = [1, 0];
    b_tlc = F_r/Parameter.m + 2*b/Parameter.delta(ex)^2 + 2*b_dot/Parameter.delta(ex);
end


function [A_tlc, b_tlc] = get_HOCBF(F_r, cur_distance, cur_vel, ex) 
    b = cur_distance - Parameter.sigma;
    b_dot = Parameter.v_0 - cur_vel;
    A_tlc = [1, 0];
    b_tlc = F_r/Parameter.m + 1*b + 2*b_dot;
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

function u = compute_SDP_rTLC(cur_state,H,solver,ex)
    
    cur_vel = cur_state(1);
    cur_distance = cur_state(2);
    F_r = get_Fr(cur_vel);
    b_val = cur_distance - Parameter.sigma;
    b_dot = Parameter.v_0 - cur_vel;
    remainder = (Parameter.f_1*(-Parameter.g * Parameter.c_d-F_r/Parameter.m)*Parameter.delta(ex)-Parameter.f_2*cur_vel^2)/Parameter.m;
    mu_bound = (-b_val - b_dot*Parameter.delta(ex))/Parameter.delta(ex)^2-remainder;
    p = [mu_bound+remainder, b_dot, b_val];
    [min_val,~] = find_poly_min(p, 0, Parameter.delta(ex));
    [A_tls, b_tls] = get_TLS(cur_vel, F_r);
    A_tlc = [1, 0];
    b_tlc = F_r/Parameter.m + 2*b_val/Parameter.delta(ex)^2 + 2*b_dot/Parameter.delta(ex) + remainder;
    F = -2 .* [F_r/Parameter.m; 0];
    A = [A_tls; A_tlc; 1, 0;-1, 0];
    b = [b_tls; b_tlc; Parameter.g * Parameter.c_a; Parameter.g * Parameter.c_d];

    if min_val>=0
        [solution,~,~,~,~] = quadprog(H,F,A,b,[],[],[],[],[],Parameter.options);
        u = solution(1)*Parameter.m;
    else
        [solution, ~] = solver([cur_vel; cur_vel^2; cur_distance; F_r; Parameter.delta(ex); Parameter.delta(ex)^2]);
        u = solution * Parameter.m;
    end
  
end