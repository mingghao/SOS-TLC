clear;
clc;
set(0,'DefaultTextInterpreter','latex');

%% simulation time and state to be shown


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
    H = 2 .* [1 ./ Parameter.m .^ 2, 0; 0, Parameter.eps];

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
        else
            tic;
            u(i) = compute_QP(state_log(:,i),H,@get_TLC,ex);
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
        subplot(3,1,1)
        plot([0 30], [Parameter.v_0,Parameter.v_0], 'k--',"LineWidth", 1.5);
        hold on;
        plot([0 30], [Parameter.v_d,Parameter.v_d], 'k--',"LineWidth", 1.5);
        plot(duration, state_log(1, :),"LineWidth", 1.5,'Color',Parameter.color(ex+1));
        xlabel('$t/s$','Fontsize',20);
        ylabel('$v/(m/s)$','Fontsize',20)
        xlim([0, 30]);
        ylim([12, 25]);
        grid on
        hold off;

        subplot(3,1,2)
        plot([0 30], [-0.7,-0.7], 'k--',"LineWidth", 1.5);
        hold on;
        plot([0 30], [0.4,0.4], 'k--',"LineWidth", 1.5);
        plot(duration(1:end-1), u./Parameter.m./Parameter.g,"LineWidth", 1.5,'Color',Parameter.color(ex+1));
        xlabel('$t/s$','Fontsize',20);
        ylabel('$u/(Mg)$','Fontsize',20);
        xlim([0, 30]);
        % ylim([-10, 100]);
        grid on;
        hold off;

        subplot(3,1,3)
        plot([0 30], [0.0, 0.0], 'k--',"LineWidth", 1.5);
        hold on;
        plot(duration, state_log(2,:)-Parameter.sigma,"LineWidth", 1.5,'Color',Parameter.color(ex+1));
        xlabel('$t/s$','Fontsize',20);
        ylabel('$h(x)$','Fontsize',20);
        xlim([0, 30]);
        grid on;
        hold off;
    else
        figure(1)
        subplot(3,1,1)
        hold on;
        plot(duration, state_log(1, :),"LineWidth", 1.5,'Color',Parameter.color(ex+1));
        hold off;

        subplot(3,1,2)
        hold on;
        plot(duration(1:end-1), u./Parameter.m./Parameter.g,"LineWidth", 1.5,'Color',Parameter.color(ex+1));
        hold off;

        subplot(3,1,3)
        hold on;
        plot(duration, state_log(2,:)-Parameter.sigma,"LineWidth", 1.5,'Color',Parameter.color(ex+1));
        hold off;
        
    end
    figure(1)
    subplot(3,1,1)
    legend('','','$HOCBF$','$ZOH-TLC \quad \Delta t=1$','$ZOH-TLC \quad \Delta t=3$','$rTLC \quad \Delta t=1$','$rTLC \quad \Delta t=3$','Interpreter', 'latex','Fontsize',18)
end


%% HOCBF














%% functions
function u = compute_QP(cur_state, H, safe_func, ex)
    cur_vel = cur_state(1);
    cur_distance = cur_state(2);
    F_r = get_Fr(cur_vel);
    [A_tls, b_tls] = get_TLS(cur_vel, F_r);
    [A_tlc, b_tlc] = safe_func(F_r, cur_distance, cur_vel,ex);
    F = -2 .* [F_r / Parameter.m ^ 2; 0];
    A = [A_tls; A_tlc; 1, 0;-1, 0];
    b = [b_tls; b_tlc; Parameter.m * Parameter.g * Parameter.c_a; Parameter.m * Parameter.g * Parameter.c_d];
    [solution,~,~,~,~] = quadprog(H,F,A,b,[],[],[],[],[],Parameter.options);
    u = solution(1);
end

function F_r = get_Fr(cur_vel)
    F_r = Parameter.f_0.*sign(cur_vel) + Parameter.f_1 .* cur_vel + ...
            Parameter.f_2 .* cur_vel .^ 2;
end

function [A_tls, b_tls] = get_TLS(cur_vel, F_r)
    v_bias = cur_vel - Parameter.v_d;
    A_tls = [2*v_bias*0.1/Parameter.m, -1];
    b_tls = 2*v_bias*F_r*0.1/Parameter.m - v_bias^2;
end

function [A_tlc, b_tlc] = get_TLC(F_r, cur_distance, cur_vel,ex)  
    b = cur_distance - Parameter.sigma;
    b_dot = Parameter.v_0 - cur_vel;

    remainder = (Parameter.f_1*(-Parameter.g * Parameter.c_d-F_r/Parameter.m)*Parameter.delta(ex)-Parameter.f_2*cur_vel^2)/Parameter.m;

    if b_dot>=0
        A_tlc = [1/Parameter.m, 0];
        b_tlc = F_r/Parameter.m + 2*b/Parameter.delta(ex)^2 + 2*b_dot/Parameter.delta(ex) + remainder;
    else
        mu_bound = -2*b/Parameter.delta(ex)^2 - 2*b_dot/Parameter.delta(ex)-remainder;
        if mu_bound > 0
            if -b_dot/mu_bound > Parameter.delta(ex)
                A_tlc = [1/Parameter.m, 0];
                b_tlc = F_r/Parameter.m + 2*b/Parameter.delta(ex)^2 + 2*b_dot/Parameter.delta(ex) + remainder;
            else
                A_tlc = [(2*b + Parameter.delta(ex)^2*remainder) / Parameter.m, 0];
                b_tlc = -b_dot^2+(2*b + Parameter.delta(ex)^2*remainder)*F_r / Parameter.m;
            end
        else
            A_tlc = [1/Parameter.m, 0];
            b_tlc = F_r/Parameter.m + 2*b/Parameter.delta(ex)^2 + 2*b_dot/Parameter.delta(ex) + remainder;
        end
    end
end


function [A_tlc, b_tlc] = get_ZOHTLC(F_r, cur_distance, cur_vel,ex) 
    b = cur_distance - Parameter.sigma;
    b_dot = Parameter.v_0 - cur_vel;
    A_tlc = [1/Parameter.m, 0];
    b_tlc = F_r/Parameter.m + 2*b/Parameter.delta(ex)^2 + 2*b_dot/Parameter.delta(ex);
end


function [A_tlc, b_tlc] = get_HOCBF(F_r, cur_distance, cur_vel,ex) 
    b = cur_distance - Parameter.sigma;
    b_dot = Parameter.v_0 - cur_vel;
    A_tlc = [1/Parameter.m, 0];
    b_tlc = F_r/Parameter.m + 1*b + 2*b_dot;
end