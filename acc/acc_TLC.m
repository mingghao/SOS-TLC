clear;
clc;
set(0,'DefaultTextInterpreter','latex');

%% simulation time and state to be shown
timestep = Parameter.dt;
duration = 0:timestep:30;
timesteps = length(duration);

u = zeros([1, timesteps - 1]);
state_log = zeros([2, timesteps]);
state_log(:, 1) = [20; 90];

%% paramters for QPs (not changed with time)
H = 2 .* [1 ./ Parameter.m .^ 2, 0; 0, Parameter.eps];

%% main
for i = 1:(timesteps - 1)
    u(i) = compute_QP(state_log(:,i),H);
    [~, y] = ode45(@(t, y) acc(t, y, u(i)),...
                    [0, timestep], state_log(:, i));
    y = y';
    state_log(:, i+1) = y(:, end);
end



%% draw 
figure(1);
b1 = state_log(1, :);
y1 = plot(duration, b1,"LineWidth", 1);
hold on;
plot([0 30], [Parameter.v_0,Parameter.v_0], '--',"LineWidth", 0.7);
% legend([y1, y2],'Linear(p=1)', 'Quadratic(p=0.02)');
xlabel('$t/s$');
ylabel('$v(x(t))$')
title('v-t figure');
% xlim([0, 30]);
% ylim([-10, 100]);



figure(2);
y1 = plot(duration(1:end-1), u./Parameter.m./Parameter.g,"LineWidth", 1);
hold on;
plot([0 30], [-0.7,-0.7], '--',"LineWidth", 0.7);
plot([0 30], [0.4,0.4], '--',"LineWidth", 0.4);
% legend([y1, y2],'Linear(p=1)', 'Quadratic(p=0.02)');
xlabel('$t/s$');
ylabel('$u(x(t))$')
title('u-t figure');
% xlim([0, 30]);
% ylim([-10, 100]);

figure(3)
y1 = plot(duration, state_log(2,:)-Parameter.sigma,"LineWidth", 1);
hold on;

%% functions
function u = compute_QP(cur_state, H)
    cur_vel = cur_state(1);
    cur_distance = cur_state(2);
    F_r = get_Fr(cur_vel);
    [A_tls, b_tls] = get_TLS(cur_vel, F_r);
    [A_tlc, b_tlc] = get_TLC(F_r, cur_distance, cur_vel);
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
    A_tls = [2*v_bias*1/Parameter.m, -1];
    b_tls = 2*v_bias*F_r*1/Parameter.m - v_bias^2;
end

function [A_tlc, b_tlc] = get_TLC(F_r, cur_distance, cur_vel)  
    b = cur_distance - Parameter.sigma;
    b_dot = Parameter.v_0 - cur_vel;

    remainder = (Parameter.f_1*(-Parameter.g * Parameter.c_d-F_r/Parameter.m)*Parameter.delta-Parameter.f_2*cur_vel^2)/Parameter.m;
    if b_dot>=0
        A_tlc = [1/Parameter.m, 0];
        b_tlc = F_r/Parameter.m + 2*b/Parameter.delta^2 + 2*b_dot/Parameter.delta + remainder;
    else
        mu_bound = -2*b/Parameter.delta^2 - 2*b_dot/Parameter.delta ;
        if mu_bound > 0
            A_tlc = [(2*b + Parameter.delta^2*remainder) / Parameter.m, 0];
            b_tlc = -b_dot^2-(2*b + Parameter.delta^2*remainder)*F_r / Parameter.m;
        else
            A_tlc = [1/Parameter.m, 0];
            b_tlc = F_r/Parameter.m + 2*b/Parameter.delta^2 + 2*b_dot/Parameter.delta + remainder;
        end
    end
end
