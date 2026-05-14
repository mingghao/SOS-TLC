% 清理环境
clear; clc; close all;

% ==================== 1. 生成模拟数据 ====================
% 生成一条穿过障碍物的基准轨迹 (紫线)
t = linspace(0, 1, 100);
x_mag = linspace(-0.5, 7.5, 100);
y_mag = -0.5 * x_mag.^2 + 4 * x_mag + 1; % 抛物线模拟
y_mag = 5 * sin(x_mag*0.2) + 1; % 或者使用正弦波模拟

% 使用样条插值生成避障轨迹 (蓝线)
ctrl_x = [-0.5, 0, 1, 3, 5, 6.5, 7, 6];
ctrl_y = [ 0.5, 0,-1, 1.2, 3, 4.5, 5, 6];
t_ctrl = linspace(0, 1, length(ctrl_x));
x_blue = spline(t_ctrl, ctrl_x, t);
y_blue = spline(t_ctrl, ctrl_y, t);

% 障碍物参数 (红圈)
obs_x = 3.3; 
obs_y = 4.0;
obs_r = 2.0;

% ==================== 2. 开始绘图 ====================
figure('Position', [100, 100, 800, 600]);
hold on; grid on;

% --- 绘制障碍物 ---
theta_circle = linspace(0, 2*pi, 100);
circle_x = obs_r * cos(theta_circle) + obs_x;
circle_y = obs_r * sin(theta_circle) + obs_y;
h_obs = plot(circle_x, circle_y, 'r-', 'LineWidth', 4); % 粗红线

% --- 绘制不带CBF的轨迹 (紫线) ---
h_mag = plot(x_mag, y_mag, 'm-', 'LineWidth', 2);

% --- 绘制带CBF的避障轨迹 (蓝线) ---
h_blue = plot(x_blue, y_blue, 'b-', 'LineWidth', 3.5);

% ==================== 3. 绘制自定义位姿标记 ====================
% 在蓝线上均匀选取一些点来画标记
num_markers = 8; 
idx_markers = round(linspace(5, 95, num_markers)); % 选取索引

% 标记的几何尺寸 (根据你的图，形状类似 [—] )
L = 0.6; % 整体长度
H = 0.2; % 两端垂直小翼的长度
marker_color = [0.5, 0.5, 0]; % 橄榄绿色

for i = 1:num_markers
    idx = idx_markers(i);
    xc = x_blue(idx);
    yc = y_blue(idx);
    
    % 计算该点的切线方向 (朝向角)
    if idx < length(x_blue)
        dx = x_blue(idx+1) - x_blue(idx);
        dy = y_blue(idx+1) - y_blue(idx);
    else
        dx = x_blue(idx) - x_blue(idx-1);
        dy = y_blue(idx) - y_blue(idx-1);
    end
    theta = atan2(dy, dx);
    
    % 定义标准形状 (中心在0,0，朝向x轴正方向)
    % 分为三段线：主轴、左翼、右翼
    shape_x = [-L/2, L/2, NaN, -L/2, -L/2, NaN, L/2, L/2];
    shape_y = [0, 0, NaN, -H/2, H/2, NaN, -H/2, H/2];
    
    % 旋转矩阵
    R = [cos(theta), -sin(theta); 
         sin(theta),  cos(theta)];
     
    % 应用旋转和平移
    pts = R * [shape_x; shape_y];
    pts_x = pts(1, :) + xc;
    pts_y = pts(2, :) + yc;
    
    % 绘制形状外壳 (注意 HandleVisibility='off' 防止它们进入图例)
    plot(pts_x, pts_y, 'Color', marker_color, 'LineWidth', 3.5, 'HandleVisibility', 'off');
    
    % 绘制中心黑点
    plot(xc, yc, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'HandleVisibility', 'off');
end

% ==================== 4. 图表美化与格式设置 ====================
% 设置坐标轴标签 (图中y轴字体特别大)
xlabel('x (m)', 'FontSize', 18, 'FontName', 'Times New Roman');
ylabel('y (m)', 'FontSize', 28, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

% 调整坐标轴刻度和比例
axis equal; % 保证圆是圆的，比例协调
set(gca, 'FontSize', 14, 'LineWidth', 1);
xlim([-1, 8]);
ylim([-1, 6.5]);

% 添加图例 (严格匹配图中文字)
legend([h_mag, h_blue, h_obs], ...
    {'Sequential-QP Controller without CBF', ...
     'Sequential-CBF-QP Controller', ...
     'Obstacle'}, ...
    'Location', 'southeast', ...
    'FontSize', 16, 'FontName', 'Arial');

hold off;