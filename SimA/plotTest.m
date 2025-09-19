
clc; clear; close all;

%% ========================================================================

r = 0.9;

% ball data
theta = linspace(0, 2*pi, 300);
x_circle = r * cos(theta);
y_circle = r * sin(theta);


% ellipsoid data
M = [3 1; 1 2];

[V, D] = eig(M);
a = 1 / sqrt(D(1,1));
b = 1 / sqrt(D(2,2));
angle = atan2(V(2,1), V(1,1));

t = linspace(0, 2*pi, 300);
x_ellipse = a * cos(t);
y_ellipse = b * sin(t);
R = [cos(angle), -sin(angle); sin(angle), cos(angle)];
ellipse_pts = R * [x_ellipse; y_ellipse];

% plot

figure('color', [1 1 1]);
set(gcf, 'Position', [100, 100, 400, 380]);

xlabel('$x_1$', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x_2$', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 12);
axis equal;
axis([-2, 2, -2, 2]);
set(gca,'fontsize',12,'fontname','Times');

hold on; axis equal;
plot(x_circle, y_circle, 'r', 'LineWidth', 2); 
plot(ellipse_pts(1,:), ellipse_pts(2,:), 'k', 'LineWidth', 2);

% 
% grid on;
box on;
% legend('$[x_1^{+}\!\!-\!x_1 \ \ x_2^{+}\!\!-\!x_2]$', ...
%        'region of attraction',...
%        'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 12);

hold off;



%% ========================================================================

clc; clear; close all;

% 定义系统矩阵 A
A = [-1 2; -2 -1]; % 你可以更改这个矩阵来修改系统的动态特性

% 生成相空间网格
[x1, x2] = meshgrid(linspace(-2, 2, 20), linspace(-2, 2, 20)); % 20x20 的网格点
dx1 = A(1,1) * x1 + A(1,2) * x2; % 计算 dx1/dt
dx2 = A(2,1) * x1 + A(2,2) * x2; % 计算 dx2/dt

% 画相平面方向场
figure; hold on;
quiver(x1, x2, dx1, dx2, 'b', 'LineWidth', 1, 'AutoScaleFactor', 0.8);

% 设置初始条件
initial_conditions = [-1.5, -1.5; 1.5, 1.5; -1.5, 1.5; 1.5, -1.5; 0, 1; 1, 0]; 

% 时间范围
tspan = [0, 5]; 

% 画出不同初始条件的轨迹
for i = 1:size(initial_conditions, 1)
    [t, x] = ode45(@(t, x) A * x, tspan, initial_conditions(i, :)');
    plot(x(:,1), x(:,2), 'r', 'LineWidth', 1.5);
    plot(x(1,1), x(1,2), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r'); % 标记起点
end

% 设定坐标轴范围
axis([-2, 2, -2, 2]);
axis equal;
box on;
grid on;

% 添加标签
xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 14);
title('Phase Portrait', 'FontSize', 14);
legend('Vector Field', 'Trajectories', 'Location', 'best');

hold off;











