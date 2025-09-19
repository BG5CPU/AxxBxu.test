clear; clc;

rr = 0.83;

% syms x
% eqn = exp(x) + 2*x == 0;
% solv = solve(eqn,x);
% x1so = -lambertw(0, 1/2);
% checkso = exp(x1so)+2*x1so
% x2sq = rr^2 - x1so^2;
% op1 = exp(x1so)-x2sq;
% op2 = exp(-1);


% plot figure
figure;

% x1 = linspace(-rr, rr, 100);
% x2 = linspace(-rr, rr, 100);
% [XX1, XX2] = meshgrid(x1, x2);
% Z = 1 + 0.45*exp(XX1) - 0.45*XX2.^2;
% surf(XX1, XX2, Z);
% hold on;
% xlabel('X1 axis');
% ylabel('X2 axis');
% zlabel('Z axis');
% colormap jet; 
% colorbar;
% shading interp;


t = linspace(0, 2*pi, 1000);
Xc1 = rr*sin(t);
Xc2 = rr*cos(t);
Zc  = 1 + 0.45*exp(Xc1) - 0.45*Xc2.^2;
plot3(Xc1, Xc2, Zc, 'LineWidth', 2);
hold on;


xdot11 = -rr;
xdot12 = 0;
Zdot1 = 1 + 0.45*exp(xdot11) - 0.45*xdot12.^2;
plot3(xdot11, xdot12, Zdot1, 'o'); 
hold on;

% this is the minimum point
xdot21 = -lambertw(0, 1/2);
xdot22 = -sqrt(rr^2 - xdot21^2);
Zdot2 = 1 + 0.45*exp(xdot21) - 0.45*xdot22.^2;
plot3(xdot21, xdot22, Zdot2, '*'); 
hold on;




a = 6; % 示例：a > pi
x = linspace(-a, a, 1000);
x(x == 0) = eps; % 避免除以0
f = sin(x) ./ x;
f(x == 0) = 1; % 补充定义

disp(['最小值: ', num2str(min(f))]);
disp(['最大值: ', num2str(max(f))]);

% 绘图
figure;
plot(x, f, 'b-', 'LineWidth', 1.5);
xlabel('x'); ylabel('sin(x)/x');
title('函数 sin(x)/x 在 [-a, a] 上的行为');
grid on;











