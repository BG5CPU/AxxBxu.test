clear; close all; clc;

%% system settings ========================================================
% Time step [s]
h = 0.002;  

% mass
mass = 1.500;
% the inertia matrix with respect to the body-fixed frame
Jx = 0.029125; 
Jy = 0.029125; 
Jz = 0.055225;
Jinertia = blkdiag(Jx,Jy,Jz);


% system dimension
dim_x = 6;
dim_u = 3;
     


%% model-based controller 1 ===============================================

coeff = 50;

Ad0 = [eye(3), h*eye(3); zeros(3), eye(3)];
Bd0 = [zeros(3); eye(3)/Jinertia]*h* coeff;


pp1 = 0.8 + 0.2i;
pp2 = 0.8 - 0.2i;
pp3 = 0.9 + 0.1i;
pp4 = 0.9 - 0.1i;
pp5 = 0.95 + 0.1i;
pp6 = 0.95 - 0.1i;
% poles = [pp1, pp2, pp3, pp4, pp5, pp6];
poles = linspace(0.8,0.9,dim_x);
vK = -place(Ad0, Bd0, poles);

disp(vK);





%% controller implementation ==============================================

% state0_imp = ones(dim_x,1)*1;
state0_imp = [0.1 0.1 0.1 -0.1 -0.1 -0.1]';
tspan_imp = [0 10];

% Torque function that, state feedback
tau_fk = @(t,state) vK*state;
[t_imp, x_imp] = discrete_attitude_equation(Jinertia, state0_imp, tau_fk, 0, h, tspan_imp);
x_imp = x_imp';


figure('color', [1 1 1]);
set(gcf, 'Position', [300, 300, 450, 350]);

subplot(2,1,1); % 2 rows, 1 column, 1st subplot
hold on;
for i = 1:3
    plot(t_imp', x_imp(i,:), 'LineWidth', 1.5);
end
hold off;
legend('$\phi$', '$\theta$', '$\psi$', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('time [s]', 'Fontname', 'Times New Roman', 'FontSize', 11);
ylabel('Euler angles [rad]', 'Fontname', 'Times New Roman', 'FontSize', 11);
set(gca,'fontsize',11,'fontname','Times');
grid on;
set(gca,'position',[0.13 0.63 0.84 0.32]);

% Bottom subplot - last 3 curves
subplot(2,1,2); % 2 rows, 1 column, 2nd subplot
hold on;
for i = 4:6
    plot(t_imp', x_imp(i,:), 'LineWidth', 1.5);
end
hold off;
legend('$\omega_1$', '$\omega_2$', '$\omega_3$', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('time [s]', 'Fontname', 'Times New Roman', 'FontSize', 11);
ylabel('angular velocity [rad/s]', 'Fontname', 'Times New Roman', 'FontSize', 11);
set(gca,'fontsize',11,'fontname','Times');
grid on;
set(gca,'position',[0.13 0.15 0.84 0.32]);






%% useful functions =======================================================

function [t_discrete, state_discrete] = discrete_attitude_equation(I, state0, tau_fun, dis_fun, h, tspan)
    % Discrete-time Attitude Dynamics Solver (Forward Euler Method)
    % Solves coupled Euler-angle kinematics and Euler's rotation equations with disturbance inputs
    %
    % Inputs:
    %   I: 3x3 inertia matrix [kg*m^2] (diagonal or full, must be positive definite)
    %   state0: Initial state [phi; theta; psi; omega_x; omega_y; omega_z] (6x1)
    %           Angles in radians, angular rates in rad/s
    %   tau_fun: Torque input (function handle @(t,state) or constant 3x1 vector) [N*m]
    %   dis_fun: Disturbance input (function handle @(t,state) or constant 6x1 vector)
    %            Applied as additive noise to state derivatives [rad/s for angles, rad/s^2 for omega]
    %   h: Time step size [s] (scalar)
    %   tspan: Time range [t_start, t_end] [s]
    %
    % Outputs:
    %   t_discrete: Time vector (Nx1) [s]
    %   state_discrete: State history [phi,theta,psi,omega_x,omega_y,omega_z] (Nx6)
    %                   Angles in radians, angular rates in rad/s

    % Initialize time vector (column vector)
    t_discrete = (tspan(1):h:tspan(2))';
    N = length(t_discrete);
    state_discrete = zeros(N, 6);
    state_discrete(1,:) = state0(:)';  % Ensure row vector storage

    % Time-stepping loop
    for k = 1:N-1
        t_k = t_discrete(k);
        state_k = state_discrete(k,:)';  % Convert to column vector
        
        % Extract Euler angles and angular velocity
        angles_k = state_k(1:3);  % [phi (roll); theta (pitch); psi (yaw)] [rad]
        omega_k = state_k(4:6);   % Body-frame angular velocity [rad/s]
        
        % Process torque input (function or constant)
        if isa(tau_fun, 'function_handle')
            tau_k = tau_fun(t_k, state_k);  % Time/state-dependent torque
        else
            tau_k = tau_fun;                % Constant torque
        end

        % Process disturbance input (function or constant)
        if isa(dis_fun, 'function_handle')
            dis_k = dis_fun(t_k, state_k);  % Time/state-dependent disturbance
        else
            dis_k = dis_fun;                % Constant disturbance
        end
        
        % --- Kinematic Equations (ZYX Euler Angles) ---
        phi = angles_k(1);    % Roll angle [rad]
        theta = angles_k(2);  % Pitch angle [rad]
        
        % Singularity check (gimbal lock at ±90° pitch)
        if abs(cos(theta)) < 1e-2
            error('Gimbal lock: Pitch=±90°. Switch to quaternion representation.');
        end
        
        % Compute Euler angle derivatives
        phi_dot = omega_k(1) + sin(phi)*tan(theta)*omega_k(2) + cos(phi)*tan(theta)*omega_k(3);
        theta_dot = cos(phi)*omega_k(2) - sin(phi)*omega_k(3);
        psi_dot = (sin(phi)/cos(theta))*omega_k(2) + (cos(phi)/cos(theta))*omega_k(3);
        euler_dot = [phi_dot; theta_dot; psi_dot];  % [rad/s]
        
        % --- Dynamic Equations (Euler's Rotation Equations) ---
        omega_dot = I \ (tau_k - cross(omega_k, I*omega_k));  % [rad/s^2]
        
        % Forward Euler integration with disturbance
        state_discrete(k+1,:) = (state_k + h * [euler_dot; omega_dot] + dis_k)';
    end
end





function BarQ = funcBarQ(a1, a2, a3, ...
                         a14, a24_u, a24_l, ...
                         a15, a25_u, a25_l, a35_u, a35_l, a45_u, a45_l, a55_u, a55_l, ...
                         a16, a26_u, a26_l, a36_u, a36_l, a46_u, a46_l, a56_u, a56_l, ...
                         b1, b2, b3)
    % Generate all 512 (2^9) combinations of upper/lower variants
    % Inputs: 
    %   a* : Component matrices
    %   *_u : Upper bound variants
    %   *_l : Lower bound variants
    % Output:
    %   BarQ : 512x1 cell array of blkdiag matrices
    
    % Define all possible combinations (binary combinations for 9 parameters)
    combinations = dec2bin(0:511) - '0'; % 512x9 binary matrix
    
    % Preallocate output
    BarQ = cell(512, 1);
    
    % Generate each combination
    for i = 1:512
        % Select variants based on binary flags
        idx = combinations(i,:);
        
        % Choose a24 variant
        a24 = idx(1)*a24_u + (1-idx(1))*a24_l;
        
        % Choose a25/a35/a45/a55 variants
        a25 = idx(2)*a25_u + (1-idx(2))*a25_l;
        a35 = idx(3)*a35_u + (1-idx(3))*a35_l;
        a45 = idx(4)*a45_u + (1-idx(4))*a45_l;
        a55 = idx(5)*a55_u + (1-idx(5))*a55_l;
        
        % Choose a26/a36/a46/a56 variants
        a26 = idx(6)*a26_u + (1-idx(6))*a26_l;
        a36 = idx(7)*a36_u + (1-idx(7))*a36_l;
        a46 = idx(8)*a46_u + (1-idx(8))*a46_l;
        a56 = idx(9)*a56_u + (1-idx(9))*a56_l;
        
        % Build the block diagonal matrix
        BarQ{i} = blkdiag(a1, a2, a3, ...
                          [a14; a24], ...
                          [a15; a25; a35; a45; a55], ...
                          [a16; a26; a36; a46; a56], ...
                          b1, b2, b3);
    end
end


function e_max = maxEig( MatrixSym )
    e_max = max(eig(MatrixSym));
end

function e_min = minEig( MatrixSym )
    e_min = min(eig(MatrixSym));
end

function y = roundUpSignificant(x, n)
    magnitude = 10.^floor(log10(abs(x)) - n + 1);
    y = ceil(x ./ magnitude) .* magnitude;
    y(x == 0) = 0; 
end

function r = uniRand(a, b)
    if b <= a
        error('a < b is required.');
    end
    r = a + (b - a) * rand();
end

