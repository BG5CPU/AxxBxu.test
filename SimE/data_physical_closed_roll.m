clear; clc; close all;

%% import data ============================================================

% file path A D F
file1 = 'physical_data_csv/physical_bad_controller_to_vK_closed_1007/imu_data.csv';
file2 = 'physical_data_csv/physical_bad_controller_to_vK_closed_1007/attitude_r_p_y.csv';
file3 = 'physical_data_csv/physical_bad_controller_to_vK_closed_1007/actuator_control.csv';
file4 = 'physical_data_csv/physical_bad_controller_to_vK_closed_1007/input_x_y_z.csv';
file5 = 'physical_data_csv/physical_bad_controller_to_vK_closed_1007/pose.csv';


% read (imu_data.csv)
opts = detectImportOptions(file1);
opts.SelectedVariableNames = [1,2,3, 5,6,7,8, 18,19,20, 30,31,32]; 
opts.VariableNamingRule = 'preserve';
imu_data = readtable(file1, opts);
matrix_imu_data = table2array(imu_data); 
check_data_continuity(matrix_imu_data(:,2), 'matrix_imu_data');


% read (attitude_r_p_y.csv)
opts = detectImportOptions(file2);
opts.SelectedVariableNames = [1,2,3, 5,6,7]; 
opts.VariableNamingRule = 'preserve';
attitude_r_p_y = readtable(file2, opts);
matrix_attitude_r_p_y = table2array(attitude_r_p_y); 
check_data_continuity(matrix_attitude_r_p_y(:,2), 'matrix_attitude_r_p_y');


% read (actuator_control.csv)
opts = detectImportOptions(file3);
opts.SelectedVariableNames = [1,2,3, 6,7,8,9]; 
opts.VariableNamingRule = 'preserve';
actuator_control = readtable(file3, opts);
matrix_actuator_control = table2array(actuator_control); 
check_data_continuity(matrix_actuator_control(:,2), 'matrix_actuator_control');


% read (input_x_y_z.csv)
opts = detectImportOptions(file4);
opts.SelectedVariableNames = [1,2,3, 5,6,7]; 
opts.VariableNamingRule = 'preserve';
input_x_y_z = readtable(file4, opts);
matrix_input_x_y_z = table2array(input_x_y_z); 
check_data_continuity(matrix_input_x_y_z(:,2), 'matrix_input_x_y_z');


% read (pose.csv)
opts = detectImportOptions(file5);
opts.SelectedVariableNames = [1,2,3, 5,6,7, 8,9,10,11]; 
opts.VariableNamingRule = 'preserve';
pose = readtable(file5, opts);
matrix_pose = table2array(pose); 
check_data_continuity(matrix_pose(:,2), 'matrix_pose');








%% data alignment =========================================================

figure(10); hold on;
plot(matrix_attitude_r_p_y(:,4), 'r'); 
plot(matrix_attitude_r_p_y(:,5), 'g');
plot(matrix_attitude_r_p_y(:,6), 'b');
hold off;



[~, idx_matrix_attitude_r_p_y_1] = unique(matrix_attitude_r_p_y(:,1), 'stable');
[~, idx_matrix_attitude_r_p_y_3] = unique(matrix_attitude_r_p_y(:,3), 'stable');

matrix_attitude_r_p_y = matrix_attitude_r_p_y(idx_matrix_attitude_r_p_y_3, :);



matrix_imu_data = matrix_imu_data(32:end,:);

figure(1); hold on;
plot(flipud(matrix_imu_data(:,3)), 'r.'); % flipud
plot(flipud(matrix_attitude_r_p_y(:,3)), 'g.');
hold off;

time = matrix_attitude_r_p_y(:,3);

velocity_angle_x = interpolate_to_match_A(time, matrix_imu_data(:,[1,8]));
velocity_angle_y = interpolate_to_match_A(time, matrix_imu_data(:,[1,9]));
velocity_angle_z = interpolate_to_match_A(time, matrix_imu_data(:,[1,10]));

velocity_angle = [velocity_angle_x(:,2), velocity_angle_y(:,2), velocity_angle_z(:,2)];

pose_angle = matrix_attitude_r_p_y(:,[4,5,6]);
pose_angle(:,1) = pose_angle(:,1)-0.059; % sensor bias
pose_angle(:,3) = pose_angle(:,3)+0.032; % sensor bias

Pstar = 6500;
Pend = 25000;

time = time(Pstar:Pend, :);
pose_angle = pose_angle(Pstar:Pend, :);

velocity_angle = velocity_angle(Pstar:Pend, :);

time = (time-time(1))/1e9;

downs = 20;

time = downsample_matrix(time, downs);
pose_angle = downsample_matrix(pose_angle, downs);
velocity_angle = downsample_matrix(velocity_angle, downs);









%% plot data ==============================================================



lm=0.095; rm=0.025; bm=0.13; tm=0.04; hs=0.065; vs=0.08;
w=(1-lm-rm-2*hs)/3;
h=(1-tm-bm-vs)/2;


figure('color', [1 1 1]);
set(gcf, 'Position', [200, 200, 600, 300]);

subplot(2,3,1)
plot(time, pose_angle(:,1), 'k', 'LineWidth', 0.5);
grid on;
xlim([0 38]);
ylim([-0.3 0.2]);
set(gca,'fontsize',10,'fontname','Times');
legend('$\phi$', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 10);
% xlabel('time [s]', 'Fontname', 'Times New Roman', 'FontSize', 10);
% ylabel('$\phi$  [rad]', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('[rad]', 'Fontname', 'Times New Roman', 'FontSize', 10);
set(gca,'Position',[lm, bm+h+vs, w, h]);


subplot(2,3,2)
plot(time, pose_angle(:,2), 'k', 'LineWidth', 0.5);
grid on;
xlim([0 38]);
ylim([-0.3 0.2]);
set(gca,'fontsize',10,'fontname','Times');
legend('$\theta$', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 10);
% xlabel('time [s]', 'Fontname', 'Times New Roman', 'FontSize', 10);
% ylabel('$\theta$  [rad]', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 10);
set(gca,'Position',[lm+w+hs, bm+h+vs, w, h]);


subplot(2,3,3)
plot(time, pose_angle(:,3), 'k', 'LineWidth', 0.5);
grid on;
xlim([0 38]);
ylim([-0.3 0.2]);
set(gca,'fontsize',10,'fontname','Times');
legend('$\psi$', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 10);
% xlabel('time [s]', 'Fontname', 'Times New Roman', 'FontSize', 10);
% ylabel('$\psi$  [rad]', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 10);
set(gca,'Position',[lm+2*(w+hs), bm+h+vs, w, h]);




subplot(2,3,4)
plot(time, velocity_angle(:,1), 'k', 'LineWidth', 0.5);
grid on;
xlim([0 38]);
ylim([-1 1]);
set(gca,'fontsize',10,'fontname','Times');
legend('$\omega_1$', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 10);
xlabel('time [s]', 'Fontname', 'Times New Roman', 'FontSize', 10);
ylabel('[rad/s]', 'Fontname', 'Times New Roman', 'FontSize', 10);
set(gca,'Position',[lm, bm, w, h]);


subplot(2,3,5)
plot(time, velocity_angle(:,2), 'k', 'LineWidth', 0.5);
grid on;
xlim([0 38]);
ylim([-1 1]);
set(gca,'fontsize',10,'fontname','Times');
legend('$\omega_2$', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 10);
xlabel('time [s]', 'Fontname', 'Times New Roman', 'FontSize', 10);
set(gca,'Position',[lm+w+hs, bm, w, h]);


subplot(2,3,6)
plot(time, velocity_angle(:,3), 'k', 'LineWidth', 0.5);
grid on;
xlim([0 38]);
ylim([-1 1]);
set(gca,'fontsize',10,'fontname','Times');
legend('$\omega_3$', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 10);
xlabel('time [s]', 'Fontname', 'Times New Roman', 'FontSize', 10);
set(gca,'Position',[lm+2*(w+hs), bm, w, h]);

















%% useful functions =======================================================
function check_data_continuity(data_matrix_column, matrix_name)

    fprintf('check "%s" continuity...\n', matrix_name);

    differences = diff(data_matrix_column);
    
    discontinuity_indices = find(differences ~= 1);
    
    if isempty(discontinuity_indices)
        fprintf('data packages complete. \n');
    else
        fprintf('data packages incomplete: %d. \n', length(discontinuity_indices));
        
        for i = 1:length(discontinuity_indices)
            idx = discontinuity_indices(i);
            fprintf('at %d -> %d: diff = %d \n', idx, idx+1, differences(idx));
        end
        
        total_missing = sum(differences(discontinuity_indices));
        fprintf('%d packages lost. \n', total_missing);
    end
    
    fprintf('\n');
end



function interpolated_data = interpolate_to_match_A(A_timestamps, B_data)
   
    B_time = B_data(:, 1);
    B_values = B_data(:, 2);
    
    if any(diff(B_time) <= 0)
        error('B matrix timestamp must be strictly incremented');
    end
    
    if any(diff(A_timestamps) <= 0)
        error('A matrix timestamp must be strictly incremented');
    end
    

    if min(A_timestamps) < min(B_time) && max(A_timestamps) > max(B_time)
        disp('timestamp of A exceeds B, and extrapolation should be used');
    end
    
    % 'linear', 'nearest', 'next', 'previous', 
    % 'pchip', 'cubic', 'v5cubic', 'makima', 'spline'
    interpolated_values = interp1(B_time, B_values, A_timestamps, 'linear', 'extrap');
    
    interpolated_data = [A_timestamps, interpolated_values];
end



function Y = downsample_matrix(X,d)
% DOWNSAMPLE_MATRIX  Downsample matrix X along rows
    [N,~] = size(X);
    N_valid = floor(N/d) * d;
    Y = X(1:d:N_valid, :);
end



