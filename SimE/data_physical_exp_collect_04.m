clear; clc; close all;

%% import data ============================================================

% file path
file1 = 'physical_data_csv/physical_mpc_sin_exp_collection/imu_data.csv';
file2 = 'physical_data_csv/physical_mpc_sin_exp_collection/attitude_r_p_y.csv';
file3 = 'physical_data_csv/physical_mpc_sin_exp_collection/actuator_control.csv';
file4 = 'physical_data_csv/physical_mpc_sin_exp_collection/input_x_y_z.csv';
file5 = 'physical_data_csv/physical_mpc_sin_exp_collection/pose.csv';


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


%% time difference ========================================================

figure(1); hold on;
plot(diff(matrix_input_x_y_z(:,1)), 'r.');
plot(diff(matrix_input_x_y_z(:,3)), 'g.');
hold off;

[~, idx_matrix_input_x_y_z_1] = unique(matrix_input_x_y_z(:,1), 'stable');
[~, idx_matrix_input_x_y_z_3] = unique(matrix_input_x_y_z(:,3), 'stable');


figure(2); hold on;
plot(diff(matrix_actuator_control(:,1)), 'r.');
plot(diff(matrix_actuator_control(:,3)), 'g.');
hold off;

[~, idx_matrix_actuator_control_1] = unique(matrix_actuator_control(:,1), 'stable');
[~, idx_matrix_actuator_control_3] = unique(matrix_actuator_control(:,3), 'stable');

matrix_actuator_control = matrix_actuator_control(idx_matrix_actuator_control_3, :);

% matrix_input_x_y_z(14531,:) = [];
% matrix_actuator_control = matrix_actuator_control(3:end,:);

figure(3); hold on;
plot(matrix_actuator_control(:,3), 'r.');
plot(matrix_input_x_y_z(:,3), 'g.');
hold off;

figure(4)
plot(matrix_actuator_control(:,6)-matrix_input_x_y_z(:,6),'r.');


matrix_imu_data = matrix_imu_data(29:end,:);
figure(5); hold on;
plot(flipud(matrix_actuator_control(:,3)), 'r.'); % flipud
plot(flipud(matrix_imu_data(:,3)), 'g.');
hold off;


% matrix_attitude_r_p_y(14531,:) = [];
figure(6); hold on;
plot((matrix_input_x_y_z(:,3)), 'r.');
plot((matrix_attitude_r_p_y(:,3)), 'g.');
hold off;


Tlen_matrix_input_x_y_z = -matrix_input_x_y_z(1,3)+matrix_input_x_y_z(end,3);
Nlen_matrix_input_x_y_z = round(Tlen_matrix_input_x_y_z/2e6);
time = matrix_input_x_y_z(1,3) + (0:Nlen_matrix_input_x_y_z-1)*2e6;
time = time';
figure(7); hold on;
plot((time), 'r.');
plot((matrix_input_x_y_z(:,3)), 'g.');
hold off;

% time = matrix_actuator_control(:,3);


pose_angle_x = interpolate_to_match_A(time, matrix_attitude_r_p_y(:,[3,4]));
pose_angle_y = interpolate_to_match_A(time, matrix_attitude_r_p_y(:,[3,5]));
pose_angle_z = interpolate_to_match_A(time, matrix_attitude_r_p_y(:,[3,6]));

% pose_angle_x = matrix_attitude_r_p_y(:,[3,4]);
% pose_angle_y = matrix_attitude_r_p_y(:,[3,5]);
% pose_angle_z = matrix_attitude_r_p_y(:,[3,6]);

velocity_angle_x = interpolate_to_match_A(time, matrix_imu_data(:,[3,8]));
velocity_angle_y = interpolate_to_match_A(time, matrix_imu_data(:,[3,9]));
velocity_angle_z = interpolate_to_match_A(time, matrix_imu_data(:,[3,10]));


figure(8); hold on;
plot(pose_angle_x(:,2), 'r');
plot(pose_angle_y(:,2), 'g');
plot(pose_angle_z(:,2), 'b');
hold off;

figure(9); hold on;
plot(velocity_angle_x(:,2), 'r');
plot(velocity_angle_y(:,2), 'g');
plot(velocity_angle_z(:,2), 'b');
hold off;


actuator_control_x = interpolate_to_match_A(time, matrix_input_x_y_z(:,[3,4]));
actuator_control_y = interpolate_to_match_A(time, matrix_input_x_y_z(:,[3,5]));
actuator_control_z = interpolate_to_match_A(time, matrix_input_x_y_z(:,[3,6]));
% actuator_control_x = matrix_input_x_y_z(:,[3,4]);
% actuator_control_y = matrix_input_x_y_z(:,[3,5]);
% actuator_control_z = matrix_input_x_y_z(:,[3,6]);


figure(10); hold on;
plot(actuator_control_x(:,2), 'r');
plot(actuator_control_y(:,2), 'g');
plot(actuator_control_z(:,2), 'b');
hold off;



data_exp_09A = [actuator_control_x, actuator_control_y, actuator_control_z, ...
                pose_angle_x, pose_angle_y, pose_angle_z, ...
                velocity_angle_x, velocity_angle_y, velocity_angle_z];


Pstart01 = 3000;
Pend01 = 18500;
% data_exp_09A = data_exp_09A(Pstart01:Pend01 , :);
% save('./physical_data_mat/data_exp_09A.mat', 'data_exp_09A');






%% plot data ==============================================================



time = time(Pstart01:Pend01 , :);

time = (time-time(1))/1e9;

actuator_control_x = actuator_control_x(Pstart01:Pend01 , :);
actuator_control_y = actuator_control_y(Pstart01:Pend01 , :);
actuator_control_z = actuator_control_z(Pstart01:Pend01 , :);
pose_angle_x = pose_angle_x(Pstart01:Pend01 , :);
pose_angle_y = pose_angle_y(Pstart01:Pend01 , :);
pose_angle_z = pose_angle_z(Pstart01:Pend01 , :);
velocity_angle_x = velocity_angle_x(Pstart01:Pend01 , :);
velocity_angle_y = velocity_angle_y(Pstart01:Pend01 , :);
velocity_angle_z = velocity_angle_z(Pstart01:Pend01 , :);


down = 10;


time = downsample_matrix(time, down);

actuator_control_x = downsample_matrix(actuator_control_x, down);
actuator_control_y = downsample_matrix(actuator_control_y, down);
actuator_control_z = downsample_matrix(actuator_control_z, down);
pose_angle_x = downsample_matrix(pose_angle_x, down);
pose_angle_y = downsample_matrix(pose_angle_y, down);
pose_angle_z = downsample_matrix(pose_angle_z, down);
velocity_angle_x = downsample_matrix(velocity_angle_x, down);
velocity_angle_y = downsample_matrix(velocity_angle_y, down);
velocity_angle_z = downsample_matrix(velocity_angle_z, down);


figure('color', [1 1 1]);
set(gcf, 'Position', [200, 200, 600, 400]);

subplot(3,1,1)
hold on;
plot(time, actuator_control_x(:,2), 'r');
plot(time, actuator_control_y(:,2), 'g');
plot(time, actuator_control_z(:,2), 'b');
hold off;
set(gca,'fontsize',10,'fontname','Times');
legend('$u_1$', '$u_2$', '$u_3$',...
       'Fontname', 'Times New Roman','Interpreter', 'latex','FontSize', 11);
ylabel('input', 'Fontname', 'Times New Roman', 'FontSize', 11);
set(gca,'Position',[0.085 0.74 0.90 0.24]);
xlim([0, max(time)]);


subplot(3,1,2)
hold on;
plot(time, pose_angle_x(:,2), 'r');
plot(time, pose_angle_y(:,2), 'g');
plot(time, pose_angle_z(:,2), 'b');
hold off;
set(gca,'fontsize',10,'fontname','Times');
legend('$\phi$', '$\theta$', '$\psi$', ...
       'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('Euler angles [rad]', 'Fontname', 'Times New Roman', 'FontSize', 11);
set(gca,'Position',[0.085 0.42 0.90 0.24]);
xlim([0, max(time)]);


subplot(3,1,3)
hold on;
plot(time, velocity_angle_x(:,2), 'r');
plot(time, velocity_angle_y(:,2), 'g');
plot(time, velocity_angle_z(:,2), 'b');
hold off;
set(gca,'fontsize',10,'fontname','Times');
legend('$\omega_1$', '$\omega_2$', '$\omega_3$', ...
       'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 11);
xlabel('time [s]', 'Fontname', 'Times New Roman', 'FontSize', 11);
ylabel('angular vel. [rad/s]', 'Fontname', 'Times New Roman', 'FontSize', 11);
set(gca,'Position',[0.085 0.10 0.90 0.24]);
xlim([0, max(time)]);



























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





