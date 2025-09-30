clear; clc; close all;

%% import data ============================================================

% file path
file1 = 'Gazebo_data_csv/sim_mpc_sin_0930C/imu_data.csv';
file2 = 'Gazebo_data_csv/sim_mpc_sin_0930C/attitude_r_p_y.csv';
file3 = 'Gazebo_data_csv/sim_mpc_sin_0930C/actuator_control.csv';
file4 = 'Gazebo_data_csv/sim_mpc_sin_0930C/input_x_y_z.csv';
file5 = 'Gazebo_data_csv/sim_mpc_sin_0930C/pose.csv';


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

matrix_input_x_y_z(14531,:) = [];
% matrix_actuator_control = matrix_actuator_control(3:end,:);

figure(3); hold on;
plot(matrix_actuator_control(:,3), 'r.');
plot(matrix_input_x_y_z(:,3), 'g.');
hold off;

figure(4)
plot(matrix_actuator_control(:,4)-matrix_input_x_y_z(:,4),'r.');


matrix_imu_data = matrix_imu_data(46:end,:);
figure(5); hold on;
plot((matrix_actuator_control(:,1)), 'r.'); % flipud
plot((matrix_imu_data(:,1)), 'g.');
hold off;


matrix_attitude_r_p_y(14531,:) = [];
figure(6); hold on;
plot((matrix_input_x_y_z(:,3)), 'r.');
plot((matrix_attitude_r_p_y(:,3)), 'g.');
hold off;


% Tlen_matrix_input_x_y_z = -matrix_input_x_y_z(1,3)+matrix_input_x_y_z(end,3);
% Nlen_matrix_input_x_y_z = round(Tlen_matrix_input_x_y_z/2e6);
% time = matrix_input_x_y_z(1,3) + (0:Nlen_matrix_input_x_y_z-1)*2e6;
% time = time';
% figure(7); hold on;
% plot((time), 'r.');
% plot((matrix_input_x_y_z(:,3)), 'g.');
% hold off;

time = matrix_actuator_control(:,3);


% pose_angle_x = interpolate_to_match_A(time, matrix_attitude_r_p_y(:,[3,4]));
% pose_angle_y = interpolate_to_match_A(time, matrix_attitude_r_p_y(:,[3,5]));
% pose_angle_z = interpolate_to_match_A(time, matrix_attitude_r_p_y(:,[3,6]));

pose_angle_x = matrix_attitude_r_p_y(:,[3,4]);
pose_angle_y = matrix_attitude_r_p_y(:,[3,5]);
pose_angle_z = matrix_attitude_r_p_y(:,[3,6]);

velocity_angle_x = interpolate_to_match_A(time, matrix_imu_data(:,[1,8]));
velocity_angle_y = interpolate_to_match_A(time, matrix_imu_data(:,[1,9]));
velocity_angle_z = interpolate_to_match_A(time, matrix_imu_data(:,[1,10]));


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


% actuator_control_x = interpolate_to_match_A(time, matrix_input_x_y_z(:,[3,4]));
% actuator_control_y = interpolate_to_match_A(time, matrix_input_x_y_z(:,[3,5]));
% actuator_control_z = interpolate_to_match_A(time, matrix_input_x_y_z(:,[3,6]));
actuator_control_x = matrix_input_x_y_z(:,[3,4]);
actuator_control_y = matrix_input_x_y_z(:,[3,5]);
actuator_control_z = matrix_input_x_y_z(:,[3,6]);


figure(10); hold on;
plot(actuator_control_x(:,2), 'r');
plot(actuator_control_y(:,2), 'g');
plot(actuator_control_z(:,2), 'b');
hold off;



data_input_output_08B = [actuator_control_x, actuator_control_y, actuator_control_z, ...
                        pose_angle_x, pose_angle_y, pose_angle_z, ...
                        velocity_angle_x, velocity_angle_y, velocity_angle_z];


Pstart01 = 1000;
Pend01 = 16000;
% data_input_output_08B = data_input_output_08B(Pstart01:Pend01 , :);
% save('./Gazebo_data_mat/data08B.mat', 'data_input_output_08B');




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
    interpolated_values = interp1(B_time, B_values, A_timestamps, 'spline', 'extrap');
    
    interpolated_data = [A_timestamps, interpolated_values];
end