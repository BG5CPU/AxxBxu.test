clear; clc; close all;

%% import data ============================================================

% file path A D F
file1 = 'Gazebo_data_csv/sim_vK_closed_1007roll/imu_data.csv';
file2 = 'Gazebo_data_csv/sim_vK_closed_1007roll/attitude_r_p_y.csv';
file3 = 'Gazebo_data_csv/sim_vK_closed_1007roll/actuator_control.csv';
file4 = 'Gazebo_data_csv/sim_vK_closed_1007roll/input_x_y_z.csv';
file5 = 'Gazebo_data_csv/sim_vK_closed_1007roll/pose.csv';


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

% figure(10); hold on;
% plot(matrix_attitude_r_p_y(:,4), 'r'); 
% plot(matrix_attitude_r_p_y(:,5), 'g');
% plot(matrix_attitude_r_p_y(:,6), 'b');
% hold off;



[~, idx_matrix_attitude_r_p_y_1] = unique(matrix_attitude_r_p_y(:,1), 'stable');
[~, idx_matrix_attitude_r_p_y_3] = unique(matrix_attitude_r_p_y(:,3), 'stable');

matrix_attitude_r_p_y = matrix_attitude_r_p_y(idx_matrix_attitude_r_p_y_3, :);



matrix_imu_data = matrix_imu_data(32:end,:);

figure(1); hold on;
plot(flipud(matrix_imu_data(:,3)), 'r.'); % flipud
plot(flipud(matrix_attitude_r_p_y(:,3)), 'g.');
hold off;

time = matrix_attitude_r_p_y(:,3);

% velocity_angle_x = interpolate_to_match_A(time, matrix_imu_data(:,[1,8]));
% velocity_angle_y = interpolate_to_match_A(time, matrix_imu_data(:,[1,9]));
% velocity_angle_z = interpolate_to_match_A(time, matrix_imu_data(:,[1,10]));
% 
% velocity_angle = [velocity_angle_x(:,2), velocity_angle_y(:,2), velocity_angle_z(:,2)];
% 
% pose_angle = matrix_attitude_r_p_y(:,[4,5,6]);
% 
% time = (time-time(1))/1e9;
% 
% downs = 20;
% 
% time = downsample_matrix(time, downs);
% pose_angle = downsample_matrix(pose_angle, downs);
% velocity_angle = downsample_matrix(velocity_angle, downs);






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