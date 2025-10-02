clear; clc; close all;

%% import data ============================================================

% file path
file1 = 'Gazebo_data_csv/sim_mpc_to_vK_1001B/imu_data.csv';
file2 = 'Gazebo_data_csv/sim_mpc_to_vK_1001B/attitude_r_p_y.csv';
file3 = 'Gazebo_data_csv/sim_mpc_to_vK_1001B/actuator_control.csv';
file4 = 'Gazebo_data_csv/sim_mpc_to_vK_1001B/input_x_y_z.csv';
file5 = 'Gazebo_data_csv/sim_mpc_to_vK_1001B/pose.csv';


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




%% plot data ==============================================================

time = matrix_attitude_r_p_y(:,1);
time = (time-time(1))/1e9;

attitude_r_p_y = matrix_attitude_r_p_y(:,[4,5,6]);
attitude_r_p_y = rad2deg(attitude_r_p_y);

down = 20;
time = downsample_matrix(time, down);
attitude_r_p_y = downsample_matrix(attitude_r_p_y, down);

figure(1); hold on;
plot(time, attitude_r_p_y(:,1), 'r');
plot(time, attitude_r_p_y(:,2), 'g');
plot(time, attitude_r_p_y(:,3), 'b');
hold off;
% ylim([-10 10]);

















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



function Y = downsample_matrix(X,d)
% DOWNSAMPLE_MATRIX  Downsample matrix X along rows
    [N,~] = size(X);
    N_valid = floor(N/d) * d;
    Y = X(1:d:N_valid, :);
end
