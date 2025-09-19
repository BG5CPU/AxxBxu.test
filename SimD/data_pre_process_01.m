clear; clc; close all;

%% import data ============================================================

% file path
file1 = 'Gazebo_data_csv/data_csv_10_04/velocity_body.csv';
file2 = 'Gazebo_data_csv/data_csv_10_04/velocity_local.csv';
file3 = 'Gazebo_data_csv/data_csv_10_04/actuator_control.csv';
file4 = 'Gazebo_data_csv/data_csv_10_04/input_x_y_z.csv';
file5 = 'Gazebo_data_csv/data_csv_10_04/pose.csv';

% read (velocity_body.csv)
opts = detectImportOptions(file1);
opts.SelectedVariableNames = [1, 2, 3, 5, 6, 7, 8, 9, 10]; 
opts.VariableNamingRule = 'preserve';
data_body = readtable(file1, opts);
matrix_velocity_body = table2array(data_body); 
check_data_continuity(matrix_velocity_body(:,2), 'Matrix_Velocity_Body');



% read (velocity_local.csv)
opts = detectImportOptions(file2);
opts.SelectedVariableNames = [1, 2, 3, 5, 6, 7, 8, 9, 10]; 
opts.VariableNamingRule = 'preserve';
data_local = readtable(file2, opts);
matrix_velocity_local = table2array(data_local); 
check_data_continuity(matrix_velocity_local(:,2), 'Matrix_Velocity_Local');


% plot(matrix_velocity_body(:,4), 'b-');
% hold on;
% plot(matrix_velocity_local(:,4), 'r:');


% read (actuator_control.csv)
opts = detectImportOptions(file3);
opts.SelectedVariableNames = [1, 2, 6, 7, 8, 9]; 
opts.VariableNamingRule = 'preserve';
actuator_control = readtable(file3, opts);
matrix_actuator_control = table2array(actuator_control); 
check_data_continuity(matrix_actuator_control(:,2), 'Matrix_Actuator_Control');



% read (input_x_y_z.csv)
opts = detectImportOptions(file4);
opts.SelectedVariableNames = [1, 2, 3, 5, 6, 7, 8, 9, 10, 11]; 
opts.VariableNamingRule = 'preserve';
input_x_y_z = readtable(file4, opts);
matrix_input_x_y_z = table2array(input_x_y_z); 
check_data_continuity(matrix_input_x_y_z(:,2), 'Matrix_Input_X_Y_Z');


dif_time_actuator_control_input_x_y_z = matrix_actuator_control(:,1) - matrix_input_x_y_z(:,1);
dif_x__actuator_control_input_x_y_z = matrix_actuator_control(:,3) - matrix_input_x_y_z(:,4);
dif_y__actuator_control_input_x_y_z = matrix_actuator_control(:,4) - matrix_input_x_y_z(:,5);
dif_z__actuator_control_input_x_y_z = matrix_actuator_control(:,5) - matrix_input_x_y_z(:,6);



% read (pose.csv)
opts = detectImportOptions(file5);
opts.SelectedVariableNames = [1, 2, 3, 5, 6, 7, 8, 9, 10, 11]; 
opts.VariableNamingRule = 'preserve';  
pose = readtable(file5, opts);
matrix_pose = table2array(pose); 
check_data_continuity(matrix_pose(:,2), 'Matrix_Pose');



% time align
figure(23);
align_pose_velocity = 12;
plot(matrix_pose(align_pose_velocity:end,3), '*');
hold on;
plot(matrix_velocity_body(:,3), 'o');

matrix_pose_A = matrix_pose(align_pose_velocity:end,:);



figure(24);
align_input_velocity = 19;
plot(matrix_actuator_control(1:end,1), '*');
hold on;
plot(matrix_velocity_body(align_input_velocity:end,1), 'o');

matrix_velocity_body_A = matrix_velocity_body(align_input_velocity:end, :);
matrix_pose_A = matrix_pose_A(align_input_velocity:end, :);

figure(25);
plot(diff(matrix_actuator_control(1:end,1)), '*');

matrix_actuator_control = remove_time_duplicates(matrix_actuator_control);

figure(26);
plot(diff(matrix_actuator_control(1:end,1)), '*');


% max(matrix_actuator_control(:,1)) - max(matrix_pose_A(:,1))
% min(matrix_actuator_control(:,1)) - min(matrix_pose_A(:,1))
pose_linear_x = interpolate_to_match_A(matrix_actuator_control(:,1), matrix_pose_A(:,[1,4]));
pose_linear_y = interpolate_to_match_A(matrix_actuator_control(:,1), matrix_pose_A(:,[1,5]));
pose_linear_z = interpolate_to_match_A(matrix_actuator_control(:,1), matrix_pose_A(:,[1,6]));
pose_angle_x = interpolate_to_match_A(matrix_actuator_control(:,1), matrix_pose_A(:,[1,7]));
pose_angle_y = interpolate_to_match_A(matrix_actuator_control(:,1), matrix_pose_A(:,[1,8]));
pose_angle_z = interpolate_to_match_A(matrix_actuator_control(:,1), matrix_pose_A(:,[1,9]));
pose_angle_w = interpolate_to_match_A(matrix_actuator_control(:,1), matrix_pose_A(:,[1,10]));


velocity_linear_x = interpolate_to_match_A(matrix_actuator_control(:,1), matrix_velocity_body_A(:,[1,4]));
velocity_linear_y = interpolate_to_match_A(matrix_actuator_control(:,1), matrix_velocity_body_A(:,[1,5]));
velocity_linear_z = interpolate_to_match_A(matrix_actuator_control(:,1), matrix_velocity_body_A(:,[1,6]));
velocity_angle_x = interpolate_to_match_A(matrix_actuator_control(:,1), matrix_velocity_body_A(:,[1,7]));
velocity_angle_y = interpolate_to_match_A(matrix_actuator_control(:,1), matrix_velocity_body_A(:,[1,8]));
velocity_angle_z = interpolate_to_match_A(matrix_actuator_control(:,1), matrix_velocity_body_A(:,[1,9]));


actuator_control_x = matrix_actuator_control(:,[1,3]);
actuator_control_y = matrix_actuator_control(:,[1,4]);
actuator_control_z = matrix_actuator_control(:,[1,5]);


data_input_output_01 = [actuator_control_x, actuator_control_y, actuator_control_z, ...
                        pose_linear_x, pose_linear_y, pose_linear_z, ...
                        pose_angle_x, pose_angle_y, pose_angle_z, pose_angle_w, ...
                        velocity_linear_x, velocity_linear_y, velocity_linear_z, ...
                        velocity_angle_x, velocity_angle_y, velocity_angle_z];

Pstart01 = 5;
Pend01 = 1466;
data_input_output_01 = data_input_output_01(Pstart01:Pend01 , :);
save('./Gazebo_data_mat/data01.mat', 'data_input_output_01');





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


function B = remove_time_duplicates(A)

    if isempty(A)
        B = A;
        return;
    end
    time = A(:,1);                 
    dupIdx = [diff(time)==0; false]; 
    B = A(~dupIdx,:);              
end



