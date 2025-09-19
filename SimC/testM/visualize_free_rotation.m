function visualize_free_rotation()
    % Parameters setup
    I = diag([1.0, 2.0, 3.0]);      % Inertia tensor (principal axes)
    omega0 = [0.5; 1.0; 0.3];        % Initial angular velocity (not aligned with principal axes)
    tspan = [0 20];                  % Simulation time range
    
    % Euler equations function (no external torque)
    function domega_dt = euler_eq(t, omega, I)
        % Compute cross product term: ω × (Iω)
        cross_term = cross(omega, I*omega);
        % Solve for dω/dt: I*dω/dt = -ω × (Iω)
        domega_dt = I \ (-cross_term);
    end

    % Rotation matrix differential equation: dR/dt = R * [ω]×
    function dRdt = rotation_matrix_eq(t, R_flat, omega)
        % Reshape flattened R matrix (9x1 vector) back to 3x3
        R = reshape(R_flat, [3 3]);
        
        % Skew-symmetric matrix of angular velocity
        omega_skew = [0, -omega(3), omega(2);
                      omega(3), 0, -omega(1);
                      -omega(2), omega(1), 0];
        
        % Compute derivative
        dRdt = R * omega_skew;
        dRdt = dRdt(:);              % Flatten to column vector
    end

    % Step 1: Solve Euler equations for ω(t)
    [t, omega] = ode45(@(t,y) euler_eq(t,y,I), tspan, omega0);
    
    % Step 2: Solve for rotation matrix R(t)
    R0 = eye(3);                     % Initial rotation matrix (body frame aligned with inertial)
    R_flat0 = R0(:);                 % Flatten initial condition
    
    % Preallocate storage for R(t)
    R = zeros(3, 3, length(t));
    R(:,:,1) = R0;
    
    % Solve R(t) in segments since it depends on ω(t)
    for i = 2:length(t)
        [~, R_flat] = ode45(@(t,y) rotation_matrix_eq(t,y,omega(i-1,:)'), ...
                           [t(i-1) t(i)], R_flat0);
        R_flat0 = R_flat(end,:)';
        R(:,:,i) = reshape(R_flat0, [3 3]);
    end

    % Step 3: Visualization
    visualize_rotation(t, R);
end

function visualize_rotation(t, R)
    % Create a simple rigid body shape (cube vertices)
    [X,Y,Z] = meshgrid([-0.5 0.5], [-0.5 0.5], [-0.5 0.5]);
    X = X(:); Y = Y(:); Z = Z(:);
    cube = [X Y Z];  % 8 vertices of the cube
    
    % Setup figure
    figure;
    axis equal;
    axis([-1 1 -1 1 -1 1]*2);  % Set axis limits
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Free Rotation of Rigid Body');
    grid on;
    view(3);  % 3D view
    
    % Draw inertial frame axes
    hold on;
    quiver3(0,0,0,1,0,0,'r','LineWidth',2); % X-axis (red)
    quiver3(0,0,0,0,1,0,'g','LineWidth',2); % Y-axis (green)
    quiver3(0,0,0,0,0,1,'b','LineWidth',2); % Z-axis (blue)
    
    % Initialize cube plot
    h = plot3(cube(:,1), cube(:,2), cube(:,3), 'ko', ...
             'MarkerSize',10, 'MarkerFaceColor','k');
    
    % Animation loop
    for i = 1:5:length(t)
        % Transform cube vertices to current orientation
        rotated_cube = (R(:,:,i) * cube')';
        
        % Update cube vertices
        set(h, 'XData', rotated_cube(:,1), ...
               'YData', rotated_cube(:,2), ...
               'ZData', rotated_cube(:,3));
        
        % Draw body frame axes at specific times
        if ismember(i, [1 round(length(t)/2) length(t)])
            R_current = R(:,:,i);
            % Delete previous body frame arrows
            delete(findobj(gca,'Type','quiver','UserData','body_frame'));
            % X-axis of body frame (red)
            quiver3(0,0,0,R_current(1,1),R_current(2,1),R_current(3,1),'r',...
                   'LineWidth',1,'UserData','body_frame');
            % Y-axis of body frame (green)
            quiver3(0,0,0,R_current(1,2),R_current(2,2),R_current(3,2),'g',...
                   'LineWidth',1,'UserData','body_frame');
            % Z-axis of body frame (blue)
            quiver3(0,0,0,R_current(1,3),R_current(2,3),R_current(3,3),'b',...
                   'LineWidth',1,'UserData','body_frame');
        end
        
        drawnow;
        pause(0.01);  % Control animation speed
    end
end







%% visualize_free_rotation_coupled ========================================
function visualize_free_rotation_coupled()
    % Parameters
    I = diag([1.0, 2.0, 3.0]);      % Inertia tensor
    omega0 = [0.5; 1.0; 0.3];        % Initial angular velocity
    R0 = eye(3);                     % Initial rotation matrix
    tspan = [0 20];                  % Simulation time
    
    % Combined state vector: [ω; R_flat] where R_flat is R(:)
    y0 = [omega0; R0(:)];
    
    % Coupled ODE function
    function dy = coupled_ode(t, y, I)
        % Extract ω and R from state vector
        omega = y(1:3);
        R = reshape(y(4:end), 3, 3);
        
        % Euler equation for dω/dt
        cross_term = cross(omega, I*omega);
        domega_dt = I \ (-cross_term);
        
        % Rotation matrix equation dR/dt = R*[ω]×
        omega_skew = [0, -omega(3), omega(2);
                     omega(3), 0, -omega(1);
                     -omega(2), omega(1), 0];
        dRdt = R * omega_skew;
        
        % Combine derivatives
        dy = [domega_dt; dRdt(:)];
    end

    % Solve coupled system
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    [t, y] = ode45(@(t,y) coupled_ode(t,y,I), tspan, y0, options);
    
    % Extract results
    omega = y(:,1:3)';
    R = reshape(y(:,4:end)', 3, 3, []);
    
    % Visualization
    visualize_rotation(t, R);
end