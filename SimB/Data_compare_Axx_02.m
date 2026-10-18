clear; close all; clc;


%% my_method test experiment =============================================

% My method: Theorem 4 in paper_2510.16451v1.pdf.
%
% The continuous-time parent system is identical to the SOS baseline:
%   xdot = (delta*I + S_n)*x
%          -(alpha/sqrt(n))*B_en*sum_i(x_i^3)+B_en*u.
%
% Theorem 4 is applied to the exact forward-Euler benchmark
%   x+ = A_d(x)*x+B_d*u+d,
% where
%   A_d(x) = I+h*(delta*I+S_n)
%            -(h*alpha/sqrt(n))*B_en*[x_1^2,...,x_n^2],
%   B_d = h*B_en.

dim_x_list = 1:13;
num_Exp = 30;

system_setting.delta = 0.1;
system_setting.alpha = 1;

% Euler sampling period h and the holding time of each input sample.
data_setting.sample_period = 0.01;
% Lower bound on the total number T of concatenated samples.
data_setting.minimum_sample_num = 100;
% Requested number of independently initialized trajectories.
data_setting.trajectory_num = 50;
% Number of consecutive samples collected within each trajectory.
data_setting.in_trajectory_samp_num = 10;
% The required total length is
% max(minimum_sample_num,sample_num_factor*(2*dim_x+1)). If the requested
% trajectory_num is insufficient, it is increased to the smallest integer
% whose product with in_trajectory_samp_num reaches that requirement.
data_setting.sample_num_factor = 10;
% Excitation input is independently sampled from [-input_bound,input_bound].
data_setting.input_bound = 10;
% Initial state is sampled uniformly by volume from the ball of this radius.
data_setting.initial_radius = 0.5;
% Each entry of the unscaled base noise is sampled from this symmetric bound.
data_setting.base_noise_bound = 1;
% Scales base noise; discrete noise is h*noise_gain*base_noise.
data_setting.noise_gain = 0;
% Defines Theta=margin*(D0*D0') for the noise-energy upper bound.
data_setting.noise_bound_margin = 1.00;
% Numerical tolerance for the known-model data-equation residual.
data_setting.data_equation_tolerance = 1e-10;
% Tolerance for lambda_max(D0*D0'-Theta) <= 0.
data_setting.noise_bound_tolerance = 1e-10;
% Reproducible master seed used to derive input, initial-state and noise seeds.
data_setting.base_seed = 20260828;

% These settings match Data_compare_SOS_02.m. Both experiments use the
% same YALMIP modeling layer and the same MOSEK solver.
% SDP solver selected through YALMIP.
lmi_setting.solver = 'mosek';
% YALMIP/MOSEK console verbosity: 0 suppresses routine solver output.
lmi_setting.verbose = 0;
% Enforces Gamma_bar >= matrix_margin*I in Theorem 4.
lmi_setting.matrix_margin = 1e-9;
% Enforces epsilon_Gamma_bar >= epsilon_margin.
lmi_setting.epsilon_margin = 1e-9;
% MOSEK primal and dual feasibility tolerances.
lmi_setting.feasibility_tolerance = 1e-4;
% MOSEK relative primal-dual objective-gap tolerance.
lmi_setting.relative_gap_tolerance = 1e-4;
% Post-solution tolerance for margins and recovered vertex LMIs.
lmi_setting.validation_tolerance = 1e-9;
% Configured stability requires rho(A_d0+B_d*K) <= 1-stability_margin.
lmi_setting.stability_margin = 1e-9;
% Maximum solver time in seconds for each trial.
lmi_setting.solve_timeout = 500;
% Run every LMI solve in an isolated MATLAB process and monitor total system
% memory, so an oversized solve cannot terminate the main MATLAB process.
% This protection requires Parallel Computing Toolbox on Linux.
lmi_setting.memory_guard_enabled = true;
% Stop a trial when 1-MemAvailable/MemTotal reaches this fraction.
lmi_setting.maximum_memory_usage_ratio = 0.90;
% Also stop when Linux MemAvailable falls below this many GiB.
lmi_setting.minimum_available_memory_gb = 1.5;
% Interval in seconds between memory checks while a worker is running.
lmi_setting.memory_poll_period = 0.25;
% Seconds allowed after SIGTERM before the isolated worker receives SIGKILL.
lmi_setting.memory_termination_grace_period = 2.0;
% An owned worker pool closes itself after this many idle minutes. The script
% does not synchronously delete the pool at the end, so it returns promptly.
lmi_setting.memory_pool_idle_timeout_minutes = 1;
% Radius r of B_r; it also sets every vertex bound x_i^2 in [0,r^2].
lmi_setting.domain_radius = 0.3;
% Prints !/? marks and all recovered-controller validation metrics.
lmi_setting.print_validation = true;

output_directory = fullfile(fileparts(mfilename('fullpath')),'output');
output_mat = fullfile(output_directory,'my_method_test_results.mat');
output_run_csv = fullfile(output_directory,'my_method_test_runs.csv');
output_summary_csv = ...
    fullfile(output_directory,'my_method_test_summary.csv');

% Output contents:
%   MAT file    - all run/summary structs, controllers, data and settings.
%   runs CSV    - one row per (dim_x, experiment_id) trial.
%   summary CSV - one aggregated row per dim_x (linear controller order 0).
% The terminal first prints the same common-summary columns as the SOS code:
%   trajectory_num              - effective number of independently
%                                 initialized trajectories after any automatic
%                                 increase needed to meet the data-length rule;
%   in_trajectory_samp_num      - number of samples in each trajectory;
%   sample_num                  - total concatenated data length, equal to the
%                                 product of the preceding two quantities;
%   required_sample_num         - lower bound imposed by minimum_sample_num
%                                 and sample_num_factor*(2*dim_x+1);
%   attempted_num               - total trials in this (dim_x,order_K) group;
%   rank_valid_num              - trials whose data regressor has full row rank;
%   solver_success_num          - trials reported solved by YALMIP/MOSEK;
%   configured/strict_feasible_num and *_rate
%                               - controller-validation counts and fractions;
%   timeout_num                 - trials that reached the solver time limit;
%   solve_time                  - YALMIP solver-call wall time (seconds) of
%                                 configured-feasible trials;
%   synthesis_time              - model build + solve + controller-recovery
%                                 wall time of configured-feasible trials;
%   attempt_*_time              - corresponding time over every solver attempt,
%                                 including failed and timed-out attempts;
%   data_condition              - cond([X0_lifted;U0]); smaller is better.
% Every time and condition quantity is summarized by both mean and median.
% The second terminal table contains these my-method-only diagnostics:
%   dim_x                       - state dimension of the tested system;
%   data_valid_num              - trials satisfying both the data-equation
%                                 residual and noise-bound eigenvalue tests;
%                                 the rank test is counted separately;
%   solver_candidate_num        - non-timeout trials for which YALMIP/MOSEK
%                                 returned finite variables and a recoverable K;
%   soft_feasible_num/rate      - candidate count/fraction additionally having
%                                 rho(A_d(0)+B_d*K)<1; certificate residuals
%                                 and standard solver-success status are not
%                                 required by this soft test;
%   memory_limit_num            - isolated solves stopped by the memory guard;
%   soft_solve_time             - solver-call time of soft-feasible trials;
%   soft_synthesis_time         - synthesis time of soft-feasible trials;
%   recovered_roa_radius        - Theorem-4 guaranteed Euclidean-ball radius,
%                                 capped by domain_radius and summarized only
%                                 over configured-feasible trials;
%   soft_spectral_radius        - rho(A_d(0)+B_d*K), summarized only over
%                                 soft-feasible trials; values below 1 mean
%                                 local discrete-time stability at the origin;
%   vertex_num                  - 2^dim_x vertices of x_i^2 in [0,r^2];
%   lmi_dimension               - row/column size 4*dim_x+1 of each vertex LMI;
%   lmi_constraint_num          - number 2^dim_x of vertex LMIs (not including
%                                 the separate Gamma_bar/epsilon constraints);
%   scalar_decision_variable_num - independent scalar entries in symmetric
%                                 Gamma_bar, scalar epsilon_Gamma_bar and Y_bar;
%                                 for dim_u=1 it is (dim_x+1)*(dim_x+2)/2.
% Every soft time, ROA and spectral-radius quantity has mean/median columns.

fprintf('My-method Theorem 4 test experiment\n');
fprintf('  dim_x          = %s\n',mat2str(dim_x_list));
fprintf('  num_Exp        = %d\n',num_Exp);
fprintf('  noise_gain     = %.3e\n',data_setting.noise_gain);
fprintf('  sample time    = %.3e s\n',data_setting.sample_period);
fprintf('  domain radius  = %.3e\n',lmi_setting.domain_radius);
fprintf('  solver          = YALMIP + %s\n',upper(lmi_setting.solver));
fprintf('  solve timeout   = %.1f s\n',lmi_setting.solve_timeout);
fprintf('  memory guard    = %.1f%% used or %.2f GiB available\n\n', ...
        100*lmi_setting.maximum_memory_usage_ratio, ...
        lmi_setting.minimum_available_memory_gb);

funcSetupYalmip();

memory_pool = [];
memory_pool_owned = false;

num_dimension = numel(dim_x_list);
num_run = num_dimension*num_Exp;
run_result = repmat(funcEmptyRunResult(),num_run,1);
controller_result = cell(num_run,1);
soft_controller_result = cell(num_run,1);
candidate_controller_result = cell(num_run,1);
experiment_data = cell(num_dimension,num_Exp);
system_result = cell(num_dimension,1);
run_index = 0;

for idim = 1:num_dimension
    dim_x = dim_x_list(idim);
    system = funcBuildSystem(dim_x,system_setting,data_setting);
    system_result{idim} = system;

    for iexp = 1:num_Exp
        experiment_data{idim,iexp} = ...
            funcGenerateExperimentData(system,data_setting,iexp);
    end

    fprintf(['\n============================================================', ...
             '================\n']);
    fprintf('dim_x = %d, vertices = %d, LMI dimension = %d\n', ...
            dim_x,2^dim_x,4*dim_x+1);

    for iexp = 1:num_Exp
        run_index = run_index+1;
        data = experiment_data{idim,iexp};
        result = funcInitializeRunResult(system,data,lmi_setting,iexp);

        fprintf('\nRun %d/%d: experiment_id = %d\n', ...
                run_index,num_run,iexp);
        fprintf(['  trajectories = %d (requested %d), samples/trajectory', ...
                 ' = %d, total = %d, required >= %d\n'], ...
                data.trajectory_num,data.requested_trajectory_num, ...
                data.in_trajectory_samp_num,data.sample_num, ...
                data.required_sample_num);
        fprintf('  rank([X0_lifted; U0]) = %d/%d, condition = %.3e\n', ...
                data.data_rank,data.expected_rank,data.data_condition);
        fprintf('  %s data equation residual = %.3e (tol %.1e)\n', ...
                funcValidationMark(data.data_equation_passed), ...
                data.data_equation_residual, ...
                data_setting.data_equation_tolerance);
        fprintf('  %s max noise-bound eig = %.3e (tol %.1e)\n', ...
                funcValidationMark(data.noise_bound_passed), ...
                data.maximum_noise_bound_eigenvalue, ...
                data_setting.noise_bound_tolerance);

        if ~data.rank_passed
            result.status = 'rank deficient';
            fprintf('  ? skipped: the data matrix is not full row rank.\n');
        elseif ~data.data_valid
            result.status = 'invalid data';
            fprintf('  ? skipped: data equation or noise bound failed.\n');
        else
            [controller,solver,validation,timing, ...
             memory_pool,memory_pool_owned] = ...
                funcRunTheorem4WithMemoryGuard( ...
                    system,data,lmi_setting, ...
                    memory_pool,memory_pool_owned);

            if solver.has_candidate
                candidate_controller_result{run_index} = controller;
            end
            if validation.soft_feasible
                soft_controller_result{run_index} = controller;
            end
            if validation.configured_feasible
                controller_result{run_index} = controller;
            end

            result = funcStoreSolveResult( ...
                result,solver,validation,timing);
        end

        run_result(run_index) = result;
    end
end

summary_result = funcSummarizeResults(run_result,dim_x_list,num_Exp);
run_table = struct2table(run_result);
summary_table = struct2table(summary_result);

common_summary_column = { ...
    'dim_x','order_K','trajectory_num','in_trajectory_samp_num', ...
    'sample_num','required_sample_num','attempted_num','rank_valid_num', ...
    'solver_success_num','configured_feasible_num', ...
    'strict_feasible_num','timeout_num','configured_solvability_rate', ...
    'strict_solvability_rate','mean_solve_time','median_solve_time', ...
    'mean_synthesis_time','median_synthesis_time', ...
    'mean_attempt_solve_time','median_attempt_solve_time', ...
    'mean_attempt_synthesis_time','median_attempt_synthesis_time', ...
    'mean_data_condition','median_data_condition'};
method_summary_column = { ...
    'dim_x','data_valid_num','memory_limit_num', ...
    'solver_candidate_num', ...
    'soft_feasible_num','soft_solvability_rate', ...
    'mean_soft_solve_time','median_soft_solve_time', ...
    'mean_soft_synthesis_time','median_soft_synthesis_time', ...
    'mean_recovered_roa_radius','median_recovered_roa_radius', ...
    'mean_soft_spectral_radius','median_soft_spectral_radius', ...
    'vertex_num','lmi_dimension','lmi_constraint_num', ...
    'scalar_decision_variable_num'};

fprintf('\n\nCommon comparison summary: my_method\n');
disp(summary_table(:,common_summary_column));
fprintf('\nMy-method-specific diagnostics\n');
disp(summary_table(:,method_summary_column));

if ~isfolder(output_directory)
    mkdir(output_directory);
end
save(output_mat,'run_result','summary_result','run_table','summary_table', ...
     'controller_result','soft_controller_result', ...
     'candidate_controller_result','experiment_data','system_result', ...
     'dim_x_list','num_Exp','system_setting','data_setting','lmi_setting');
writetable(run_table,output_run_csv);
writetable(summary_table,output_summary_csv);

fprintf('Results saved to:\n  %s\n  %s\n  %s\n', ...
        output_mat,output_run_csv,output_summary_csv);

if memory_pool_owned && ~isempty(memory_pool)
    fprintf(['Memory-guard worker pool will close automatically after ', ...
             '%.1f idle minute(s).\n'], ...
            lmi_setting.memory_pool_idle_timeout_minutes);
end








%% useful functions =======================================================

function system = funcBuildSystem(dim_x,system_setting,data_setting)
    S = zeros(dim_x);
    if dim_x > 1
        S(1:dim_x-1,2:dim_x) = eye(dim_x-1);
    end

    B_en = zeros(dim_x,1);
    B_en(end) = 1;

    h = data_setting.sample_period;
    cubic_coefficient = ones(dim_x,1)/sqrt(dim_x);
    A_continuous = system_setting.delta*eye(dim_x)+S;
    A_d0 = eye(dim_x)+h*A_continuous;
    B_d = h*B_en;

    E_A = zeros(dim_x,2*dim_x);
    for ix = 1:dim_x
        E_A(:,2*ix-1) = A_d0(:,ix);
        E_A(:,2*ix) = ...
            -h*system_setting.alpha*cubic_coefficient(ix)*B_en;
    end

    system.dim_x = dim_x;
    system.dim_u = 1;
    system.delta = system_setting.delta;
    system.alpha = system_setting.alpha;
    system.sample_period = h;
    system.S = S;
    system.B_en = B_en;
    system.cubic_coefficient = cubic_coefficient;
    system.A_continuous = A_continuous;
    system.A_d0 = A_d0;
    system.B_d = B_d;
    system.E_A = E_A;
    system.E_B = B_d;
    system.lifted_state_dimension = 2*dim_x;
    system.regressor_dimension = 2*dim_x+1;
end


function data = funcGenerateExperimentData(system,setting,experiment_id)
    dim_x = system.dim_x;
    [trajectory_num,in_trajectory_samp_num,sample_num, ...
     required_sample_num] = ...
        funcResolveTrajectoryLayout(system,setting);

    input_seed = setting.base_seed+experiment_id;
    initial_seed = setting.base_seed+1000*dim_x+experiment_id;
    noise_seed = setting.base_seed+100000+1000*dim_x+experiment_id;

    rng(input_seed,'twister');
    U0 = setting.input_bound*(2*rand(1,sample_num)-1);

    rng(initial_seed,'twister');
    x_initial = zeros(dim_x,trajectory_num);
    for itrajectory = 1:trajectory_num
        x_initial(:,itrajectory) = ...
            funcRandomInitialState(dim_x,setting.initial_radius);
    end

    rng(noise_seed,'twister');
    base_noise = setting.base_noise_bound* ...
                 (2*rand(dim_x,sample_num)-1);
    D0_continuous = setting.noise_gain*base_noise;
    D0 = setting.sample_period*D0_continuous;

    % Each block is a separate trajectory. Only valid within-trajectory
    % sample pairs are retained; no artificial transition is introduced
    % between the last state of one block and the first state of the next.
    X_state = zeros(dim_x,sample_num);
    X1 = zeros(dim_x,sample_num);
    for itrajectory = 1:trajectory_num
        column = (itrajectory-1)*in_trajectory_samp_num+ ...
                 (1:in_trajectory_samp_num);
        [X_state(:,column),X1(:,column)] = ...
            funcCollectEulerData( ...
                system,x_initial(:,itrajectory),U0(:,column),D0(:,column));
    end
    X0_lifted = funcLiftedStateData(X_state);
    U0_lifted = U0;
    M0 = [X0_lifted;U0_lifted];
    Theta = setting.noise_bound_margin*(D0*D0');

    singular_value = svd(M0);
    if isempty(singular_value) || singular_value(end) <= 0
        data_condition = Inf;
    else
        data_condition = singular_value(1)/singular_value(end);
    end

    % Concatenated-data check in the lifted Theorem-4 notation:
    % X1 = E_A*X0_lifted+E_B*U0_lifted+D0.
    data_equation_residual = norm( ...
        X1-system.E_A*X0_lifted-system.E_B*U0_lifted-D0,'fro');
    noise_bound_residual = D0*D0'-Theta;
    noise_bound_residual = 0.5*(noise_bound_residual+ ...
                                      noise_bound_residual');
    maximum_noise_bound_eigenvalue = ...
        max(real(eig(noise_bound_residual)));

    data.dim_x = dim_x;
    data.requested_trajectory_num = setting.trajectory_num;
    data.trajectory_num = trajectory_num;
    data.in_trajectory_samp_num = in_trajectory_samp_num;
    data.required_sample_num = required_sample_num;
    data.sample_num = sample_num;
    data.sample_period = setting.sample_period;
    data.input_seed = input_seed;
    data.initial_seed = initial_seed;
    data.noise_seed = noise_seed;
    data.x_initial = x_initial;
    data.trajectory_id = repelem(1:trajectory_num, ...
                                  in_trajectory_samp_num);
    data.U0 = U0;
    data.X_state = X_state;
    data.X0_lifted = X0_lifted;
    data.U0_lifted = U0_lifted;
    data.X1 = X1;
    data.base_noise = base_noise;
    data.D0_continuous = D0_continuous;
    data.D0 = D0;
    data.Theta = Theta;
    data.M0 = M0;
    data.data_A = M0*M0';
    data.data_B = -X1*M0';
    data.data_C = X1*X1'-Theta;
    data.data_rank = rank(M0);
    data.expected_rank = system.regressor_dimension;
    if size(M0,1) ~= data.expected_rank
        error(['The lifted data regressor has %d rows, but Theorem 4 ', ...
               'requires 2*dim_x+1 = %d rows.'], ...
              size(M0,1),data.expected_rank);
    end
    data.rank_passed = data.data_rank == data.expected_rank;
    data.data_condition = data_condition;
    data.data_equation_residual = data_equation_residual;
    data.maximum_noise_bound_eigenvalue = ...
        maximum_noise_bound_eigenvalue;
    data.data_equation_passed = ...
        data_equation_residual <= setting.data_equation_tolerance;
    data.noise_bound_passed = ...
        maximum_noise_bound_eigenvalue <= setting.noise_bound_tolerance;
    data.data_valid = ...
        data.data_equation_passed && data.noise_bound_passed;
end


function [trajectory_num,in_trajectory_samp_num,sample_num, ...
          required_sample_num] = ...
    funcResolveTrajectoryLayout(system,setting)

    requested_trajectory_num = setting.trajectory_num;
    in_trajectory_samp_num = setting.in_trajectory_samp_num;
    if ~isscalar(requested_trajectory_num) || ...
            ~isfinite(requested_trajectory_num) || ...
            requested_trajectory_num < 1 || ...
            requested_trajectory_num ~= fix(requested_trajectory_num)
        error('data_setting.trajectory_num must be a positive integer.');
    end
    if ~isscalar(in_trajectory_samp_num) || ...
            ~isfinite(in_trajectory_samp_num) || ...
            in_trajectory_samp_num < 1 || ...
            in_trajectory_samp_num ~= fix(in_trajectory_samp_num)
        error(['data_setting.in_trajectory_samp_num must be a ', ...
               'positive integer.']);
    end
    if ~isscalar(setting.minimum_sample_num) || ...
            ~isfinite(setting.minimum_sample_num) || ...
            setting.minimum_sample_num < 0
        error('data_setting.minimum_sample_num must be nonnegative.');
    end
    if ~isscalar(setting.sample_num_factor) || ...
            ~isfinite(setting.sample_num_factor) || ...
            setting.sample_num_factor < 0
        error('data_setting.sample_num_factor must be nonnegative.');
    end

    required_sample_num = ceil(max( ...
        setting.minimum_sample_num, ...
        setting.sample_num_factor*system.regressor_dimension));
    minimum_trajectory_num = ...
        ceil(required_sample_num/in_trajectory_samp_num);
    trajectory_num = max(requested_trajectory_num,minimum_trajectory_num);
    sample_num = trajectory_num*in_trajectory_samp_num;
end


function x_initial = funcRandomInitialState(dim_x,initial_radius)
    direction = randn(dim_x,1);
    if norm(direction) == 0
        direction(1) = 1;
    end
    radius = initial_radius*rand()^(1/dim_x);
    x_initial = radius*direction/norm(direction);
end


function [X_state,X1] = ...
    funcCollectEulerData(system,x_initial,U0,D0)

    sample_num = numel(U0);
    dim_x = system.dim_x;
    X_state = zeros(dim_x,sample_num);
    X1 = zeros(dim_x,sample_num);
    state = x_initial;

    for it = 1:sample_num
        X_state(:,it) = state;
        state = funcDiscreteSystem(state,U0(it),D0(:,it),system);
        X1(:,it) = state;
    end
end


function state_next = funcDiscreteSystem(x,u,d,system)
    dense_cubic = sum(system.cubic_coefficient.*(x(:).^3));
    state_next = system.A_d0*x ...
                 -system.sample_period*system.alpha* ...
                    system.B_en*dense_cubic ...
                 +system.B_d*u+d;
end


function X0_lifted = funcLiftedStateData(X_state)
    dim_x = size(X_state,1);
    sample_num = size(X_state,2);
    X0_lifted = zeros(2*dim_x,sample_num);
    X0_lifted(1:2:end,:) = X_state;
    X0_lifted(2:2:end,:) = X_state.^3;
end


function [controller,solver,validation,timing, ...
          memory_pool,memory_pool_owned] = ...
    funcRunTheorem4WithMemoryGuard( ...
        system,data,setting,memory_pool,memory_pool_owned)

    if ~setting.memory_guard_enabled
        [controller,solver,validation,timing] = ...
            funcRunTheorem4(system,data,setting);
        return;
    end

    [memory_usage_ratio,available_memory_gb] = ...
        funcReadLinuxMemoryUsage();
    if funcMemoryLimitExceeded( ...
            memory_usage_ratio,available_memory_gb,setting)
        reason = funcMemoryLimitMessage( ...
            memory_usage_ratio,available_memory_gb,setting);
        if memory_pool_owned
            funcDeleteMemoryWorkerPool(memory_pool);
            memory_pool = [];
        end
        [controller,solver,validation,timing] = ...
            funcFailedIsolatedSolve(-997,reason,NaN,true, ...
                                    memory_usage_ratio,available_memory_gb);
        fprintf('  ? memory guard: solve skipped; %s\n',reason);
        return;
    end

    [memory_pool,memory_pool_owned] = ...
        funcEnsureMemoryWorkerPool( ...
            memory_pool,memory_pool_owned,setting);
    pid_queue = parallel.pool.PollableDataQueue;
    future = parfeval(memory_pool,@funcRunTheorem4Worker,5, ...
                      system,data,setting,pid_queue);
    guard_timer = tic;
    worker_pid = NaN;

    while ~strcmpi(future.State,'finished')
        received_pid = poll(pid_queue,0);
        if ~isempty(received_pid)
            worker_pid = double(received_pid);
        end

        [memory_usage_ratio,available_memory_gb] = ...
            funcReadLinuxMemoryUsage();
        if strcmpi(future.State,'finished')
            break;
        end
        if funcMemoryLimitExceeded( ...
                memory_usage_ratio,available_memory_gb,setting)
            elapsed_time = toc(guard_timer);
            reason = funcMemoryLimitMessage( ...
                memory_usage_ratio,available_memory_gb,setting);
            try
                cancel(future);
            catch
            end
            if isfinite(worker_pid)
                funcTerminateWorkerProcess( ...
                    worker_pid,setting.memory_termination_grace_period);
            end
            funcDeleteMemoryWorkerPool(memory_pool);
            memory_pool = [];
            % A replacement pool created on the next trial belongs to this
            % script even if the terminated pool existed before the run.
            memory_pool_owned = true;
            [controller,solver,validation,timing] = ...
                funcFailedIsolatedSolve( ...
                    -997,reason,elapsed_time,true, ...
                    memory_usage_ratio,available_memory_gb);
            fprintf('  ? memory guard: worker terminated; %s\n',reason);
            return;
        end

        pause(setting.memory_poll_period);
    end

    try
        [controller,solver,validation,timing,worker_output] = ...
            fetchOutputs(future);
        fprintf('%s',worker_output);
    catch worker_error
        elapsed_time = toc(guard_timer);
        reason = sprintf('isolated LMI worker failed: %s', ...
                         worker_error.message);
        [controller,solver,validation,timing] = ...
            funcFailedIsolatedSolve(-996,reason,elapsed_time,false, ...
                                    memory_usage_ratio,available_memory_gb);
        fprintf('  ? %s\n',reason);
    end
end


function [controller,solver,validation,timing,worker_output] = ...
    funcRunTheorem4Worker(system,data,setting,pid_queue)

    send(pid_queue,feature('getpid'));
    funcSetupYalmip();
    worker_argument = {system,data,setting}; %#ok<NASGU>
    controller = [];
    solver = funcEmptySolverResult();
    validation = funcEmptyValidation();
    timing = funcEmptyTimingResult();
    worker_output = evalc([ ...
        '[controller,solver,validation,timing] = ', ...
        'funcRunTheorem4(worker_argument{1},', ...
        'worker_argument{2},worker_argument{3});']);
end


function [controller,solver,validation,timing] = ...
    funcRunTheorem4(system,data,setting)

    controller = [];
    solver = funcEmptySolverResult();
    validation = funcEmptyValidation();
    timing = funcEmptyTimingResult();
    synthesis_timer = tic;

    yalmip('clear');

    vertex_timer = tic;
    Q_vertex = funcBuildVertices(system.dim_x,setting.domain_radius);
    timing.vertex_generation_time = toc(vertex_timer);

    model_timer = tic;
    dim_x = system.dim_x;
    Gamma_bar = sdpvar(dim_x,dim_x,'symmetric');
    epsilon_Gamma_bar = sdpvar(1,1);
    Y_bar = sdpvar(system.dim_u,dim_x,'full');
    Gamma_Y = [Gamma_bar;Y_bar];

    constraint = [ ...
        Gamma_bar >= setting.matrix_margin*eye(dim_x), ...
        epsilon_Gamma_bar >= setting.epsilon_margin];

    for iv = 1:numel(Q_vertex)
        Qv = Q_vertex{iv};
        vertex_lmi = [ ...
            -Gamma_bar-data.data_C, zeros(dim_x), data.data_B; ...
            zeros(dim_x), ...
                -Gamma_bar+epsilon_Gamma_bar*eye(dim_x), ...
                -Gamma_Y'*Qv'; ...
            data.data_B',-Qv*Gamma_Y,-data.data_A];
        vertex_lmi = 0.5*(vertex_lmi+vertex_lmi');
        constraint = [constraint,vertex_lmi <= 0]; %#ok<AGROW>
    end

    timing.model_build_time = toc(model_timer);

    solve_timer = tic;
    try
        diagnostic = optimize( ...
            constraint,[],funcLmiOptions(setting));
        timing.attempt_solve_time = toc(solve_timer);
        solver.problem = diagnostic.problem;
        solver.info = diagnostic.info;
        solver.timed_out = funcSolverTimedOut( ...
            diagnostic.info,timing.attempt_solve_time, ...
            setting.solve_timeout);
        solver.solved = diagnostic.problem == 0 && ~solver.timed_out;
        solver.native_solve_time = diagnostic.solvertime;

        raw_controller.Gamma_bar = value(Gamma_bar);
        raw_controller.epsilon_Gamma_bar = value(epsilon_Gamma_bar);
        raw_controller.Y_bar = value(Y_bar);
        raw_value = [raw_controller.Gamma_bar(:); ...
                     raw_controller.epsilon_Gamma_bar(:); ...
                     raw_controller.Y_bar(:)];
        finite_raw_value = ...
            ~isempty(raw_value) && all(isfinite(raw_value));
        candidate_status = funcSolverReturnedCandidate( ...
            diagnostic.problem,diagnostic.info);

        if ~solver.timed_out && candidate_status && finite_raw_value
            try
                raw_controller.K = ...
                    raw_controller.Y_bar/raw_controller.Gamma_bar;
                finite_gain = all(isfinite(raw_controller.K(:)));
            catch
                finite_gain = false;
            end

            if finite_gain
                solver.has_candidate = true;
                controller.name = 'Theorem 4';
                controller.Gamma_bar = raw_controller.Gamma_bar;
                controller.epsilon_Gamma_bar = ...
                    raw_controller.epsilon_Gamma_bar;
                controller.Y_bar = raw_controller.Y_bar;
                controller.K = raw_controller.K;
                controller.Q_vertex = Q_vertex;
            end
        end

        if solver.timed_out
            solver.info = sprintf( ...
                'Timeout after %.3f s (limit %.3f s). %s', ...
                timing.attempt_solve_time,setting.solve_timeout, ...
                diagnostic.info);
        end
    catch solve_error
        timing.attempt_solve_time = toc(solve_timer);
        solver.problem = -999;
        solver.info = solve_error.message;
        solver.timed_out = ...
            timing.attempt_solve_time >= setting.solve_timeout;
        solver.solved = false;
        solver.has_candidate = false;
    end

    if ~solver.has_candidate
        timing.attempt_synthesis_time = toc(synthesis_timer);
        if solver.timed_out
            fprintf('  ? timeout: no controller was accepted.\n');
        else
            fprintf('  ? no recoverable solver candidate: %s\n',solver.info);
        end
        return;
    end

    % End synthesis timing before post-solution validation, matching the SOS
    % baseline's timing boundary. Validation has its own validation_time.
    timing.attempt_synthesis_time = toc(synthesis_timer);
    validation_timer = tic;
    validation = funcValidateController( ...
        controller,solver,system,data,setting);
    timing.validation_time = toc(validation_timer);
    controller.validation = validation;
end


function Q_vertex = funcBuildVertices(dim_x,domain_radius)
    vertex_num = 2^dim_x;
    Q_vertex = cell(vertex_num,1);

    for iv = 0:vertex_num-1
        active = bitget(iv,1:dim_x);
        theta = domain_radius^2*double(active(:));
        Xi_A = zeros(2*dim_x,dim_x);
        for ix = 1:dim_x
            Xi_A(2*ix-1,ix) = 1;
            Xi_A(2*ix,ix) = theta(ix);
        end
        Q_vertex{iv+1} = blkdiag(Xi_A,1);
    end
end


function validation = funcValidateController( ...
    controller,solver,system,data,setting)

    validation = funcEmptyValidation();
    tolerance = setting.validation_tolerance;

    Gamma_symmetric = 0.5*(controller.Gamma_bar+ ...
                                  controller.Gamma_bar');
    validation.Gamma_margin = min(eig( ...
        Gamma_symmetric-setting.matrix_margin*eye(system.dim_x)));
    validation.epsilon_margin_value = ...
        controller.epsilon_Gamma_bar-setting.epsilon_margin;
    validation.Gamma_condition = cond(controller.Gamma_bar);

    validation.maximum_vertex_lmi_eigenvalue = -Inf;
    Gamma_Y = [controller.Gamma_bar;controller.Y_bar];
    for iv = 1:numel(controller.Q_vertex)
        Qv = controller.Q_vertex{iv};
        vertex_lmi = [ ...
            -controller.Gamma_bar-data.data_C, ...
                zeros(system.dim_x),data.data_B; ...
            zeros(system.dim_x), ...
                -controller.Gamma_bar+ ...
                    controller.epsilon_Gamma_bar*eye(system.dim_x), ...
                -Gamma_Y'*Qv'; ...
            data.data_B',-Qv*Gamma_Y,-data.data_A];
        vertex_lmi = 0.5*(vertex_lmi+vertex_lmi');
        validation.maximum_vertex_lmi_eigenvalue = max( ...
            validation.maximum_vertex_lmi_eigenvalue, ...
            max(real(eig(vertex_lmi))));
    end

    closed_loop_jacobian = ...
        system.A_d0+system.B_d*controller.K;
    validation.closed_loop_eigenvalue = eig(closed_loop_jacobian);
    validation.closed_loop_spectral_radius = ...
        max(abs(validation.closed_loop_eigenvalue));

    lambda_Gamma = real(eig(Gamma_symmetric));
    lambda_min_Gamma = min(lambda_Gamma);
    lambda_max_Gamma = max(lambda_Gamma);
    validation.roa_denominator = lambda_max_Gamma^2- ...
        controller.epsilon_Gamma_bar*lambda_min_Gamma;
    if validation.roa_denominator > 0 && ...
            lambda_min_Gamma > 0 && lambda_max_Gamma > 0
        roa_scale = sqrt(lambda_max_Gamma*lambda_min_Gamma/ ...
                         validation.roa_denominator);
        validation.recovered_roa_radius = ...
            min(setting.domain_radius,roa_scale*setting.domain_radius);
    end

    validation.soft_solver_solved = solver.has_candidate;
    validation.soft_local_stability_passed = ...
        isfinite(validation.closed_loop_spectral_radius) && ...
        validation.closed_loop_spectral_radius < 1;
    validation.soft_feasible = ...
        validation.soft_solver_solved && ...
        validation.soft_local_stability_passed;

    validation.Gamma_passed = ...
        validation.Gamma_margin >= -tolerance;
    validation.epsilon_passed = ...
        validation.epsilon_margin_value >= -tolerance;
    validation.vertex_lmi_passed = ...
        isfinite(validation.maximum_vertex_lmi_eigenvalue) && ...
        validation.maximum_vertex_lmi_eigenvalue <= tolerance;
    validation.local_stability_passed = ...
        isfinite(validation.closed_loop_spectral_radius) && ...
        validation.closed_loop_spectral_radius <= ...
            1-setting.stability_margin;
    validation.roa_radius_passed = ...
        isfinite(validation.recovered_roa_radius) && ...
        validation.recovered_roa_radius > 0 && ...
        validation.recovered_roa_radius <= ...
            setting.domain_radius+tolerance;

    validation.configured_feasible = ...
        solver.solved && validation.Gamma_passed && ...
        validation.epsilon_passed && ...
        validation.vertex_lmi_passed && ...
        validation.local_stability_passed && ...
        validation.roa_radius_passed;
    validation.strict_feasible = validation.configured_feasible;

    if setting.print_validation
        funcPrintValidation(validation,solver);
    end
end


function funcPrintValidation(validation,solver)
    fprintf('  Controller feasibility validation:\n');
    fprintf('    %s standard solver status accepted      = %d\n', ...
        funcValidationMark(solver.solved),solver.solved);
    fprintf('    %s recoverable solver candidate         = %d\n', ...
        funcValidationMark(solver.has_candidate),solver.has_candidate);
    fprintf('    %s minimum Gamma_bar-margin eigenvalue  = %.3e\n', ...
        funcValidationMark(validation.Gamma_passed), ...
        validation.Gamma_margin);
    fprintf('    %s epsilon_Gamma_bar margin             = %.3e\n', ...
        funcValidationMark(validation.epsilon_passed), ...
        validation.epsilon_margin_value);
    fprintf('    %s maximum vertex-LMI eigenvalue        = %.3e\n', ...
        funcValidationMark(validation.vertex_lmi_passed), ...
        validation.maximum_vertex_lmi_eigenvalue);
    fprintf('    %s closed-loop spectral radius          = %.3e\n', ...
        funcValidationMark(validation.local_stability_passed), ...
        validation.closed_loop_spectral_radius);
    fprintf('    %s soft unit-circle test                = %d\n', ...
        funcValidationMark(validation.soft_local_stability_passed), ...
        validation.soft_local_stability_passed);
    fprintf('    %s recovered ROA radius                 = %.3e\n', ...
        funcValidationMark(validation.roa_radius_passed), ...
        validation.recovered_roa_radius);
    fprintf('      condition number of Gamma_bar         = %.3e\n', ...
        validation.Gamma_condition);
    fprintf('    %s soft validation passed               = %d\n', ...
        funcValidationMark(validation.soft_feasible), ...
        validation.soft_feasible);
    fprintf('    %s configured validation passed         = %d\n', ...
        funcValidationMark(validation.configured_feasible), ...
        validation.configured_feasible);
    fprintf('    %s strict validation passed             = %d\n', ...
        funcValidationMark(validation.strict_feasible), ...
        validation.strict_feasible);
end


function result = funcInitializeRunResult( ...
    system,data,setting,experiment_id)

    result = funcEmptyRunResult();
    result.dim_x = system.dim_x;
    result.order_K = 0;
    result.experiment_id = experiment_id;
    result.input_seed = data.input_seed;
    result.initial_seed = data.initial_seed;
    result.noise_seed = data.noise_seed;
    result.requested_trajectory_num = data.requested_trajectory_num;
    result.trajectory_num = data.trajectory_num;
    result.in_trajectory_samp_num = data.in_trajectory_samp_num;
    result.required_sample_num = data.required_sample_num;
    result.sample_num = data.sample_num;
    result.data_rank = data.data_rank;
    result.expected_rank = data.expected_rank;
    result.rank_passed = data.rank_passed;
    result.data_equation_passed = data.data_equation_passed;
    result.noise_bound_passed = data.noise_bound_passed;
    result.data_valid = data.data_valid;
    result.data_condition = data.data_condition;
    result.data_equation_residual = data.data_equation_residual;
    result.maximum_noise_bound_eigenvalue = ...
        data.maximum_noise_bound_eigenvalue;
    result.domain_radius = setting.domain_radius;
    result.vertex_num = 2^system.dim_x;
    result.lmi_dimension = 4*system.dim_x+1;
    result.lmi_constraint_num = result.vertex_num;
    result.scalar_decision_variable_num = ...
        (system.dim_x+1)*(system.dim_x+2)/2;
end


function result = funcStoreSolveResult(result,solver,validation,timing)
    result.solver_problem = solver.problem;
    result.solver_info = solver.info;
    result.solver_solved = solver.solved;
    result.solver_candidate = solver.has_candidate;
    result.solver_native_time = solver.native_solve_time;
    result.timed_out = solver.timed_out;
    result.memory_exceeded = solver.memory_exceeded;
    result.memory_usage_ratio_at_stop = ...
        solver.memory_usage_ratio_at_stop;
    result.available_memory_gb_at_stop = ...
        solver.available_memory_gb_at_stop;
    result.vertex_generation_time = timing.vertex_generation_time;
    result.model_build_time = timing.model_build_time;
    result.attempt_solve_time = timing.attempt_solve_time;
    result.validation_time = timing.validation_time;
    result.attempt_synthesis_time = timing.attempt_synthesis_time;

    if solver.memory_exceeded
        result.status = 'memory limit';
        return;
    elseif solver.timed_out
        result.status = 'timeout';
        return;
    end

    if ~solver.has_candidate
        if solver.solved
            result.status = 'controller recovery failed';
        else
            result.status = 'solver failed';
        end
        return;
    end

    result.Gamma_margin = validation.Gamma_margin;
    result.epsilon_margin_value = validation.epsilon_margin_value;
    result.maximum_vertex_lmi_eigenvalue = ...
        validation.maximum_vertex_lmi_eigenvalue;
    result.closed_loop_spectral_radius = ...
        validation.closed_loop_spectral_radius;
    result.recovered_roa_radius = validation.recovered_roa_radius;
    result.Gamma_condition = validation.Gamma_condition;
    result.soft_solver_solved = validation.soft_solver_solved;
    result.soft_local_stability_passed = ...
        validation.soft_local_stability_passed;
    result.soft_feasible = validation.soft_feasible;
    result.configured_feasible = validation.configured_feasible;
    result.strict_feasible = validation.strict_feasible;
    result.feasible = validation.configured_feasible;

    if validation.soft_feasible
        result.soft_solve_time = timing.attempt_solve_time;
        result.soft_synthesis_time = timing.attempt_synthesis_time;
    end
    if validation.configured_feasible
        result.solve_time = timing.attempt_solve_time;
        result.synthesis_time = timing.attempt_synthesis_time;
        result.status = 'validated feasible';
    elseif validation.soft_feasible
        result.status = 'soft feasible only';
    else
        result.status = 'numerical validation failed';
    end
end


function summary = funcSummarizeResults(result,dim_x_list,num_Exp)
    summary = repmat(funcEmptySummaryResult(),numel(dim_x_list),1);

    for idim = 1:numel(dim_x_list)
        dim_x = dim_x_list(idim);
        current = result([result.dim_x] == dim_x);
        configured_select = [current.configured_feasible];
        soft_select = [current.soft_feasible];

        summary(idim).dim_x = dim_x;
        summary(idim).order_K = 0;
        summary(idim).requested_trajectory_num = ...
            current(1).requested_trajectory_num;
        summary(idim).trajectory_num = current(1).trajectory_num;
        summary(idim).in_trajectory_samp_num = ...
            current(1).in_trajectory_samp_num;
        summary(idim).required_sample_num = ...
            current(1).required_sample_num;
        summary(idim).sample_num = current(1).sample_num;
        summary(idim).attempted_num = num_Exp;
        summary(idim).rank_valid_num = sum([current.rank_passed]);
        summary(idim).data_valid_num = sum([current.data_valid]);
        summary(idim).solver_success_num = sum([current.solver_solved]);
        summary(idim).solver_candidate_num = sum([current.solver_candidate]);
        summary(idim).soft_feasible_num = sum(soft_select);
        summary(idim).configured_feasible_num = sum(configured_select);
        summary(idim).strict_feasible_num = ...
            sum([current.strict_feasible]);
        summary(idim).timeout_num = sum([current.timed_out]);
        summary(idim).memory_limit_num = ...
            sum([current.memory_exceeded]);
        summary(idim).soft_solvability_rate = ...
            summary(idim).soft_feasible_num/num_Exp;
        summary(idim).configured_solvability_rate = ...
            summary(idim).configured_feasible_num/num_Exp;
        summary(idim).strict_solvability_rate = ...
            summary(idim).strict_feasible_num/num_Exp;
        summary(idim).solvability_rate = ...
            summary(idim).configured_solvability_rate;
        summary(idim).mean_soft_solve_time = ...
            funcFiniteMean([current.soft_solve_time]);
        summary(idim).median_soft_solve_time = ...
            funcFiniteMedian([current.soft_solve_time]);
        summary(idim).mean_soft_synthesis_time = ...
            funcFiniteMean([current.soft_synthesis_time]);
        summary(idim).median_soft_synthesis_time = ...
            funcFiniteMedian([current.soft_synthesis_time]);
        summary(idim).mean_solve_time = ...
            funcFiniteMean([current.solve_time]);
        summary(idim).median_solve_time = ...
            funcFiniteMedian([current.solve_time]);
        summary(idim).mean_synthesis_time = ...
            funcFiniteMean([current.synthesis_time]);
        summary(idim).median_synthesis_time = ...
            funcFiniteMedian([current.synthesis_time]);
        summary(idim).mean_attempt_solve_time = ...
            funcFiniteMean([current.attempt_solve_time]);
        summary(idim).median_attempt_solve_time = ...
            funcFiniteMedian([current.attempt_solve_time]);
        summary(idim).mean_attempt_synthesis_time = ...
            funcFiniteMean([current.attempt_synthesis_time]);
        summary(idim).median_attempt_synthesis_time = ...
            funcFiniteMedian([current.attempt_synthesis_time]);
        summary(idim).mean_data_condition = ...
            funcFiniteMean([current.data_condition]);
        summary(idim).median_data_condition = ...
            funcFiniteMedian([current.data_condition]);
        summary(idim).mean_recovered_roa_radius = ...
            funcFiniteMean([current(configured_select).recovered_roa_radius]);
        summary(idim).median_recovered_roa_radius = ...
            funcFiniteMedian([current(configured_select).recovered_roa_radius]);
        summary(idim).mean_soft_spectral_radius = ...
            funcFiniteMean([current(soft_select).closed_loop_spectral_radius]);
        summary(idim).median_soft_spectral_radius = ...
            funcFiniteMedian([current(soft_select).closed_loop_spectral_radius]);
        summary(idim).vertex_num = 2^dim_x;
        summary(idim).lmi_dimension = 4*dim_x+1;
        summary(idim).lmi_constraint_num = 2^dim_x;
        summary(idim).scalar_decision_variable_num = ...
            (dim_x+1)*(dim_x+2)/2;
    end
end


function value = funcFiniteMean(candidate)
    candidate = candidate(isfinite(candidate));
    if isempty(candidate)
        value = NaN;
    else
        value = mean(candidate);
    end
end


function value = funcFiniteMedian(candidate)
    candidate = candidate(isfinite(candidate));
    if isempty(candidate)
        value = NaN;
    else
        value = median(candidate);
    end
end


function result = funcEmptyRunResult()
    result.dim_x = NaN;
    result.order_K = NaN;
    result.experiment_id = NaN;
    result.input_seed = NaN;
    result.initial_seed = NaN;
    result.noise_seed = NaN;
    result.requested_trajectory_num = NaN;
    result.trajectory_num = NaN;
    result.in_trajectory_samp_num = NaN;
    result.required_sample_num = NaN;
    result.sample_num = NaN;
    result.data_rank = NaN;
    result.expected_rank = NaN;
    result.rank_passed = false;
    result.data_equation_passed = false;
    result.noise_bound_passed = false;
    result.data_valid = false;
    result.data_condition = NaN;
    result.data_equation_residual = NaN;
    result.maximum_noise_bound_eigenvalue = NaN;
    result.domain_radius = NaN;
    result.vertex_num = NaN;
    result.lmi_dimension = NaN;
    result.lmi_constraint_num = NaN;
    result.scalar_decision_variable_num = NaN;
    result.solver_problem = NaN;
    result.solver_info = '';
    result.solver_solved = false;
    result.solver_candidate = false;
    result.solver_native_time = NaN;
    result.timed_out = false;
    result.memory_exceeded = false;
    result.memory_usage_ratio_at_stop = NaN;
    result.available_memory_gb_at_stop = NaN;
    result.status = 'not run';
    result.vertex_generation_time = NaN;
    result.model_build_time = NaN;
    result.attempt_solve_time = NaN;
    result.validation_time = NaN;
    result.attempt_synthesis_time = NaN;
    result.soft_solve_time = NaN;
    result.soft_synthesis_time = NaN;
    result.solve_time = NaN;
    result.synthesis_time = NaN;
    result.Gamma_margin = NaN;
    result.epsilon_margin_value = NaN;
    result.maximum_vertex_lmi_eigenvalue = NaN;
    result.closed_loop_spectral_radius = NaN;
    result.recovered_roa_radius = NaN;
    result.Gamma_condition = NaN;
    result.soft_solver_solved = false;
    result.soft_local_stability_passed = false;
    result.soft_feasible = false;
    result.configured_feasible = false;
    result.strict_feasible = false;
    result.feasible = false;
end


function summary = funcEmptySummaryResult()
    summary.dim_x = NaN;
    summary.order_K = NaN;
    summary.requested_trajectory_num = NaN;
    summary.trajectory_num = NaN;
    summary.in_trajectory_samp_num = NaN;
    summary.required_sample_num = NaN;
    summary.sample_num = NaN;
    summary.attempted_num = NaN;
    summary.rank_valid_num = NaN;
    summary.data_valid_num = NaN;
    summary.solver_success_num = NaN;
    summary.solver_candidate_num = NaN;
    summary.soft_feasible_num = NaN;
    summary.configured_feasible_num = NaN;
    summary.strict_feasible_num = NaN;
    summary.timeout_num = NaN;
    summary.memory_limit_num = NaN;
    summary.soft_solvability_rate = NaN;
    summary.configured_solvability_rate = NaN;
    summary.strict_solvability_rate = NaN;
    summary.solvability_rate = NaN;
    summary.mean_soft_solve_time = NaN;
    summary.median_soft_solve_time = NaN;
    summary.mean_soft_synthesis_time = NaN;
    summary.median_soft_synthesis_time = NaN;
    summary.mean_solve_time = NaN;
    summary.median_solve_time = NaN;
    summary.mean_synthesis_time = NaN;
    summary.median_synthesis_time = NaN;
    summary.mean_attempt_solve_time = NaN;
    summary.median_attempt_solve_time = NaN;
    summary.mean_attempt_synthesis_time = NaN;
    summary.median_attempt_synthesis_time = NaN;
    summary.mean_data_condition = NaN;
    summary.median_data_condition = NaN;
    summary.mean_recovered_roa_radius = NaN;
    summary.median_recovered_roa_radius = NaN;
    summary.mean_soft_spectral_radius = NaN;
    summary.median_soft_spectral_radius = NaN;
    summary.vertex_num = NaN;
    summary.lmi_dimension = NaN;
    summary.lmi_constraint_num = NaN;
    summary.scalar_decision_variable_num = NaN;
end


function solver = funcEmptySolverResult()
    solver.problem = NaN;
    solver.info = '';
    solver.solved = false;
    solver.has_candidate = false;
    solver.timed_out = false;
    solver.memory_exceeded = false;
    solver.memory_usage_ratio_at_stop = NaN;
    solver.available_memory_gb_at_stop = NaN;
    solver.native_solve_time = NaN;
end


function [controller,solver,validation,timing] = ...
    funcFailedIsolatedSolve( ...
        problem,reason,elapsed_time,memory_exceeded, ...
        memory_usage_ratio,available_memory_gb)

    controller = [];
    solver = funcEmptySolverResult();
    validation = funcEmptyValidation();
    timing = funcEmptyTimingResult();
    timing.attempt_solve_time = elapsed_time;
    timing.attempt_synthesis_time = elapsed_time;
    solver.problem = problem;
    solver.info = reason;
    solver.memory_exceeded = memory_exceeded;
    solver.memory_usage_ratio_at_stop = memory_usage_ratio;
    solver.available_memory_gb_at_stop = available_memory_gb;
end


function [memory_pool,memory_pool_owned] = ...
    funcEnsureMemoryWorkerPool(memory_pool,memory_pool_owned,setting)

    if ~license('test','Distrib_Computing_Toolbox') || ...
            exist('parfeval','file') ~= 2
        error(['The memory guard requires Parallel Computing Toolbox ', ...
               'so each solve can run in an isolated process.']);
    end

    [memory_usage_ratio,available_memory_gb] = ...
        funcReadLinuxMemoryUsage();
    if ~isfinite(memory_usage_ratio) || ~isfinite(available_memory_gb)
        error(['The memory guard could not read /proc/meminfo. ', ...
               'This implementation requires Linux memory information.']);
    end

    if ~isempty(memory_pool)
        try
            if isvalid(memory_pool) && memory_pool.NumWorkers > 0
                if ~isa(memory_pool,'parallel.ProcessPool')
                    error(['The memory guard cannot use a thread pool ', ...
                           'because terminating it would also terminate ', ...
                           'the main MATLAB process.']);
                end
                return;
            end
        catch pool_error
            if contains(pool_error.message,'thread pool')
                rethrow(pool_error);
            end
        end
    end

    existing_pool = gcp('nocreate');
    if ~isempty(existing_pool) && existing_pool.NumWorkers > 0
        if ~isa(existing_pool,'parallel.ProcessPool')
            error(['Close the existing thread pool before running this ', ...
                   'experiment; the memory guard requires process workers.']);
        end
        memory_pool = existing_pool;
        return;
    end

    if ~isempty(existing_pool)
        funcDeleteMemoryWorkerPool(existing_pool);
    end
    memory_pool = parpool('Processes',1,'IdleTimeout', ...
        setting.memory_pool_idle_timeout_minutes);
    memory_pool_owned = true;
end


function funcDeleteMemoryWorkerPool(memory_pool)
    if isempty(memory_pool)
        return;
    end
    try
        delete(memory_pool);
    catch
    end
end


function [memory_usage_ratio,available_memory_gb] = ...
    funcReadLinuxMemoryUsage()

    memory_usage_ratio = NaN;
    available_memory_gb = NaN;
    try
        memory_information = fileread('/proc/meminfo');
        total_token = regexp(memory_information, ...
            'MemTotal:\s+(\d+)\s+kB','tokens','once');
        available_token = regexp(memory_information, ...
            'MemAvailable:\s+(\d+)\s+kB','tokens','once');
        if isempty(total_token) || isempty(available_token)
            return;
        end
        total_memory_kb = str2double(total_token{1});
        available_memory_kb = str2double(available_token{1});
        memory_usage_ratio = ...
            1-available_memory_kb/total_memory_kb;
        available_memory_gb = available_memory_kb/1024^2;
    catch
    end
end


function exceeded = funcMemoryLimitExceeded( ...
    memory_usage_ratio,available_memory_gb,setting)

    exceeded = ...
        isfinite(memory_usage_ratio) && ...
        memory_usage_ratio >= setting.maximum_memory_usage_ratio;
    exceeded = exceeded || ...
        (isfinite(available_memory_gb) && ...
         available_memory_gb <= setting.minimum_available_memory_gb);
end


function message = funcMemoryLimitMessage( ...
    memory_usage_ratio,available_memory_gb,setting)

    message = sprintf([ ...
        'system memory %.1f%% used, %.2f GiB available; ', ...
        'limits are %.1f%% and %.2f GiB'], ...
        100*memory_usage_ratio,available_memory_gb, ...
        100*setting.maximum_memory_usage_ratio, ...
        setting.minimum_available_memory_gb);
end


function funcTerminateWorkerProcess(worker_pid,grace_period)
    worker_pid = round(double(worker_pid));
    if ~isfinite(worker_pid) || worker_pid <= 1 || ...
            worker_pid == feature('getpid')
        return;
    end

    process_directory = sprintf('/proc/%d',worker_pid);
    system(sprintf('/bin/kill -TERM %d',worker_pid));
    termination_timer = tic;
    while isfolder(process_directory) && ...
            toc(termination_timer) < grace_period
        pause(0.1);
    end
    if isfolder(process_directory)
        system(sprintf('/bin/kill -KILL %d',worker_pid));
    end
end


function validation = funcEmptyValidation()
    validation.Gamma_margin = NaN;
    validation.epsilon_margin_value = NaN;
    validation.maximum_vertex_lmi_eigenvalue = NaN;
    validation.closed_loop_eigenvalue = [];
    validation.closed_loop_spectral_radius = NaN;
    validation.roa_denominator = NaN;
    validation.recovered_roa_radius = NaN;
    validation.Gamma_condition = NaN;
    validation.soft_solver_solved = false;
    validation.soft_local_stability_passed = false;
    validation.soft_feasible = false;
    validation.Gamma_passed = false;
    validation.epsilon_passed = false;
    validation.vertex_lmi_passed = false;
    validation.local_stability_passed = false;
    validation.roa_radius_passed = false;
    validation.configured_feasible = false;
    validation.strict_feasible = false;
end


function timing = funcEmptyTimingResult()
    timing.vertex_generation_time = NaN;
    timing.model_build_time = NaN;
    timing.attempt_solve_time = NaN;
    timing.validation_time = NaN;
    timing.attempt_synthesis_time = NaN;
end


function mark = funcValidationMark(is_passed)
    if is_passed
        mark = '!';
    else
        mark = '?';
    end
end


function candidate_status = funcSolverReturnedCandidate(problem,info)
    info_lower = lower(char(info));
    hard_failure = contains(info_lower,'infeasible') || ...
                   contains(info_lower,'unbounded') || ...
                   contains(info_lower,'no solution');
    candidate_code = ismember(problem,[0,3,4,5]);
    candidate_message = contains(info_lower,'solved') || ...
                        contains(info_lower,'optimal') || ...
                        contains(info_lower,'near optimal') || ...
                        contains(info_lower,'inaccurate');
    candidate_status = ...
        (candidate_code || candidate_message) && ~hard_failure;
end


function timed_out = funcSolverTimedOut(info,solve_time,time_limit)
    info_lower = lower(char(info));
    timeout_message = contains(info_lower,'time limit') || ...
                      contains(info_lower,'maximum time') || ...
                      contains(info_lower,'optimizer_max_time');
    timed_out = solve_time >= time_limit || timeout_message;
end


function option = funcLmiOptions(setting)
    option = sdpsettings( ...
        'solver',setting.solver, ...
        'verbose',setting.verbose, ...
        'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS', ...
            setting.feasibility_tolerance, ...
        'mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS', ...
            setting.feasibility_tolerance, ...
        'mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP', ...
            setting.relative_gap_tolerance, ...
        'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME', ...
            setting.solve_timeout);
end


function funcSetupYalmip()
    if isempty(which('sdpvar'))
        yalmip_path = '/home/lidong/Documents/MATLAB/YALMIP-master';
        if ~isfolder(yalmip_path)
            error('YALMIP was not found. Set yalmip_path in funcSetupYalmip().');
        end
        addpath(genpath(yalmip_path));
    end

    if isempty(which('mosekopt'))
        mosek_path = ['/home/lidong/Documents/MATLAB/mosektoolslinux64x86/', ...
                      'mosek/11.0/toolbox/r2019b'];
        if isfolder(mosek_path)
            addpath(mosek_path);
        end
    end

    rehash toolboxcache;

    if isempty(which('sdpvar'))
        error('YALMIP function sdpvar was not found on the MATLAB path.');
    end
    if isempty(which('mosekopt'))
        error('MOSEK function mosekopt was not found on the MATLAB path.');
    end
end
