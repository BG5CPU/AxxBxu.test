clear; close all; clc;


%% baseline test experiment ==============================================

% Baseline: Corollary 1 in
% M. Guo, C. De Persis, and P. Tesi, "Data-Driven Stabilization of
% Nonlinear Polynomial Systems With Noisy Data", IEEE TAC, 2022.
%
% Dense continuous-time benchmark:
%   xdot = (delta*I + S_n)*x
%          -(alpha/sqrt(n))*B_en*sum_i(x_i^3)+B_en*u,
% where B_en = [0,...,0,1]' and S_n is the shift matrix. The factor
% 1/sqrt(n) keeps the Euclidean norm of the dense cubic coefficient
% vector independent of the state dimension.

dim_x_list = 1:9;
order_K_list = 0:4;
num_Exp = 30;

system_setting.delta = 0.1;
system_setting.alpha = 1;

% Continuous-time sampling interval and holding time of each input sample.
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
% Each entry of the unscaled derivative noise lies in this symmetric bound.
data_setting.base_noise_bound = 1;
% Scales the additive derivative-data noise: D0=noise_gain*base_noise.
data_setting.noise_gain = 0;
% Defines RD2=margin*(D0*D0') for the noise-energy upper bound.
data_setting.noise_bound_margin = 1.00;
% Numerical tolerance for the known-model data-equation residual.
data_setting.data_equation_tolerance = 1e-10;
% Tolerance for lambda_max(D0*D0'-RD2) <= 0.
data_setting.noise_bound_tolerance = 1e-10;
% Reproducible master seed used to derive input, initial-state and noise seeds.
data_setting.base_seed = 20260828;
% Relative and absolute tolerances used by ode45 during data generation.
data_setting.ode = odeset('RelTol',1e-9,'AbsTol',1e-11);

% SDP solver selected through YALMIP's SOS module.
sos_setting.solver = 'mosek';
% YALMIP/MOSEK console verbosity: 0 suppresses routine solver output.
sos_setting.verbose = 0;
% Enforces P >= matrix_margin*I in Corollary 1.
sos_setting.matrix_margin = 1e-9;
% Enforces the epsilon_1 Gram matrix >= epsilon_margin*I.
sos_setting.epsilon_margin = 1e-9;
% MOSEK primal and dual feasibility tolerances.
sos_setting.feasibility_tolerance = 1e-4;
% MOSEK relative primal-dual objective-gap tolerance.
sos_setting.relative_gap_tolerance = 1e-4;
% Post-solution tolerance for matrix and SOS Gram eigenvalue checks.
sos_setting.validation_tolerance = 1e-4;
% Maximum accepted SOS coefficient-matching residual for strict validation.
sos_setting.sos_coefficient_tolerance = 1e-4;
% Requires max(real(eig(J_closed_loop(0)))) <= -stability_margin.
sos_setting.stability_margin = 1e-6;
% Central finite-difference step used to compute the closed-loop Jacobian.
sos_setting.jacobian_step = 1e-6;
% Maximum solver time in seconds for each trial.
sos_setting.solve_timeout = 500;
% Run every SOS solve in an isolated MATLAB process and monitor total system
% memory, so an oversized solve cannot terminate the main MATLAB process.
% This protection requires Parallel Computing Toolbox on Linux.
sos_setting.memory_guard_enabled = true;
% Stop a trial when 1-MemAvailable/MemTotal reaches this fraction.
sos_setting.maximum_memory_usage_ratio = 0.90;
% Also stop when Linux MemAvailable falls below this many GiB.
sos_setting.minimum_available_memory_gb = 1.5;
% Interval in seconds between memory checks while a worker is running.
sos_setting.memory_poll_period = 0.25;
% Seconds allowed after SIGTERM before the isolated worker receives SIGKILL.
sos_setting.memory_termination_grace_period = 2.0;
% An owned worker pool closes itself after this many idle minutes. The script
% does not synchronously delete the pool at the end, so it returns promptly.
sos_setting.memory_pool_idle_timeout_minutes = 1;
% False: configured feasibility records but does not require coefficient fit.
sos_setting.use_sos_coefficient_as_hard_metric = false;
% Prints !/? marks and all recovered-controller validation metrics.
sos_setting.print_validation = true;

output_directory = fullfile(fileparts(mfilename('fullpath')),'output');
output_mat = fullfile(output_directory,'baseline_sos_test_results.mat');
output_run_csv = fullfile(output_directory,'baseline_sos_test_runs.csv');
output_summary_csv = ...
    fullfile(output_directory,'baseline_sos_test_summary.csv');

% Output contents:
%   MAT file    - all run/summary structs, controllers, data and settings.
%   runs CSV    - one row per (dim_x, order_K, experiment_id) trial.
%   summary CSV - one aggregated row per (dim_x, order_K) group.
% The terminal first prints the same common-summary columns as the Axx code:
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
%   data_condition              - cond([Ubar0;Z0]); smaller is better.
% Every time and condition quantity is summarized by both mean and median.
% The second terminal table contains these SOS-only diagnostics:
%   dim_x                       - state dimension of the tested system;
%   order_K                     - maximum total degree allowed in K(x);
%   epsilon_degree              - maximum (even) degree used for epsilon_1(x)
%                                 and epsilon_2(x), equal to
%                                 max(2,2*ceil(order_K/2));
%   data_valid_num              - trials satisfying both the concatenated
%                                 data-equation residual and noise-bound
%                                 eigenvalue tests; rank is counted separately;
%   memory_limit_num            - isolated solves stopped by the memory guard;
%   solver_candidate_num        - non-timeout trials for which YALMIP/MOSEK
%                                 returned finite variables and a recoverable
%                                 Corollary-1 controller;
%   soft_feasible_num/rate      - candidate count/fraction additionally having
%                                 max(real(eig(J_cl(0))))<0; certificate
%                                 residuals and standard solver-success status
%                                 are not required by this soft test;
%   soft_solve_time             - solver-call time of soft-feasible trials;
%   soft_synthesis_time         - synthesis time of soft-feasible trials;
%   soft_max_real_part          - max(real(eig(J_cl(0)))), summarized only over
%                                 soft-feasible trials; negative means local
%                                 continuous-time stability at the origin;
%   sos_monomial_num            - size of the scalar state-monomial basis used
%                                 by the largest matrix-SOS Gram block, inferred
%                                 as gram_dimension/(3*dim_x+1);
%   sos_gram_dimension          - row/column dimension of the largest Gram
%                                 matrix block returned by YALMIP's sosd;
%   sos_gram_block_num          - number of Gram matrix blocks returned by sosd.
% The SOS size quantities are available for recoverable candidate trials whose
% Gram decomposition was obtained; they do not require configured feasibility.
% Each size quantity is summarized by both mean and median. Larger monomial,
% Gram-dimension and Gram-block values generally indicate a larger SDP.

fprintf('Baseline SOS test experiment\n');
fprintf('  dim_x       = %s\n',mat2str(dim_x_list));
fprintf('  order_K     = %s\n',mat2str(order_K_list));
fprintf('  num_Exp     = %d\n',num_Exp);
fprintf('  noise_gain  = %.3e\n',data_setting.noise_gain);
fprintf('  sample time = %.3e s\n',data_setting.sample_period);
fprintf('  solve timeout = %.1f s\n',sos_setting.solve_timeout);
fprintf('  memory guard = %.1f%% used or %.2f GiB available\n', ...
        100*sos_setting.maximum_memory_usage_ratio, ...
        sos_setting.minimum_available_memory_gb);
fprintf(['  SOS coefficient residual is recorded but is not a hard ', ...
         'metric in this test.\n\n']);

funcSetupYalmip();

memory_pool = [];
memory_pool_owned = false;

% record the results of each run and the summary of all runs
num_dimension = numel(dim_x_list);
num_order = numel(order_K_list);
num_run = num_dimension*num_order*num_Exp;
run_result = repmat(funcEmptyRunResult(),num_run,1);
controller_result = cell(num_run,1);
soft_controller_result = cell(num_run,1);
candidate_controller_result = cell(num_run,1);
experiment_data = cell(num_dimension,num_Exp);
system_result = cell(num_dimension,1);
run_index = 0;

for idim = 1:num_dimension
    dim_x = dim_x_list(idim);
    system = funcBuildSystem(dim_x,system_setting);
    system_result{idim} = system;

    % Generate each dataset once. All controller orders use the exact same
    % dataset for a fixed (dim_x, experiment_id).
    for iexp = 1:num_Exp
        experiment_data{idim,iexp} = ...
            funcGenerateExperimentData(system,data_setting,iexp);
    end

    for iorder = 1:num_order
        order_K = order_K_list(iorder);
        epsilon_degree = max(2,2*ceil(order_K/2));

        setting = sos_setting;
        setting.controller_degree = order_K;
        setting.epsilon1_degree = epsilon_degree;
        setting.epsilon2_degree = epsilon_degree;

        fprintf(['\n============================================================', ...
                 '================\n']);
        fprintf('dim_x = %d, order_K = %d, epsilon degree = %d\n', ...
                dim_x,order_K,epsilon_degree);

        for iexp = 1:num_Exp
            run_index = run_index+1;
            data = experiment_data{idim,iexp};
            result = funcInitializeRunResult( ...
                system,data,order_K,epsilon_degree,iexp);

            fprintf('\nRun %d/%d: experiment_id = %d\n', ...
                    run_index,num_run,iexp);
            fprintf(['  trajectories = %d (requested %d), ', ...
                     'samples/trajectory = %d, total = %d, ', ...
                     'required >= %d\n'], ...
                    data.trajectory_num,data.requested_trajectory_num, ...
                    data.in_trajectory_samp_num,data.sample_num, ...
                    data.required_sample_num);
            fprintf('  rank([Ubar0; Z0]) = %d/%d, condition = %.3e\n', ...
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
                    funcRunCorollary1WithMemoryGuard( ...
                        system,data,setting,memory_pool,memory_pool_owned);
                if solver.has_candidate
                    candidate_controller_result{run_index} = controller;
                end
                if validation.soft_feasible
                    soft_controller_result{run_index} = controller;
                end
                if validation.passed
                    controller_result{run_index} = controller;
                end
                result = funcStoreSolveResult( ...
                    result,solver,validation,timing);
            end

            run_result(run_index) = result;
        end
    end
end

summary_result = ...
    funcSummarizeResults(run_result,dim_x_list,order_K_list,num_Exp);
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
    'dim_x','order_K','epsilon_degree', ...
    'data_valid_num','memory_limit_num','solver_candidate_num', ...
    'soft_feasible_num','soft_solvability_rate', ...
    'mean_soft_solve_time','median_soft_solve_time', ...
    'mean_soft_synthesis_time','median_soft_synthesis_time', ...
    'mean_soft_max_real_part','median_soft_max_real_part', ...
    'mean_sos_monomial_num','median_sos_monomial_num', ...
    'mean_sos_gram_dimension','median_sos_gram_dimension', ...
    'mean_sos_gram_block_num','median_sos_gram_block_num'};

fprintf('\n\nCommon comparison summary: baseline SOS\n');
disp(summary_table(:,common_summary_column));
fprintf('\nSOS-specific diagnostics\n');
disp(summary_table(:,method_summary_column));

if ~isfolder(output_directory)
    mkdir(output_directory);
end
save(output_mat,'run_result','summary_result','run_table','summary_table', ...
     'controller_result','soft_controller_result', ...
     'candidate_controller_result', ...
     'experiment_data','system_result', ...
     'dim_x_list','order_K_list','num_Exp', ...
     'system_setting','data_setting','sos_setting');
writetable(run_table,output_run_csv);
writetable(summary_table,output_summary_csv);

fprintf('Results saved to:\n  %s\n  %s\n  %s\n', ...
        output_mat,output_run_csv,output_summary_csv);

if memory_pool_owned && ~isempty(memory_pool)
    fprintf(['Memory-guard worker pool will close automatically after ', ...
             '%.1f idle minute(s).\n'], ...
            sos_setting.memory_pool_idle_timeout_minutes);
end








%% useful functions =======================================================

function system = funcBuildSystem(dim_x,setting)
    S = zeros(dim_x);
    if dim_x > 1
        S(1:dim_x-1,2:dim_x) = eye(dim_x-1);
    end

    B_en = zeros(dim_x,1);
    B_en(end) = 1;

    system.dim_x = dim_x;
    system.delta = setting.delta;
    system.alpha = setting.alpha;
    system.S = S;
    system.B_en = B_en;
    system.A_linear = setting.delta*eye(dim_x)+S;
    system.cubic_coefficient = ones(dim_x,1)/sqrt(dim_x);
    system.A = [system.A_linear, ...
        -setting.alpha*B_en*system.cubic_coefficient'];
    system.B = B_en;
    system.feature_dimension = 2*dim_x;
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
    D0 = setting.noise_gain*base_noise;

    % Each block is a separate trajectory. Only valid within-trajectory
    % samples are concatenated, so the ODE integrator is restarted from the
    % independently sampled initial state at every block boundary.
    X0 = zeros(dim_x,sample_num);
    Xdot_true = zeros(dim_x,sample_num);
    for itrajectory = 1:trajectory_num
        column = (itrajectory-1)*in_trajectory_samp_num+ ...
                 (1:in_trajectory_samp_num);
        [X0(:,column),Xdot_true(:,column)] = ...
            funcCollectTrajectoryData( ...
                system,x_initial(:,itrajectory),U0(:,column),setting);
    end
    X1 = Xdot_true+D0;
    RD2 = setting.noise_bound_margin*(D0*D0');

    Z0 = funcPowerVector(X0);
    Ubar0 = funcLiftInputData(X0,U0);
    Wbar0 = [Ubar0;Z0];
    singular_value = svd(Wbar0);
    if isempty(singular_value) || singular_value(end) <= 0
        data_condition = Inf;
    else
        data_condition = singular_value(1)/singular_value(end);
    end

    % Concatenated-data check in Corollary-1 notation. Here the polynomial
    % state data called X0 in the generic equation are denoted by Z0, while
    % the lifted input data are Ubar0:
    % X1 = A*Z0+B*Ubar0+D0.
    data_equation_residual = norm( ...
        X1-system.A*Z0-system.B*Ubar0-D0,'fro');
    maximum_noise_bound_eigenvalue = ...
        max(real(eig(0.5*((D0*D0'-RD2)+(D0*D0'-RD2)'))));

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
    data.X0 = X0;
    data.Z0 = Z0;
    data.Ubar0 = Ubar0;
    data.Xdot_true = Xdot_true;
    data.base_noise = base_noise;
    data.D0 = D0;
    data.X1 = X1;
    data.RD2 = RD2;
    data.data_rank = rank(Wbar0);
    data.expected_rank = system.regressor_dimension;
    if size(Wbar0,1) ~= data.expected_rank
        error(['The dense data regressor has %d rows, but Corollary 1 ', ...
               'requires 2*dim_x+1 = %d rows.'], ...
              size(Wbar0,1),data.expected_rank);
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


function [X0,Xdot_true] = ...
    funcCollectTrajectoryData(system,x_initial,U0,setting)

    sample_num = numel(U0);
    dim_x = system.dim_x;
    X0 = zeros(dim_x,sample_num);
    Xdot_true = zeros(dim_x,sample_num);
    state = x_initial;

    for it = 1:sample_num
        input = U0(it);
        X0(:,it) = state;
        Xdot_true(:,it) = funcPolynomialSystem(state,input,system);

        [~,state_interval] = ode45( ...
            @(~,x) funcPolynomialSystem(x,input,system), ...
            [0,setting.sample_period],state,setting.ode);
        state = state_interval(end,:)';
    end
end


function dx = funcPolynomialSystem(x,u,system)
    dense_cubic = sum(system.cubic_coefficient.*(x(:).^3));
    dx = system.A_linear*x ...
         -system.alpha*system.B_en*dense_cubic ...
         +system.B_en*u;
end


function [controller,solver,validation,timing, ...
          memory_pool,memory_pool_owned] = ...
    funcRunCorollary1WithMemoryGuard( ...
        system,data,setting,memory_pool,memory_pool_owned)

    if ~setting.memory_guard_enabled
        [controller,solver,validation,timing] = ...
            funcRunCorollary1(system,data,setting);
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
    future = parfeval(memory_pool,@funcRunCorollary1Worker,5, ...
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
        reason = sprintf('isolated SOS worker failed: %s', ...
                         worker_error.message);
        [controller,solver,validation,timing] = ...
            funcFailedIsolatedSolve(-996,reason,elapsed_time,false, ...
                                    memory_usage_ratio,available_memory_gb);
        fprintf('  ? %s\n',reason);
    end
end


function [controller,solver,validation,timing,worker_output] = ...
    funcRunCorollary1Worker(system,data,setting,pid_queue)

    send(pid_queue,feature('getpid'));
    funcSetupYalmip();
    worker_argument = {system,data,setting}; %#ok<NASGU>
    controller = [];
    solver = funcEmptySolverResult();
    validation = funcEmptyValidation();
    timing.solve_time = NaN;
    timing.synthesis_time = NaN;
    timing.validation_time = NaN;
    worker_output = evalc([ ...
        '[controller,solver,validation,timing] = ', ...
        'funcRunCorollary1(worker_argument{1},', ...
        'worker_argument{2},worker_argument{3});']);
end


function [controller,solver,validation,timing] = ...
    funcRunCorollary1(system,data,setting)

    controller = [];
    solver = funcEmptySolverResult();
    validation = funcEmptyValidation();
    timing.solve_time = NaN;
    timing.synthesis_time = NaN;
    timing.validation_time = NaN;

    synthesis_timer = tic;
    yalmip('clear');
    x = sdpvar(system.dim_x,1);
    Zhat = funcZhat(x);
    dZhat_dx = jacobian(Zhat,x);
    H = funcHmatrix(x);
    W = funcWmatrix(x);
    Wbar0 = [data.Ubar0;data.Z0];

    [K,coeff_K,monomial_power] = ...
        funcPolynomialMatrix(x,1,system.dim_x, ...
                             setting.controller_degree);
    P = sdpvar(system.dim_x,system.dim_x,'symmetric');
    [epsilon1,Qepsilon1] = ...
        funcSosPolynomial(x,setting.epsilon1_degree,false);
    [epsilon2,Qepsilon2] = ...
        funcSosPolynomial(x,setting.epsilon2_degree,true);

    Upsilon_D = epsilon2*dZhat_dx* ...
        (data.X1*data.X1'-data.RD2)*dZhat_dx';
    W0KP = [W*K;H*P];
    block_12 = -W0KP' ...
               -epsilon2*dZhat_dx*data.X1*Wbar0';
    block_22 = epsilon2*(Wbar0*Wbar0');
    block_sos = [Upsilon_D-epsilon1*eye(system.dim_x),block_12;
                 block_12',block_22];
    block_sos = 0.5*(block_sos+block_sos');

    constraint = [ ...
        P >= setting.matrix_margin*eye(system.dim_x), ...
        Qepsilon1 >= setting.epsilon_margin*eye(size(Qepsilon1,1)), ...
        Qepsilon2 >= 0, ...
        sos(block_sos)];

    parameter = funcDecisionVariables( ...
        [coeff_K(:);P(:);Qepsilon1(:);Qepsilon2(:)],x);

    try
        solve_timer = tic;
        diagnostic = solvesos( ...
            constraint,[],funcSosOptions(setting),parameter);
        timing.solve_time = toc(solve_timer);
        solver.problem = diagnostic.problem;
        solver.info = diagnostic.info;
        solver.timed_out = funcSolverTimedOut( ...
            diagnostic.info,timing.solve_time,setting.solve_timeout);
        solver.solved = diagnostic.problem == 0 && ~solver.timed_out;

        candidate_status = funcSolverReturnedCandidate( ...
            diagnostic.problem,diagnostic.info);
        if ~solver.timed_out && candidate_status
            try
                controller.name = 'Corollary 1';
                controller.P = value(P);
                controller.Qepsilon1 = value(Qepsilon1);
                controller.Qepsilon2 = value(Qepsilon2);
                controller.gain_coefficient = value(coeff_K);
                controller.monomial_power = monomial_power;

                controller_value = [controller.P(:); ...
                                    controller.Qepsilon1(:); ...
                                    controller.Qepsilon2(:); ...
                                    controller.gain_coefficient(:)];
                solver.has_candidate = ...
                    ~isempty(controller_value) && ...
                    all(isfinite(controller_value));
            catch recovery_error
                solver.info = recovery_error.message;
                solver.has_candidate = false;
            end
        end

        if solver.timed_out
            solver.info = sprintf( ...
                'Timeout after %.3f s (limit %.3f s). %s', ...
                timing.solve_time,setting.solve_timeout,diagnostic.info);
        end
    catch solve_error
        timing.solve_time = toc(solve_timer);
        solver.problem = -999;
        solver.info = solve_error.message;
        solver.timed_out = timing.solve_time >= setting.solve_timeout;
        solver.solved = false;
        solver.has_candidate = false;
    end

    if ~solver.has_candidate
        timing.synthesis_time = toc(synthesis_timer);
        if solver.timed_out
            fprintf('  ? timeout: no controller was accepted.\n');
        else
            fprintf('  ? no recoverable solver candidate: %s\n', ...
                    solver.info);
        end
        return;
    end

    timing.synthesis_time = toc(synthesis_timer);
    validation_timer = tic;
    validation = funcValidateController( ...
        controller,constraint,system,setting,solver);
    timing.validation_time = toc(validation_timer);
    controller.validation = validation;
end


function validation = ...
    funcValidateController(controller,constraint,system,setting,solver)

    validation = funcEmptyValidation();
    tolerance = setting.validation_tolerance;

    validation.P_margin = min(eig(0.5*(controller.P+controller.P') ...
        -setting.matrix_margin*eye(size(controller.P))));
    validation.epsilon1_margin = min(eig( ...
        0.5*(controller.Qepsilon1+controller.Qepsilon1') ...
        -setting.epsilon_margin*eye(size(controller.Qepsilon1))));
    validation.epsilon2_margin = min(eig( ...
        0.5*(controller.Qepsilon2+controller.Qepsilon2')));

    try
        coefficient_check = check(constraint(end));
        validation.sos_coefficient_residual = ...
            max(abs(coefficient_check(:)));
    catch
        validation.sos_coefficient_residual = Inf;
    end

    try
        [~,~,sos_gram] = sosd(constraint(end));
        validation.sos_gram_margin = ...
            funcMinimumGramEigenvalue(sos_gram);
        [validation.sos_monomial_num, ...
         validation.sos_gram_dimension, ...
         validation.sos_gram_block_num] = ...
            funcSosDecompositionSize(sos_gram,system.dim_x);
    catch
        validation.sos_gram_margin = -Inf;
    end

    try
        closed_loop_jacobian = funcClosedLoopJacobian( ...
            controller,system,setting.jacobian_step);
        validation.closed_loop_eigenvalue = eig(closed_loop_jacobian);
        validation.max_closed_loop_real_part = ...
            max(real(validation.closed_loop_eigenvalue));
    catch
        validation.max_closed_loop_real_part = Inf;
    end

    validation.P_passed = validation.P_margin >= -tolerance;
    validation.epsilon1_passed = ...
        validation.epsilon1_margin >= -tolerance;
    validation.epsilon2_passed = ...
        validation.epsilon2_margin >= -tolerance;
    validation.sos_gram_passed = ...
        validation.sos_gram_margin >= -tolerance;
    validation.sos_coefficient_passed = ...
        isfinite(validation.sos_coefficient_residual) && ...
        validation.sos_coefficient_residual <= ...
            setting.sos_coefficient_tolerance;
    validation.local_stability_passed = ...
        isfinite(validation.max_closed_loop_real_part) && ...
        validation.max_closed_loop_real_part <= -setting.stability_margin;
    validation.soft_solver_solved = solver.has_candidate;
    validation.soft_local_stability_passed = ...
        isfinite(validation.max_closed_loop_real_part) && ...
        validation.max_closed_loop_real_part < 0;
    validation.soft_feasible = ...
        validation.soft_solver_solved && ...
        validation.soft_local_stability_passed;

    validation.candidate_certificate_passed = ...
        validation.P_passed && validation.epsilon1_passed && ...
        validation.epsilon2_passed && validation.sos_gram_passed;
    validation.strict_certificate_passed = ...
        validation.candidate_certificate_passed && ...
        validation.sos_coefficient_passed;
    validation.candidate_passed = ...
        solver.solved && validation.candidate_certificate_passed && ...
        validation.local_stability_passed;
    validation.strict_passed = ...
        solver.solved && validation.strict_certificate_passed && ...
        validation.local_stability_passed;

    if setting.use_sos_coefficient_as_hard_metric
        validation.passed = validation.strict_passed;
    else
        validation.passed = validation.candidate_passed;
    end

    if setting.print_validation
        funcPrintValidation(validation,setting,solver);
    end
end


function [monomial_num,gram_dimension,gram_block_num] = ...
    funcSosDecompositionSize(sos_gram,dim_x)

    gram_cell = funcCellify(sos_gram);
    gram_block_num = numel(gram_cell);
    gram_dimension = 0;

    for iblock = 1:gram_block_num
        gram_dimension = max(gram_dimension,size(gram_cell{iblock},1));
    end

    % For this matrix-SOS constraint, the Gram dimension is the scalar
    % state-monomial basis size times the (3*n+1) SOS-matrix dimension.
    sos_matrix_dimension = 3*dim_x+1;
    if mod(gram_dimension,sos_matrix_dimension) == 0
        monomial_num = gram_dimension/sos_matrix_dimension;
    else
        monomial_num = NaN;
    end
end


function output = funcCellify(input)
    if iscell(input)
        output = input;
    else
        output = {input};
    end
end


function minimum_eigenvalue = funcMinimumGramEigenvalue(sos_gram)
    minimum_eigenvalue = Inf;
    if iscell(sos_gram)
        for igram = 1:numel(sos_gram)
            gram = double(sos_gram{igram});
            if ~isempty(gram)
                gram = 0.5*(gram+gram');
                minimum_eigenvalue = ...
                    min(minimum_eigenvalue,min(real(eig(gram))));
            end
        end
    else
        gram = double(sos_gram);
        if ~isempty(gram)
            gram = 0.5*(gram+gram');
            minimum_eigenvalue = min(real(eig(gram)));
        end
    end

    if isinf(minimum_eigenvalue) && minimum_eigenvalue > 0
        minimum_eigenvalue = NaN;
    end
end


function funcPrintValidation(validation,setting,solver)
    fprintf('  Controller feasibility validation:\n');
    fprintf('    %s standard solver status accepted      = %d\n', ...
        funcValidationMark(solver.solved),solver.solved);
    fprintf('    %s recoverable solver candidate         = %d\n', ...
        funcValidationMark(solver.has_candidate),solver.has_candidate);
    fprintf('    %s minimum P-margin eigenvalue          = %.3e\n', ...
        funcValidationMark(validation.P_passed),validation.P_margin);
    fprintf('    %s minimum epsilon_1-margin eigenvalue  = %.3e\n', ...
        funcValidationMark(validation.epsilon1_passed), ...
        validation.epsilon1_margin);
    fprintf('    %s minimum epsilon_2 eigenvalue         = %.3e\n', ...
        funcValidationMark(validation.epsilon2_passed), ...
        validation.epsilon2_margin);
    fprintf('    %s minimum SOS Gram eigenvalue          = %.3e\n', ...
        funcValidationMark(validation.sos_gram_passed), ...
        validation.sos_gram_margin);
    fprintf(['      SOS monomial number / Gram dimension ', ...
             '/ Gram blocks = %g / %g / %g\n'], ...
        validation.sos_monomial_num,validation.sos_gram_dimension, ...
        validation.sos_gram_block_num);
    fprintf('    %s SOS coefficient residual             = %.3e\n', ...
        funcValidationMark(validation.sos_coefficient_passed), ...
        validation.sos_coefficient_residual);
    fprintf('    %s maximum closed-loop real eigenvalue  = %.3e\n', ...
        funcValidationMark(validation.local_stability_passed), ...
        validation.max_closed_loop_real_part);
    fprintf('    %s soft left-half-plane test            = %d\n', ...
        funcValidationMark(validation.soft_local_stability_passed), ...
        validation.soft_local_stability_passed);
    fprintf('    %s soft validation passed               = %d\n', ...
        funcValidationMark(validation.soft_feasible), ...
        validation.soft_feasible);
    fprintf('    %s configured validation passed         = %d\n', ...
        funcValidationMark(validation.passed),validation.passed);
    fprintf('    %s strict validation passed             = %d\n', ...
        funcValidationMark(validation.strict_passed), ...
        validation.strict_passed);

    if ~setting.use_sos_coefficient_as_hard_metric
        fprintf(['    Note: the configured result records but does not use ', ...
                 'the SOS coefficient residual as a hard metric.\n']);
    end
end


function result = funcInitializeRunResult( ...
    system,data,order_K,epsilon_degree,experiment_id)

    result = funcEmptyRunResult();
    result.dim_x = system.dim_x;
    result.order_K = order_K;
    result.epsilon_degree = epsilon_degree;
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
end


function result = funcStoreSolveResult(result,solver,validation,timing)
    result.solver_problem = solver.problem;
    result.solver_info = solver.info;
    result.solver_solved = solver.solved;
    result.solver_candidate = solver.has_candidate;
    result.timed_out = solver.timed_out;
    result.memory_exceeded = solver.memory_exceeded;
    result.memory_usage_ratio_at_stop = ...
        solver.memory_usage_ratio_at_stop;
    result.available_memory_gb_at_stop = ...
        solver.available_memory_gb_at_stop;
    result.attempt_solve_time = timing.solve_time;
    result.attempt_synthesis_time = timing.synthesis_time;
    result.validation_time = timing.validation_time;

    if solver.memory_exceeded
        result.status = 'memory limit';
    elseif solver.timed_out
        result.status = 'timeout';
    elseif solver.has_candidate
        result.P_margin = validation.P_margin;
        result.epsilon1_margin = validation.epsilon1_margin;
        result.epsilon2_margin = validation.epsilon2_margin;
        result.sos_gram_margin = validation.sos_gram_margin;
        result.sos_monomial_num = validation.sos_monomial_num;
        result.sos_gram_dimension = validation.sos_gram_dimension;
        result.sos_gram_block_num = validation.sos_gram_block_num;
        result.sos_coefficient_residual = ...
            validation.sos_coefficient_residual;
        result.max_closed_loop_real_part = ...
            validation.max_closed_loop_real_part;
        result.soft_solver_solved = validation.soft_solver_solved;
        result.soft_local_stability_passed = ...
            validation.soft_local_stability_passed;
        result.soft_feasible = validation.soft_feasible;
        result.candidate_feasible = validation.candidate_passed;
        result.strict_feasible = validation.strict_passed;
        result.configured_feasible = validation.passed;
        % Backward-compatible alias for older result-processing scripts.
        result.feasible = validation.passed;

        if validation.soft_feasible
            result.soft_solve_time = timing.solve_time;
            result.soft_synthesis_time = timing.synthesis_time;
        end
        if validation.passed
            result.status = 'validated feasible';
            result.solve_time = timing.solve_time;
            result.synthesis_time = timing.synthesis_time;
        elseif validation.soft_feasible
            result.status = 'soft feasible only';
        else
            result.status = 'numerical validation failed';
        end
    else
        result.status = 'solver failed';
    end
end


function summary = ...
    funcSummarizeResults(result,dim_x_list,order_K_list,num_Exp)

    summary = repmat(funcEmptySummaryResult(), ...
                     numel(dim_x_list)*numel(order_K_list),1);
    summary_index = 0;

    for idim = 1:numel(dim_x_list)
        for iorder = 1:numel(order_K_list)
            summary_index = summary_index+1;
            dim_x = dim_x_list(idim);
            order_K = order_K_list(iorder);
            select = [result.dim_x] == dim_x & ...
                     [result.order_K] == order_K;
            current = result(select);
            soft_select = [current.soft_feasible];
            solve_time = [current.solve_time];
            synthesis_time = [current.synthesis_time];
            attempt_solve_time = [current.attempt_solve_time];
            attempt_synthesis_time = [current.attempt_synthesis_time];
            data_condition = [current.data_condition];
            sos_monomial_num = [current.sos_monomial_num];
            sos_gram_dimension = [current.sos_gram_dimension];
            sos_gram_block_num = [current.sos_gram_block_num];

            summary(summary_index).dim_x = dim_x;
            summary(summary_index).order_K = order_K;
            summary(summary_index).epsilon_degree = ...
                max(2,2*ceil(order_K/2));
            summary(summary_index).requested_trajectory_num = ...
                current(1).requested_trajectory_num;
            summary(summary_index).trajectory_num = ...
                current(1).trajectory_num;
            summary(summary_index).in_trajectory_samp_num = ...
                current(1).in_trajectory_samp_num;
            summary(summary_index).required_sample_num = ...
                current(1).required_sample_num;
            summary(summary_index).sample_num = current(1).sample_num;
            summary(summary_index).attempted_num = num_Exp;
            summary(summary_index).rank_valid_num = ...
                sum([current.rank_passed]);
            summary(summary_index).data_valid_num = ...
                sum([current.data_valid]);
            summary(summary_index).solver_success_num = ...
                sum([current.solver_solved]);
            summary(summary_index).solver_candidate_num = ...
                sum([current.solver_candidate]);
            summary(summary_index).soft_feasible_num = sum(soft_select);
            summary(summary_index).configured_feasible_num = ...
                sum([current.configured_feasible]);
            % Backward-compatible alias for older result-processing scripts.
            summary(summary_index).feasible_num = ...
                summary(summary_index).configured_feasible_num;
            summary(summary_index).strict_feasible_num = ...
                sum([current.strict_feasible]);
            summary(summary_index).timeout_num = ...
                sum([current.timed_out]);
            summary(summary_index).memory_limit_num = ...
                sum([current.memory_exceeded]);
            summary(summary_index).soft_solvability_rate = ...
                summary(summary_index).soft_feasible_num/num_Exp;
            summary(summary_index).configured_solvability_rate = ...
                summary(summary_index).configured_feasible_num/num_Exp;
            summary(summary_index).solvability_rate = ...
                summary(summary_index).configured_solvability_rate;
            summary(summary_index).strict_solvability_rate = ...
                summary(summary_index).strict_feasible_num/num_Exp;
            summary(summary_index).mean_soft_solve_time = ...
                funcFiniteMean([current.soft_solve_time]);
            summary(summary_index).median_soft_solve_time = ...
                funcFiniteMedian([current.soft_solve_time]);
            summary(summary_index).mean_soft_synthesis_time = ...
                funcFiniteMean([current.soft_synthesis_time]);
            summary(summary_index).median_soft_synthesis_time = ...
                funcFiniteMedian([current.soft_synthesis_time]);
            summary(summary_index).mean_soft_max_real_part = ...
                funcFiniteMean( ...
                    [current(soft_select).max_closed_loop_real_part]);
            summary(summary_index).median_soft_max_real_part = ...
                funcFiniteMedian( ...
                    [current(soft_select).max_closed_loop_real_part]);
            summary(summary_index).mean_solve_time = ...
                funcFiniteMean(solve_time);
            summary(summary_index).median_solve_time = ...
                funcFiniteMedian(solve_time);
            summary(summary_index).mean_synthesis_time = ...
                funcFiniteMean(synthesis_time);
            summary(summary_index).median_synthesis_time = ...
                funcFiniteMedian(synthesis_time);
            summary(summary_index).mean_attempt_solve_time = ...
                funcFiniteMean(attempt_solve_time);
            summary(summary_index).median_attempt_solve_time = ...
                funcFiniteMedian(attempt_solve_time);
            summary(summary_index).mean_attempt_synthesis_time = ...
                funcFiniteMean(attempt_synthesis_time);
            summary(summary_index).median_attempt_synthesis_time = ...
                funcFiniteMedian(attempt_synthesis_time);
            summary(summary_index).mean_data_condition = ...
                funcFiniteMean(data_condition);
            summary(summary_index).median_data_condition = ...
                funcFiniteMedian(data_condition);
            summary(summary_index).mean_sos_monomial_num = ...
                funcFiniteMean(sos_monomial_num);
            summary(summary_index).median_sos_monomial_num = ...
                funcFiniteMedian(sos_monomial_num);
            summary(summary_index).mean_sos_gram_dimension = ...
                funcFiniteMean(sos_gram_dimension);
            summary(summary_index).median_sos_gram_dimension = ...
                funcFiniteMedian(sos_gram_dimension);
            summary(summary_index).mean_sos_gram_block_num = ...
                funcFiniteMean(sos_gram_block_num);
            summary(summary_index).median_sos_gram_block_num = ...
                funcFiniteMedian(sos_gram_block_num);
        end
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
    result.epsilon_degree = NaN;
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
    result.solver_problem = NaN;
    result.solver_info = '';
    result.solver_solved = false;
    result.solver_candidate = false;
    result.timed_out = false;
    result.memory_exceeded = false;
    result.memory_usage_ratio_at_stop = NaN;
    result.available_memory_gb_at_stop = NaN;
    result.status = 'not run';
    result.attempt_solve_time = NaN;
    result.attempt_synthesis_time = NaN;
    result.soft_solve_time = NaN;
    result.soft_synthesis_time = NaN;
    result.solve_time = NaN;
    result.synthesis_time = NaN;
    result.validation_time = NaN;
    result.P_margin = NaN;
    result.epsilon1_margin = NaN;
    result.epsilon2_margin = NaN;
    result.sos_gram_margin = NaN;
    result.sos_monomial_num = NaN;
    result.sos_gram_dimension = NaN;
    result.sos_gram_block_num = NaN;
    result.sos_coefficient_residual = NaN;
    result.max_closed_loop_real_part = NaN;
    result.soft_solver_solved = false;
    result.soft_local_stability_passed = false;
    result.soft_feasible = false;
    result.candidate_feasible = false;
    result.strict_feasible = false;
    result.configured_feasible = false;
    result.feasible = false;
end


function summary = funcEmptySummaryResult()
    summary.dim_x = NaN;
    summary.order_K = NaN;
    summary.epsilon_degree = NaN;
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
    summary.feasible_num = NaN;
    summary.strict_feasible_num = NaN;
    summary.timeout_num = NaN;
    summary.memory_limit_num = NaN;
    summary.soft_solvability_rate = NaN;
    summary.configured_solvability_rate = NaN;
    summary.solvability_rate = NaN;
    summary.strict_solvability_rate = NaN;
    summary.mean_soft_solve_time = NaN;
    summary.median_soft_solve_time = NaN;
    summary.mean_soft_synthesis_time = NaN;
    summary.median_soft_synthesis_time = NaN;
    summary.mean_soft_max_real_part = NaN;
    summary.median_soft_max_real_part = NaN;
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
    summary.mean_sos_monomial_num = NaN;
    summary.median_sos_monomial_num = NaN;
    summary.mean_sos_gram_dimension = NaN;
    summary.median_sos_gram_dimension = NaN;
    summary.mean_sos_gram_block_num = NaN;
    summary.median_sos_gram_block_num = NaN;
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
end


function [controller,solver,validation,timing] = ...
    funcFailedIsolatedSolve( ...
        problem,reason,elapsed_time,memory_exceeded, ...
        memory_usage_ratio,available_memory_gb)

    controller = [];
    solver = funcEmptySolverResult();
    validation = funcEmptyValidation();
    timing.solve_time = elapsed_time;
    timing.synthesis_time = elapsed_time;
    timing.validation_time = NaN;
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


function validation = funcEmptyValidation()
    validation.P_margin = NaN;
    validation.epsilon1_margin = NaN;
    validation.epsilon2_margin = NaN;
    validation.sos_gram_margin = NaN;
    validation.sos_monomial_num = NaN;
    validation.sos_gram_dimension = NaN;
    validation.sos_gram_block_num = NaN;
    validation.sos_coefficient_residual = NaN;
    validation.max_closed_loop_real_part = NaN;
    validation.closed_loop_eigenvalue = [];
    validation.soft_solver_solved = false;
    validation.soft_local_stability_passed = false;
    validation.soft_feasible = false;
    validation.P_passed = false;
    validation.epsilon1_passed = false;
    validation.epsilon2_passed = false;
    validation.sos_gram_passed = false;
    validation.sos_coefficient_passed = false;
    validation.local_stability_passed = false;
    validation.candidate_certificate_passed = false;
    validation.strict_certificate_passed = false;
    validation.candidate_passed = false;
    validation.strict_passed = false;
    validation.passed = false;
end


function mark = funcValidationMark(is_passed)
    if is_passed
        mark = '!';
    else
        mark = '?';
    end
end


function [sos_polynomial,gram_matrix] = ...
    funcSosPolynomial(x,degree_num,include_constant)

    if degree_num < 2 || mod(degree_num,2) ~= 0
        error(['The maximum SOS polynomial degree must be a positive ', ...
               'even integer.']);
    end

    monomial = monolist(x,degree_num/2);
    if ~include_constant
        monomial = monomial(2:end);
    end

    gram_matrix = ...
        sdpvar(length(monomial),length(monomial),'symmetric');
    sos_polynomial = monomial'*gram_matrix*monomial;
end


function [poly_matrix,coefficient,monomial_power] = ...
    funcPolynomialMatrix(x,row_num,column_num,degree_num)

    monomial = monolist(x,degree_num);
    coefficient = ...
        sdpvar(row_num*column_num,length(monomial),'full');
    poly_matrix = reshape(coefficient*monomial,row_num,column_num);

    monomial_power = zeros(length(monomial),length(x));
    for im = 1:length(monomial)
        for ix = 1:length(x)
            monomial_power(im,ix) = degree(monomial(im),x(ix));
        end
    end
end


function parameter = funcDecisionVariables(candidate,x)
    parameter_id = setdiff(depends(candidate),depends(x));
    parameter = recover(parameter_id);
end


function option = funcSosOptions(setting)
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
            setting.solve_timeout, ...
        'sos.model',1, ...
        'sos.newton',1, ...
        'sos.congruence',1, ...
        'sos.scale',1, ...
        'sos.clean',1e-12);
end


function closed_loop_jacobian = ...
    funcClosedLoopJacobian(controller,system,step)

    dim_x = system.dim_x;
    state_origin = zeros(dim_x,1);
    closed_loop_jacobian = zeros(dim_x,dim_x);

    for ix = 1:dim_x
        state_step = zeros(dim_x,1);
        state_step(ix) = step;
        derivative_plus = funcClosedLoopSystem( ...
            state_origin+state_step,controller,system);
        derivative_minus = funcClosedLoopSystem( ...
            state_origin-state_step,controller,system);
        closed_loop_jacobian(:,ix) = ...
            (derivative_plus-derivative_minus)/(2*step);
    end
end


function dx = funcClosedLoopSystem(x,controller,system)
    input = funcControllerValue(x,controller);
    dx = funcPolynomialSystem(x,input,system);
end


function input = funcControllerValue(x,controller)
    monomial = prod(x(:)'.^controller.monomial_power,2);
    gain = (controller.gain_coefficient*monomial)';
    input = gain*(controller.P\funcZhat(x));
end


function Z = funcPowerVector(x)
    Z = [x;x.^3];
end


function Zhat = funcZhat(x)
    Zhat = x(:);
end


function Ubar0 = funcLiftInputData(~,U0)
    % The benchmark uses W(x) = 1, hence Ubar0 = U0.
    Ubar0 = U0;
end


function W = funcWmatrix(~)
    W = 1;
end


function H = funcHmatrix(x)
    dim_x = length(x);
    H = [eye(dim_x);diag(x.^2)];
end


function funcSetupYalmip()
    if exist('sdpvar','file') ~= 2
        yalmip_path = '/home/lidong/Documents/MATLAB/YALMIP-master';
        if ~isfolder(yalmip_path)
            error('YALMIP was not found. Set yalmip_path in funcSetupYalmip().');
        end
        addpath(genpath(yalmip_path));
    end

    if exist('mosekopt','file') ~= 2
        mosek_path = ['/home/lidong/Documents/MATLAB/mosektoolslinux64x86/', ...
                      'mosek/11.0/toolbox/r2019b'];
        if isfolder(mosek_path)
            addpath(mosek_path);
        end
    end
end
