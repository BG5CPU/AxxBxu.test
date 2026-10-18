clear; close all; clc;


%% experiment settings ====================================================
dim_x = 2;
dim_u = 1;
ph = 0.1;                         % discrete time step

rxtt = linspace(0.5, 1.3, 21);    % tested values of r
r0tt = nan(size(rxtt));
vKtt = cell(size(rxtt));
cvxStatus = cell(size(rxtt));


%% load one fixed data set =================================================
% Load once before the r sweep. Every controller synthesis below uses the
% same dataA, dataB, and dataC from this saved experiment.
savedData = load("Mdata/data_stabilization_09.mat");
requiredFields = {'dataA', 'dataB', 'dataC'};
assert(all(isfield(savedData, requiredFields)), ...
    'The saved data set must contain dataA, dataB, and dataC.');

fixedData.dataA = savedData.dataA;
fixedData.dataB = savedData.dataB;
fixedData.dataC = savedData.dataC;

fprintf('Minimum eigenvalue of dataA: %.6g\n', ...
    min(eig((fixedData.dataA + fixedData.dataA') / 2)));


%% solve one data-based controller for every r ============================
for ir = 1:numel(rxtt)
    rx = rxtt(ir);

    [vKtt{ir}, r0tt(ir), cvxStatus{ir}] = ...
        solveDataControllerAndR0(rx, fixedData, dim_x, dim_u);

    fprintf('r = %.4f, r0 = %.6f, CVX status = %s\n', ...
        rx, r0tt(ir), cvxStatus{ir});
end


%% collect and plot the results ===========================================
experimentResults = table(rxtt(:), r0tt(:), cvxStatus(:), ...
    'VariableNames', {'r', 'r0', 'cvx_status'});
disp(experimentResults);

plotRAndR0(rxtt, r0tt);


%% useful functions =======================================================
function [vK, r0, solverStatus] = solveDataControllerAndR0( ...
        rx, fixedData, dim_x, dim_u)
    % Only barQ depends on rx; fixedData is unchanged between solves.
    barQ = buildBarQ(rx);
    dataA = fixedData.dataA;
    dataB = fixedData.dataB;
    dataC = fixedData.dataC;

    cvx_begin sdp quiet
        variable epGa
        variable lyGa(dim_x, dim_x) symmetric
        variable kY(dim_u, dim_x)
        variable mu_l nonnegative
        variable mu_u nonnegative

        minimize(mu_u - mu_l - 0.1*epGa)

        subject to
            for iv = 1:numel(barQ)
                blockS = [ -lyGa-dataC,   zeros(dim_x),          dataB;
                           zeros(dim_x), -lyGa+epGa*eye(dim_x), -[lyGa; kY]'*barQ{iv}';
                           dataB',       -barQ{iv}*[lyGa; kY],   -dataA ];
                blockS <= 0;
            end

            epGa >= 0.00001;
            lyGa >= 0.0001 * eye(dim_x);
            lyGa >= mu_l * eye(dim_x);
            lyGa <= mu_u * eye(dim_x);
    cvx_end

    solverStatus = cvx_status;
    if contains(solverStatus, 'Solved')
        vK = kY / lyGa;
        r0 = calculateR0(rx, lyGa, epGa);
    else
        vK = nan(dim_u, dim_x);
        r0 = NaN;
        warning('No data-based controller was found for r = %.4f (%s).', ...
            rx, solverStatus);
    end
end


function barQ = buildBarQ(rx)
    % Bounds for the seven independently varying basis-function entries.
    entryBounds = {
        [sin(rx)/rx, 1], ...          % sin(x1)/x1
        [-rx, rx], ...                % x1
        [0, rx^2], ...                % x1^2
        [0, rx], ...                  % |x1|
        [0, rx], ...                  % |x2|
        [exp(-rx), exp(rx)], ...      % exp(x1)
        [exp(-rx), exp(rx)] ...       % exp(x2)
    };

    % Seven binary choices give 2^7 = 128 polytope vertices.
    vertexBits = dec2bin(0:127, 7) - '0';
    barQ = cell(size(vertexBits, 1), 1);

    for iv = 1:size(vertexBits, 1)
        entryValue = zeros(1, 7);
        for ib = 1:7
            entryValue(ib) = entryBounds{ib}(vertexBits(iv, ib) + 1);
        end

        xiA1 = [1; entryValue(1)];
        xiA2 = [1; entryValue(2); entryValue(3)];
        xiB1 = [1; entryValue(4); entryValue(5); ...
                   entryValue(6); entryValue(7)];
        barQ{iv} = blkdiag(xiA1, xiA2, xiB1);
    end
end


function r0 = calculateR0(rx, lyGa, epGa)
    lyGa = (lyGa + lyGa') / 2;
    eigLyGa = eig(lyGa);
    eigMax = max(eigLyGa);
    eigMin = min(eigLyGa);

    numerator = eigMax * eigMin;
    denominator = eigMax^2 - epGa*eigMin;
    numericalTolerance = 1e-10 * max(1, eigMax^2);

    if denominator > numericalTolerance
        r0Candidate = sqrt(max(0, numerator / denominator));
    elseif denominator >= -numericalTolerance
        r0Candidate = Inf;
    else
        warning('The r0 expression has a negative denominator for r = %.4f.', rx);
        r0 = NaN;
        return;
    end

    r0 = min(1, r0Candidate) * rx;
end


function plotRAndR0(rxtt, r0tt)
    validResult = isfinite(r0tt);

    figure('color', [1 1 1]);
    set(gcf, 'Position', [100, 100, 400, 380], 'Renderer', 'Painters');

    scatter(rxtt(validResult), r0tt(validResult), 40, 'k', 'filled');
    axis equal;

    xlabel('$r$', 'Fontname', 'Times New Roman', ...
        'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$\bar{r}_0$', 'Fontname', 'Times New Roman', ...
        'Interpreter', 'latex', 'FontSize', 12);
    set(gca, 'fontsize', 12, 'fontname', 'Times');
end
