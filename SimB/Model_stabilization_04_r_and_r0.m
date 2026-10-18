clear; close all; clc;


%% experiment settings ====================================================
dim_x = 2;
dim_u = 1;
ph = 0.1;                         % discrete time step

rxtt = linspace(0.5, 1.5, 21);    % tested values of rx
r0tt = nan(size(rxtt));
vKtt = cell(size(rxtt));
cvxStatus = cell(size(rxtt));


%% solve one controller for every rx ======================================
for ir = 1:numel(rxtt)
    rx = rxtt(ir);

    [vKtt{ir}, r0tt(ir), cvxStatus{ir}] = ...
        solveControllerAndR0(rx, ph, dim_x, dim_u);

    fprintf('rx = %.4f, r0 = %.6f, CVX status = %s\n', ...
        rx, r0tt(ir), cvxStatus{ir});
end


%% collect and plot the results ===========================================
experimentResults = table(rxtt(:), r0tt(:), cvxStatus(:), ...
    'VariableNames', {'rx', 'r0', 'cvx_status'});
disp(experimentResults);

plotRxAndR0(rxtt, r0tt);


%% useful functions =======================================================
function [vK, r0, solverStatus] = solveControllerAndR0(rx, ph, dim_x, dim_u)
    % Construct the uncertain-system vertices for the current rx.
    barG = buildBarG(rx, ph);

    % Solve the controller synthesis problem for the current rx.
    cvx_begin sdp quiet
        variable epGa
        variable lyGa(dim_x, dim_x) symmetric
        variable kY(dim_u, dim_x)
        variables mu_l mu_u

        minimize(mu_u - mu_l)

        subject to
            for iv = 1:numel(barG)
                blockS = [ -lyGa,                  barG{iv}*[lyGa; kY];
                           [lyGa; kY]'*barG{iv}', -lyGa + epGa*eye(dim_x) ];
                blockS <= 0;
            end

            epGa >= 0.001;
            lyGa >= 0.01 * eye(dim_x);
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
        warning('No controller was found for rx = %.4f (%s).', ...
            rx, solverStatus);
    end
end


function barG = buildBarG(rx, ph)
    % Bounds of the entries in A(x) and B(x) over |x| <= rx.
    a11Bounds = [1 + ph*sin(rx)/rx, 1 + ph];
    a22Bounds = [1 - ph, 1 - ph + ph*rx^2];
    b11Bounds = [ph, (1 + rx)*ph];
    b12Bounds = [ph*exp(-rx), ph*exp(rx)];

    a12 = 2*ph;
    a21 = 2*ph;

    % Four independently bounded entries give 2^4 = 16 vertices.
    vertexBits = dec2bin(0:15, 4) - '0';
    barG = cell(size(vertexBits, 1), 1);

    for iv = 1:size(vertexBits, 1)
        a11 = a11Bounds(vertexBits(iv, 1) + 1);
        a22 = a22Bounds(vertexBits(iv, 2) + 1);
        b11 = b11Bounds(vertexBits(iv, 3) + 1);
        b12 = b12Bounds(vertexBits(iv, 4) + 1);

        barG{iv} = [a11, a12, b11;
                    a21, a22, b12];
    end
end


function r0 = calculateR0(rx, lyGa, epGa)
    % Symmetrization suppresses tiny asymmetry introduced by the solver.
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
        warning('The r0 expression has a negative denominator for rx = %.4f.', rx);
        r0 = NaN;
        return;
    end

    r0 = min(1, r0Candidate) * rx;
end


function plotRxAndR0(rxtt, r0tt)
    validResult = isfinite(r0tt);

    figure('color', [1 1 1]);
    set(gcf, 'Position', [100, 100, 400, 380], 'Renderer', 'Painters');

    scatter(rxtt(validResult), r0tt(validResult), 40, 'k', 'filled');
    axis equal;

    xlabel('$r$', 'Fontname', 'Times New Roman', ...
        'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$r_0$', 'Fontname', 'Times New Roman', ...
        'Interpreter', 'latex', 'FontSize', 12);
    set(gca, 'fontsize', 12, 'fontname', 'Times');
end
