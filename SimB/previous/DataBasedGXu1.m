close all; 
clear; clc;

dimN = 2;
dimM = 1;



%% system settings ========================================================
% dt = 0.01;

Ad = [ 0.8, 0.5;
       0.4, 1.2 ];

% qx = [1, x2^2, exp(x1)]';
Ed = [ 1, -0.45,  0.45;
       2, -0.30, -0.30 ];


dimQ = 6;
% qx_library = [1, x1^2, x1*x2, x2^2, exp(x1), exp(x2)]';
ED = [1, 0, 0, -0.45,  0.45, 0;
      2, 0, 0, -0.30, -0.30, 0 ];




%% data collection ========================================================
NT =  4;
trajt   =  10;

engyDis = 0.01;
engyInp = 1;

QU = []; X0 = [];
X1 = []; D0 = [];

for it = 1:trajt
    uu = (rand(dimM, NT)-0.5)*2* engyInp;
    xx = zeros(dimN, NT+1); xx(:,1) = randn(dimN,1) *0.001;
    dd = (rand(dimN, NT)-0.5)*2* engyDis;    
    for ii = 1:NT
        xx(:,ii+1) =    Ad*xx(:,ii) ...
                      + Ed*[1; xx(2,ii)^2; exp(xx(1,ii))]*uu(:,ii) ...
                      + dd(:,ii);
    end
    % creat data matrices
    QU_ = [ ones(1,NT)             .*uu;
            xx(1,1:NT).^2          .*uu;
            xx(1,1:NT).*xx(2,1:NT) .*uu;
            xx(2,1:NT).^2          .*uu;
            exp(xx(1,1:NT))        .*uu;
            exp(xx(2,1:NT))        .*uu ];
    X0_ = xx(:, 1:NT);
    X1_ = xx(:, 2:NT+1);
    D0_ = dd;
   
    X0  = [X0, X0_];
    X1  = [X1, X1_];
    QU  = [QU, QU_];
    D0  = [D0, D0_];
end

plot(X1(1,:));

dataTest = Ad*X0+ED*QU+D0-X1;

de =  ceil( max(svd(D0*D0'))* 1e4 )* 1e-4;
De =  de * eye(dimN);

W0 = [X0; QU];

% data consistent set
Adata =  W0*W0';
Bdata = -X1*W0';
Cdata =  X1*X1'-De;

eig_W0  = eig(W0*W0');

disp( "min_eig_Adata = " + min(eig_W0) );
disp( "max_eig_Adata = " + max(eig_W0) );
disp( "noise energy = " + max(svd(D0*D0')) );
disp( "energy bound = " + de );


% load("data\data_based_data.mat");






%% data based controllers =================================================
% set the domain: |x| <= rx
rx = 0.8;
rho11   = 1;                           % 1
rho12_u = rx^2;    rho12_l = 0;        % x1^2
rho13_u = rx^2;    rho13_l = -rx^2;    % x1*x2
rho14_u = rx^2;    rho14_l = 0;        % x2^2
rho15_u = exp(rx); rho15_l = exp(-rx); % exp(x1)
rho16_u = exp(rx); rho16_l = exp(-rx); % exp(x2)

Vrho = { [rho11, rho12_u, rho13_u, rho14_u, rho15_u, rho16_u]';
         [rho11, rho12_u, rho13_u, rho14_u, rho15_u, rho16_l]'; 
         [rho11, rho12_u, rho13_u, rho14_u, rho15_l, rho16_u]'; 
         [rho11, rho12_u, rho13_u, rho14_u, rho15_l, rho16_l]';
         [rho11, rho12_u, rho13_u, rho14_l, rho15_u, rho16_u]';
         [rho11, rho12_u, rho13_u, rho14_l, rho15_u, rho16_l]'; 
         [rho11, rho12_u, rho13_u, rho14_l, rho15_l, rho16_u]'; 
         [rho11, rho12_u, rho13_u, rho14_l, rho15_l, rho16_l]';
         [rho11, rho12_u, rho13_l, rho14_u, rho15_u, rho16_u]';
         [rho11, rho12_u, rho13_l, rho14_u, rho15_u, rho16_l]'; 
         [rho11, rho12_u, rho13_l, rho14_u, rho15_l, rho16_u]'; 
         [rho11, rho12_u, rho13_l, rho14_u, rho15_l, rho16_l]';
         [rho11, rho12_u, rho13_l, rho14_l, rho15_u, rho16_u]';
         [rho11, rho12_u, rho13_l, rho14_l, rho15_u, rho16_l]'; 
         [rho11, rho12_u, rho13_l, rho14_l, rho15_l, rho16_u]'; 
         [rho11, rho12_u, rho13_l, rho14_l, rho15_l, rho16_l]';
         [rho11, rho12_l, rho13_u, rho14_u, rho15_u, rho16_u]';
         [rho11, rho12_l, rho13_u, rho14_u, rho15_u, rho16_l]'; 
         [rho11, rho12_l, rho13_u, rho14_u, rho15_l, rho16_u]'; 
         [rho11, rho12_l, rho13_u, rho14_u, rho15_l, rho16_l]';
         [rho11, rho12_l, rho13_u, rho14_l, rho15_u, rho16_u]';
         [rho11, rho12_l, rho13_u, rho14_l, rho15_u, rho16_l]'; 
         [rho11, rho12_l, rho13_u, rho14_l, rho15_l, rho16_u]'; 
         [rho11, rho12_l, rho13_u, rho14_l, rho15_l, rho16_l]';
         [rho11, rho12_l, rho13_l, rho14_u, rho15_u, rho16_u]';
         [rho11, rho12_l, rho13_l, rho14_u, rho15_u, rho16_l]'; 
         [rho11, rho12_l, rho13_l, rho14_u, rho15_l, rho16_u]'; 
         [rho11, rho12_l, rho13_l, rho14_u, rho15_l, rho16_l]';
         [rho11, rho12_l, rho13_l, rho14_l, rho15_u, rho16_u]';
         [rho11, rho12_l, rho13_l, rho14_l, rho15_u, rho16_l]'; 
         [rho11, rho12_l, rho13_l, rho14_l, rho15_l, rho16_u]'; 
         [rho11, rho12_l, rho13_l, rho14_l, rho15_l, rho16_l]' };


cvx_begin sdp % quiet
    variable epGa semidefinite;
    variable lyGa(dimN,dimN) semidefinite;
    variable kY(dimM,dimN);

    minimize( lambda_max(lyGa) );

    subject to   
        for iv = 1:length(Vrho)
            blockS = [ -lyGa-Cdata,   zeros(dimN),          Bdata;
                        zeros(dimN), -lyGa+epGa*eye(dimN), -[lyGa; Vrho{iv}*kY]';    
                        Bdata',      -[lyGa; Vrho{iv}*kY], -Adata ];
            blockS <= 0;
        end
        epGa >= 0.001;
        lyGa >= 0.001 * eye(dimN);
cvx_end

vK = kY/lyGa;

% region of attraction
roac = sqrt( maxEig(lyGa)*minEig(lyGa) / ( maxEig(lyGa)^2 - epGa*minEig(lyGa) ) );
roa  = min( 1, roac ) * rx;
disp("radius of ROA is " + roa);










%% controller implementation ==============================================
ed = 0.01;
tn = 100; timeS = 0:1:tn;
xt = zeros(dimN,tn+1); 
xt(:,1) = ones(dimN,1) *0.1;
disp("initial value is " + xt(1,1) + ", " + xt(2,1));
for ik  = 1:tn
    dt = (rand(dimN,1)-0.5)*2*ed;
    ut = vK*xt(:,ik);
    xt(:,ik+1) = Ad*xt(:,ik) + Ed*[1; xt(2,ik)^2; exp(xt(1,ik))]*ut + dt;
end

figure('color', [1 1 1]);
set(gcf, 'Position', [100, 100, 500, 200]);
subplot(2,1,1);
plot(timeS, xt(1,:), 'K', 'LineWidth', 1); grid on;
ylabel('$x_1$', 'Fontname', 'Times New Roman','Interpreter', 'latex','FontSize',12);
set(gca,'fontsize',12,'fontname','Times');
set(gca,'position',[0.11 0.65 0.86 0.30]);
subplot(2,1,2);
plot(timeS, xt(2,:), 'K', 'LineWidth', 1); grid on;
set(gca,'fontsize',12,'fontname','Times');
xlabel('steps', 'Fontname', 'Times New Roman','FontSize',12);
ylabel('$x_2$', 'Fontname', 'Times New Roman','Interpreter', 'latex','FontSize',12);
set(gca,'position',[0.11 0.22 0.86 0.30]);






%% phase portrait =========================================================
xrr = 2.5; dx = xrr/10;
xbis = 0.8; ybis = -0.1; 
xleft = -xrr-xbis; xright = xrr-xbis;
yleft = -xrr-ybis; yright = xrr-ybis;
[x1, x2] = meshgrid(xleft:dx:xright, yleft:dx:yright);
Ux = 0.8*x1 + 0.5*x2 + (1-0.45*x2.^2+0.45*exp(x1)).*(vK(1)*x1+vK(2)*x2) - x1;
Vx = 0.4*x1 + 1.2*x2 + (2-0.30*x2.^2-0.30*exp(x1)).*(vK(1)*x1+vK(2)*x2) - x2;
figure('color', [1 1 1]);
set(gcf, 'Position', [100, 100, 400, 380]);
quiver(x1, x2 ,Ux ,Vx, 3); hold on;
xlabel('$x_1$', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x_2$', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 12);
axis equal;
axis([xleft xright yleft yright]);
set(gca,'fontsize',12,'fontname','Times');
fill(cos(0:0.01*pi:2*pi)*roa, sin(0:0.01*pi:2*pi)*roa, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
legend('$[x_1^{+}\!\!-\!x_1 \ \ x_2^{+}\!\!-\!x_2]$', ...
       'region of attraction',...
       'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 12);








%% useful fuction =========================================================
function e_max = maxEig( MatrixSym )
    e_max = max(eig(MatrixSym));
end

function e_min = minEig( MatrixSym )
    e_min = min(eig(MatrixSym));
end

















