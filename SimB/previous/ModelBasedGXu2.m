close all; 
clear; clc;

dimN = 2;
dimM = 1;

%% system settings ========================================================
% dt = 0.01;

Ad = [ 0.8, 0.5;
       0.4, 1.2 ];
Ed = [ 1, -0.45,  0.45;
       2, -0.30, -0.30 ];
% qx = [1, x2^2, exp(x1)]';



%% model based controllers ================================================
% set the domain: |x| <= rx
rx = 0.96;                         
rho11_u = 1 - 0.45*0^2 + 0.45*exp(rx);  

% rho11_l = 1 - 0.45*0^2 + 0.45*exp(-rx);
x1so = -lambertw(0, 1/2);
x2sq = rx^2 - x1so^2;
rho11_l = 1 - 0.45*x2sq + 0.45*exp(x1so);

rho12_u = 2 - 0.30*0^2 - 0.30*exp(-rx); 
rho12_l = 2 - 0.30*0^2 - 0.30*exp(rx); 

Vrho = { [rho11_u, rho12_u]';
         [rho11_u, rho12_l]'; 
         [rho11_l, rho12_u]'; 
         [rho11_l, rho12_l]' };



cvx_begin sdp % quiet
    variable epGa semidefinite;
    variable lyGa(dimN,dimN) semidefinite;
    variable kY(dimM,dimN);

    % maximize( lambda_min(lyGa) - lambda_max(lyGa) + epGa );
    minimize( lambda_max(lyGa) );

    subject to   
        for iv = 1:length(Vrho)
            blockS = [ -lyGa,                       Ad*lyGa + Vrho{iv}*kY;
                        lyGa*Ad' + kY'*Vrho{iv}',  -lyGa + epGa*eye(dimN) ];
            blockS <= 0;
        end
        epGa >= 1;
        lyGa >= 1 * eye(dimN);
cvx_end

vK = kY/lyGa;


% region of attraction
roac = sqrt( maxEig(lyGa)*minEig(lyGa) / ( maxEig(lyGa)^2 - epGa*minEig(lyGa) ) );
roa  = min( 1, roac ) * rx;
disp("radius of ROA is " + roa);


% load("data\model_based_data.mat");


%% controller implementation ==============================================
ed = 0.1;
tn = 100; timeS = 0:1:tn;
xt = zeros(dimN,tn+1); 
xt(:,1) = ones(dimN,1) *0.5;
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
set(gca,'position',[0.095 0.65 0.88 0.30]);
subplot(2,1,2);
plot(timeS, xt(2,:), 'K', 'LineWidth', 1); grid on;
set(gca,'fontsize',12,'fontname','Times');
xlabel('steps', 'Fontname', 'Times New Roman','FontSize',12);
ylabel('$x_2$', 'Fontname', 'Times New Roman','Interpreter', 'latex','FontSize',12);
set(gca,'position',[0.095 0.22 0.88 0.30]);







%% phase portrait =========================================================
xrr = 2.5; dx = xrr/10;
xbis = 0.8; ybis = -0.2; 
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

































