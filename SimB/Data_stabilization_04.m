clear; close all; clc;

dim_x = 2;
dim_u = 1;

ph = 0.1; % discrete time



%% original system generate data ==========================================

% xiA1 = [ 1; sin(x1)/x1 ];
EA1 = [1, ph; 2*ph, 0];
dim_xiA1 = 2;
% x1A2 = [1; x1; x1^2];
EA2 = [2*ph, 0, 0; 1-ph, 0, ph];
dim_xiA2 = 3;
dim_xiA = dim_xiA1+dim_xiA2;
% xiB1 = [1; abs(x1); abs(x2); exp(x1); exp(x2)];
EB1 = [ph, 0, ph, 0, 0; 0, 0, 0, ph, 0];
dim_xiB1 = 5;
dim_xiB = dim_xiB1;



trajt = 10;
n_s = 13;

euin = 1.3;  % bound of input
edis = 0.0064;  % bound of output noise

U0 = []; X0 = [];
X1 = []; D0 = [];
Dstate0 = zeros(dim_x, n_s, trajt); 
Dstate1 = Dstate0; 
Dinput0 = zeros(dim_u, n_s, trajt);
Ddistr0 = Dstate0;

for ij = 1:trajt
    us = (rand(dim_u, n_s)-0.5)*2* euin;
    xs = zeros(dim_x, n_s+1); xs(:,1) = 0.05*randn(2,1);
    ds = (rand(dim_x, n_s)-0.5)*2* edis;
    
    XIA = zeros(dim_xiA, n_s);
    XIB = zeros(dim_xiB, n_s);
    for it = 1:n_s
        Ax = [ 1+ph*sin(xs(1,it))/xs(1,it), ph*2; ph*2, 1-ph+ph*xs(1,it)^2];
        Bx = [ ph*(1+abs(xs(2,it))); ph*exp(xs(1,it)) ];
        xs(:,it+1) = Ax*xs(:,it) + Bx*us(:,it) + ds(:,it);
        
        x1t = xs(1,it); x2t = xs(2,it);
        xiA1 = [ 1; sin(x1t)/x1t ]; xiA2 = [ 1; x1t; x1t^2 ];
        XIA(:,it) = blkdiag(xiA1, xiA2)*xs(:,it);
        
        xiB1 = [1; abs(x1t); abs(x2t); exp(x1t); exp(x2t)];
        XIB(:,it) = xiB1*us(:,it);
    end
    X0 = [X0, XIA];
    U0 = [U0, XIB];
    X1 = [X1, xs(:,2:n_s+1)];
    D0 = [D0, ds(:,1:n_s)];
    Dstate0(:,:,ij) = xs(:,1:n_s);
    Dstate1(:,:,ij) = xs(:,2:n_s+1);
    Dinput0(:,:,ij) = us;
    Ddistr0(:,:,ij) = ds;
end

DATATEST = [EA1, EA2]*X0 +EB1*U0 + D0 - X1;
figure(1)
plot(xs(2,:));
% figure(2)
% plot(DATATEST(2,:));


dataA = [X0; U0]*[X0; U0]';
dataB = -X1*[X0; U0]';
% dataC = X1*X1'- edis^2*dim_x*length(D0);

THETAnoise = roundUpSignificant(max(svd(D0*D0')), 2);
disp("THETAnoise eig is " + THETAnoise);

dataC = X1*X1'- THETAnoise*eye(dim_x);

disp("min eig is " + min(eig(dataA)));





%% data from previous experiment ==========================================
load("Mdata/data_stabilization_09.mat");




%% data based controllers =================================================
close;

% set the domain: |x| <= rx
rx = 0.92; 

% XIA
xiA11 = 1;
xiA21_u = 1; xiA21_l = sin(rx)/rx;
xiA12 = 1;
xiA22_u = rx; xiA22_l = -rx;
xiA32_u = rx^2; xiA32_l = 0;
% XIB
xiB11 = 1;
xiB21_u = rx; xiB21_l = 0;
xiB31_u = rx; xiB31_l = 0;
xiB41_u = exp(rx); xiB41_l = exp(-rx);
xiB51_u = exp(rx); xiB51_l = exp(-rx);

% the set of barQ
barQ = funcBarQ(xiA11, xiA21_u,... 
                xiA21_l, xiA12, xiA22_u, xiA22_l, xiA32_u, xiA32_l,... 
                xiB11, xiB21_u, xiB21_l, xiB31_u, xiB31_l, xiB41_u, xiB41_l, xiB51_u, xiB51_l);




% solve sdp LMI
cvx_begin sdp % quiet
    variable epGa semidefinite;
    variable lyGa(dim_x,dim_x) semidefinite;
    variable kY(dim_u,dim_x);
    % variable mu_l semidefinite;
    % variable mu_u semidefinite;
    
    % minimize( mu_u - 0*mu_l + 1*epGa );
    minimize( lambda_max(lyGa) - lambda_min(lyGa) - 0.1*epGa - 0.0*trace(lyGa) );

    subject to   
        for iv = 1:length(barQ)
            blockS = [ -lyGa-dataC,   zeros(dim_x),          dataB;
                       zeros(dim_x), -lyGa+epGa*eye(dim_x), -[lyGa;kY]'*barQ{iv}'; 
                       dataB',       -barQ{iv}*[lyGa;kY],   -dataA ];
            blockS <= 0;
        end
        epGa >= 0.00001;
        lyGa >= 0.0001 * eye(dim_x);
        % lyGa >= mu_l * eye(dim_x);
        % lyGa <= (mu_l+mu_u) * eye(dim_x);
        
cvx_end

vK = kY/lyGa;

% region of attraction
roac = sqrt( maxEig(lyGa)*minEig(lyGa) / ( maxEig(lyGa)^2 - epGa*minEig(lyGa) ) );
roa  = min( 1, roac ) * rx;
disp("radius of ROA is " + roa);






%% phase portrait =========================================================
xrr = 2; dx = xrr/10;
xbis = 0.8; ybis = -0.8; 
xleft = -xrr-xbis; xright = xrr-xbis;
yleft = -xrr-ybis; yright = xrr-ybis;
[x1, x2] = meshgrid(xleft:dx:xright, yleft:dx:yright);

Px1 = x1 + ph*sin(x1) + ph*2*x2 + ph*(1+abs(x2)).*(vK(1)*x1+vK(2)*x2) - x1;
Px2 = ph*2*x1 + (1-ph+ph*x1.^2).*x2 + ph*exp(x1).*(vK(1)*x1+vK(2)*x2) - x2;

figure('color', [1 1 1]);
set(gcf, 'Position', [100, 100, 400, 380], 'Renderer', 'Painters');

fill(cos(0:0.01*pi:2*pi)*roa, sin(0:0.01*pi:2*pi)*roa, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); 

hold on;

quiver(x1, x2 ,Px1 ,Px2, 4); 

xlabel('$x_1$', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x_2$', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 12);
axis equal;
axis([xleft xright yleft yright]);
set(gca,'fontsize',12,'fontname','Times');

legend('$\mathcal{B}_{\bar{r}_0}$', ...
       '$(x_1^{+}\!\!-\!x_1 , x_2^{+}\!\!-\!x_2)$', ...
       'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 12);




   

%% controller implementation ==============================================
ed = 0.1;
tn = 300; timeS = 0:1:tn;
xt = zeros(dim_x,tn+1); 
xt(:,1) = [-0.4, -0.4];
% disp("initial value is " + xt(1,1) + ", " + xt(2,1));
for ik  = 1:tn
    dt = (rand(dim_x,1)-0.5)*2*ed;
     
    Ax = [ 1+ph*sin(xt(1,ik))/xt(1,ik), ph*2; ph*2, 1-ph+ph*xt(1,ik)^2];
    Bx = [ ph*(1+abs(xt(2,ik))); ph*exp(xt(1,ik)) ];
    
    ut = vK*xt(:,ik);
    xt(:,ik+1) = Ax*xt(:,ik) + Bx*ut + dt;
end

figure('color', [1 1 1]);
set(gcf, 'Position', [100, 100, 500, 200]);
subplot(2,1,1);
plot(timeS, xt(1,:), 'K', 'LineWidth', 1); grid on;
ylabel('$x_1$', 'Fontname', 'Times New Roman','Interpreter', 'latex','FontSize',12);
set(gca,'fontsize',12,'fontname','Times');
set(gca,'position',[0.11 0.67 0.86 0.30]);
subplot(2,1,2);
plot(timeS, xt(2,:), 'K', 'LineWidth', 1); grid on;
set(gca,'fontsize',12,'fontname','Times');
xlabel('steps', 'Fontname', 'Times New Roman','FontSize',12);
ylabel('$x_2$', 'Fontname', 'Times New Roman','Interpreter', 'latex','FontSize',12);
set(gca,'position',[0.11 0.24 0.86 0.30]);









%% useful fuction =========================================================
function e_max = maxEig( MatrixSym )
    e_max = max(eig(MatrixSym));
end

function e_min = minEig( MatrixSym )
    e_min = min(eig(MatrixSym));
end

function y = roundUpSignificant(x, n)
    magnitude = 10.^floor(log10(abs(x)) - n + 1);
    y = ceil(x ./ magnitude) .* magnitude;
    y(x == 0) = 0; 
end

function BarQ = funcBarQ(xiA11, xiA21_u,... 
                         xiA21_l, xiA12, xiA22_u, xiA22_l, xiA32_u, xiA32_l,... 
                         xiB11, xiB21_u, xiB21_l, xiB31_u, xiB31_l, xiB41_u, xiB41_l, xiB51_u, xiB51_l)
BarQ = {
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_u; xiB31_u; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_u; xiB31_u; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_u; xiB31_u; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_u; xiB31_u; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_u; xiB31_u; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_u; xiB31_u; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_u; xiB31_u; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_u; xiB31_u; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_l; xiB31_u; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_l; xiB31_u; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_l; xiB31_u; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_l; xiB31_u; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_l; xiB31_u; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_l; xiB31_u; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_l; xiB31_u; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_l; xiB31_u; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_u; xiB31_l; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_u; xiB31_l; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_u; xiB31_l; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_u; xiB31_l; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_u; xiB31_l; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_u; xiB31_l; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_u; xiB31_l; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_u; xiB31_l; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_l; xiB31_l; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_l; xiB31_l; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_l; xiB31_l; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_l; xiB31_l; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_l; xiB31_l; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_l; xiB31_l; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_l; xiB31_l; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_l; xiB31_l; xiB41_u; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_u; xiB31_u; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_u; xiB31_u; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_u; xiB31_u; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_u; xiB31_u; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_u; xiB31_u; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_u; xiB31_u; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_u; xiB31_u; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_u; xiB31_u; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_l; xiB31_u; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_l; xiB31_u; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_l; xiB31_u; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_l; xiB31_u; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_l; xiB31_u; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_l; xiB31_u; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_l; xiB31_u; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_l; xiB31_u; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_u; xiB31_l; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_u; xiB31_l; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_u; xiB31_l; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_u; xiB31_l; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_u; xiB31_l; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_u; xiB31_l; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_u; xiB31_l; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_u; xiB31_l; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_l; xiB31_l; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_l; xiB31_l; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_l; xiB31_l; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_l; xiB31_l; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_l; xiB31_l; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_l; xiB31_l; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_l; xiB31_l; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_l; xiB31_l; xiB41_l; xiB51_u] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_u; xiB31_u; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_u; xiB31_u; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_u; xiB31_u; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_u; xiB31_u; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_u; xiB31_u; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_u; xiB31_u; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_u; xiB31_u; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_u; xiB31_u; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_l; xiB31_u; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_l; xiB31_u; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_l; xiB31_u; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_l; xiB31_u; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_l; xiB31_u; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_l; xiB31_u; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_l; xiB31_u; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_l; xiB31_u; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_u; xiB31_l; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_u; xiB31_l; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_u; xiB31_l; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_u; xiB31_l; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_u; xiB31_l; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_u; xiB31_l; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_u; xiB31_l; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_u; xiB31_l; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_l; xiB31_l; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_l; xiB31_l; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_l; xiB31_l; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_l; xiB31_l; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_l; xiB31_l; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_l; xiB31_l; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_l; xiB31_l; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_l; xiB31_l; xiB41_u; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_u; xiB31_u; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_u; xiB31_u; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_u; xiB31_u; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_u; xiB31_u; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_u; xiB31_u; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_u; xiB31_u; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_u; xiB31_u; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_u; xiB31_u; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_l; xiB31_u; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_l; xiB31_u; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_l; xiB31_u; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_l; xiB31_u; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_l; xiB31_u; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_l; xiB31_u; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_l; xiB31_u; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_l; xiB31_u; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_u; xiB31_l; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_u; xiB31_l; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_u; xiB31_l; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_u; xiB31_l; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_u; xiB31_l; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_u; xiB31_l; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_u; xiB31_l; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_u; xiB31_l; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_l; xiB31_l; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_u], [xiB11; xiB21_l; xiB31_l; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_l; xiB31_l; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_u], [xiB11; xiB21_l; xiB31_l; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_l; xiB31_l; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_u; xiA32_l], [xiB11; xiB21_l; xiB31_l; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_u], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_l; xiB31_l; xiB41_l; xiB51_l] );
    blkdiag( [xiA11; xiA21_l], [xiA12; xiA22_l; xiA32_l], [xiB11; xiB21_l; xiB31_l; xiB41_l; xiB51_l] )
};
end









