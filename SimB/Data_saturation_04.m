clear; close all; clc;

dim_x = 2;
dim_u = 1;

ph = 0.1; % discrete time

%% original system =================================================

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
edis = 0.001;  % bound of output noise

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



%% data from previous experiment ===================================
load("Mdata/data_stabilization_09.mat");




%% data based controllers ==========================================
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



barU = 0.25; % actuator saturation

% solve sdp LMI
cvx_begin sdp % quiet
    variable epGa semidefinite;
    variable lyGa(dim_x,dim_x) semidefinite;
    variable kY(dim_u,dim_x);
    % variable mu_l semidefinite;
    % variable mu_u semidefinite;
    variable kW(dim_u,dim_x);
    variable kS(dim_u,dim_u) diagonal semidefinite;
    
    % minimize( 0.01*mu_u - 0*mu_l + 0*epGa - 1*trace(epGa) );
    minimize( lambda_max(lyGa) - lambda_min(lyGa) - 0*epGa - 2*trace(lyGa) );

    subject to   
        for iu = 1:dim_u
            blockW = [ lyGa,     kW(iu,:)';
                       kW(iu,:), barU^2 ];
            blockW >= 0;
        end
        
        for iv = 1:length(barQ)
            blockS = [ -lyGa+epGa*eye(dim_x), -(kY+kW)',                          zeros(dim_x),       [lyGa;kY]'*barQ{iv}';
                       -(kY+kW),              -2*kS,                              zeros(dim_u,dim_x), [zeros(dim_x,dim_u);kS]'*barQ{iv}';
                        zeros(dim_x),          zeros(dim_x,dim_u),               -lyGa-dataC,        -dataB;
                        barQ{iv}*[lyGa;kY],    barQ{iv}*[zeros(dim_x,dim_u);kS], -dataB',            -dataA ];
            blockS <= 0;
        end
        epGa >= 0.00001;
        lyGa >= 0.0001 * eye(dim_x);
        % lyGa >= mu_l * eye(dim_x);
        % lyGa <= (mu_l+mu_u) * eye(dim_x);
cvx_end

vK = kY/lyGa;
lyP = inv(lyGa);

disp(lyGa);

% region of attraction
roac = sqrt( maxEig(lyGa)*minEig(lyGa) / ( maxEig(lyGa)^2 - epGa*minEig(lyGa) ) );
roa  = min( 1, roac ) * rx;
disp("radius of ROA is " + roa);






%% phase portrait ==================================================

% for the Br, circle
theta = linspace(0, 2*pi, 100); 
x1br = roa * cos(theta);
x2br = roa * sin(theta);


% for ellipsoid
[Q, Lambda] = eig(lyP);
a = sqrt(1 / Lambda(1,1));  % axis1
b = sqrt(1 / Lambda(2,2));  % axis2
% points for circle
theta = linspace(0, 2*pi, 100);
x_circle = a * cos(theta);
y_circle = b * sin(theta);
% points for ellipsoid
xy_rotated = Q * [x_circle; y_circle];
x1Ep = xy_rotated(1,:);
x2Ep = xy_rotated(2,:);



xrr = 0.5; dx = xrr/8;
xbis = 0.0; ybis = -0.0; 
xleft = -xrr-xbis; xright = xrr-xbis;
yleft = -xrr-ybis; yright = xrr-ybis;
[x1, x2] = meshgrid(xleft:dx:xright, yleft:dx:yright);

Px1 = x1 + ph*sin(x1) + ph*2*x2 + ph*(1+abs(x2)).*saturU(vK(1)*x1+vK(2)*x2,barU) - x1;
Px2 = ph*2*x1 + (1-ph+ph*x1.^2).*x2 + ph*exp(x1).*saturU(vK(1)*x1+vK(2)*x2,barU) - x2;

figure('color', [1 1 1]);
set(gcf, 'Position', [100, 100, 400, 380]);

quiver(x1, x2 ,Px1 ,Px2, 2); hold on;
plot(x1br, x2br, 'r-', 'LineWidth', 2); hold on;
plot(x1Ep, x2Ep, 'k:', 'LineWidth', 2);

xlabel('$x_1$', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$x_2$', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 18);

title(sprintf('$\\bar{u} = %.2f$', barU), 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 18);

axis equal;
axis([xleft xright yleft yright]);
set(gca,'fontsize',18,'fontname','Times');




  








%% useful fuction =========================================================
function e_max = maxEig( MatrixSym )
    e_max = max(eig(MatrixSym));
end

function e_min = minEig( MatrixSym )
    e_min = min(eig(MatrixSym));
end

function uSat = saturU(u, bar_u)
    uSat = min(max(u, -bar_u), bar_u);
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









