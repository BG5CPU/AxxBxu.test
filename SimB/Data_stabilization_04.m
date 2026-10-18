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
    variable mu_l nonnegative;
    variable mu_u nonnegative;
    
    % epigraph form of  lambda_max(lyGa) - lambda_min(lyGa) - 0.1*epGa,
    % with mu_l <= lambda_min(lyGa) and mu_u >= lambda_max(lyGa)
    minimize( mu_u - mu_l - 0.1*epGa );

    subject to   
        for iv = 1:length(barQ)
            blockS = [ -lyGa-dataC,   zeros(dim_x),          dataB;
                       zeros(dim_x), -lyGa+epGa*eye(dim_x), -[lyGa;kY]'*barQ{iv}'; 
                       dataB',       -barQ{iv}*[lyGa;kY],   -dataA ];
            blockS <= 0;
        end
        epGa >= 0.00001; %0.00001;
        lyGa >= 0.0001 * eye(dim_x);
        lyGa >= mu_l * eye(dim_x);
        lyGa <= mu_u * eye(dim_x);
        
cvx_end

vK = kY/lyGa;

% region of attraction
roac = sqrt( maxEig(lyGa)*minEig(lyGa) / ( maxEig(lyGa)^2 - epGa*minEig(lyGa) ) );
roa  = min( 1, roac ) * rx;
disp("radius of ROA is " + roa);


% robustness parameters of (58): gammaAcl, bar gammaAcl, delta_x0, delta_w
rbst = funcRobustPara( rx, 0.01, 0.01*pi, [EA1, EA2], EB1, vK, barQ, ...
                       dataA, dataB, dataC, lyGa, epGa );






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









%% robustness parameters ==================================================
% Computes the quantities of (58) and Remark 1 by sampling x over B_r:
%   gAcl    : gamma_{Acl}    = max_{x in B_r} || A(x) + B(x)K ||      (58d)
%   gBarMat : bar gamma_{Acl} with max || [XiA(x); XiB(x)K] ||, i.e. the induced
%             norm of the map, since A(x)+B(x)K = [EA EB]*[XiA(x); XiB(x)K]
%   gBarVtx : max_{Q in barQ} || Q*[I; K] ||, a sampling-free upper bound of gBarMat
%             ([XiA(x); XiB(x)K] = Q(x)*[I; K] is linear in Q(x), the spectral norm
%             is convex, and Q(x) lies in conv(barQ) by (49)). Reported for reference
%             only: barQ over-approximates Q, so this value is valid but loose.
%   dx0     : delta_x0 of (58b)
%   dw      : delta_w of (58c), evaluated at gAcl
%   dwB     : delta_w of (58c), evaluated at gBarMat
function res = funcRobustPara( rx, dr, dth, EA, EB, vK, barQ, ...
                               dataA, dataB, dataC, lyGa, epGa )
% data-computable constant || Zc || + || Q^{1/2} || || A^{-1/2} ||  of (44b)
EA_EB = [EA, EB];
Zc = -dataB/dataA;
Qc = dataB/dataA*dataB' - dataC; Qc = (Qc+Qc')/2;
cst = norm(Zc) + sqrt(max(eig(Qc))) / sqrt(min(eig(dataA)));

gAcl = 0; gMat = 0;
for rs = dr:dr:rx
    for ts = 0:dth:2*pi-dth
        x1t = rs*cos(ts); x2t = rs*sin(ts);
        if abs(x1t) < 1e-12, sc1 = 1; else, sc1 = sin(x1t)/x1t; end
        XiA = blkdiag( [1; sc1], [1; x1t; x1t^2] );
        XiB = [1; abs(x1t); abs(x2t); exp(x1t); exp(x2t)];
        Acl = EA*XiA + EB*XiB*vK;
        gAcl = max( gAcl, norm(Acl) );
        gMat = max( gMat, norm( [XiA; XiB*vK] ) );
    end
end

% vertex bound over barQ, no sampling involved
IK = [eye(size(vK,2)); vK];
gVtx = 0;
for iv = 1:length(barQ)
    gVtx = max( gVtx, norm( barQ{iv}*IK ) );
end

lmax = max(eig(lyGa)); lmin = min(eig(lyGa));
dx0  = ( 2*lmax^2 - epGa*lmin ) / ( 2*lmin*lmax );
dwf  = @(g) ( 4*g^2*lmax^5 + 2*epGa*lmin*lmax^3 ) / ( epGa^2*lmin^3 );

res = struct( 'gAcl', gAcl, 'gBarMat', cst*gMat, 'gBarVtx', cst*gVtx, ...
              'cst', cst, ...
              'dx0', dx0, 'dw', dwf(gAcl), 'dwB', dwf(cst*gMat) );

disp("norm of [EA EB] is " + norm(EA_EB) + ", data-driven constant is " + res.cst);
disp("gamma_Acl is " + res.gAcl);
disp("bar gamma_Acl (sampled, induced norm) is " + res.gBarMat);
disp("bar gamma_Acl (vertices of barQ, loose) is " + res.gBarVtx);
disp("delta_x0 is " + res.dx0);
disp("delta_w is " + res.dw + " with gamma_Acl, " + res.dwB + " with bar gamma_Acl");
end




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









