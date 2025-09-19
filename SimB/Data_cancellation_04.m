clear; close all; clc;



%% data from previous experiment ==========================================
load("Mdata/data_stabilization_09.mat");



%% data matrices ==========================================================
dim_Q = 4;
dim_Z = 7;

XI1 = []; ZXI0 = []; V0 = [];

% for it = 1:trajt
%     u0 = Dinput0(:,:,it);
%     x0 = Dstate0(:,:,it);
%     v0 = [Dinput0(:,2:end,it), zeros(dim_u,1)];
%     V0 = [V0, v0];
%     XI1 = [XI1, [Dstate1(:,:,it); v0]];
%     zxi0 = [ x0;
%              u0;
%              sin(x0(1,:));
%              x0(1,:).^2.*x0(2,:);
%              abs(x0(2,:)).*u0;
%              exp(x0(1,:)).*u0 ];
%     ZXI0 = [ZXI0, zxi0];
% end
% 
% Acoi = [ 1,   0.2, 0.1, 0.1, 0,   0.1, 0;
%          0.2, 0.9, 0,   0,   0.1, 0,   0.1;
%          0,   0,   0,   0,   0,   0,   0 ];
% Bcoi = [0; 0; 1];
% Ecoi = [eye(2); [0, 0]];
% TestDataCancel = Acoi*ZXI0 + Bcoi*V0 + Ecoi*D0 - XI1;

for it = 1:trajt
    u0 = Dinput0(:,:,it);
    x0 = Dstate0(:,:,it);
    v0 = [Dinput0(:,2:end,it), zeros(dim_u,1)];
    V0 = [V0, v0];
    XI1 = [XI1, [Dstate1(:,:,it); v0]];
    zxi0 = [ x0;
             u0;
             sin(x0(1,:)) - x0(1,:);
             x0(1,:).^2.*x0(2,:);
             abs(x0(2,:)).*u0;
             exp(x0(1,:)).*u0 - u0];
    ZXI0 = [ZXI0, zxi0];
end

Acoi = [ 1.1, 0.2, 0.1, 0.1, 0,   0.1, 0;
         0.2, 0.9, 0.1, 0,   0.1, 0,   0.1;
         0,   0,   0,   0,   0,   0,   0 ];
Bcoi = [0; 0; 1];
Ecoi = [eye(2); [0, 0]];
TestDataCancel = Acoi*ZXI0 + Bcoi*V0 + Ecoi*D0 - XI1;

fQx = @(qx) [ sin(qx(1)) - qx(1); 
              qx(1)^2*qx(2);
              abs(qx(2))*qx(3);
              exp(qx(1))*qx(3) - qx(3) ];


%% solve controller =======================================================

Tlen = length(XI1);
dim_w = dim_x+dim_u;

OMega = 1*eye(dim_w);
EDelt = Ecoi*THETAnoise*eye(dim_x)*Ecoi';

% cvx_begin sdp % quiet
%     variable cY1(Tlen,dim_w);
%     variable cG2(Tlen,dim_Q);
%     variable cP1(dim_w,dim_w) semidefinite
% 
%     minimize( norm(XI1*cG2) );
% 
%     subject to
%         ZXI0*cY1 == [cP1; zeros(dim_Q,dim_w)];
%         [ cP1, (XI1*cY1)'; XI1*cY1, cP1 ] >= 0;
%         ZXI0*cG2 == [zeros(dim_w,dim_Q); eye(dim_Q)];
%         cP1 >= 0.0001 * eye(dim_w);
% cvx_end
% 
% cK = V0*[cY1/cP1, cG2];


cvx_begin sdp
    variable cY1(Tlen,dim_w);
    variable cG2(Tlen,dim_Q);
    variable cP1(dim_w,dim_w) semidefinite
    variable cep(1,1) semidefinite

    minimize( norm(XI1*cG2) );
    
    subject to 
        ZXI0*cY1 == [cP1; zeros(dim_Q,dim_w)];

        [ cP1-OMega,  (XI1*cY1)',        cY1';
          XI1*cY1,    cP1-cep*EDelt,     zeros(dim_w,Tlen);
          cY1,        zeros(Tlen,dim_w), cep*eye(Tlen) ] >= 0;

        ZXI0*cG2 == [zeros(dim_w,dim_Q); eye(dim_Q)];
cvx_end

% this is the controller
cK = V0*[cY1/cP1, cG2];





%% region of attraction ===================================================

cG1 = cY1/cP1;
hdlt = sqrt(THETAnoise);

dtsx = 0.005;
tsx1 = -1:dtsx:1;   % x1-axis (horizontal)
tsx2 = -1:dtsx:1;   % x2-axis (vertical)

nCol = length(tsx1);        % x1 -> columns
nRow = length(tsx2);        % x2 -> rows

% Initialize mask
mask = false(nRow, nCol);   % rows = x2, cols = x1

% for the set L
% Evaluate isNegative at each (x1, x2)
for i = 1:nRow
    for j = 1:nCol
        ellval1 = funcValue_ell([tsx1(j); tsx2(i); 0.0], cP1, OMega, XI1, cG1, cG2, fQx, hdlt, Ecoi);
        mask(i, j) = ellval1<0;
    end
end


% for the sublevel set
cPx = inv(cP1);
cPx = cPx(1:2,1:2);
cr = 0.0018;
[cQ, cLambda] = eig(cPx);
ca = sqrt(cr / cLambda(1,1));  % axis1
cb = sqrt(cr / cLambda(2,2));  % axis2
% points for circle
ctheta = linspace(0, 2*pi, 100);
cx_circle = ca * cos(ctheta);
cy_circle = cb * sin(ctheta);
% points for ellipsoid
cxy_rotated = cQ * [cx_circle; cy_circle];
cx1Ep = cxy_rotated(1,:);
cx2Ep = cxy_rotated(2,:);



% for the Br, circle
theta = linspace(0, 2*pi, 100); 
x1br = roa * cos(theta);
x2br = roa * sin(theta);


% Plot
figure('Color', 'w');
set(gcf, 'Position', [100, 100, 400, 380], 'Renderer', 'Painters');

imagesc(tsx1, tsx2, mask);           % x-axis: x1, y-axis: x2
colormap([1 1 1; 0.6 0.6 0.6]);      % false white，true gray
axis xy;                             % Make y-axis increase upward

hold on;

plot(cx1Ep, cx2Ep, 'k:', 'LineWidth', 2);

plot(x1br, x2br, 'r-', 'LineWidth', 2);

xlabel('$x_1$', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x_2$', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 12);
axis equal;
grid on;
axis([-1 1 -1 1]);
set(gca,'fontsize',12,'fontname','Times');

legend('$\mathcal{R}_{\gamma}$',...
       '$\mathcal{B}_{\bar{r}_0}$', ...
       'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 12);











%% test ===================================================================
ellval1 = funcValue_ell([0; 0; 0], cP1, OMega, XI1, cG1, cG2, fQx, hdlt, Ecoi);
disp(ellval1);





%% useful functions =======================================================
function ellval = funcValue_ell( x, P, Omega, XI1, G1, G2, Q, Del, E)
    Psi_ = P\Omega/P;
    ell0 = -x'*Psi_*x;
    ell1 = ( 2*XI1*G1*x + XI1*G2*Q(x) )' /P * XI1*G2*Q(x);
    ell2 = Del * norm( ( 2*XI1*G1*x + XI1*G2*Q(x) )' /P*E ) * norm( G2*Q(x) );
    ell3 = Del * norm( 2*G1*x + G2*Q(x) ) * norm( E'/P*XI1*G2*Q(x) );
    ell4 = Del^2 * norm( E'/P*E ) * norm( 2*G1*x + G2*Q(x) ) * norm( G2*Q(x) );

    ellval = ell0+ell1+ell2+ell3+ell4;
end














