clear; close all; clc;

dim_x = 2;
dim_u = 1;

ph = 0.1; % discrete time



%% original system =================================================

sample_s = 272;
Xs = zeros(dim_x, sample_s+1); Xs(:,1) = 0.01*randn(2,1);
Us = zeros(dim_u, sample_s );

for it = 1:sample_s
    Ax = [ 1+ph*sin(Xs(1,it))/Xs(1,it), ph*2; ph*2, 1-ph+ph*Xs(1,it)^2];
    Bx = [ ph*(1+abs(Xs(2,it))); ph*exp(Xs(1,it)) ];
    Xs(:,it+1) = Ax*Xs(:,it) + Bx*Us(:,it);
end
% figure(1)
% plot(Xs(2,:));

% phase portrait of the original system
xrr = 1; dx = xrr/10;
xbis = 0.0; ybis = -0.0; 
xleft = -xrr-xbis; xright = xrr-xbis;
yleft = -xrr-ybis; yright = xrr-ybis;
[x1, x2] = meshgrid(xleft:dx:xright, yleft:dx:yright);
Px1 = x1 + ph*sin(x1) + ph*2*x2 + ph*(1+abs(x2)).*0 - x1;
Px2 = ph*2*x1 + (1-ph+ph*x1.^2).*x2 + ph*exp(x1).*0 - x2;
figure('color', [1 1 1]);
set(gcf, 'Position', [100, 100, 400, 380]);
quiver(x1, x2 ,Px1 ,Px2, 3); hold on;
xlabel('$x_1$', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x_2$', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 12);
axis equal;
axis([xleft xright yleft yright]);
set(gca,'fontsize',12,'fontname','Times');








%% model based controllers =========================================
close;

% set the domain: |x| <= rx
rx = 1.1; 

% A(x)
a11_u = 1+ph;
a11_l = 1+ph*sin(rx)/rx;
a12 = 2*ph;
a21 = 2*ph;
a22_u = 1-ph+ph*rx^2;
a22_l = 1-ph;
% B(x)
b11_u = (1+rx)*ph;
b11_l = 1*ph;
b12_u = (0+exp(rx))*ph*1; 
b12_l = (0+exp(-rx))*ph*1; 

% the set of barG
barG = {
    [ a11_u, a12,   b11_u;
      a21,   a22_u, b12_u ];
    [ a11_l, a12,   b11_u;
      a21,   a22_u, b12_u ];
    [ a11_u, a12,   b11_u;
      a21,   a22_l, b12_u ];
    [ a11_l, a12,   b11_u;
      a21,   a22_l, b12_u ];
    [ a11_u, a12,   b11_l;
      a21,   a22_u, b12_u ];
    [ a11_l, a12,   b11_l;
      a21,   a22_u, b12_u ];
    [ a11_u, a12,   b11_l;
      a21,   a22_l, b12_u ];
    [ a11_l, a12,   b11_l;
      a21,   a22_l, b12_u ];
    [ a11_u, a12,   b11_u;
      a21,   a22_u, b12_l ];
    [ a11_l, a12,   b11_u;
      a21,   a22_u, b12_l ];
    [ a11_u, a12,   b11_u;
      a21,   a22_l, b12_l ];
    [ a11_l, a12,   b11_u;
      a21,   a22_l, b12_l ];
    [ a11_u, a12,   b11_l;
      a21,   a22_u, b12_l ];
    [ a11_l, a12,   b11_l;
      a21,   a22_u, b12_l ];
    [ a11_u, a12,   b11_l;
      a21,   a22_l, b12_l ];
    [ a11_l, a12,   b11_l;
      a21,   a22_l, b12_l ];
};


barU = 0.5; % actuator saturation

% solve sdp LMI
cvx_begin sdp % quiet
    variable epGa semidefinite;
    variable lyGa(dim_x,dim_x) semidefinite;
    variable kY(dim_u,dim_x);
    % variable mu_l semidefinite;
    % variable mu_u semidefinite;
    variable kW(dim_u,dim_x);
    variable kS(dim_u,dim_u) diagonal semidefinite;
    
    % minimize( mu_u - 0*mu_l + 1*epGa - 0.5*trace(lyGa) );
    minimize( lambda_max(lyGa) - lambda_min(lyGa) - 1*epGa - 0.1*trace(lyGa) );


    subject to   
        for iu = 1:dim_u
            blockW = [ lyGa,     kW(iu,:)';
                       kW(iu,:), barU^2 ];
            blockW >= 0;
        end

        for iv = 1:length(barG)
            blockS = [ -lyGa+epGa*eye(dim_x), -(kY+kW)',                          [lyGa;kY]'*barG{iv}';
                       -(kY+kW),              -2*kS,                              [zeros(dim_x,dim_u);kS]'*barG{iv}';
                       barG{iv}*[lyGa;kY],     barG{iv}*[zeros(dim_x,dim_u);kS], -lyGa ];
            blockS <= 0;
        end
        epGa >= 0.001;
        lyGa >= 0.01 * eye(dim_x);
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



xrr = 1.2; dx = xrr/8;
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











