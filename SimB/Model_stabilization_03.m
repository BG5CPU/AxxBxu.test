clear; close all; clc;

dim_x = 2;
dim_u = 1;

ph = 0.5; % discrete time
pe = 2; % van der pol parameter


%% original system =================================================
sample_s = 100;
Xs = zeros(dim_x, sample_s+1); Xs(:,1) = 0.00001*randn(2,1);
Us = zeros(dim_u, sample_s );

for it = 1:sample_s
    Ax = [ 1+1*ph, -ph; ph, 1-ph*pe+ph*pe*Xs(1,it)^2 ];
    Bx = ph*[ 1-0.45*Xs(2,it)^2+0.45*exp(Xs(1,it));
              2-0.30*Xs(2,it)^2-0.30*exp(Xs(1,it)) ];
    Xs(:,it+1) = Ax*Xs(:,it) + Bx*Us(:,it);
end
% plot(Xs(2,:));

% phase portrait of the original system
xrr = 1; dx = xrr/10;
xbis = 0.0; ybis = -0.0; 
xleft = -xrr-xbis; xright = xrr-xbis;
yleft = -xrr-ybis; yright = xrr-ybis;
[x1, x2] = meshgrid(xleft:dx:xright, yleft:dx:yright);
Px1 = (1+1*ph)*x1 - ph*x2 + ph*(1-0.45*x2.^2+0.45*exp(x1)).*0 - x1;
Px2 = ph*x1 + (1-ph*pe+ph*pe*x1.^2).*x2 + ph*(2-0.30*x2.^2-0.30*exp(x1)).*0 - x2;
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
rx = 0.27; 

% A(x)
a11 = 1+1*ph; 
a12 = -1*ph;
a21 = 1*ph; 
% a22 = 1-ph*pe+ph*pe*x1^2
a22_u = 1-ph*pe+ph*pe*rx^2;
a22_l = 1-ph*pe;
% B(x)
b11_u = ( 1 - 0.45*0^2 + 0.45*exp(rx) )*ph;
x1so = -lambertw(0, 1/2);
x2sq = rx^2 - x1so^2;
b11_l = ( 1 - 0.45*x2sq + 0.45*exp(x1so) )*ph;
b12_u = ( 2 - 0.30*0^2 - 0.30*exp(-rx) )*ph; 
b12_l = ( 2 - 0.30*0^2 - 0.30*exp(rx) ) *ph; 

% the set of barG
barG = {
    [ a11, a12,   b11_u;
      a21, a22_u, b12_u ];
    [ a11, a12,   b11_l;
      a21, a22_u, b12_u ];
    [ a11, a12,   b11_u;
      a21, a22_l, b12_u ];
    [ a11, a12,   b11_l;
      a21, a22_l, b12_u ];
    [ a11, a12,   b11_u;
      a21, a22_u, b12_l ];
    [ a11, a12,   b11_l;
      a21, a22_u, b12_l ];
    [ a11, a12,   b11_u;
      a21, a22_l, b12_l ];
    [ a11, a12,   b11_l;
      a21, a22_l, b12_l ];
};


% solve sdp LMI
cvx_begin sdp % quiet
    variable epGa semidefinite;
    variable lyGa(dim_x,dim_x) semidefinite;
    variable kY(dim_u,dim_x);
    variable mu_l semidefinite;
    variable mu_u semidefinite;
    
    % minimize( lambda_max(lyGa) - 1.0*lambda_min(lyGa) - 0.1*epGa );
    % minimize( lambda_max(lyGa) );
    % minimize( lambda_max(lyGa) - 1*epGa );
    % minimize( lambda_max(lyGa) - 1*trace(lyGa) - 0.1*epGa );
    minimize( mu_u - 0.0*mu_l + 1*epGa );
    % minimize( trace(lyGa) - 1*epGa );
    % minimize( lambda_max(lyGa) - lambda_min(lyGa) + epGa - 0.1*trace(lyGa) );

    subject to   
        for iv = 1:length(barG)
            blockS = [ -lyGa,                   barG{iv}*[lyGa; kY];
                        [lyGa; kY]'*barG{iv}', -lyGa + epGa*eye(dim_x) ];
            blockS <= 0;
        end
        epGa >= 0.01;
        % lyGa >= 0.1 * eye(dim_x);
        lyGa >= mu_l * eye(dim_x);
        lyGa <= (mu_l+mu_u) * eye(dim_x);
        
cvx_end

vK = kY/lyGa;

% region of attraction
roac = sqrt( maxEig(lyGa)*minEig(lyGa) / ( maxEig(lyGa)^2 - epGa*minEig(lyGa) ) );
roa  = min( 1, roac ) * rx;
disp("radius of ROA is " + roa);






%% phase portrait ==================================================
xrr = 0.5; dx = xrr/10;
xbis = 0.15; ybis = -0.15; 
xleft = -xrr-xbis; xright = xrr-xbis;
yleft = -xrr-ybis; yright = xrr-ybis;
[x1, x2] = meshgrid(xleft:dx:xright, yleft:dx:yright);
Px1 = (1+1*ph)*x1 - ph*x2 + ph*(1-0.45*x2.^2+0.45*exp(x1)).*(vK(1)*x1+vK(2)*x2) - x1;
Px2 = ph*x1 + (1-ph*pe+ph*pe*x1.^2).*x2 + ph*(2-0.30*x2.^2-0.30*exp(x1)).*(vK(1)*x1+vK(2)*x2) - x2;
figure('color', [1 1 1]);
set(gcf, 'Position', [100, 100, 400, 380]);
quiver(x1, x2 ,Px1 ,Px2, 2); hold on;
xlabel('$x_1$', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x_2$', 'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 12);
axis equal;
axis([xleft xright yleft yright]);
set(gca,'fontsize',12,'fontname','Times');
fill(cos(0:0.01*pi:2*pi)*roa, sin(0:0.01*pi:2*pi)*roa, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
legend('$(x_1^{+}\!\!-\!x_1 , x_2^{+}\!\!-\!x_2)$', ...
       'region of attraction',...
       'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 12);




   

%% controller implementation ==============================================
ed = 0.01;
tn = 100; timeS = 0:1:tn;
xt = zeros(dim_x,tn+1); 
xt(:,1) = ones(dim_x,1) *0.15;
disp("initial value is " + xt(1,1) + ", " + xt(2,1));
for ik  = 1:tn
    dt = (rand(dim_x,1)-0.5)*2*ed;
    
    Ax = [ 1+1*ph, -ph; ph, 1-ph*pe+ph*pe*xt(1,ik)^2 ];
    Bx = ph*[ 1-0.45*xt(2,ik)^2+0.45*exp(xt(1,ik));
              2-0.30*xt(2,ik)^2-0.30*exp(xt(1,ik)) ];
          
    ut = vK*xt(:,ik);
    xt(:,ik+1) = Ax*xt(:,ik) + Bx*ut + dt;
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









%% useful fuction =========================================================
function e_max = maxEig( MatrixSym )
    e_max = max(eig(MatrixSym));
end

function e_min = minEig( MatrixSym )
    e_min = min(eig(MatrixSym));
end












