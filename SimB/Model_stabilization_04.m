clear; close all; clc;

dim_x = 2;
dim_u = 1;

ph = 0.1; % discrete time



%% original system ========================================================

sample_s = 272;
Xs = zeros(dim_x, sample_s+1); Xs(:,1) = 0.01*randn(2,1);
Us = zeros(dim_u, sample_s );

for it = 1:sample_s
    Ax = [ 1+ph*sin(Xs(1,it))/Xs(1,it), ph*2; ph*2, 1-ph+ph*Xs(1,it)^2];
    Bx = [ ph*(1+abs(Xs(2,it))); ph*exp(Xs(1,it)) ];
    Xs(:,it+1) = Ax*Xs(:,it) + Bx*Us(:,it);
end
figure(1)
plot(Xs(2,:));

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







%% model based controllers ================================================
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


% solve sdp LMI
cvx_begin sdp % quiet
    variable epGa semidefinite;
    variable lyGa(dim_x,dim_x) semidefinite;
    variable kY(dim_u,dim_x);
    % variable mu_l semidefinite;
    % variable mu_u semidefinite;
    
    % minimize( lambda_max(lyGa) - 1.0*lambda_min(lyGa) - 0.1*epGa );
    % minimize( lambda_max(lyGa) );
    % minimize( lambda_max(lyGa) - 1*epGa );
    % minimize( lambda_max(lyGa) - 1*trace(lyGa) - 0.1*epGa );
    % minimize( mu_u - 0*mu_l + 1*epGa );
    % minimize( trace(lyGa) - 1*epGa );
    minimize( lambda_max(lyGa) - lambda_min(lyGa) - 0.0*epGa - 0.0*trace(lyGa) );

    subject to   
        for iv = 1:length(barG)
            blockS = [ -lyGa,                   barG{iv}*[lyGa; kY];
                        [lyGa; kY]'*barG{iv}', -lyGa + epGa*eye(dim_x) ];
            blockS <= 0;
        end
        epGa >= 0.001;
        lyGa >= 0.01 * eye(dim_x);
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

legend('$\mathcal{B}_{r_0}$', ...
       '$(x_1^{+}\!\!-\!x_1 , x_2^{+}\!\!-\!x_2)$', ...
       'Fontname', 'Times New Roman', 'Interpreter', 'latex', 'FontSize', 12);




   

%% controller implementation ==============================================
ed = 0.1;
tn = 300; timeS = 0:1:tn;
xt = zeros(dim_x,tn+1); 
xt(:,1) = [-0.5, -0.5];
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
set(gca,'position',[0.11 0.65 0.86 0.30]);
subplot(2,1,2);
plot(timeS, xt(2,:), 'K', 'LineWidth', 1); grid on;
set(gca,'fontsize',12,'fontname','Times');
xlabel('steps', 'Fontname', 'Times New Roman','FontSize',12);
ylabel('$x_2$', 'Fontname', 'Times New Roman','Interpreter', 'latex','FontSize',12);
set(gca,'position',[0.11 0.22 0.86 0.30]);









%% useful fuction =========================================================
function e_max = maxEig( MatrixSym )
    e_max = max(eig(MatrixSym));
end

function e_min = minEig( MatrixSym )
    e_min = min(eig(MatrixSym));
end












