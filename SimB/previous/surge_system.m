%% NUMERICAL EXAMPLE ON CONTRACTIVITY
% system: surge system
%  REMARKS
%  In this example we pick Q(x)=[x1^2 x1^3]' or more nonlinearities Q(x)=[x1^2 x1^3 x2^2 x2^3]' 

%% Initialization
clear,clc
rng(3);

%% System parameters
global K

n = 2; % dimension of state
m = 1; % dimension of input

%% Controllability test
Abar = [0 -1; 0 0];
B = [0 1]';
rank(ctrb(Abar, B))

%% Data acquisition phase via Simulink
T    = 10;  % number of samples
Ts   = 0.1; % sampling interval
Tsim = T*Ts; % duration of simulation

mag = 1; % magnitude of initial conditions
x0  = (2*mag).*rand(n,1)-mag; % initial state

sim('data_collection_surge_system');

x  = state.signals.values'; 
xd = state_deriv.signals.values'; 
u  = input.signals.values'; 

X0  = x(:,1:T);
U0  = u(:,1:T);
X1  = xd(:,1:T);

s = 4; % dimension of Z(x) ----rng(3)
Z0  = [X0;X0(1,:).^2;X0(1,:).^3];
w = 1;
RQ = [sqrt(4*w^2+9*w^4) 0; 0 0]; % Bound for Jacobian of Q(x)

% s = 6; % dimension of Z(x) ----rng(1)
% Z0  = [X0;X0(1,:).^2;X0(1,:).^3;X0(2,:).^2;X0(2,:).^3];
% w1 = 1;
% w2 = 0.1;
% RQ = [sqrt(4*w1^2+9*w1^4) 0; 0 sqrt(4*w2^2+9*w2^4)]; % Bound for Jacobian of Q(x)

rank([U0; Z0]) % full row rank check

%% Controller design (CONTRACTIVITY)
r = n; % column number of RQ

cvx_begin sdp
    variable P1(n,n) symmetric
    variable Y1(T,n)
    variable G2(T,s-n)
    variable a 
    P1 >= 0*eye(n);
    a >= 0;
    Z0*Y1 == [P1;zeros(s-n,n)];
    [X1*Y1+transpose(X1*Y1)+a*eye(n) X1*G2 P1*RQ;
        transpose(X1*G2) -eye(s-n) zeros(s-n,r);
        transpose(P1*RQ) zeros(r,s-n) -eye(r)] <= 0;
    Z0*G2 == [zeros(n,s-n);eye(s-n)];
cvx_end

G1 = Y1/P1;
G  = [G1 G2];
K  = U0*G; 
closed = X1*G;

%% Estimate of ROA for closed-loop system
range = 200;
v = -range:0.1:range;  
[x1 x2] = meshgrid(v);
[e,f] = size(x1);

P1inv=inv(P1);

% Vector Q(\xi)
Qx1 = x1.^2;
Qx2 = x1.^3;
Qx3 = x2.^2;
Qx4 = x2.^3;

% Vector \dot \xi
M = X1*G1;
N = X1*G2;
xdot1 = M(1,1)*x1+M(1,2)*x2+N(1,1)*Qx1+N(1,2)*Qx2; 
xdot2 = M(2,1)*x1+M(2,2)*x2+N(2,1)*Qx1+N(2,2)*Qx2; 
% xdot1 = M(1,1)*x1+M(1,2)*x2+N(1,1)*Qx1+N(1,2)*Qx2+N(1,3)*Qx3+N(1,4)*Qx4; 
% xdot2 = M(2,1)*x1+M(2,2)*x2+N(2,1)*Qx1+N(2,2)*Qx2+N(2,3)*Qx3+N(2,4)*Qx4; 

% H1 = \dot \xi'*P1^{-1}*\xi 
H1 =     P1inv(1,1)*x1.*xdot1 + P1inv(1,2)*x1.*xdot2...
       + P1inv(2,1)*x2.*xdot1 + P1inv(2,2)*x2.*xdot2;

% H2 = \xi'*P1^{-1}*\xi
H2 =     P1inv(1,1)*x1.*x1 + P1inv(1,2)*x2.*x1...
       + P1inv(2,1)*x1.*x2 + P1inv(2,2)*x2.*x2;

H3 = zeros(e,f); % set of the points for which the Lyapunov derivative is positive inside the Lyapunov sublevel set H2
index   = []; % record the number of points in H3

gamma = 95; % range of Lyapunov sublevel set H2; modify gamma until H3 is empty 300

for i= 1:e
    for j= 1:f
        if H2(i,j) <= gamma
            if H1(i,j) > 0
                index = [index [i j]'];
            end
        else
            continue;
        end
    end
end

number_index = size(index);
for i=1 : number_index(2)
    H3(index(1,i),index(2,i)) = 1;
end

cond1 = H1 < 0; % Lyapunov derivative is negative
cond1 = double(cond1);  
cond1(cond1 == 0) = NaN;
colormap([0.5 0.5 0.5;0 0 1;1 0 0])% gray, blue and red  
surf(x1,x2,cond1,ones(length(v))); % gray
shading interp

hold on
cond2 = H2 <= gamma; % Lyapunov sublevel set
cond2 = double(cond2);
cond2(cond2 == 0) = NaN;
surf(x1,x2,cond2,ones(length(v))+1); % blue
shading interp

hold on
cond3 = H3 > 0; 
cond3 = double(cond3);
cond3(cond3 == 0) = NaN;
surf(x1,x2,cond3,ones(length(v))+2); % red
shading interp

axis([-30 30 -200 200])
view(0,90) 
xlabel('x_1') 
ylabel('x_2') 
% set(gca,'xtick',-range:0.2:range) 
% set(gca,'ytick',-range:2:range) 

%% Evaluation of obtained controller via ODE45 function
mx0 = 1;   % magnitude of initial conditions
x0  = (2*mx0).*rand(n,1)-mx0;
tspan = [0,10];  % duration of simulation

[t,x] = ode45(@inverted,tspan,x0);

figure;
plot(t,x(:,1),'r');
hold on;
plot(t,x(:,2),'b');
xlabel('t');
legend('x(1) ','x(2)');

function dxdt = inverted(t,x)  
    global K
    dxdt = zeros(2,1);
    dxdt(1) = -x(2) - 1.5*x(1)^2 - 0.5*x(1)^3;
    dxdt(2) = K*[x;x(1)^2;x(1)^3];
    % dxdt(2) = K*[x;x(1)^2;x(1)^3;x(2)^2;x(2)^3]; % more nonlinearities
    % dxdt(2) = 0; % check whether open-loop system is unstable
end
