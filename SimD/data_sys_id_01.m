clear; clc; close all;


%% load data ==============================================================
% load("Gazebo_data_mat\data02.mat");
load("Gazebo_data_mat/data04.mat");


data_input_output = data_input_output_04(4205:8705,:);

[nRow, nColumn] = size(data_input_output);
num_steps = nRow;

cl_time = 1;
cl_u1 = 2; cl_u2 = 4; cl_u3 = 6;
% cl_q1 = 14; cl_q2 = 16; cl_q3 = 18; cl_q0 = 20;
% cl_x4 = 28; cl_x5 = 30; cl_x6 = 32;
cl_x1 = 08; cl_x2 = 10; cl_x3 = 12;
cl_x4 = 14; cl_x5 = 16; cl_x6 = 18;

% dif_time = data_input_output(:,1) - data_input_output(:,3);

sample_time_diff = diff(data_input_output(:,cl_time));
% plot(sample_time, '*');
h = mode(sample_time_diff)/1e9;


% quat2eul (Robotics System Toolbox)
% [w, x, y, z] -> quat2eul -> [Z, Y, X] (Yaw, Pitch, Roll)
% Euler_rad = quat2eul(data_input_output(:,[cl_q0, cl_q1, cl_q2, cl_q3]), 'ZYX'); % [Yaw, Pitch, Roll]
% Euler_RPY_rad = [Euler_rad(:, 3), Euler_rad(:, 2), Euler_rad(:, 1)]; % [Roll, Pitch, Yaw]

% x_exp = [Euler_RPY_rad, data_input_output(:,[cl_x4, cl_x5, cl_x6])]';
x_exp = data_input_output(:,[cl_x1, cl_x2, cl_x3, cl_x4, cl_x5, cl_x6])';
u_exp = data_input_output(:,[cl_u1, cl_u2, cl_u3])';



%% arrange data ===========================================================

% system dimension
dim_x = 6;
dim_u = 3;

% the inertia matrix with respect to the body-fixed frame
Jx = 0.029125; Jy = 0.029125; Jz = 0.055225;
Jinertia = blkdiag(Jx,Jy,Jz);
coeUnk = 1;
BinM = [zeros(3,3); eye(3)/Jinertia]*h* coeUnk;
      

% xiA1 = [ 1 ];
f_xiA1 = @(x) 1;
EA1 = [1; 0; 0; 0; 0; 0];
dim_xiA1 = 1;
% xiA2 = [ 1 ];
f_xiA2 = @(x) 1;
EA2 = [0; 1; 0; 0; 0; 0];
dim_xiA2 = 1;
% xiA3 = [ 1 ];
f_xiA3 = @(x) 1;
EA3 = [0; 0; 1; 0; 0; 0];
dim_xiA3 = 1;
% xiA4 = [ 1; x6 ];
f_xiA4 = @(x) [1; x(6)];
EA4 = [ h,            0; 
        0,            0; 
        0,            0; 
        1,            0; 
        0, h*(Jz-Jx)/Jy; 
        0,            0 ];
dim_xiA4 = 2;
% xiA5 = [ 1; sin(x1)*tan(x2); cos(x1); sin(x1)/cos(x2); x4 ];
f_xiA5 = @(x) [1; sin(x(1))*tan(x(2)); cos(x(1)); sin(x(1))/cos(x(2)); x(4)];
EA5 = [ 0, h, 0, 0, 0; 
        0, 0, h, 0, 0;
        0, 0, 0, h, 0;
        0, 0, 0, 0, 0;
        1, 0, 0, 0, 0;
        0, 0, 0, 0, h*(Jx-Jy)/Jz ];
dim_xiA5 = 5;
% xiA6 = [ 1; cos(x1)*tan(x2); sin(x1); cos(x1)/cos(x2); x5 ];
f_xiA6 = @(x) [ 1; cos(x(1))*tan(x(2)); sin(x(1)); cos(x(1))/cos(x(2)); x(5) ];
EA6 = [ 0, h,  0, 0, 0; 
        0, 0, -h, 0, 0;
        0, 0,  0, h, 0;
        0, 0,  0, 0, h*(Jy-Jz)/Jx;
        0, 0,  0, 0, 0;
        1, 0,  0, 0, 0 ];
dim_xiA6 = 5;
dim_xiA = dim_xiA1+dim_xiA2+dim_xiA3+dim_xiA4+dim_xiA5+dim_xiA6;

% xiB1 = [ 1 ];
f_xiB1 = @(x) 1;
EB1 = BinM(:,1);
dim_xiB1 = 1;
% xiB2 = [ 1 ];
f_xiB2 = @(x) 1;
EB2 = BinM(:,2);
dim_xiB2 = 1;
% xiB3 = [ 1 ];
f_xiB3 = @(x) 1;
EB3 = BinM(:,3);
dim_xiB3 = 1;
dim_xiB = dim_xiB1+dim_xiB2+dim_xiB3;


bXI1 = x_exp(:,2:num_steps);
bXI0 = zeros(dim_xiA, num_steps-1); 
bUI0 = zeros(dim_xiB, num_steps-1);
for ik = 1:num_steps-1
    bXI0(:,ik) = blkdiag( f_xiA1(x_exp(:,ik)), ...
                          f_xiA2(x_exp(:,ik)), ...
                          f_xiA3(x_exp(:,ik)), ...
                          f_xiA4(x_exp(:,ik)), ...
                          f_xiA5(x_exp(:,ik)), ...
                          f_xiA6(x_exp(:,ik)) ...
                         ) * x_exp(:,ik); 
    
    bUI0(:,ik) = blkdiag( f_xiB1(x_exp(:,ik)), ...
                          f_xiB2(x_exp(:,ik)), ...
                          f_xiB3(x_exp(:,ik)) ...
                         ) * u_exp(:,ik); 
end

% DATATEST = [EA1, EA2, EA3, EA4, EA5, EA6]*bXI0 + [EB1, EB2, EB3, EB4]*bUI0 - bXI1;

bWI0 = [bXI0; bUI0];
EAEB = bXI1*bWI0'/(bWI0*bWI0');
EA = EAEB(:,1:dim_xiA);
EB = EAEB(:,dim_xiA+1:dim_xiA+dim_xiB);



dataA = bWI0*bWI0';
dataB = -bXI1*bWI0';
dataC = bXI1*bXI1' - eye(dim_x)*0.0001;

disp("min eig is " + min(eig(dataA)));



%% solve the controller ===================================================

% set the domain: |x| <= rx
rx = 0.1; 

% XIA
xiA11 = 1;

xiA12 = 1;

xiA13 = 1;

xiA14 = 1;
xiA24_u = rx; xiA24_l = -rx; 

xiA15 = 1;
xiA25_u = sin(rx)*tan(rx) ; xiA25_l = -sin(rx)*tan(rx);
xiA35_u = 1; xiA35_l = cos(rx);
xiA45_u = sin(rx)/cos(rx); xiA45_l = -sin(rx)/cos(rx);
xiA55_u = rx; xiA55_l = -rx;

xiA16 = 1;
xiA26_u = tan(rx) ; xiA26_l = -tan(rx);
xiA36_u = sin(rx) ; xiA36_l = -sin(rx);
xiA46_u = 1/cos(rx) ; xiA46_l = cos(rx);
xiA56_u = rx; xiA56_l = -rx;

% XIB
xiB11 = 1;
xiB12 = 1;
xiB13 = 1;

% the set of barQ
barQ = funcBarQ(xiA11, xiA12, xiA13, ...
                xiA14, xiA24_u, xiA24_l, ...
                xiA15, xiA25_u, xiA25_l, xiA35_u, xiA35_l, xiA45_u, xiA45_l, xiA55_u, xiA55_l, ...
                xiA16, xiA26_u, xiA26_l, xiA36_u, xiA36_l, xiA46_u, xiA46_l, xiA56_u, xiA56_l, ...
                xiB11, xiB12, xiB13);


% solve sdp LMI
cvx_begin sdp % quiet
    variable epGa semidefinite;
    variable lyGa(dim_x,dim_x) semidefinite;
    variable kY(dim_u,dim_x);
    
    % minimize( 0.3*lambda_max(lyGa) - lambda_min(lyGa) - 0.1*epGa - 0*trace(lyGa) );
    minimize( 1.0*lambda_max(lyGa) - 0.0*lambda_min(lyGa) - 0.0*epGa - 0.1*trace(lyGa) );

    subject to   
        for iv = 1:length(barQ)
            blockS = [ -lyGa-dataC,   zeros(dim_x),          dataB;
                       zeros(dim_x), -lyGa+epGa*eye(dim_x), -[lyGa;kY]'*barQ{iv}'; 
                       dataB',       -barQ{iv}*[lyGa;kY],   -dataA ];
            blockS <= 0;
        end
        epGa >= 1e-12;
        lyGa >= 1e-12 * eye(dim_x);
        
cvx_end

vK = kY/lyGa;

disp(vK);

% region of attraction
roac = sqrt( maxEig(lyGa)*minEig(lyGa) / ( maxEig(lyGa)^2 - epGa*minEig(lyGa) ) );
roa  = min( 1, roac ) * rx;
disp("radius of ROA is " + roa);




















%% useful functions =======================================================

function BarQ = funcBarQ(a1, a2, a3, ...
                         a14, a24_u, a24_l, ...
                         a15, a25_u, a25_l, a35_u, a35_l, a45_u, a45_l, a55_u, a55_l, ...
                         a16, a26_u, a26_l, a36_u, a36_l, a46_u, a46_l, a56_u, a56_l, ...
                         b1, b2, b3)
    % Generate all 512 (2^9) combinations of upper/lower variants
    % Inputs: 
    %   a* : Component matrices
    %   *_u : Upper bound variants
    %   *_l : Lower bound variants
    % Output:
    %   BarQ : 512x1 cell array of blkdiag matrices
    
    % Define all possible combinations (binary combinations for 9 parameters)
    combinations = dec2bin(0:511) - '0'; % 512x9 binary matrix
    
    % Preallocate output
    BarQ = cell(512, 1);
    
    % Generate each combination
    for i = 1:512
        % Select variants based on binary flags
        idx = combinations(i,:);
        
        % Choose a24 variant
        a24 = idx(1)*a24_u + (1-idx(1))*a24_l;
        
        % Choose a25/a35/a45/a55 variants
        a25 = idx(2)*a25_u + (1-idx(2))*a25_l;
        a35 = idx(3)*a35_u + (1-idx(3))*a35_l;
        a45 = idx(4)*a45_u + (1-idx(4))*a45_l;
        a55 = idx(5)*a55_u + (1-idx(5))*a55_l;
        
        % Choose a26/a36/a46/a56 variants
        a26 = idx(6)*a26_u + (1-idx(6))*a26_l;
        a36 = idx(7)*a36_u + (1-idx(7))*a36_l;
        a46 = idx(8)*a46_u + (1-idx(8))*a46_l;
        a56 = idx(9)*a56_u + (1-idx(9))*a56_l;
        
        % Build the block diagonal matrix
        BarQ{i} = blkdiag(a1, a2, a3, ...
                          [a14; a24], ...
                          [a15; a25; a35; a45; a55], ...
                          [a16; a26; a36; a46; a56], ...
                          b1, b2, b3);
    end
end


function e_max = maxEig( MatrixSym )
    e_max = max(eig(MatrixSym));
end

function e_min = minEig( MatrixSym )
    e_min = min(eig(MatrixSym));
end