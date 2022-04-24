close all
clear all
clc
% this script runs simulations to show the difference in perfomances for
% the regular Fisher method and the periodic Fisher method
%--------------------------------------------------------------------------
%% Initialization

M = 100;                          
P = 25;
K_1 = 15;
K_3 = 1;
K = 3 * K_3 + K_1;

alpha = pi/4 + 0.1;
v_0 = 9600;                     % P waves velocities can be from 6 Km/sec to 11 Km/sec
f = 3;                        
w = f * ones(M,1);              % the omega_m is constant at this point for all the frquencies

SNR = 10;
sigma_source = 10;
sigma_noise = 1;                % for white and colored noise

data = load("data.mat");
rm = data.r_m;                  % r_m is loaded form the data matrix we extrcted from the getData python script

acc = 1e-5;                    % this is the accuracy for the Fisher methods 
iters = 10000;                 % the number of iterations for the Fisher's scoring (second method)
step_size_white = 1;
step_size_colored = 1;
gamma = 0.95;
theta_0 = pi/2;                 % starting estimate at 45 deg

[a_model, da_model] = model(rm, K_1, K_3, w);

%% Estimation 

Tests = 400;

RMSPE_F_1_colored_colored_b = [];
CyclicErr_F_1_colored_colored_b = [];

RMSPE_F_1_colored_white_b = [];
CyclicErr_F_1_colored_white_b = [];

RMSPE_F_1_white_white_b = [];
CyclicErrF_1_white_white_b = [];

RMSPE_F_1_white_colored_b = [];
CyclicErr_F_1_white_colored_b = [];

RMSPE_F_2_colored_colored_b = [];
CyclicErr_F_2_colored_colored_b = [];

RMSPE_F_2_colored_white_b = [];
CyclicErr_F_2_colored_white_b = [];

RMSPE_F_2_white_white_b = [];
CyclicErrF_2_white_white_b = [];

RMSPE_F_2_white_colored_b = [];
CyclicErr_F_2_white_colored_b = [];

RMSPE_F_3_colored_colored_b = [];
CyclicErr_F_3_colored_colored_b = [];

RMSPE_F_3_colored_white_b = [];
CyclicErr_F_3_colored_white_b = [];

RMSPE_F_3_white_white_b = [];
CyclicErrF_3_white_white_b = [];

RMSPE_F_3_white_colored_b = [];
CyclicErr_F_3_white_colored_b = [];
%--------------------------------------------------------------------------

CRB_white_reg_b = [];
CRB_colored_reg_b = [];
CRB_white_cyc1_b = [];
CRB_colored_cyc1_b = [];
CRB_white_cyc2_b = [];
CRB_colored_cyc2_b = [];

%f = waitbar(0,'Please wait...');
J = 400;
for j=1:J
%waitbar(j/(J+1), f, append('iter: ', string(j), ' of ', string(J)));

RMSPE_F_1_colored_colored = [];
CyclicErr_F_1_colored_colored = [];
ThetaEst_F_1_colored_colored = [];

RMSPE_F_1_colored_white = [];
CyclicErr_F_1_colored_white = [];
ThetaEst_F_1_colored_white = [];

RMSPE_F_1_white_white = [];
CyclicErr_F_1_white_white = [];
ThetaEst_F_1_white_white = [];

RMSPE_F_1_white_colored = [];
CyclicErr_F_1_white_colored = [];
ThetaEst_F_1_white_colored = [];

RMSPE_F_2_colored_colored = [];
CyclicErr_F_2_colored_colored = [];
ThetaEst_F_2_colored_colored = [];

RMSPE_F_2_colored_white = [];
CyclicErr_F_2_colored_white = [];
ThetaEst_F_2_colored_white = [];

RMSPE_F_2_white_white = [];
CyclicErr_F_2_white_white = [];
ThetaEst_F_2_white_white = [];

RMSPE_F_2_white_colored = [];
CyclicErr_F_2_white_colored = [];
ThetaEst_F_2_white_colored = [];

RMSPE_F_3_colored_colored = [];
CyclicErr_F_3_colored_colored = [];
ThetaEst_F_3_colored_colored = [];

RMSPE_F_3_colored_white = [];
CyclicErr_F_3_colored_white = [];
ThetaEst_F_3_colored_white = [];

RMSPE_F_3_white_white = [];
CyclicErr_F_3_white_white = [];
ThetaEst_F_3_white_white = [];

RMSPE_F_3_white_colored = [];
CyclicErr_F_3_white_colored = [];
ThetaEst_F_3_white_colored = [];

%--------------------------------------------------------------------------

CRB_white_reg = [];
CRB_colored_reg = [];
CRB_white_cyc1 = [];
CRB_colored_cyc1 = [];
CRB_white_cyc2 = [];
CRB_colored_cyc2 = [];

ThetaEst_F_1_colored_colored_b = [];
ThetaEst_F_1_colored_white_b = [];
ThetaEst_F_1_white_white_b = [];
ThetaEst_F_1_white_colored_b = [];

ThetaEst_F_2_colored_colored_b = [];
ThetaEst_F_2_colored_white_b = [];
ThetaEst_F_2_white_white_b = [];
ThetaEst_F_2_white_colored_b = [];

ThetaEst_F_3_colored_colored_b = [];
ThetaEst_F_3_colored_white_b = [];
ThetaEst_F_3_white_white_b = [];
ThetaEst_F_3_white_colored_b = [];

% theta_og = -pi +2*pi/Tests : 2*pi/Tests : pi;
epsilon = 0.1;

SNR = logspace(-1000, -10, Tests);
theta_og = pi - epsilon;
tic;
parfor i = 1 : Tests

    theta = pi / 5 + epsilon;         % theta is in [-pi,pi]
    sigma_source = SNR(i) * sigma_noise;

    [X_colored,~,Rv_colored,~] = synData(rm, theta, alpha, v_0, sigma_source, sigma_noise, M, 'colored', w, K_1, K_3, P);
    [X_white,~,Rv_white,~] = synData(rm, theta, alpha, v_0, sigma_source, sigma_noise, M, 'white', w, K_1, K_3, P);


    % Fisher version 1 for colored noise-colored
    [~,theta_colored_F_1] = Fisher_1('syn',theta_0,Rv_colored,v_0,alpha,K,X_colored,iters,step_size_colored,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_F_1_colored_colored = [RMSPE_F_1_colored_colored, MSPE(theta, theta_colored_F_1, 'MSPE')];
    CyclicErr_F_1_colored_colored = [CyclicErr_F_1_colored_colored, MSPE(theta, theta_colored_F_1, 'cyclic')];
    ThetaEst_F_1_colored_colored = [ThetaEst_F_1_colored_colored, theta_colored_F_1];


    [~,theta_colored_F_1] = Fisher_1('syn',theta_0,Rv_white,v_0,alpha,K,X_colored,iters,step_size_colored,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_F_1_colored_white = [RMSPE_F_1_colored_white, MSPE(theta, theta_colored_F_1, 'MSPE')];
    CyclicErr_F_1_colored_white = [CyclicErr_F_1_colored_white, MSPE(theta, theta_colored_F_1, 'cyclic')];
    ThetaEst_F_1_colored_white = [ThetaEst_F_1_colored_white, theta_colored_F_1];

    % Fisher version 1 for white noise-white
    [~,theta_white_F_1] = Fisher_1('syn',theta_0,Rv_white,v_0,alpha,K,X_white,iters,step_size_white,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_F_1_white_white = [RMSPE_F_1_white_white, MSPE(theta, theta_white_F_1, 'MSPE')];
    CyclicErr_F_1_white_white = [CyclicErr_F_1_white_white, MSPE(theta, theta_white_F_1, 'cyclic')];
    ThetaEst_F_1_white_white = [ThetaEst_F_1_white_white, theta_white_F_1];

    [~,theta_white_F_1] = Fisher_1('syn',theta_0,Rv_colored,v_0,alpha,K,X_white,iters,step_size_white,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_F_1_white_colored = [RMSPE_F_1_white_colored, MSPE(theta, theta_white_F_1, 'MSPE')];
    CyclicErr_F_1_white_colored = [CyclicErr_F_1_white_colored, MSPE(theta, theta_white_F_1, 'cyclic')];
    ThetaEst_F_1_white_colored = [ThetaEst_F_1_white_colored, theta_white_F_1];

    %----------------------------------------------------------------------

    % Fisher version 2 for colored noise-colored
    [~,theta_colored_F_2] = Fisher_2('syn',theta_0,Rv_colored,v_0,alpha,K,X_colored,iters,step_size_colored,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_F_2_colored_colored = [RMSPE_F_2_colored_colored, MSPE(theta, theta_colored_F_2, 'MSPE')];
    CyclicErr_F_2_colored_colored = [CyclicErr_F_2_colored_colored, MSPE(theta, theta_colored_F_2, 'cyclic')];
    ThetaEst_F_2_colored_colored = [ThetaEst_F_2_colored_colored, theta_colored_F_2];


    [~,theta_colored_F_2] = Fisher_2('syn',theta_0,Rv_white,v_0,alpha,K,X_colored,iters,step_size_colored,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_F_2_colored_white = [RMSPE_F_2_colored_white, MSPE(theta, theta_colored_F_2, 'MSPE')];
    CyclicErr_F_2_colored_white = [CyclicErr_F_2_colored_white, MSPE(theta, theta_colored_F_2, 'cyclic')];
    ThetaEst_F_2_colored_white = [ThetaEst_F_2_colored_white, theta_colored_F_2];

    % Fisher version 2 for white noise-white
    [~,theta_white_F_2] = Fisher_2('syn',theta_0,Rv_white,v_0,alpha,K,X_white,iters,step_size_white,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_F_2_white_white = [RMSPE_F_2_white_white, MSPE(theta, theta_white_F_2, 'MSPE')];
    CyclicErr_F_2_white_white = [CyclicErr_F_2_white_white, MSPE(theta, theta_white_F_2, 'cyclic')];
    ThetaEst_F_2_white_white = [ThetaEst_F_2_white_white, theta_white_F_2];

    [~,theta_white_F_2] = Fisher_2('syn',theta_0,Rv_colored,v_0,alpha,K,X_white,iters,step_size_white,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_F_2_white_colored = [RMSPE_F_2_white_colored, MSPE(theta, theta_white_F_2, 'MSPE')];
    CyclicErr_F_2_white_colored = [CyclicErr_F_2_white_colored, MSPE(theta, theta_white_F_2, 'cyclic')];
    ThetaEst_F_2_white_colored = [ThetaEst_F_2_white_colored, theta_white_F_2];

    %----------------------------------------------------------------------

    % Fisher version 3 for colored noise-colored
    [~,theta_colored_F_3] = Fisher_3('syn',theta_0,Rv_colored,v_0,alpha,K,X_colored,iters,step_size_colored,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_F_3_colored_colored = [RMSPE_F_3_colored_colored, MSPE(theta, theta_colored_F_3, 'MSPE')];
    CyclicErr_F_3_colored_colored = [CyclicErr_F_3_colored_colored, MSPE(theta, theta_colored_F_3, 'cyclic')];
    ThetaEst_F_3_colored_colored = [ThetaEst_F_3_colored_colored, theta_colored_F_3];


    [~,theta_colored_F_3] = Fisher_3('syn',theta_0,Rv_white,v_0,alpha,K,X_colored,iters,step_size_colored,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_F_3_colored_white = [RMSPE_F_3_colored_white, MSPE(theta, theta_colored_F_3, 'MSPE')];
    CyclicErr_F_3_colored_white = [CyclicErr_F_3_colored_white, MSPE(theta, theta_colored_F_3, 'cyclic')];
    ThetaEst_F_3_colored_white = [ThetaEst_F_3_colored_white, theta_colored_F_3];

    % Fisher version 3 for white noise-white
    [~,theta_white_F_3] = Fisher_3('syn',theta_0,Rv_white,v_0,alpha,K,X_white,iters,step_size_white,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_F_3_white_white = [RMSPE_F_3_white_white, MSPE(theta, theta_white_F_3, 'MSPE')];
    CyclicErr_F_3_white_white = [CyclicErr_F_3_white_white, MSPE(theta, theta_white_F_3, 'cyclic')];
    ThetaEst_F_3_white_white = [ThetaEst_F_3_white_white, theta_white_F_3];

    [~,theta_white_F_3] = Fisher_3('syn',theta_0,Rv_colored,v_0,alpha,K,X_white,iters,step_size_white,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_F_3_white_colored = [RMSPE_F_3_white_colored, MSPE(theta, theta_white_F_3, 'MSPE')];
    CyclicErr_F_3_white_colored = [CyclicErr_F_3_white_colored, MSPE(theta, theta_white_F_3, 'cyclic')];
    ThetaEst_F_3_white_colored = [ThetaEst_F_3_white_colored, theta_white_F_3];

    %----------------------------------------------------------------------

    % Calculating the CRB (calculating all kinds)

    CRB_white_reg = [CRB_white_reg, CRB('regular',v_0,alpha,theta,Rv_white,M,X_white,a_model,da_model,P)];
    CRB_colored_reg = [CRB_colored_reg, CRB('regular',v_0,alpha,theta,Rv_colored,M,X_colored,a_model,da_model,P)];

    CRB_white_cyc1 = [CRB_white_cyc1, CRB('cyclic 1',v_0,alpha,theta,Rv_white,M,X_white,a_model,da_model,P)];
    CRB_colored_cyc1 = [CRB_colored_cyc1, CRB('cyclic 1',v_0,alpha,theta,Rv_colored,M,X_colored,a_model,da_model,P)];

    CRB_white_cyc2 = [CRB_white_cyc2, CRB('cyclic 2',v_0,alpha,theta,Rv_white,M,X_white,a_model,da_model,P)];
    CRB_colored_cyc2 = [CRB_colored_cyc2, CRB('cyclic 2',v_0,alpha,theta,Rv_colored,M,X_colored,a_model,da_model,P)];
end
toc;

RMSPE_F_1_colored_colored_b = [RMSPE_F_1_colored_colored_b ; RMSPE_F_1_colored_colored];
CyclicErr_F_1_colored_colored_b = [CyclicErr_F_1_colored_colored_b; CyclicErr_F_1_colored_colored];

RMSPE_F_1_colored_white_b = [RMSPE_F_1_colored_white_b; RMSPE_F_1_colored_white];
CyclicErr_F_1_colored_white_b = [CyclicErr_F_1_colored_white_b; CyclicErr_F_1_colored_white];

RMSPE_F_1_white_white_b = [RMSPE_F_1_white_white_b; RMSPE_F_1_white_white];
CyclicErrF_1_white_white_b = [CyclicErrF_1_white_white_b; CyclicErr_F_1_white_white];

RMSPE_F_1_white_colored_b = [RMSPE_F_1_white_colored_b; RMSPE_F_1_white_colored];
CyclicErr_F_1_white_colored_b = [CyclicErr_F_1_white_colored_b; CyclicErr_F_1_white_colored];

RMSPE_F_2_colored_colored_b = [RMSPE_F_2_colored_colored_b ; RMSPE_F_2_colored_colored];
CyclicErr_F_2_colored_colored_b = [CyclicErr_F_2_colored_colored_b; CyclicErr_F_2_colored_colored];

RMSPE_F_2_colored_white_b = [RMSPE_F_2_colored_white_b; RMSPE_F_2_colored_white];
CyclicErr_F_2_colored_white_b = [CyclicErr_F_2_colored_white_b; CyclicErr_F_2_colored_white];

RMSPE_F_2_white_white_b = [RMSPE_F_2_white_white_b; RMSPE_F_2_white_white];
CyclicErrF_2_white_white_b = [CyclicErrF_2_white_white_b; CyclicErr_F_2_white_white];

RMSPE_F_2_white_colored_b = [RMSPE_F_2_white_colored_b; RMSPE_F_2_white_colored];
CyclicErr_F_2_white_colored_b = [CyclicErr_F_2_white_colored_b; CyclicErr_F_2_white_colored];

RMSPE_F_3_colored_colored_b = [RMSPE_F_3_colored_colored_b ; RMSPE_F_3_colored_colored];
CyclicErr_F_3_colored_colored_b = [CyclicErr_F_3_colored_colored_b; CyclicErr_F_3_colored_colored];

RMSPE_F_3_colored_white_b = [RMSPE_F_3_colored_white_b; RMSPE_F_3_colored_white];
CyclicErr_F_3_colored_white_b = [CyclicErr_F_3_colored_white_b; CyclicErr_F_3_colored_white];

RMSPE_F_3_white_white_b = [RMSPE_F_3_white_white_b; RMSPE_F_3_white_white];
CyclicErrF_3_white_white_b = [CyclicErrF_3_white_white_b; CyclicErr_F_3_white_white];

RMSPE_F_3_white_colored_b = [RMSPE_F_3_white_colored_b; RMSPE_F_3_white_colored];
CyclicErr_F_3_white_colored_b = [CyclicErr_F_3_white_colored_b; CyclicErr_F_3_white_colored];

%--------------------------------------------------------------------------

CRB_white_reg_b = [CRB_white_reg_b; CRB_white_reg];
CRB_colored_reg_b = [CRB_colored_reg_b; CRB_colored_reg];
CRB_white_cyc1_b = [CRB_white_cyc1_b; CRB_white_cyc1];
CRB_colored_cyc1_b = [CRB_colored_cyc1_b; CRB_colored_cyc1];
CRB_white_cyc2_b = [CRB_white_cyc2_b; CRB_white_cyc2];
CRB_colored_cyc2_b = [CRB_colored_cyc2_b; CRB_colored_cyc2];

ThetaEst_F_1_colored_colored_b = [ThetaEst_F_1_colored_colored_b; ThetaEst_F_1_colored_colored];
ThetaEst_F_1_colored_white_b = [ThetaEst_F_1_colored_white_b; ThetaEst_F_1_colored_white];
ThetaEst_F_1_white_white_b = [ThetaEst_F_1_white_white_b; ThetaEst_F_1_white_white];
ThetaEst_F_1_white_colored_b = [ThetaEst_F_1_white_colored_b; ThetaEst_F_1_white_colored];

ThetaEst_F_2_colored_colored_b = [ThetaEst_F_2_colored_colored_b; ThetaEst_F_2_colored_colored];
ThetaEst_F_2_colored_white_b = [ThetaEst_F_2_colored_white_b; ThetaEst_F_2_colored_white];
ThetaEst_F_2_white_white_b = [ThetaEst_F_2_white_white_b; ThetaEst_F_2_white_white];
ThetaEst_F_2_white_colored_b = [ThetaEst_F_2_white_colored_b; ThetaEst_F_2_white_colored];

ThetaEst_F_3_colored_colored_b = [ThetaEst_F_3_colored_colored_b; ThetaEst_F_3_colored_colored];
ThetaEst_F_3_colored_white_b = [ThetaEst_F_3_colored_white_b; ThetaEst_F_3_colored_white];
ThetaEst_F_3_white_white_b = [ThetaEst_F_3_white_white_b; ThetaEst_F_3_white_white];
ThetaEst_F_3_white_colored_b = [ThetaEst_F_3_white_colored_b; ThetaEst_F_3_white_colored];

end

RMSPE_F_1_colored_colored = mean(RMSPE_F_1_colored_colored_b, 1);
CyclicErr_F_1_colored_colored = mean(CyclicErr_F_1_colored_colored_b, 1);

RMSPE_F_1_colored_white = mean(RMSPE_F_1_colored_white_b, 1);
CyclicErr_F_1_colored_white = mean(CyclicErr_F_1_colored_white_b, 1);

RMSPE_F_1_white_white = mean(RMSPE_F_1_white_white_b, 1);
CyclicErr_F_1_white_white = mean(CyclicErrF_1_white_white_b, 1);

RMSPE_F_1_white_colored = mean(RMSPE_F_1_white_colored_b, 1);
CyclicErr_F_1_white_colored = mean(CyclicErr_F_1_white_colored_b, 1);

RMSPE_F_2_colored_colored = mean(RMSPE_F_2_colored_colored_b, 1);
CyclicErr_F_2_colored_colored = mean(CyclicErr_F_2_colored_colored_b, 1);

RMSPE_F_2_colored_white = mean(RMSPE_F_2_colored_white_b, 1);
CyclicErr_F_2_colored_white = mean(CyclicErr_F_2_colored_white_b, 1);

RMSPE_F_2_white_white = mean(RMSPE_F_2_white_white_b, 1);
CyclicErr_F_2_white_white = mean(CyclicErrF_2_white_white_b, 1);

RMSPE_F_2_white_colored = mean(RMSPE_F_2_white_colored_b, 1);
CyclicErr_F_2_white_colored = mean(CyclicErr_F_2_white_colored_b, 1);

RMSPE_F_3_colored_colored = mean(RMSPE_F_3_colored_colored_b, 1);
CyclicErr_F_3_colored_colored = mean(CyclicErr_F_3_colored_colored_b, 1);

RMSPE_F_3_colored_white = mean(RMSPE_F_3_colored_white_b, 1);
CyclicErr_F_3_colored_white = mean(CyclicErr_F_3_colored_white_b, 1);

RMSPE_F_3_white_white = mean(RMSPE_F_3_white_white_b, 1);
CyclicErr_F_3_white_white = mean(CyclicErrF_3_white_white_b, 1);

RMSPE_F_3_white_colored = mean(RMSPE_F_3_white_colored_b, 1);
CyclicErr_F_3_white_colored = mean(CyclicErr_F_3_white_colored_b, 1);

%--------------------------------------------------------------------------
CRB_white_reg = mean(CRB_white_reg_b ,1);
CRB_colored_reg = mean(CRB_colored_reg_b ,1);
CRB_white_cyc1 = mean(CRB_white_cyc1_b ,1);
CRB_colored_cyc1 = mean(CRB_colored_cyc1_b ,1);
CRB_white_cyc2 = mean(CRB_white_cyc2_b ,1);
CRB_colored_cyc2 = mean(CRB_colored_cyc2_b ,1);

ThetaEst_F_1_colored_colored_mean = mean(wrapToPi(ThetaEst_F_1_colored_colored_b - theta_og), 1);
ThetaEst_F_1_colored_white_mean = mean(wrapToPi(ThetaEst_F_1_colored_white_b - theta_og), 1);
ThetaEst_F_1_white_white_mean = mean(wrapToPi(ThetaEst_F_1_white_white_b - theta_og), 1);
ThetaEst_F_1_white_colored_mean = mean(wrapToPi(ThetaEst_F_1_white_colored_b - theta_og), 1);

ThetaEst_F_2_colored_colored_mean = mean(wrapToPi(ThetaEst_F_2_colored_colored_b - theta_og), 1);
ThetaEst_F_2_colored_white_mean = mean(wrapToPi(ThetaEst_F_2_colored_white_b - theta_og), 1);
ThetaEst_F_2_white_white_mean = mean(wrapToPi(ThetaEst_F_2_white_white_b - theta_og), 1);
ThetaEst_F_2_white_colored_mean = mean(wrapToPi(ThetaEst_F_2_white_colored_b - theta_og), 1);

ThetaEst_F_3_colored_colored_mean = mean(wrapToPi(ThetaEst_F_3_colored_colored_b - theta_og), 1);
ThetaEst_F_3_colored_white_mean = mean(wrapToPi(ThetaEst_F_3_colored_white_b - theta_og), 1);
ThetaEst_F_3_white_white_mean = mean(wrapToPi(ThetaEst_F_3_white_white_b - theta_og), 1);
ThetaEst_F_3_white_colored_mean = mean(wrapToPi(ThetaEst_F_3_white_colored_b - theta_og), 1);

%--------------------------------------------------------------------------

res = dir('./res/F_regVSper_results_*.mat');
save(append('./res/F_regVSper_results_', string(length(res)+1)));
delete(gcp);
clear all;

%% Circle
% angles = linspace(0, 2*pi, 500);
% radius = 20;
% CenterX = 50;
% CenterY = 40;
% x = radius * cos(angles) + CenterX;
% y = radius * sin(angles) + CenterY;
% 
% x_og = radius * cos(theta_og) + CenterX;
% x_MLE_white = radius * cos(ThetaEst_MLE_white_white) + CenterX;
% x_MLE_colored = radius * cos(ThetaEst_MLE_colored_colored) + CenterX;
% x_Fisher_white = radius * cos(ThetaEst_fisher_white_white) + CenterX;
% x_Fisher_colored = radius * cos(ThetaEst_fisher_colored_colored) + CenterX;
% 
% y_og = radius * sin(theta_og) + CenterY;
% y_MLE_white = radius * sin(ThetaEst_MLE_white_white) + CenterY;
% y_MLE_colored = radius * sin(ThetaEst_MLE_colored_colored) + CenterY;
% y_Fisher_white = radius * sin(ThetaEst_fisher_white_white) + CenterY;
% y_Fisher_colored = radius * sin(ThetaEst_fisher_colored_colored) + CenterY;
% 
% 
% figure
% hold on; grid on ; axis equal;
% plot(x, y, 'k-', 'LineWidth', 1);
% plot(x_og, y_og, 'r*','LineWidth', 2.5)
% plot(x_MLE_white, y_MLE_white, 'o','LineWidth', 1.5)
% plot(x_MLE_colored, y_MLE_colored, 'o','LineWidth', 1.5)
% plot(x_Fisher_white, y_Fisher_white, 'x','LineWidth', 1.5)
% plot(x_Fisher_colored, y_Fisher_colored, 'x','LineWidth', 1.5)
% xlabel('X', 'FontSize', 10); ylabel('Y', 'FontSize', 10);
% legend('','Original', 'MLE W-W', 'MLE C-C', 'Fisher W-W', 'Fisher C-C')
% title('Comparison of the Estimations on Circle')


%--------------------------------------------------------------------------
%% Functions
function [theta_estimate,theta_final] = Fisher_1(type_sim,theta_0,Rv,v_0,alpha,K,X,iters,step_size,gamma,M,P,a_f,da_f,acc)
% This function estimates the DOA (theta) using Fisher's scoring
% This version uses the regular CRB and doesn't perform mod 2pi every
% iteration. We limit the angle only after the convergance

if (strcmp(type_sim , 'realData'))
    mcdRv = mcdcov(X.','cor', 1, 'plots', 0);
    Rv = mcdRv.cov;
    X_w = zeros(P, K, M); %X_w(p,:,m) to address the pth segment and mth frequency
    for i = 0:P-1
        X_w(i+1,:,:) = fft(X(:,i*L+1:min(N, (i+1)*L)), L, 2);
    end
    X = X_w;
end

theta_estimate = zeros(1,iters);
theta_estimate(1) = theta_0;

pdf = zeros(1,iters);
s = zeros(1,M);

Rv_inv = pinv(Rv);
X_w = squeeze(sum(X,1));

%--------------------------------------------------------------------------
for i = 1 : iters - 1
    
    FIM = 0;
    df_theta = 0;

    for j = 1 : M
           
        x_i = X_w(:,j);
        % calculate the a and a derivative 
        a = a_f(j, theta_estimate(i), alpha, v_0);
        da = da_f(j, theta_estimate(i), alpha, v_0);

        % estimating s(m) using the MLE of the deterministic signal
        down = a' * Rv_inv * a;   
        s_j = a' * Rv_inv * x_i;   % per the formula only x is a variable of p
        s(j) = s_j / down;

        mue = a * s(j);
        dmue = da * s(j);
        

        df_theta = df_theta + (x_i - mue)' * Rv_inv * dmue + (dmue' * Rv_inv * (x_i - mue)).';
        FIM = FIM + real(dmue' * Rv_inv * dmue);

    end

    FIM = FIM * 2;
    crb = FIM ^ (-1);
    pdf(i+1) = log_likelihood(X_w,s,Rv_inv,M,a_f,v_0,alpha,theta_estimate(i+1));

    theta_estimate(i+1) = (real(theta_estimate(i) + step_size * crb * df_theta));

    if (pdf(i+1) < pdf(i))
        step_size = step_size * gamma;
    end

    if (mod(abs(theta_estimate(i+1) - theta_estimate(i)),2*pi) < acc)
        theta_estimate(end) = theta_estimate(i+1);
        break;
    end
end
theta_final = wrapToPi(theta_estimate(end));
end

%--------------------------------------------------------------------------
function [theta_estimate,theta_final] = Fisher_2(type_sim,theta_0,Rv,v_0,alpha,K,X,iters,step_size,gamma,M,P,a_f,da_f,acc)
% This function estimates the DOA (theta) using Fisher's scoring
% This version uses the regular CRB and performs mod 2pi every iteration
if (strcmp(type_sim , 'realData'))
    mcdRv = mcdcov(X.','cor', 1, 'plots', 0);
    Rv = mcdRv.cov;
    X_w = zeros(P, K, M); %X_w(p,:,m) to address the pth segment and mth frequency
    for i = 0:P-1
        X_w(i+1,:,:) = fft(X(:,i*L+1:min(N, (i+1)*L)), L, 2);
    end
    X = X_w;
end

theta_estimate = zeros(1,iters);
theta_estimate(1) = theta_0;

pdf = zeros(1,iters);
s = zeros(1,M);

Rv_inv = pinv(Rv);
X_w = squeeze(sum(X,1));

%--------------------------------------------------------------------------
for i = 1 : iters - 1
    
    FIM = 0;
    df_theta = 0;

    for j = 1 : M
           
        x_i = X_w(:,j);
        % calculate the a and a derivative 
        a = a_f(j, theta_estimate(i), alpha, v_0);
        da = da_f(j, theta_estimate(i), alpha, v_0);

        % estimating s(m) using the MLE of the deterministic signal
        down = a' * Rv_inv * a;   
        s_j = a' * Rv_inv * x_i;   % per the formula only x is a variable of p
        s(j) = s_j / down;

        mue = a * s(j);
        dmue = da * s(j);
        

        df_theta = df_theta + (x_i - mue)' * Rv_inv * dmue + (dmue' * Rv_inv * (x_i - mue)).';
        FIM = FIM + real(dmue' * Rv_inv * dmue);

    end

    FIM = FIM * 2;
    crb = FIM ^ (-1);
    pdf(i+1) = log_likelihood(X_w,s,Rv_inv,M,a_f,v_0,alpha,theta_estimate(i+1));

    theta_estimate(i+1) = wrapToPi(real(theta_estimate(i) + step_size * crb * df_theta));

    if (pdf(i+1) < pdf(i))
        step_size = step_size * gamma;
    end

    if (mod(abs(theta_estimate(i+1) - theta_estimate(i)),2*pi) < acc)
        theta_estimate(end) = theta_estimate(i+1);
        break;
    end
end
theta_final = theta_estimate(end);
end

%--------------------------------------------------------------------------
function [theta_estimate,theta_final] = Fisher_3(type_sim,theta_0,Rv,v_0,alpha,K,X,iters,step_size,gamma,M,P,a_f,da_f,acc)
% This function estimates the DOA (theta) using Fisher's scoring
% This version usese the cyclic CRB and performs mod 2pi every iteration

if (strcmp(type_sim , 'realData'))
    mcdRv = mcdcov(X.','cor', 1, 'plots', 0);
    Rv = mcdRv.cov;
    X_w = zeros(P, K, M); %X_w(p,:,m) to address the pth segment and mth frequency
    for i = 0:P-1
        X_w(i+1,:,:) = fft(X(:,i*L+1:min(N, (i+1)*L)), L, 2);
    end
    X = X_w;
end

theta_estimate = zeros(1,iters);
theta_estimate(1) = theta_0;

pdf = zeros(1,iters);
s = zeros(1,M);

Rv_inv = pinv(Rv);
X_w = squeeze(sum(X,1));

%--------------------------------------------------------------------------
for i = 1 : iters - 1
    
    FIM = 0;
    df_theta = 0;

    for j = 1 : M
           
        x_i = X_w(:,j);
        % calculate the a and a derivative 
        a = a_f(j, theta_estimate(i), alpha, v_0);
        da = da_f(j, theta_estimate(i), alpha, v_0);

        % estimating s(m) using the MLE of the deterministic signal
        down = a' * Rv_inv * a;   
        s_j = a' * Rv_inv * x_i;   % per the formula only x is a variable of p
        s(j) = s_j / down;

        mue = a * s(j);
        dmue = da * s(j);
        

        df_theta = df_theta + (x_i - mue)' * Rv_inv * dmue + (dmue' * Rv_inv * (x_i - mue)).';
        FIM = FIM + real(dmue' * Rv_inv * dmue);

    end

    FIM = FIM * 2;
    crb = FIM ^ (-1);
    ccrb = 2 - 2/(sqrt(crb + 1));
    pdf(i+1) = log_likelihood(X_w,s,Rv_inv,M,a_f,v_0,alpha,theta_estimate(i+1));

    theta_estimate(i+1) = wrapToPi(real(theta_estimate(i) + step_size * ccrb * df_theta));

    if (pdf(i+1) < pdf(i))
        step_size = step_size * gamma;
    end

    if (mod(abs(theta_estimate(i+1) - theta_estimate(i)),2*pi) < acc)
        theta_estimate(end) = theta_estimate(i+1);
        break;
    end
end
theta_final = theta_estimate(end);
end
