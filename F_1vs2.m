close all
clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: This script runs simulations to show the difference in
% perfomance for the regular Fisher method and the periodic Fisher method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
iters = 100;                    % the number of iterations for the Fisher's scoring (second method)
step_size_white = 1;
step_size_colored = 1;
gamma = 0.95;
theta_0 = pi/3;                 % starting estimate at 45 deg

[a_model, da_model] = model(rm, K_1, K_3, w);

%% Estimation 

Tests = 50;

% alphas = (1:Tests)/Tests * pi/2 - 1/(2*Tests);
% alpha_0 = alpha;

RMSPE_F_1_colored_colored_b = [];
CyclicErr_F_1_colored_colored_b = [];
MCE_F_1_colored_colored_b = [];


RMSPE_F_1_colored_white_b = [];
CyclicErr_F_1_colored_white_b = [];
MCE_F_1_colored_white_b = [];


RMSPE_F_1_white_white_b = [];
CyclicErr_F_1_white_white_b = [];
MCE_F_1_white_white_b = [];


RMSPE_F_1_white_colored_b = [];
CyclicErr_F_1_white_colored_b = [];
MCE_F_1_white_colored_b = [];


RMSPE_F_2_colored_colored_b = [];
CyclicErr_F_2_colored_colored_b = [];
MCE_F_2_colored_colored_b = [];

RMSPE_F_2_colored_white_b = [];
CyclicErr_F_2_colored_white_b = [];
MCE_F_2_colored_white_b = [];

RMSPE_F_2_white_white_b = [];
CyclicErr_F_2_white_white_b = [];
MCE_F_2_white_white_b = [];

RMSPE_F_2_white_colored_b = [];
CyclicErr_F_2_white_colored_b = [];
MCE_F_2_white_colored_b = [];

RMSPE_F_3_colored_colored_b = [];
CyclicErr_F_3_colored_colored_b = [];
MCE_F_3_colored_colored_b = [];

RMSPE_F_3_colored_white_b = [];
CyclicErr_F_3_colored_white_b = [];
MCE_F_3_colored_white_b = [];

RMSPE_F_3_white_white_b = [];
CyclicErr_F_3_white_white_b = [];
MCE_F_3_white_white_b = [];

RMSPE_F_3_white_colored_b = [];
CyclicErr_F_3_white_colored_b = [];
MCE_F_3_white_colored_b = [];

%--------------------------------------------------------------------------

CRB_white_reg_b = [];
CRB_colored_reg_b = [];
CRB_white_cyc1_b = [];
CRB_colored_cyc1_b = [];
CRB_white_cyc2_b = [];
CRB_colored_cyc2_b = [];

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

Iter_F_1_colored_colored_b = [];
Iter_F_1_colored_white_b = [];
Iter_F_1_white_white_b = [];
Iter_F_1_white_colored_b = [];

Iter_F_2_colored_colored_b = [];
Iter_F_2_colored_white_b = [];
Iter_F_2_white_white_b = [];
Iter_F_2_white_colored_b = [];

Iter_F_3_colored_colored_b = [];
Iter_F_3_colored_white_b = [];
Iter_F_3_white_white_b = [];
Iter_F_3_white_colored_b = [];

J = 5;

for j=1:J

RMSPE_F_1_colored_colored = [];
CyclicErr_F_1_colored_colored = [];
MCE_F_1_colored_colored = [];
ThetaEst_F_1_colored_colored = [];
Iter_F_1_colored_colored = [];


RMSPE_F_1_colored_white = [];
CyclicErr_F_1_colored_white = [];
MCE_F_1_colored_white = [];
ThetaEst_F_1_colored_white = [];
Iter_F_1_colored_white = [];

RMSPE_F_1_white_white = [];
CyclicErr_F_1_white_white = [];
MCE_F_1_white_white = [];
ThetaEst_F_1_white_white = [];
Iter_F_1_white_white = [];

RMSPE_F_1_white_colored = [];
CyclicErr_F_1_white_colored = [];
MCE_F_1_white_colored = [];
ThetaEst_F_1_white_colored = [];
Iter_F_1_white_colored = [];


RMSPE_F_2_colored_colored = [];
CyclicErr_F_2_colored_colored = [];
MCE_F_2_colored_colored = [];
ThetaEst_F_2_colored_colored = [];
Iter_F_2_colored_colored = [];


RMSPE_F_2_colored_white = [];
CyclicErr_F_2_colored_white = [];
MCE_F_2_colored_white = [];
ThetaEst_F_2_colored_white = [];
Iter_F_2_colored_white = [];

RMSPE_F_2_white_white = [];
CyclicErr_F_2_white_white = [];
MCE_F_2_white_white = [];
ThetaEst_F_2_white_white = [];
Iter_F_2_white_white = [];


RMSPE_F_2_white_colored = [];
CyclicErr_F_2_white_colored = [];
MCE_F_2_white_colored = [];
ThetaEst_F_2_white_colored = [];
Iter_F_2_white_colored = [];


RMSPE_F_3_colored_colored = [];
CyclicErr_F_3_colored_colored = [];
MCE_F_3_colored_colored = [];
ThetaEst_F_3_colored_colored = [];
Iter_F_3_colored_colored = [];

RMSPE_F_3_colored_white = [];
CyclicErr_F_3_colored_white = [];
MCE_F_3_colored_white = [];
ThetaEst_F_3_colored_white = [];
Iter_F_3_colored_white = [];

RMSPE_F_3_white_white = [];
CyclicErr_F_3_white_white = [];
MCE_F_3_white_white = [];
ThetaEst_F_3_white_white = [];
Iter_F_3_white_white = [];

RMSPE_F_3_white_colored = [];
CyclicErr_F_3_white_colored = [];
MCE_F_3_white_colored = [];
ThetaEst_F_3_white_colored = [];
Iter_F_3_white_colored = [];

%--------------------------------------------------------------------------

CRB_white_reg = [];
CRB_colored_reg = [];
CRB_white_cyc1 = [];
CRB_colored_cyc1 = [];
CRB_white_cyc2 = [];
CRB_colored_cyc2 = [];


% theta_og = -pi +2*pi/Tests : 2*pi/Tests : pi;
epsilon = 0.01;

SNR = logspace(-2,2, Tests);
theta_og = pi / 6;

tic;
parfor i = 1 : Tests

    theta = theta_og ;        % theta is in [-pi,pi]
    sigma_source = SNR(i) * sigma_noise;

    [X_colored,~,Rv_colored] = synData(rm, theta, alpha, v_0, sigma_source, sigma_noise, M, 'colored', w, K_1, K_3, P);
    [X_white,~,Rv_white] = synData(rm, theta, alpha, v_0, sigma_source, sigma_noise, M, 'white', w, K_1, K_3, P);



    % Fisher version 1 for colored noise-colored
    [Iter_F_1_CC,theta_colored_F_1] = Fisher_1(theta_0,Rv_colored,v_0,alpha,K,X_colored,iters,step_size_colored,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_F_1_colored_colored = [RMSPE_F_1_colored_colored, MSPE(theta, theta_colored_F_1, 'MSPE')];
    CyclicErr_F_1_colored_colored = [CyclicErr_F_1_colored_colored, MSPE(theta, theta_colored_F_1, 'cyclic')];
    MCE_F_1_colored_colored = [MCE_F_1_colored_colored, MSPE(theta, theta_colored_F_1, 'MCE')];
    ThetaEst_F_1_colored_colored = [ThetaEst_F_1_colored_colored, theta_colored_F_1];
    Iter_F_1_colored_colored = [Iter_F_1_colored_colored, Iter_F_1_CC];


    [Iter_F_1_CW,theta_colored_F_1] = Fisher_1(theta_0,Rv_white,v_0,alpha,K,X_colored,iters,step_size_colored,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_F_1_colored_white = [RMSPE_F_1_colored_white, MSPE(theta, theta_colored_F_1, 'MSPE')];
    CyclicErr_F_1_colored_white = [CyclicErr_F_1_colored_white, MSPE(theta, theta_colored_F_1, 'cyclic')];
    MCE_F_1_colored_white = [MCE_F_1_colored_white, MSPE(theta, theta_colored_F_1, 'MCE')];
    ThetaEst_F_1_colored_white = [ThetaEst_F_1_colored_white, theta_colored_F_1];
    Iter_F_1_colored_white = [Iter_F_1_colored_white, Iter_F_1_CW];

    % Fisher version 1 for white noise-white
    [Iter_F_1_WW,theta_white_F_1] = Fisher_1(theta_0,Rv_white,v_0,alpha,K,X_white,iters,step_size_white,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_F_1_white_white = [RMSPE_F_1_white_white, MSPE(theta, theta_white_F_1, 'MSPE')];
    CyclicErr_F_1_white_white = [CyclicErr_F_1_white_white, MSPE(theta, theta_white_F_1, 'cyclic')];
    MCE_F_1_white_white = [MCE_F_1_white_white, MSPE(theta, theta_white_F_1, 'MCE')];
    ThetaEst_F_1_white_white = [ThetaEst_F_1_white_white, theta_white_F_1];
    Iter_F_1_white_white = [Iter_F_1_white_white, Iter_F_1_WW];

    [Iter_F_1_WC,theta_white_F_1] = Fisher_1(theta_0,Rv_colored,v_0,alpha,K,X_white,iters,step_size_white,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_F_1_white_colored = [RMSPE_F_1_white_colored, MSPE(theta, theta_white_F_1, 'MSPE')];
    CyclicErr_F_1_white_colored = [CyclicErr_F_1_white_colored, MSPE(theta, theta_white_F_1, 'cyclic')];
    MCE_F_1_white_colored = [MCE_F_1_white_colored, MSPE(theta, theta_white_F_1, 'MCE')];
    ThetaEst_F_1_white_colored = [ThetaEst_F_1_white_colored, theta_white_F_1];
    Iter_F_1_white_colored = [Iter_F_1_white_colored, Iter_F_1_WC];

    %----------------------------------------------------------------------

    % Fisher version 2 for colored noise-colored
    [Iter_F_2_CC,theta_colored_F_2] = Fisher_2(theta_0,Rv_colored,v_0,alpha,K,X_colored,iters,step_size_colored,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_F_2_colored_colored = [RMSPE_F_2_colored_colored, MSPE(theta, theta_colored_F_2, 'MSPE')];
    CyclicErr_F_2_colored_colored = [CyclicErr_F_2_colored_colored, MSPE(theta, theta_colored_F_2, 'cyclic')];
    MCE_F_2_colored_colored = [MCE_F_2_colored_colored, MSPE(theta, theta_colored_F_2, 'MCE')];
    ThetaEst_F_2_colored_colored = [ThetaEst_F_2_colored_colored, theta_colored_F_2];
    Iter_F_2_colored_colored = [Iter_F_2_colored_colored, Iter_F_2_CC];


    [Iter_F_2_CW,theta_colored_F_2] = Fisher_2(theta_0,Rv_white,v_0,alpha,K,X_colored,iters,step_size_colored,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_F_2_colored_white = [RMSPE_F_2_colored_white, MSPE(theta, theta_colored_F_2, 'MSPE')];
    CyclicErr_F_2_colored_white = [CyclicErr_F_2_colored_white, MSPE(theta, theta_colored_F_2, 'cyclic')];
    MCE_F_2_colored_white = [MCE_F_2_colored_white, MSPE(theta, theta_colored_F_2, 'MCE')];
    ThetaEst_F_2_colored_white = [ThetaEst_F_2_colored_white, theta_colored_F_2];
    Iter_F_2_colored_white = [Iter_F_2_colored_white, Iter_F_2_CW];

    % Fisher version 2 for white noise-white
    [Iter_F_2_WW,theta_white_F_2] = Fisher_2(theta_0,Rv_white,v_0,alpha,K,X_white,iters,step_size_white,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_F_2_white_white = [RMSPE_F_2_white_white, MSPE(theta, theta_white_F_2, 'MSPE')];
    CyclicErr_F_2_white_white = [CyclicErr_F_2_white_white, MSPE(theta, theta_white_F_2, 'cyclic')];
    MCE_F_2_white_white = [MCE_F_2_white_white, MSPE(theta, theta_white_F_2, 'MCE')];
    ThetaEst_F_2_white_white = [ThetaEst_F_2_white_white, theta_white_F_2];
    Iter_F_2_white_white = [Iter_F_2_white_white, Iter_F_2_WW];

    [Iter_F_2_WC,theta_white_F_2] = Fisher_2(theta_0,Rv_colored,v_0,alpha,K,X_white,iters,step_size_white,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_F_2_white_colored = [RMSPE_F_2_white_colored, MSPE(theta, theta_white_F_2, 'MSPE')];
    CyclicErr_F_2_white_colored = [CyclicErr_F_2_white_colored, MSPE(theta, theta_white_F_2, 'cyclic')];
    MCE_F_2_white_colored = [MCE_F_2_white_colored, MSPE(theta, theta_white_F_2, 'MCE')];
    ThetaEst_F_2_white_colored = [ThetaEst_F_2_white_colored, theta_white_F_2];
    Iter_F_2_white_colored = [Iter_F_2_white_colored, Iter_F_2_WC];

    %----------------------------------------------------------------------

    % Fisher version 3 for colored noise-colored
    [Iter_F_3_CC,theta_colored_F_3] = Fisher_3(theta_0,Rv_colored,v_0,alpha,K,X_colored,iters,step_size_colored,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_F_3_colored_colored = [RMSPE_F_3_colored_colored, MSPE(theta, theta_colored_F_3, 'MSPE')];
    CyclicErr_F_3_colored_colored = [CyclicErr_F_3_colored_colored, MSPE(theta, theta_colored_F_3, 'cyclic')];
    MCE_F_3_colored_colored = [MCE_F_3_colored_colored, MSPE(theta, theta_colored_F_3, 'MCE')];
    ThetaEst_F_3_colored_colored = [ThetaEst_F_3_colored_colored, theta_colored_F_3];
    Iter_F_3_colored_colored = [Iter_F_3_colored_colored, Iter_F_3_CC];


    [Iter_F_3_CW,theta_colored_F_3] = Fisher_3(theta_0,Rv_white,v_0,alpha,K,X_colored,iters,step_size_colored,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_F_3_colored_white = [RMSPE_F_3_colored_white, MSPE(theta, theta_colored_F_3, 'MSPE')];
    CyclicErr_F_3_colored_white = [CyclicErr_F_3_colored_white, MSPE(theta, theta_colored_F_3, 'cyclic')];
    MCE_F_3_colored_white = [MCE_F_3_colored_white, MSPE(theta, theta_colored_F_3, 'MCE')];
    ThetaEst_F_3_colored_white = [ThetaEst_F_3_colored_white, theta_colored_F_3];
    Iter_F_3_colored_white = [Iter_F_3_colored_white, Iter_F_3_CW];

    % Fisher version 3 for white noise-white
    [Iter_F_3_WW,theta_white_F_3] = Fisher_3(theta_0,Rv_white,v_0,alpha,K,X_white,iters,step_size_white,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_F_3_white_white = [RMSPE_F_3_white_white, MSPE(theta, theta_white_F_3, 'MSPE')];
    CyclicErr_F_3_white_white = [CyclicErr_F_3_white_white, MSPE(theta, theta_white_F_3, 'cyclic')];
    MCE_F_3_white_white = [MCE_F_3_white_white, MSPE(theta, theta_white_F_3, 'MCE')];
    ThetaEst_F_3_white_white = [ThetaEst_F_3_white_white, theta_white_F_3];
    Iter_F_3_white_white = [Iter_F_3_white_white, Iter_F_3_WW];

    [Iter_F_3_WC,theta_white_F_3] = Fisher_3(theta_0,Rv_colored,v_0,alpha,K,X_white,iters,step_size_white,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_F_3_white_colored = [RMSPE_F_3_white_colored, MSPE(theta, theta_white_F_3, 'MSPE')];
    CyclicErr_F_3_white_colored = [CyclicErr_F_3_white_colored, MSPE(theta, theta_white_F_3, 'cyclic')];
    MCE_F_3_white_colored = [MCE_F_3_white_colored, MSPE(theta, theta_white_F_3, 'MCE')];
    ThetaEst_F_3_white_colored = [ThetaEst_F_3_white_colored, theta_white_F_3];
    Iter_F_3_white_colored = [Iter_F_3_white_colored, Iter_F_3_WC];

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
MCE_F_1_colored_colored_b = [MCE_F_1_colored_colored_b; MCE_F_1_colored_colored];


RMSPE_F_1_colored_white_b = [RMSPE_F_1_colored_white_b; RMSPE_F_1_colored_white];
CyclicErr_F_1_colored_white_b = [CyclicErr_F_1_colored_white_b; CyclicErr_F_1_colored_white];
MCE_F_1_colored_white_b = [MCE_F_1_colored_white_b; MCE_F_1_colored_white];


RMSPE_F_1_white_white_b = [RMSPE_F_1_white_white_b; RMSPE_F_1_white_white];
CyclicErr_F_1_white_white_b = [CyclicErr_F_1_white_white_b; CyclicErr_F_1_white_white];
MCE_F_1_white_white_b = [MCE_F_1_white_white_b; MCE_F_1_white_white];


RMSPE_F_1_white_colored_b = [RMSPE_F_1_white_colored_b; RMSPE_F_1_white_colored];
CyclicErr_F_1_white_colored_b = [CyclicErr_F_1_white_colored_b; CyclicErr_F_1_white_colored];
MCE_F_1_white_colored_b = [MCE_F_1_white_colored_b; MCE_F_1_white_colored];


RMSPE_F_2_colored_colored_b = [RMSPE_F_2_colored_colored_b ; RMSPE_F_2_colored_colored];
CyclicErr_F_2_colored_colored_b = [CyclicErr_F_2_colored_colored_b; CyclicErr_F_2_colored_colored];
MCE_F_2_colored_colored_b = [MCE_F_2_colored_colored_b; MCE_F_2_colored_colored];

RMSPE_F_2_colored_white_b = [RMSPE_F_2_colored_white_b; RMSPE_F_2_colored_white];
CyclicErr_F_2_colored_white_b = [CyclicErr_F_2_colored_white_b; CyclicErr_F_2_colored_white];
MCE_F_2_colored_white_b = [MCE_F_2_colored_white_b; MCE_F_2_colored_white];

RMSPE_F_2_white_white_b = [RMSPE_F_2_white_white_b; RMSPE_F_2_white_white];
CyclicErr_F_2_white_white_b = [CyclicErr_F_2_white_white_b; CyclicErr_F_2_white_white];
MCE_F_2_white_white_b = [MCE_F_2_white_white_b; MCE_F_2_white_white];

RMSPE_F_2_white_colored_b = [RMSPE_F_2_white_colored_b; RMSPE_F_2_white_colored];
CyclicErr_F_2_white_colored_b = [CyclicErr_F_2_white_colored_b; CyclicErr_F_2_white_colored];
MCE_F_2_white_colored_b = [MCE_F_2_white_colored_b; MCE_F_2_white_colored];

RMSPE_F_3_colored_colored_b = [RMSPE_F_3_colored_colored_b ; RMSPE_F_3_colored_colored];
CyclicErr_F_3_colored_colored_b = [CyclicErr_F_3_colored_colored_b; CyclicErr_F_3_colored_colored];
MCE_F_3_colored_colored_b = [MCE_F_3_colored_colored_b; MCE_F_3_colored_colored];

RMSPE_F_3_colored_white_b = [RMSPE_F_3_colored_white_b; RMSPE_F_3_colored_white];
CyclicErr_F_3_colored_white_b = [CyclicErr_F_3_colored_white_b; CyclicErr_F_3_colored_white];
MCE_F_3_colored_white_b = [MCE_F_3_colored_white_b; MCE_F_3_colored_white];

RMSPE_F_3_white_white_b = [RMSPE_F_3_white_white_b; RMSPE_F_3_white_white];
CyclicErr_F_3_white_white_b = [CyclicErr_F_3_white_white_b; CyclicErr_F_3_white_white];
MCE_F_3_white_white_b = [MCE_F_3_white_white_b; MCE_F_3_white_white];

RMSPE_F_3_white_colored_b = [RMSPE_F_3_white_colored_b; RMSPE_F_3_white_colored];
CyclicErr_F_3_white_colored_b = [CyclicErr_F_3_white_colored_b; CyclicErr_F_3_white_colored];
MCE_F_3_white_colored_b = [MCE_F_3_white_colored_b; MCE_F_3_white_colored];

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

%--------------------------------------------------------------------------

Iter_F_1_colored_colored_b = [Iter_F_1_colored_colored_b; Iter_F_1_colored_colored];
Iter_F_1_colored_white_b = [Iter_F_1_colored_white_b; Iter_F_1_colored_white];
Iter_F_1_white_white_b = [Iter_F_1_white_white_b; Iter_F_1_white_white];
Iter_F_1_white_colored_b = [Iter_F_1_white_colored_b; Iter_F_1_white_colored];

Iter_F_2_colored_colored_b = [Iter_F_2_colored_colored_b; Iter_F_2_colored_colored];
Iter_F_2_colored_white_b = [Iter_F_2_colored_white_b; Iter_F_2_colored_white];
Iter_F_2_white_white_b = [Iter_F_2_white_white_b; Iter_F_2_white_white];
Iter_F_2_white_colored_b = [Iter_F_2_white_colored_b; Iter_F_2_white_colored];

Iter_F_3_colored_colored_b = [Iter_F_3_colored_colored_b; Iter_F_3_colored_colored];
Iter_F_3_colored_white_b = [Iter_F_3_colored_white_b; Iter_F_3_colored_white];
Iter_F_3_white_white_b = [Iter_F_3_white_white_b; Iter_F_3_white_white];
Iter_F_3_white_colored_b = [Iter_F_3_white_colored_b; Iter_F_3_white_colored];
end

% theta = pi/2   ;   %%% DONT FORGET TO CHANGE!!! 
%%%% THETA_OG = THETA !!!!

RMSPE_F_1_colored_colored = mean(RMSPE_F_1_colored_colored_b, 1);
CyclicErr_F_1_colored_colored = mean(CyclicErr_F_1_colored_colored_b, 1);
MCE_F_1_colored_colored = mean(MCE_F_1_colored_colored_b, 1);

RMSPE_F_1_colored_white = mean(RMSPE_F_1_colored_white_b, 1);
CyclicErr_F_1_colored_white = mean(CyclicErr_F_1_colored_white_b, 1);
MCE_F_1_colored_white = mean(MCE_F_1_colored_white_b, 1);

RMSPE_F_1_white_white = mean(RMSPE_F_1_white_white_b, 1);
CyclicErr_F_1_white_white = mean(CyclicErr_F_1_white_white_b, 1);
MCE_F_1_white_white = mean(MCE_F_1_white_white_b, 1);

RMSPE_F_1_white_colored = mean(RMSPE_F_1_white_colored_b, 1);
CyclicErr_F_1_white_colored = mean(CyclicErr_F_1_white_colored_b, 1);
MCE_F_1_white_colored = mean(MCE_F_1_white_colored_b, 1);

RMSPE_F_2_colored_colored = mean(RMSPE_F_2_colored_colored_b, 1);
CyclicErr_F_2_colored_colored = mean(CyclicErr_F_2_colored_colored_b, 1);
MCE_F_2_colored_colored = mean(MCE_F_2_colored_colored_b, 1);

RMSPE_F_2_colored_white = mean(RMSPE_F_2_colored_white_b, 1);
CyclicErr_F_2_colored_white = mean(CyclicErr_F_2_colored_white_b, 1);
MCE_F_2_colored_white = mean(MCE_F_2_colored_white_b, 1);

RMSPE_F_2_white_white = mean(RMSPE_F_2_white_white_b, 1);
CyclicErr_F_2_white_white = mean(CyclicErr_F_2_white_white_b, 1);
MCE_F_2_white_white = mean(MCE_F_2_white_white_b, 1);

RMSPE_F_2_white_colored = mean(RMSPE_F_2_white_colored_b, 1);
CyclicErr_F_2_white_colored = mean(CyclicErr_F_2_white_colored_b, 1);
MCE_F_2_white_colored = mean(MCE_F_2_white_colored_b, 1);

RMSPE_F_3_colored_colored = mean(RMSPE_F_3_colored_colored_b, 1);
CyclicErr_F_3_colored_colored = mean(CyclicErr_F_3_colored_colored_b, 1);
MCE_F_3_colored_colored = mean(MCE_F_3_colored_colored_b, 1);

RMSPE_F_3_colored_white = mean(RMSPE_F_3_colored_white_b, 1);
CyclicErr_F_3_colored_white = mean(CyclicErr_F_3_colored_white_b, 1);
MCE_F_3_colored_white = mean(MCE_F_3_colored_white_b, 1);

RMSPE_F_3_white_white = mean(RMSPE_F_3_white_white_b, 1);
CyclicErr_F_3_white_white = mean(CyclicErr_F_3_white_white_b, 1);
MCE_F_3_white_white = mean(MCE_F_3_white_white_b, 1);

RMSPE_F_3_white_colored = mean(RMSPE_F_3_white_colored_b, 1);
CyclicErr_F_3_white_colored = mean(CyclicErr_F_3_white_colored_b, 1);
MCE_F_3_white_colored = mean(MCE_F_3_white_colored_b, 1);

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

ThetaEst_F_1_colored_colored = mean(ThetaEst_F_1_colored_colored_b ,1);
ThetaEst_F_1_colored_white = mean(ThetaEst_F_1_colored_white_b ,1);
ThetaEst_F_1_white_colored = mean(ThetaEst_F_1_white_colored_b ,1);
ThetaEst_F_1_white_white = mean(ThetaEst_F_1_white_white_b ,1);

ThetaEst_F_2_colored_colored = mean(ThetaEst_F_2_colored_colored_b ,1);
ThetaEst_F_2_colored_white = mean(ThetaEst_F_2_colored_white_b ,1);
ThetaEst_F_2_white_colored = mean(ThetaEst_F_2_white_colored_b ,1);
ThetaEst_F_2_white_white = mean(ThetaEst_F_2_white_white_b ,1);

ThetaEst_F_3_colored_colored = mean(ThetaEst_F_3_colored_colored_b ,1);
ThetaEst_F_3_colored_white = mean(ThetaEst_F_3_colored_white_b ,1);
ThetaEst_F_3_white_colored = mean(ThetaEst_F_3_white_colored_b ,1);
ThetaEst_F_3_white_white = mean(ThetaEst_F_3_white_white_b ,1);
%--------------------------------------------------------------------------
Iter_F_1_colored_colored = mean(Iter_F_1_colored_colored_b ,1);
Iter_F_1_colored_white = mean(Iter_F_1_colored_white_b ,1);
Iter_F_1_white_colored = mean(Iter_F_1_white_colored_b ,1);
Iter_F_1_white_white = mean(Iter_F_1_white_white_b ,1);

Iter_F_2_colored_colored = mean(Iter_F_2_colored_colored_b ,1);
Iter_F_2_colored_white = mean(Iter_F_2_colored_white_b ,1);
Iter_F_2_white_colored = mean(Iter_F_2_white_colored_b ,1);
Iter_F_2_white_white = mean(Iter_F_2_white_white_b ,1);

Iter_F_3_colored_colored = mean(Iter_F_3_colored_colored_b ,1);
Iter_F_3_colored_white = mean(Iter_F_3_colored_white_b ,1);
Iter_F_3_white_colored = mean(Iter_F_3_white_colored_b ,1);
Iter_F_3_white_white = mean(Iter_F_3_white_white_b ,1);
%--------------------------------------------------------------------------

res = dir('./res/Fvers_*.mat');
save(append('./res/Fvers_', string(length(res)+1)));
delete(gcp);
clear all;

%--------------------------------------------------------------------------
%% Functions
function [i,theta_final] = Fisher_1(theta_0,Rv,v_0,alpha,K,X,iters,step_size,gamma,M,P,a_f,da_f,acc)
% This function estimates the DOA (theta) using Fisher's scoring
% This version uses the regular CRB and doesn't perform mod 2pi every
% iteration. We limit the angle only after the convergance


theta_estimate = zeros(1,iters);
theta_estimate(1) = theta_0;

pdf = zeros(1,iters);
% s = zeros(1,M);

Rv_inv = pinv(Rv);
X_w = squeeze(sum(X,1));

% a = @(m, theta) a_f(m,theta,alpha,v_0);
% da = @(m, theta) da_f(m,theta,alpha,v_0);

%--------------------------------------------------------------------------
for i = 1 : iters - 1
    
    theta = theta_estimate(i);
    A = cell2mat(arrayfun(a_f,1 : M,theta*ones(1,M),alpha*ones(1,M), v_0*ones(1,M),'uniformoutput',false));
    dA = cell2mat(arrayfun(da_f,1 : M,theta*ones(1,M),alpha*ones(1,M), v_0*ones(1,M),'uniformoutput',false));
    
    Up = diag(A' * Rv_inv * X_w);
    Down = diag(A' * Rv_inv * A);
    s = (Up ./ Down).' / P;
    Mue = A .* s;
    dMue = dA .* s;
    df_theta = sum(2 * real(diag((X_w - Mue)' * Rv_inv * dMue)));
    FIM = sum(2 * real(diag(dMue' * Rv_inv * dMue)));
    crb = FIM ^ (-1);
    pdf(i+1) = log_likelihood(X_w,Rv_inv,M,a_f,v_0,alpha,theta_estimate(i+1),P);

    theta_estimate(i+1) = real(theta_estimate(i) + step_size * crb * df_theta);

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
function [i,theta_final] = Fisher_2(theta_0,Rv,v_0,alpha,K,X,iters,step_size,gamma,M,P,a_f,da_f,acc)
% This function estimates the DOA (theta) using Fisher's scoring
% This version uses the regular CRB and performs mod 2pi every iteration


theta_estimate = zeros(1,iters);
theta_estimate(1) = theta_0;

pdf = zeros(1,iters);
% s = zeros(1,M);

Rv_inv = pinv(Rv);
X_w = squeeze(sum(X,1));

% a = @(m, theta) a_f(m,theta,alpha,v_0);
% da = @(m, theta) da_f(m,theta,alpha,v_0);

%--------------------------------------------------------------------------
for i = 1 : iters - 1
        
    theta = theta_estimate(i);
    A = cell2mat(arrayfun(a_f,1 : M,theta*ones(1,M),alpha*ones(1,M), v_0*ones(1,M),'uniformoutput',false));
    dA = cell2mat(arrayfun(da_f,1 : M,theta*ones(1,M),alpha*ones(1,M), v_0*ones(1,M),'uniformoutput',false));
    
    Up = diag(A' * Rv_inv * X_w);
    Down = diag(A' * Rv_inv * A);
    s = (Up ./ Down).' / P;
    Mue = A .* s;
    dMue = dA .* s;
    df_theta = sum(2 * real(diag((X_w - Mue)' * Rv_inv * dMue)),1);
    FIM = sum(2 * real(diag(dMue' * Rv_inv * dMue)));
    crb = FIM ^ (-1);
    pdf(i+1) = log_likelihood(X_w,Rv_inv,M,a_f,v_0,alpha,theta_estimate(i+1),P);

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
function [i,theta_final] = Fisher_3(theta_0,Rv,v_0,alpha,K,X,iters,step_size,gamma,M,P,a_f,da_f,acc)
% This function estimates the DOA (theta) using Fisher's scoring
% This version usese the cyclic CRB and performs mod 2pi every iteration


theta_estimate = zeros(1,iters);
theta_estimate(1) = theta_0;

pdf = zeros(1,iters);
% s = zeros(1,M);

Rv_inv = pinv(Rv);
X_w = squeeze(sum(X,1));

% a = @(m, theta) a_f(m,theta,alpha,v_0);
% da = @(m, theta) da_f(m,theta,alpha,v_0);

%--------------------------------------------------------------------------
for i = 1 : iters - 1

    theta = theta_estimate(i);
    A = cell2mat(arrayfun(a_f,1 : M,theta*ones(1,M),alpha*ones(1,M), v_0*ones(1,M),'uniformoutput',false));
    dA = cell2mat(arrayfun(da_f,1 : M,theta*ones(1,M),alpha*ones(1,M), v_0*ones(1,M),'uniformoutput',false));

    
    Up = diag(A' * Rv_inv * X_w);
    Down = diag(A' * Rv_inv * A);
    s = (Up ./ Down).' / P;
    Mue = A .* s;
    dMue = dA .* s;
    df_theta = sum(2 * real(diag((X_w - Mue)' * Rv_inv * dMue)),1);
    FIM = sum(2 * real(diag(dMue' * Rv_inv * dMue)));
    crb = FIM ^ (-1);
    ccrb = 2 - 2/(sqrt(crb + 1));

    pdf(i+1) = log_likelihood(X_w,Rv_inv,M,a_f,v_0,alpha,theta_estimate(i+1),P);

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

