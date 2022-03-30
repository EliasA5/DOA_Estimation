close all
clc
%--------------------------------------------------------------------------
% This scripts performs the comparison between the MLE and Fisher's scoring
% for white noise and colored noise
% All signal's are in the frequency domain
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

acc = 0.01;                    % this is the accuracy for the MLE sweep in method 1
iters = 1e3;                    % the number of iterations for the Fisher's scoring (second method)
step_size_white = 1;
step_size_colored = 1;
gamma = 0.95;
theta_0 = pi/4-1e-1;                 % starting estimate at 45 deg

[a_model, da_model] = model(rm, K_1, K_3, w);

%% Estimation 

Tests = 102;

RMSPE_MLE_colored_colored_b = [];
CyclicErr_MLE_colored_colored_b = [];

RMSPE_MLE_colored_white_b = [];
CyclicErr_MLE_colored_white_b = [];

RMSPE_MLE_white_white_b = [];
CyclicErr_MLE_white_white_b = [];

RMSPE_MLE_white_colored_b = [];
CyclicErr_MLE_white_colored_b = [];

RMSPE_fisher_colored_colored_b = [];
CyclicErr_fisher_colored_colored_b = [];

RMSPE_fisher_colored_white_b = [];
CyclicErr_fisher_colored_white_b = [];

RMSPE_fisher_white_white_b = [];
CyclicErr_fisher_white_white_b = [];

RMSPE_fisher_white_colored_b = [];
CyclicErr_fisher_white_colored_b = [];

CRB_white_reg_b = [];
CRB_colored_reg_b = [];
CRB_white_cyc1_b = [];
CRB_colored_cyc1_b = [];
CRB_white_cyc2_b = [];
CRB_colored_cyc2_b = [];

%f = waitbar(0,'Please wait...');
J = 20;
for j=1:J
%waitbar(j/(J+1), f, append('iter: ', string(j), ' of ', string(J)));
RMSPE_MLE_colored_colored = [];
CyclicErr_MLE_colored_colored = [];
ThetaEst_MLE_colored_colored = [];

RMSPE_MLE_colored_white = [];
CyclicErr_MLE_colored_white = [];
ThetaEst_MLE_colored_white = [];

RMSPE_MLE_white_white = [];
CyclicErr_MLE_white_white = [];
ThetaEst_MLE_white_white = [];

RMSPE_MLE_white_colored = [];
CyclicErr_MLE_white_colored = [];
ThetaEst_MLE_white_colored = [];

RMSPE_fisher_colored_colored = [];
CyclicErr_fisher_colored_colored = [];
ThetaEst_fisher_colored_colored = [];

RMSPE_fisher_colored_white = [];
CyclicErr_fisher_colored_white = [];
ThetaEst_fisher_colored_white = [];

RMSPE_fisher_white_white = [];
CyclicErr_fisher_white_white = [];
ThetaEst_fisher_white_white = [];

RMSPE_fisher_white_colored = [];
CyclicErr_fisher_white_colored = [];
ThetaEst_fisher_white_colored = [];

CRB_white_reg = [];
CRB_colored_reg = [];
CRB_white_cyc1 = [];
CRB_colored_cyc1 = [];
CRB_white_cyc2 = [];
CRB_colored_cyc2 = [];

% theta_og = -pi +2*pi/Tests : 2*pi/Tests : pi;
epsilon = 0.1;

SNR = logspace(-2, 3, Tests);
theta_og = pi / 5 + epsilon;
tic;
parfor i = 1 : Tests

    theta = pi / 5 + epsilon;         % theta is in [-pi,pi]
    sigma_source = SNR(i) * sigma_noise;

    [X_colored,~,Rv_colored,~] = synData(rm, theta, alpha, v_0, sigma_source, sigma_noise, M, 'colored', w, K_1, K_3, P);
    [X_white,~,Rv_white,~] = synData(rm, theta, alpha, v_0, sigma_source, sigma_noise, M, 'white', w, K_1, K_3, P);


    % MLE for colored noise-colored
    fun_colored = toMaximizeMLE(a_model, Rv_colored, X_colored, M, P);
    theta_colored_MLE = real(MaximizeTheta(fun_colored, alpha, v_0, acc));
    RMSPE_MLE_colored_colored = [RMSPE_MLE_colored_colored, MSPE(theta, theta_colored_MLE, 'MSPE')];
    CyclicErr_MLE_colored_colored = [CyclicErr_MLE_colored_colored, MSPE(theta, theta_colored_MLE, 'cyclic')];
    ThetaEst_MLE_colored_colored = [ThetaEst_MLE_colored_colored, theta_colored_MLE];


    fun_colored = toMaximizeMLE(a_model, Rv_white, X_colored, M, P);
    theta_colored_MLE = real(MaximizeTheta(fun_colored, alpha, v_0, acc));
    RMSPE_MLE_colored_white = [RMSPE_MLE_colored_white, MSPE(theta, theta_colored_MLE, 'MSPE')];
    CyclicErr_MLE_colored_white = [CyclicErr_MLE_colored_white, MSPE(theta, theta_colored_MLE, 'cyclic')];
    ThetaEst_MLE_colored_white = [ThetaEst_MLE_colored_white, theta_colored_MLE];

    % MLE for white noise-white
    fun_white = toMaximizeMLE(a_model, Rv_white, X_white, M, P);
    theta_white_MLE = real(MaximizeTheta(fun_white, alpha, v_0, acc));
    RMSPE_MLE_white_white = [RMSPE_MLE_white_white, MSPE(theta, theta_white_MLE, 'MSPE')];
    CyclicErr_MLE_white_white = [CyclicErr_MLE_white_white, MSPE(theta, theta_white_MLE, 'cyclic')];
    ThetaEst_MLE_white_white = [ThetaEst_MLE_white_white, theta_white_MLE];

    fun_white = toMaximizeMLE(a_model, Rv_colored, X_white, M, P);
    theta_white_MLE = real(MaximizeTheta(fun_white, alpha, v_0, acc));
    RMSPE_MLE_white_colored = [RMSPE_MLE_white_colored, MSPE(theta, theta_white_MLE, 'MSPE')];
    CyclicErr_MLE_white_colored = [CyclicErr_MLE_white_colored, MSPE(theta, theta_white_MLE, 'cyclic')];
    ThetaEst_MLE_white_colored = [ThetaEst_MLE_white_colored, theta_white_MLE];

    % Fisher's scoring for colored noise-colored
    [~,theta_colored_fisher] = Fisher_scoring('syn',theta_0,Rv_colored,v_0,alpha,K,X_colored,iters,step_size_colored,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_fisher_colored_colored = [RMSPE_fisher_colored_colored, MSPE(theta, theta_colored_fisher, 'MSPE')];
    CyclicErr_fisher_colored_colored = [CyclicErr_fisher_colored_colored, MSPE(theta, theta_colored_fisher, 'cyclic')];
    ThetaEst_fisher_colored_colored = [ThetaEst_fisher_colored_colored, theta_colored_fisher];

    [~,theta_colored_fisher] = Fisher_scoring('syn',theta_0,Rv_white,v_0,alpha,K,X_colored,iters,step_size_colored,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_fisher_colored_white = [RMSPE_fisher_colored_white, MSPE(theta, theta_colored_fisher, 'MSPE')];
    CyclicErr_fisher_colored_white = [CyclicErr_fisher_colored_white, MSPE(theta, theta_colored_fisher, 'cyclic')];
    ThetaEst_fisher_colored_white = [ThetaEst_fisher_colored_white, theta_colored_fisher];
    
    % Fisher's scoring for white noise-white
    [~,theta_white_fisher] = Fisher_scoring('syn',theta_0,Rv_white,v_0,alpha,K,X_white,iters,step_size_white,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_fisher_white_white = [RMSPE_fisher_white_white, MSPE(theta, theta_white_fisher, 'MSPE')];
    CyclicErr_fisher_white_white = [CyclicErr_fisher_white_white, MSPE(theta, theta_white_fisher, 'cyclic')];
    ThetaEst_fisher_white_white = [ThetaEst_fisher_white_white, theta_white_fisher];

    [~,theta_white_fisher] = Fisher_scoring('syn',theta_0,Rv_colored,v_0,alpha,K,X_white,iters,step_size_white,gamma,M,P,a_model,da_model,1e-5);
    RMSPE_fisher_white_colored = [RMSPE_fisher_white_colored, MSPE(theta, theta_white_fisher, 'MSPE')];
    CyclicErr_fisher_white_colored = [CyclicErr_fisher_white_colored, MSPE(theta, theta_white_fisher, 'cyclic')];
    ThetaEst_fisher_white_colored = [ThetaEst_fisher_white_colored, theta_white_fisher];

    % Calculating the CRB (calculating all kinds)

    CRB_white_reg = [CRB_white_reg, CRB('regular',v_0,alpha,theta,Rv_white,M,X_white,a_model,da_model,P)];
    CRB_colored_reg = [CRB_colored_reg, CRB('regular',v_0,alpha,theta,Rv_colored,M,X_colored,a_model,da_model,P)];

    CRB_white_cyc1 = [CRB_white_cyc1, CRB('cyclic 1',v_0,alpha,theta,Rv_white,M,X_white,a_model,da_model,P)];
    CRB_colored_cyc1 = [CRB_colored_cyc1, CRB('cyclic 1',v_0,alpha,theta,Rv_colored,M,X_colored,a_model,da_model,P)];

    CRB_white_cyc2 = [CRB_white_cyc2, CRB('cyclic 2',v_0,alpha,theta,Rv_white,M,X_white,a_model,da_model,P)];
    CRB_colored_cyc2 = [CRB_colored_cyc2, CRB('cyclic 2',v_0,alpha,theta,Rv_colored,M,X_colored,a_model,da_model,P)];
end
toc;
RMSPE_MLE_colored_colored_b = [RMSPE_MLE_colored_colored_b ; RMSPE_MLE_colored_colored];
CyclicErr_MLE_colored_colored_b = [CyclicErr_MLE_colored_colored_b; CyclicErr_MLE_colored_colored];

RMSPE_MLE_colored_white_b = [RMSPE_MLE_colored_white_b; RMSPE_MLE_colored_white];
CyclicErr_MLE_colored_white_b = [CyclicErr_MLE_colored_white_b; CyclicErr_MLE_colored_white];

RMSPE_MLE_white_white_b = [RMSPE_MLE_white_white_b; RMSPE_MLE_white_white];
CyclicErr_MLE_white_white_b = [CyclicErr_MLE_white_white_b; CyclicErr_MLE_white_white];

RMSPE_MLE_white_colored_b = [RMSPE_MLE_white_colored_b; RMSPE_MLE_white_colored];
CyclicErr_MLE_white_colored_b = [CyclicErr_MLE_white_colored_b; CyclicErr_MLE_white_colored];

RMSPE_fisher_colored_colored_b = [RMSPE_fisher_colored_colored_b; RMSPE_fisher_colored_colored];
CyclicErr_fisher_colored_colored_b = [CyclicErr_fisher_colored_colored_b; CyclicErr_fisher_colored_colored];

RMSPE_fisher_colored_white_b = [RMSPE_fisher_colored_white_b; RMSPE_fisher_colored_white];
CyclicErr_fisher_colored_white_b = [CyclicErr_fisher_colored_white_b; CyclicErr_fisher_colored_white];

RMSPE_fisher_white_white_b = [RMSPE_fisher_white_white_b; RMSPE_fisher_white_white];
CyclicErr_fisher_white_white_b = [CyclicErr_fisher_white_white_b; CyclicErr_fisher_white_white];

RMSPE_fisher_white_colored_b = [RMSPE_fisher_white_colored_b; RMSPE_fisher_white_colored];
CyclicErr_fisher_white_colored_b = [CyclicErr_fisher_white_colored_b; CyclicErr_fisher_white_colored];

CRB_white_reg_b = [CRB_white_reg_b; CRB_white_reg];
CRB_colored_reg_b = [CRB_colored_reg_b; CRB_colored_reg];
CRB_white_cyc1_b = [CRB_white_cyc1_b; CRB_white_cyc1];
CRB_colored_cyc1_b = [CRB_colored_cyc1_b; CRB_colored_cyc1];
CRB_white_cyc2_b = [CRB_white_cyc2_b; CRB_white_cyc2];
CRB_colored_cyc2_b = [CRB_colored_cyc2_b; CRB_colored_cyc2];

end

RMSPE_MLE_colored_colored = mean(RMSPE_MLE_colored_colored_b, 1);
CyclicErr_MLE_colored_colored = mean(CyclicErr_MLE_colored_colored_b, 1);

RMSPE_MLE_colored_white = mean(RMSPE_MLE_colored_white_b, 1);
CyclicErr_MLE_colored_white = mean(CyclicErr_MLE_colored_white_b, 1);

RMSPE_MLE_white_white = mean(RMSPE_MLE_white_white_b, 1);
CyclicErr_MLE_white_white = mean(CyclicErr_MLE_white_white_b, 1);

RMSPE_MLE_white_colored = mean(RMSPE_MLE_white_colored_b, 1);
CyclicErr_MLE_white_colored = mean(CyclicErr_MLE_white_colored_b, 1);

RMSPE_fisher_colored_colored = mean(RMSPE_fisher_colored_colored_b ,1);
CyclicErr_fisher_colored_colored = mean(CyclicErr_fisher_colored_colored_b ,1);

RMSPE_fisher_colored_white = mean(RMSPE_fisher_colored_white_b, 1);
CyclicErr_fisher_colored_white = mean(CyclicErr_fisher_colored_white_b, 1);

RMSPE_fisher_white_white = mean(RMSPE_fisher_white_white_b ,1);
CyclicErr_fisher_white_white = mean(CyclicErr_fisher_white_white_b ,1);

RMSPE_fisher_white_colored = mean(RMSPE_fisher_white_colored_b, 1);
CyclicErr_fisher_white_colored = mean(CyclicErr_fisher_white_colored_b, 1);

CRB_white_reg = mean(CRB_white_reg_b ,1);
CRB_colored_reg = mean(CRB_colored_reg_b ,1);
CRB_white_cyc1 = mean(CRB_white_cyc1_b ,1);
CRB_colored_cyc1 = mean(CRB_colored_cyc1_b ,1);
CRB_white_cyc2 = mean(CRB_white_cyc2_b ,1);
CRB_colored_cyc2 = mean(CRB_colored_cyc2_b ,1);
%--------------------------------------------------------------------------

res = dir('./res/ML_simulation_results_*.mat');
save(append('./res/ML_simulation_results_', string(length(res)+1)));
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



