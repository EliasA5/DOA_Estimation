close all
clear all
clc
%--------------------------------------------------------------------------
% This scripts performs the comparison between the MLE and Fisher's scoring
% for white noise and colored noise
% All signal's are in the frequency domain
%--------------------------------------------------------------------------
%% Initialization

M = 10;                          
P = 1;
K_1 = 15;
K_3 = 1;
K = 3 * K_3 + K_1;

alpha = pi/4;
v_0 = 9600;                     % P waves velocities can be from 6 Km/sec to 11 Km/sec
f = 3;                        
w = f * ones(M,1);              % the omega_m is constant at this point for all the frquencies

SNR = 10;
sigma_source = 10;
sigma_noise = 1;                % for white and colored noise

data = load("data.mat");
rm = data.r_m;                  % r_m is loaded form the data matrix we extrcted from the getData python script

acc = 0.001;                    % this is the accuracy for the MLE sweep in method 1
iters = 400;                    % the number of iterations for the Fisher's scoring (second method)
step_size_white = 1;
step_size_colored = 1;
gamma = 0.95;
theta_0 = pi/4;                 % starting estimate at 45 deg

[a_model, da_model] = model(rm, K_1, K_3, w);

%% Estimation 

Tests = 10;

RMSPE_MLE_colored = [];
CyclicErr_MLE_colored = [];
ThetaEst_MLE_colored = [];

RMSPE_MLE_white = [];
CyclicErr_MLE_white = [];
ThetaEst_MLE_white = [];

RMSPE_fisher_colored = [];
CyclicErr_fisher_colored = [];
ThetaEst_fisher_colored = [];

RMSPE_fisher_white = [];
CyclicErr_fisher_white = [];
ThetaEst_fisher_white = [];

CRB_white_reg = [];
CRB_colored_reg = [];
CRB_white_cyc1 = [];
CRB_colored_cyc1 = [];
CRB_white_cyc2 = [];
CRB_colored_cyc2 = [];

theta_og = -pi +2*pi/Tests : 2*pi/Tests : pi;

parfor i = 1 : Tests

    theta = theta_og(i);         % theta is in [-pi,pi]

    [X_colored,~,Rv_colored,~] = synData(rm, theta, alpha, v_0, sigma_source, sigma_noise, M, 'colored', w, K_1, K_3, P);
    [X_white,~,Rv_white,~] = synData(rm, theta, alpha, v_0, sigma_source, sigma_noise, M, 'white', w, K_1, K_3, P);


    % MLE for colored noise
    fun_colored = toMaximize(a_model, Rv_colored, X_colored, M, P);
    theta_colored_MLE = real(MaximizeTheta(fun_colored, alpha, v_0, acc));
    RMSPE_MLE_colored = [RMSPE_MLE_colored, MSPE(theta, theta_colored_MLE, 'MSPE')];
    CyclicErr_MLE_colored = [CyclicErr_MLE_colored, MSPE(theta, theta_colored_MLE, 'cyclic')];
    ThetaEst_MLE_colored = [ThetaEst_MLE_colored, theta_colored_MLE];

    % MLE for white noise
    fun_white = toMaximize(a_model, Rv_white, X_white, M, P);
    theta_white_MLE = real(MaximizeTheta(fun_white, alpha, v_0, acc));
    RMSPE_MLE_white = [RMSPE_MLE_white, MSPE(theta, theta_white_MLE, 'MSPE')];
    CyclicErr_MLE_white = [CyclicErr_MLE_white, MSPE(theta, theta_white_MLE, 'cyclic')];
    ThetaEst_MLE_white = [ThetaEst_MLE_white, theta_white_MLE];

    % Fisher's scoring for colored noise
    [~,theta_colored_fisher] = Fisher_scoring('syn',theta_0,Rv_colored,v_0,alpha,K,X_colored,iters,step_size_colored,gamma,M,P,a_model,da_model);
    RMSPE_fisher_colored = [RMSPE_fisher_colored, MSPE(theta, theta_colored_fisher, 'MSPE')];
    CyclicErr_fisher_colored = [CyclicErr_fisher_colored, MSPE(theta, theta_colored_fisher, 'cyclic')];
    ThetaEst_fisher_colored = [ThetaEst_fisher_colored, theta_colored_fisher];
    
    % Fisher's scoring for white noise
    [~,theta_white_fisher] = Fisher_scoring('syn',theta_0,Rv_white,v_0,alpha,K,X_white,iters,step_size_white,gamma,M,P,a_model,da_model);
    RMSPE_fisher_white = [RMSPE_fisher_white, MSPE(theta, theta_white_fisher, 'MSPE')];
    CyclicErr_fisher_white = [CyclicErr_fisher_white, MSPE(theta, theta_white_fisher, 'cyclic')];
    ThetaEst_fisher_white = [ThetaEst_fisher_white, theta_white_fisher];

    % Calculating the CRB (calculating all kinds)

    CRB_white_reg = [CRB_white_reg, CRB('regular',v_0,alpha,theta,Rv_white,M,X_white,a_model,da_model,P)];
    CRB_colored_reg = [CRB_colored_reg, CRB('regular',v_0,alpha,theta,Rv_colored,M,X_colored,a_model,da_model,P)];

    CRB_white_cyc1 = [CRB_white_cyc1, CRB('cyclic 1',v_0,alpha,theta,Rv_white,M,X_white,a_model,da_model,P)];
    CRB_colored_cyc1 = [CRB_colored_cyc1, CRB('cyclic 1',v_0,alpha,theta,Rv_colored,M,X_colored,a_model,da_model,P)];

    CRB_white_cyc2 = [CRB_white_cyc2, CRB('cyclic 2',v_0,alpha,theta,Rv_white,M,X_white,a_model,da_model,P)];
    CRB_colored_cyc2 = [CRB_colored_cyc2, CRB('cyclic 2',v_0,alpha,theta,Rv_colored,M,X_colored,a_model,da_model,P)];
end

%--------------------------------------------------------------------------
%% Graphs 

% figure;
% hold on;
% plot(theta_og,RMSPE_fisher_white)
% plot(theta_og,RMSPE_fisher_colored)
% plot(theta_og,RMSPE_MLE_white)
% plot(theta_og,RMSPE_MLE_colored)
% plot(theta_og,CRB_white)
% plot(theta_og,CRB_colored)
% hold off;
% grid on;
% legend('fisher null' , 'fisher colored' , 'MLE null' , 'MLE colored' , 'CRB white' , 'CRB colored');

% figure;
% hold on;
% plot(theta_og,RMSPE_fisher_white(1:Tests-1))
% plot(theta_og,RMSPE_fisher_colored(1:Tests-1))
% plot(theta_og,RMSPE_MLE_white(1:Tests-1))
% plot(theta_og,RMSPE_MLE_colored(1:Tests-1))
% plot(theta_og,CRB_white(1:Tests-1),'*')
% plot(theta_og,CRB_colored(1:Tests-1),'-o')
% hold off;
% grid on;
% legend('fisher null' , 'fisher colored' , 'MLE null' , 'MLE colored' , 'CRB white' , 'CRB colored');

figure;
hold on; grid on;
plot(theta_og,'*','LineWidth', 2)
plot(ThetaEst_MLE_white,'o'); plot(ThetaEst_MLE_colored,'o');
plot(ThetaEst_fisher_white,'x'); plot(ThetaEst_fisher_colored,'x');
legend('Original', 'MLE W', 'MLE C', 'Fisher W', 'Fisher C')
xlabel('Test'); ylabel('\theta'); title('Comparison of the Estimations')
hold off



%% Circle
angles = linspace(0, 2*pi, 500);
radius = 20;
CenterX = 50;
CenterY = 40;
x = radius * cos(angles) + CenterX;
y = radius * sin(angles) + CenterY;

x_og = radius * cos(theta_og) + CenterX;
x_MLE_white = radius * cos(ThetaEst_MLE_white) + CenterX;
x_MLE_colored = radius * cos(ThetaEst_MLE_colored) + CenterX;
x_Fisher_white = radius * cos(ThetaEst_fisher_white) + CenterX;
x_Fisher_colored = radius * cos(ThetaEst_fisher_colored) + CenterX;

y_og = radius * sin(theta_og) + CenterY;
y_MLE_white = radius * sin(ThetaEst_MLE_white) + CenterY;
y_MLE_colored = radius * sin(ThetaEst_MLE_colored) + CenterY;
y_Fisher_white = radius * sin(ThetaEst_fisher_white) + CenterY;
y_Fisher_colored = radius * sin(ThetaEst_fisher_colored) + CenterY;


figure
hold on; grid on ; axis equal;
plot(x, y, 'k-', 'LineWidth', 1);
plot(x_og, y_og, 'r*','LineWidth', 2.5)
plot(x_MLE_white, y_MLE_white, 'o','LineWidth', 1.5)
plot(x_MLE_colored, y_MLE_colored, 'o','LineWidth', 1.5)
plot(x_Fisher_white, y_Fisher_white, 'x','LineWidth', 1.5)
plot(x_Fisher_colored, y_Fisher_colored, 'x','LineWidth', 1.5)
xlabel('X', 'FontSize', 10); ylabel('Y', 'FontSize', 10);
legend('','Original', 'MLE W', 'MLE C', 'Fisher W', 'Fisher C')
title('Comparison of the Estimations on Circle')


%--------------------------------------------------------------------------
%% Functions
function [fun] = toMaximize(a, R, X_w, M, P)
    R_inv = pinv(R);
    X_w_p = squeeze(sum(X_w,1));
    norm = @(m, theta, alpha, v_0) a(m, theta, alpha, v_0)' * R_inv * a(m, theta, alpha, v_0);
    upper = @(m, theta, alpha, v_0) a(m, theta, alpha, v_0)' * R_inv * X_w_p(:,m);
    f = @(m, theta, alpha, v_0) 1/P * upper(m, theta, alpha, v_0)' * upper(m, theta, alpha, v_0);
    fun = @(theta, alpha, v_0) 0;
    for m = 1:M
        fun = @(theta, alpha, v_0) fun(theta, alpha, v_0) + f(m, theta, alpha, v_0)/norm(m, theta, alpha, v_0);
    end
    fun = @(theta, alpha, v_0) real(fun(theta, alpha, v_0));
end

function [val] = MaximizeTheta(fun, alpha, v_0, acc)
    t_vec = -pi:acc:pi;
    Alphas = ones(size(t_vec)) * alpha;
    V_0 = ones(size(t_vec)) * v_0;
    max_vals = arrayfun(fun, t_vec, Alphas, V_0);
    [~, I] = max(max_vals);
    val = t_vec(I);
end


