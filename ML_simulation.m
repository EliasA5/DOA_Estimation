close all
clear all
clc
%--------------------------------------------------------------------------
% This scripts performs the comparison between the MLE and Fisher's scoring
% for white noise and colored noise
% All signal's are in the frequency domain
%%
K_1 = 15;
K_3 = 1;
alpha = pi/4;
v_0 = 1;  % P waves velocities can be from 6 Km/sec to 11 Km/sec
sigma_source = 1;
sigma_noise = 1;  % for white and colored noise
data = load("data.mat");
rm = data.r_m;   % r_m is loaded form the data matrix we extrcted from the getData python script
M = 10;   % for these simlutaions we decided to use one segement
P = 1;
f = 3;      % the omega_m is constant at this point for all the frquencies
acc = 0.001;   % this is the accuracy for the MLE sweep in method 1
iters = 100;    % the number of iterations for the Fisher's scoring (second method)
[a_model, ~] = model(rm, K_1, K_3, f*ones(M,1));
%%
Tests = 11;

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

CRB_white = [];
CRB_colored = [];
theta_og = -pi +2*pi/Tests : 2*pi/Tests : pi-2*pi/Tests;
%theta_og = zeros(Tests,1);

parfor i = 1:Tests-1

    theta = theta_og(i); % theta is in [-pi,pi]
    theta_0 = 0;
    [X_colored,~,Rv_colored,~] = synData(rm, theta, alpha, v_0, sigma_source, sigma_noise, M, 'colored', f, K_1, K_3, P);
    [X_white,s,Rv_white,a] = synData(rm, theta, alpha, v_0, sigma_source, sigma_noise, M, 'white', f, K_1, K_3, P);


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
    theta_colored_fisher = real(Fisher_scoring(theta_0,s,Rv_colored,f,v_0,alpha,K_3,K_1,X_colored,iters,rm));
    RMSPE_fisher_colored = [RMSPE_fisher_colored, MSPE(theta, theta_colored_fisher, 'MSPE')];
    CyclicErr_fisher_colored = [CyclicErr_fisher_colored, MSPE(theta, theta_colored_fisher, 'cyclic')];
    ThetaEst_fisher_colored = [ThetaEst_fisher_colored, theta_colored_fisher];
    
    % Fisher's scoring for white noise
    theta_white_fisher = real(Fisher_scoring(theta_0,s,Rv_white,f,v_0,alpha,K_3,K_1,X_white,iters,rm));
    RMSPE_fisher_white = [RMSPE_fisher_white, MSPE(theta, theta_white_fisher, 'MSPE')];
    CyclicErr_fisher_white = [CyclicErr_fisher_white, MSPE(theta, theta_white_fisher, 'cyclic')];
    ThetaEst_fisher_white = [ThetaEst_fisher_white, theta_white_fisher];

    % Calculating the CRB
    CRB_white(i) = CRB(v_0,alpha,rm,K_3,K_1,theta,f,Rv_white,M,s);
    CRB_colored(i) = CRB(v_0,alpha,rm,K_3,K_1,theta,f,Rv_colored,M,s);
end

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

figure;
hold on;
plot(theta_og,RMSPE_fisher_white(1:Tests-1))
plot(theta_og,RMSPE_fisher_colored(1:Tests-1))
plot(theta_og,RMSPE_MLE_white(1:Tests-1))
plot(theta_og,RMSPE_MLE_colored(1:Tests-1))
plot(theta_og,CRB_white(1:Tests-1),'*')
plot(theta_og,CRB_colored(1:Tests-1),'-o')
hold off;
grid on;
legend('fisher null' , 'fisher colored' , 'MLE null' , 'MLE colored' , 'CRB white' , 'CRB colored');




%--------------------------------------------------------------------------
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


