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
alpha = rand * pi/2;
v_0 = 1;  % P waves velocities can be from 6 Km/sec to 11 Km/sec
sigma_source = 1;
sigma_noise = 1;  % for white and colored noise
data = load("data.mat");
rm = data.r_m;   % r_m is loaded form the data matrix we extrcted from the getData python script
M = 10;   % for these simlutaions we decided to use one segement
P = 1;
f = 3;      % the omega_m s constant at this point for all the frquencies
acc = 0.001;   % this is the accuracy for the MLE sweep in method 1
iters = 1500;    % the number of iterations for the Fisher's scoring (second method)
a_model = model(rm, K_1, K_3, f*ones(M,1));
%%
Tests = 2;

RMSPE_MLE_colored = zeros(Tests,1);
RMSPE_MLE_white = zeros(Tests,1);
RMSPE_fisher_colored = zeros(Tests,1);
RMSPE_fisher_white = zeros(Tests,1);
theta_og = zeros(Tests,1);

for i = 1:Tests

    theta = rand*2*pi - pi; % theta is in [-pi,pi)
    theta_og(i) = theta;

    [X_colored,~,Rv_colored,~] = synData(rm, theta, alpha, v_0, sigma_source, sigma_noise, M, 'colored', f, K_1, K_3, P);
    [X_white,s,Rv_white,a] = synData(rm, theta, alpha, v_0, sigma_source, sigma_noise, M, 'white', f, K_1, K_3, P);


    % MLE for colored noise
    fun_colored = toMaximize(a_model, Rv_colored, X_colored, M, P);
    theta_colored = MaximizeTheta(fun_colored, alpha, v_0, acc);
    RMSPE_MLE_colored(i) = sqrt(mean((theta-theta_colored).^2));

    % MLE for white noise
    fun_white = toMaximize(a_model, Rv_white, X_white, M, P);
    theta_white = MaximizeTheta(fun_white, alpha, v_0, acc);
    RMSPE_MLE_white(i) = sqrt(mean((theta-theta_white).^2));

    % Fisher's scoring for colored noise
    theta_colored = Fisher_scoring(theta,s,Rv_colored,f,v_0,alpha,K_3,K_1,X_colored,iters,rm);
    RMSPE_fisher_colored(i) = sqrt(mean((theta-theta_colored).^2));

    % Fisher's scoring for white noise
    theta_white = Fisher_scoring(theta,s,Rv_white,f,v_0,alpha,K_3,K_1,X_white,iters,rm);
    RMSPE_fisher_white(i) = sqrt(mean((theta-theta_white).^2));

end

figure;
hold on;
stem(RMSPE_fisher_white)
stem(RMSPE_fisher_colored)
stem(RMSPE_MLE_white)
stem(RMSPE_MLE_colored)
hold off;
grid on;
legend('fisher null' , 'fisher colored' , 'MLE null' , 'MLE colored');





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


