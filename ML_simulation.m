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
acc = 0.1;   % this is the accuracy for the MLE sweep in method 1
iters = 1000;    % the number of iterations for the Fisher's scoring (second method)
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

    [X_colored,~,Rv_colored,~] = synData(rm, theta, alpha, v_0, sigma_source, sigma_noise, M, 'colored', f, K_1, K_3);
    [X_white,s,Rv_white,a] = synData(rm, theta, alpha, v_0, sigma_source, sigma_noise, M, 'white', f, K_1, K_3);


    % MLE for colored noise
    fun_colored = toMaximize(a, Rv_colored, X_colored, M, P);
    theta_colored = MaximizeTheta(fun_colored, alpha, v_0, acc);
    RMSPE_MLE_colored(i) = sqrt(mean((theta-theta_colored).^2));

    % MLE for white noise
    fun_white = toMaximize(a, Rv_white, X_white, M, P);
    theta_white = MaximizeTheta(fun_white, alpha, v_0, acc);
    RMSPE_MLE_white(i) = sqrt(mean((theta-theta_white).^2));

    % Fisher's scoring for colored noise
    theta_colored = Fisher_scoring(theta,s,Rv_colored,f,v_0,alpha,K_3,K_1,X_colored,iters);
    RMSPE_fisher_colored = sqrt(mean((theta-theta_colored).^2));

    % Fisher's scoring for white noise
    theta_white = Fisher_scoring(theta,s,Rv_white,f,v_0,alpha,K_3,K_1,X_white,iters);
    RMSPE_fisher_white = sqrt(mean((theta-theta_white).^2));

end

figure;
hold on;
plot(RMSPE_fisher_white)
plot(RMSPE_fisher_colored)
plot(RMSPE_MLE_white)
plot(RMSPE_MLE_colored)
hold off;
legend('fisher null' , 'fisher colored' , 'MLE null' , 'MLE colored');





function [fun] = toMaximize(a, R, X_w, M, P)
    R_inv = inv(R);
    %X_w_p = squeeze(sum(X_w,1));
    X_w_p = X_w;
    norm = @(m, theta, alpha, v_0) a(m, theta, alpha, v_0)' * R_inv * a(m, theta, alpha, v_0);
    upper = @(m, theta, alpha, v_0) a(m, theta, alpha, v_0)' * R_inv * X_w_p;
    f = @(m, theta, alpha, v_0) 1/P * upper(m, theta, alpha, v_0)' * upper(m, theta, alpha, v_0);
    fun = @(theta, alpha, v_0) 0;
    for m = 1:M
        fun = @(theta, alpha, v_0) fun(theta, alpha, v_0) + f(m, theta, alpha, v_0)/norm(m, theta, alpha, v_0);
    end
    fun = @(theta, alpha, v_0) real(fun(theta, alpha, v_0));
end

function [val] = MaximizeTheta(fun, alpha, v_0, acc)
    t_vec = -pi:acc:pi;
    val = -pi;
    curr_max = 0;
    for t = t_vec
        tmp = fun(t, alpha, v_0);
        if(tmp > curr_max)
            curr_max = tmp;
            val = t;
        end
    end
end


