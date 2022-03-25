close all
clear all
clc
%--------------------------------------------------------------------------
% This scripts checks the effect of the number of iterations of the Fisher's scoring ethod
% Check the maximisation of the PDF
%%
K_1 = 15;
K_3 = 1;
alpha = pi/4;
v_0 = 1;  % P waves velocities can be from 6 Km/sec to 11 Km/sec
SNR = 10;
sigma_source = 10;
sigma_noise = 1;  % for white and colored noise
data = load("data.mat");
rm = data.r_m;   % r_m is loaded form the data matrix we extrcted from the getData python script
M = 100;   % for these simlutaions we decided to use one segement
P = 1;
f = 3;      % the omega_m is constant at this point for all the frquencies

%%
Tests = 1;

RMSPE_colored = zeros(Tests,1);
RMSPE_white= zeros(Tests,1);
PDF_colored = zeros(5000,1);
PDF_white = zeros(5000,1);

RMS_colored = zeros(5000,1);
RMS_white = zeros(5000,1);

%theta_og = -pi +2*pi/Tests : 2*pi/Tests : pi;
%theta_og = zeros(Tests,1);
theta_og = rand*2*pi - pi;

step_size = 1;
gamma = 0.8;  % tuning parameter

for i = 1:Tests

    theta = theta_og(i); % theta is in [-pi,pi]

    [X_colored,~,Rv_colored,~] = synData(rm, theta, alpha, v_0, sigma_source, sigma_noise, M, 'colored', f, K_1, K_3, P);
    [X_white,s,Rv_white,a] = synData(rm, theta, alpha, v_0, sigma_source, sigma_noise, M, 'white', f, K_1, K_3, P);

    % Fisher's scoring for colored noise
    theta_estimate = theta;

    [theta_estimate,PDF_colored(1)] = Fisher_step(theta_estimate,v_0,alpha,rm,K_3,K_1,f,s,Rv_colored,X_colored,M,step_size);
    theta_colored = real(theta_estimate);
    RMS_colored(1) = sqrt(mean((theta-theta_colored).^2));
    for j = 2:5000
        [theta_estimate,PDF_colored(j)] = Fisher_step(theta_estimate,v_0,alpha,rm,K_3,K_1,f,s,Rv_colored,X_colored,M,step_size);
        theta_colored = real(theta_estimate);
        if abs(theta_colored) > pi
              theta_colored = mod(theta_colored,pi);
        end
        RMS_colored(j) = sqrt(mean((theta-theta_colored).^2));
%         if (PDF_colored(j) >= PDF_colored(j-1))
%             step_size = gamma*step_size;
%         end
        if RMS_colored(j) >= RMS_colored(j-1)
             step_size = gamma*step_size;
        end
    end
    theta_colored = real(theta_estimate);
%     if abs(theta_colored) > pi
%           theta_colored = mod(theta_colored,pi);
%     end
    RMSPE_colored(i) = mod(sqrt(mean((theta-theta_colored).^2)),2*pi);
    fprintf('The final step size for colored: %f\n',step_size);

    %%
    step_size = 1;
    
    % Fisher's scoring for white noise
    theta_estimate = 0;

    [theta_estimate,PDF_white(1)] = Fisher_step(theta_estimate,v_0,alpha,rm,K_3,K_1,f,s,Rv_white,X_white,M,step_size);
    theta_white = real(theta_estimate);

    RMS_white(1) = sqrt(mean((theta-theta_white).^2));
    for j = 2:5000
        [theta_estimate,PDF_white(j)] = Fisher_step(theta_estimate,v_0,alpha,rm,K_3,K_1,f,s,Rv_white,X_white,M,step_size);
        theta_white = real(theta_estimate);
        RMS_white(j) = sqrt(mean((theta-theta_white).^2));

        if RMS_white(j) >= RMS_white(j-1)
            step_size = gamma*step_size;
        end
%         if abs(theta_colored) > pi
%           theta_white = mod(theta_white,pi) - pi ;
%         end
%         if (PDF_white(j) >= PDF_white(j-1))
%             step_size = gamma*step_size;
%         end
    end
    theta_white = real(theta_estimate);
%     if abs(theta_colored) > pi
%           theta_colored = mod(theta_colored,pi);
%     end
    RMSPE_white(i) = mod(sqrt(mean((theta-theta_white).^2)),2*pi);
    fprintf('The final step size for white: %f\n',step_size);
end

%--------------------------------------------------------------------------
function [theta_k_1,df_theta] = Fisher_step(theta_k,v_0,alpha,rm,K_3,K_1,f,s,Rv,X,M,step_size)

    da =  A_derivative(v_0,alpha,rm,K_3,K_1,theta_k,f);
    df_theta = PDF_derivative(K_3,K_1,v_0,alpha,f,s,rm,Rv,X,theta_k,da,M);
    F = FIM(s,M,Rv,da); %scalar for only theta
    theta_k_1 = theta_k + step_size*df_theta/F;

end







