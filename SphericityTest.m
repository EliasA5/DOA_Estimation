close all
clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: This script runs the sphericity test on mount meron data.
% Note: Inconclusive results.
% classifing every measurement as colored (= positive) if psi < c and white = false) if psi > c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K_1 = 15;
K_3 = 1;
c = 10^(-100);
sigma_noise = 0.01;  % for white and colored noise
data = load("data.mat");
rm = data.r_m;   
M = 10;   
P = 1;
K = K_1 + 3 * K_3;
TPR = zeros(4,1);
FPR = zeros(4,1);
for j = 1:10

    psi_colored = zeros(100,1);
    psi_white = zeros(100,1);
    for i = 1:100
        X_colored = synData_noise(rm,sigma_noise, M, 'colored', P);
        X_colored = squeeze(sum(X_colored,1));
        psi_colored(i) = NoiseTest(X_colored,K,M); 
    
        X_white = synData_noise(rm,sigma_noise, M, 'white', P);
        X_white = squeeze(sum(X_white,1));
        psi_white(i) = NoiseTest(X_white,K,M); 
    
    end
    TPR(j) = sum(abs(psi_colored) < c);
    FPR(j) = sum(abs(psi_white) < c);
end
TPR = TPR/100;
FPR = FPR/100;

figure
stem(FPR,TPR)
xlabel('False Positive Rate')
ylabel('True Positive Rate')
title('ROC Curve')
grid on
%%
% %all the 19 sensors 
% for i=1:1
%     %system('conda activate obspy & python getData.py'); %uncomment to run the python data getter
%     data = load("data.mat","data");
%     x = data.data;
%     %x(6,:) = [];
%     [K,N] = size(x); %k is number of sensors
%     c = 10^-10;
% 
%     [Rv,psi_] = NoiseTest(x,K,N,10^(-10));
%     figure;
%     subplot(1,2,1)
%     heatmap(db(Rv), 'Colormap', bone);
%     title('ML Estimator')
%     
%     %oldpath = addpath('./LIBRA', '-end');
%     %https://wis.kuleuven.be/stat/robust/LIBRAfiles/LIBRA-home-orig
%     mcdRv = mcdcov(x.','cor', 1, 'plots', 0);
%     mcdCov = mcdRv.cov;
%     subplot(1,2,2)
%     heatmap(db(mcdCov), 'Colormap', bone);
%     title('MCD Estimator')
% end

%useful code
%heatmap(db(Rv), 'Colormap', bone); %, 'Colorlimits', [-350, -270]
%hcn = dsp.ColoredNoise('InverseFrequencyPower',1,'SamplesPerFrame',9600,'NumChannels',19);
%%
%--------------------------------------------------------------------------
function [X_w] = synData_noise(rm,lambda_noise, num_samples, noise_type, P)
    % buildes syn data in the frequency domain with s = 0
zx = rm(:,1);zy = rm(:,2);
distance_mat = sqrt((zx - zx.').^2 + (zy - zy.').^2);
K = length(rm);
switch noise_type %Q is covariance matrix
    case 'colored'
        Q = besselj(0, 2*pi/lambda_noise * distance_mat); %positive semidefinite matrix
    case 'white'
        Q = eye(K);
end

noise = zeros(P, K, num_samples);
X_w = zeros(P, K, num_samples);
for i=1:P
    noise_all = mvnrnd(zeros(2*K, 1), kron(eye(2), 0.5 * Q), num_samples);
    noise(i,:,:) = (noise_all(:,1:K)+1i*noise_all(:,(K+1):end)).';
    X_w(i,:,:) = (noise(i,:,:));
end
end

function [psi_] = NoiseTest(x,K,N) 

  f = @(k) k.' * k;
  g = @(k) k * k.';
  sigma_squared_array = zeros(1,N);
  Rv_array = zeros(K,K,N);

  for i = 1:N
    k = x(:,i);
    sigma_squared_array(i) = f(k);
    Rv_array(:,:,i) = g(k);

  end

  %sigma_squared = 1/(K*N) * sum(sigma_squared_array);
  Rv = 1/N * sum(Rv_array,3);

  psi_ = ( det(Rv) / ((trace(Rv)/K)^K));
  %figure;
  %imagesc(db(co));
end