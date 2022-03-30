close all
clear all
clc
%--------------------------------------------------------------------------
% classifing every measurement as colored (= positive) if psi < c and white = false) if psi > c
%%
c = 10^(-100);
sigma_noise = 0.01;
P = 1;
files = dir('./matFiles/*.mat');

TPR = zeros(4,1);
FPR = zeros(4,1);
real_psi = [];
syn_psi_white = [];
syn_psi_colored = [];
j = 1;
limit = false;

for file = files'
    load(fullfile(file.folder, file.name));
    data_size = size(data);
    err_size = size(err);

    for i = 1:data_size(1)
        x = squeeze(data(i,:,:));
        rm = squeeze(distances(i,:,:));
        [M,N] = size(x); %M is number of sensors
        K_3 = 3*floor((M - min([length(unique(rm(:,1))), length(unique(rm(:,2)))]))/3);
        K_1 = M - K_3;
        K = K_1 + K_3;
        real_psi = [real_psi, NoiseTest(x,K,N)]; 


        X_colored = synData_noise(rm, sigma_noise, N, 'colored', P);
        X_colored = squeeze(sum(X_colored,1));
        syn_psi_colored = [syn_psi_colored, NoiseTest(X_colored,K,N)]; 
    
        X_white = synData_noise(rm,sigma_noise, N, 'white', P);
        X_white = squeeze(sum(X_white,1));
        syn_psi_white = [syn_psi_white, NoiseTest(X_white,K,N)]; 
    
    end
    j = j+1;
    if(limit && j == 3), break; end
end


figure;
stem(MSPE(real_psi,abs(syn_psi_colored),'MSE'));hold on; stem(MSPE(real_psi,abs(syn_psi_white),'MSE'));
xlabel('False Positive Rate');ylabel('True Positive Rate');
title('ROC Curve');grid on;
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

  psi_ = ( det(Rv) / ((trace(Rv)/K)^K))^(-K/2);
  %figure;
  %imagesc(db(co));
end