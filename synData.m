%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: This function generates synthetic sensor data,
% Inputs:
%   rm: The distances matrix of size (3,K), where K is the number of
%   Sensors.
%   theta: DOA angle.
%   alpha: Incidence angle.
%   v0: Wave velocity.
%   sigma_squared: Power of the signal.
%   lambda_noise: Power of the noise.
%   num_samples: Number of samples in each segment P.
%   noise_type: The type of generated noise, can have the values:
%   ['colored', 'white'].
%   f: The frequency of the wave.
%   K_1: The number of 1d sensors.
%   K_3: The number of 3d sensors.
%   P: The numbner of segments in the data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X_w,signal,Q,a] = synData(rm, theta, alpha, v0, sigma_squared, lambda_noise, num_samples, noise_type, f, K_1, K_3, P)
    % buildes syn data in the frequency domain,
zx = rm(:,1);zy = rm(:,2);
distance_mat = sqrt((zx - zx.').^2 + (zy - zy.').^2);
K = length(rm);
switch noise_type %Q is covariance matrix
    case 'colored'
        Q = besselj(0, 2*pi/lambda_noise * distance_mat); %positive semidefinite matrix
    case 'white'
        Q = eye(K);
end

noise = zeros(P, K, M);
X_w = zeros(P, K, M);
for i=1:P
    noise_all = mvnrnd(zeros(2*K, 1), kron(eye(2), 0.5 * Q), M);
    noise(i,:,:) = (noise_all(:,1:K)+1i*noise_all(:,(K+1):end)).';
end

signal_all = mvnrnd(zeros(2,1), 0.5*eye(2)*sigma_squared, M);
signal = (signal_all(:,1) + 1i*signal_all(:,2)).';   % signal = source
[a_func, ~] = model(rm, K_1, K_3, [f]);
% a = a_func(1, theta, alpha, v0).';
% Mue = repmat(a, M, 1).' .* repmat(signal, K, 1);
A = cell2mat(arrayfun(a_func,1 : M,theta*ones(1,M),alpha*ones(1,M), v_0*ones(1,M),'uniformoutput',false));
Mue = A .* signal;

for i=1:P
    X_w(i,:,:) = Mue + squeeze(noise(i,:,:));
end

end