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

noise = zeros(P, K, num_samples);
X_w = zeros(P, K, num_samples);
for i=1:P
    noise_all = mvnrnd(zeros(2*K, 1), kron(eye(2), 0.5 * Q), num_samples);
    noise(i,:,:) = (noise_all(:,1:K)+1i*noise_all(:,(K+1):end)).';
end

signal_all = mvnrnd(zeros(2,1), 0.5*eye(2)*sigma_squared, num_samples);
signal = (signal_all(:,1) + 1i*signal_all(:,2)).';   % signal = source
a_func = model(rm, K_1, K_3, [f]);
a = a_func(1, theta, alpha, v0).';
a_s = repmat(a, num_samples, 1).' .* repmat(signal, K, 1);
for i=1:P
    X_w(i,:,:) = a_s + squeeze(noise(i,:,:));
end

end