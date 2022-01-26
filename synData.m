function [X_w] = synData(rm, theta, alpha, v0, sigma_squared, lambda_noise, num_samples, noise_type, f, K_1, K_3)
zx = rm(:,1);zy = rm(:,2);
distance_mat = sqrt((zx - zx.').^2 + (zy - zy.').^2);
K = length(rm);
switch noise_type %Q is covariance matrix
    case 'colored'
        Q = besselj(0, 2*pi/lambda_noise * distance_mat); %positive semidefinite matrix
    case 'white'
        Q = eye(K);
end

noise_all = mvnrnd(zeros(2*K, 1), kron(eye(2), 0.5 * Q), num_samples);
noise = noise_all(:,1:K)+1i*noise_all(:,(K+1):end);
signal_all = mvnrnd(zeros(2,1), 0.5*eye(2)*sigma_squared, num_samples);
signal = signal_all(:,1) + 1i*signal_all(:,2);
a_func = model(rm, K_1, K_3, [f]);
a = a_func(1, theta, alpha, v0).';
X_w = repmat(a, num_samples, 1).* repmat(signal, 1, K) + noise;

end