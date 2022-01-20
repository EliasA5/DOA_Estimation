function [theta_estimate] = Fisher_scoring(theta,lambda,sigma_s,Rv,zx,zy,x,iter)
% This function estimates the DOA (theta) using Fisher's scoring
K = length(zx);
Z = [zx,zy];
N = length(x);
J = FIM(N,SNR,K,lambda,Z);
J_inv = inv(J);
theta_estimate = theta;
for i = 1:iter
    df_theta = PDF_derivative(lambda,zx,zy,sigma_s,Rv,x,theta_estimate);
    theta_estimate = theta_estimate + J_inv*df_theta;
end
end