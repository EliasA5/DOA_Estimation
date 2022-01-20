function [J_p] = FIM(N,SNR,K,lambda,Z)
%This function calculates the periodic FIM based on equation 3.15
% lambda - wave length. K sensors, N time samples.
% Z - sensors position matrix. size[Z] = [K,2]
% SNR = var(s)/var(v)
A = ones(k);
B = Z*Z.';
J_p = (N*SNR^2*K)/(1+K*SNR)*(2*pi/lambda)^2*(-A.'*B*A/K+trace(B));
end