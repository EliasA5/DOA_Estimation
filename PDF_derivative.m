function [df_theta] = PDF_derivative(lambda,zx,zy,sigma_s,Rv,x,theta)
% This function calculates the derivative of the log of the PDF

A = exp(-1i*2*(pi/lambda)*[zx.',zy.']*[sind(theta);cosd(theta)]);
Rx = sigma_s*A + Rv;
R_inv = inv(Rx);
dR_theta = R_derivative(lambda,zx,zy,theta,sigma_s,A);
df_theta = x'*R_inv*dR_theta*R_inv*x - trace(R_inv*dR_theta);
end