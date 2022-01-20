function [dR_theta] = R_derivative(lambda,zx,zy,theta,sigma_s,A)
% This function calculates the covarince matrix R's derivative
% R = sigma_s*a(Z,theta)*(a(z,theta)^H)+ sigma_v*I - white noise
% R = sigma_s*a(Z,theta)*(a(z,theta)^H)+ Rv - colored noise

d2 = (-1i*2*(pi/lambda))*[zx.',zy.']*[cosd(theta); -sind(theta)];
dR_theta = sigma_s*A*d2;

end