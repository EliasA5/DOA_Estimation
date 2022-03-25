function [log_pdf] = log_likelihood(x,s,R_inv,M,a_f,v_0,alpha,theta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

log_pdf = 0;
for i = 1:M
    x_i = x(:,i);
    a = a_f(i, theta, alpha, v_0);
    mue = a * s(i);
    log_pdf = log_pdf - (x_i - mue)' * R_inv * (x_i - mue);

end
end