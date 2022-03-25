function [log_pdf] = log_likelihood(x,mue,R_inv,M)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

log_pdf = 0;

for i = 1:M
    x_i = x(:,i);
    log_pdf = log_pdf - (x_i - mue)' * R_inv * (x_i - mue);

end


end