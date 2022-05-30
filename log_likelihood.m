function [log_pdf] = log_likelihood(x,R_inv,M,a_f,v_0,alpha,theta,P)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

a = @(m) a_f(m,theta,alpha,v_0);
A = cell2mat(arrayfun(a,1 : M,'uniformoutput',false));
Up = diag(A' * R_inv * x);
Down = diag(A' * R_inv * A);
s = (Up ./ Down).' / P;
Mue = A .* s;
log_pdf = - sum(diag((x - Mue)' * R_inv * (x - Mue)));

end