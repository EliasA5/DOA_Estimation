%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: This function finds the maximum of the function fun with respect
% to theta in [0,2*pi] by searching through all values with a specified
% resolution.
% Inputs:
%   fun: The function to maximize, takes 3 parameters.
%   alpha: The value of the second parameter.
%   v_0: The value of the third parameter.
%   acc: The accuracy of the search, prefed value is in 1e-3 for 2.5 second
%   runtime with our used functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [val] = MaximizeTheta(fun, alpha, v_0, acc)
    t_vec = 0:acc:2*pi;
    Alphas = ones(size(t_vec)) * alpha;
    V_0 = ones(size(t_vec)) * v_0;
    max_vals = arrayfun(fun, t_vec, Alphas, V_0);
    [~, I] = max(max_vals);
    val = t_vec(I);
end