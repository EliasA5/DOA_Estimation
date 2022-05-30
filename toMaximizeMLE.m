%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: This function returns the function that we need to maximize for
% the MLE estimator, feed it's output to MaximizeTheta.
% Inputs:
%   a: The function used to calculate the model, can be obtained from
%   model.m file.
%   R: The covariance matrix, if assuming white noise use the function
%   named toMaximizeMLEWhite.
%   X_w: The signal in the frequency domain.
%   M: The number of frequencies, in other files this is given by L.
%   P: The number of segments.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fun] = toMaximizeMLE(a, R, X_w, M, P)
    R_inv = pinv(R);
    X_w_p = squeeze(sum(X_w,1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % old method for calculating the function, kept for documentation and
    % debugging.
%     norm = @(m, theta, alpha, v_0) a(m, theta, alpha, v_0)' * R_inv * a(m, theta, alpha, v_0);
%     upper = @(m, theta, alpha, v_0) a(m, theta, alpha, v_0)' * R_inv * X_w_p(:,m);
%     f = @(m, theta, alpha, v_0) 1/P * upper(m, theta, alpha, v_0)' * upper(m, theta, alpha, v_0);
%     fun = @(theta, alpha, v_0) 0;
%     for m = 1:M
%         fun = @(theta, alpha, v_0) fun(theta, alpha, v_0) + f(m, theta, alpha, v_0)/norm(m, theta, alpha, v_0);
%     end
%     d = @(m, theta, alpha, v_0) (1/P * (a(m, theta, alpha, v_0)' * R_inv * X_w_p(:,m))' * (a(m, theta, alpha, v_0)' * R_inv * X_w_p(:,m)))/...
%                                 (a(m, theta, alpha, v_0)' * R_inv * a(m, theta, alpha, v_0));
%     fun =  @(theta, alpha, v_0) real(sum(arrayfun(d, (1:M), theta*ones(1,M), alpha*ones(1,M), v_0*ones(1,M))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fun = @(theta, alpha, v_0) d(theta, alpha, v_0, a, R_inv, X_w_p, M, P);
end

function [val] = d(theta, alpha, v_0, a, R_inv,X_w_p,M,P)
a_val = cell2mat(arrayfun(a,(1:M), theta*ones(1,M), alpha*ones(1,M), v_0*ones(1,M),'UniformOutput',false));
d = @(m, theta, alpha, v_0) (1/P * (a_val(:, m)' * R_inv * X_w_p(:,m))' * (a_val(:, m)' * R_inv * X_w_p(:,m)))/...
                                (a_val(:, m)' * R_inv * a_val(:, m));
val = real(sum(arrayfun(d, (1:M), theta*ones(1,M), alpha*ones(1,M), v_0*ones(1,M))));
end