%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: This function returns the function that we need to maximize for
% the MLE estimator, feed it's output to MaximizeTheta.
% Inputs:
%   a: The function used to calculate the model, can be obtained from
%   model.m file.
%   R: The covariance matrix, this function assumes white noise, so the
%   covariance matrix gets canceled and is unused, but kept for compliance
%   for other functions of the same purpose.
%   X_w: The signal in the frequency domain.
%   M: The number of frequencies, in other files this is given by L.
%   P: The number of segments.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fun] = toMaximizeMLEWhite(a, R, X_w, M, P)
    X_w_p = squeeze(sum(X_w,1));
    fun = @(theta, alpha, v_0) d(theta, alpha, v_0, a, [], X_w_p, M, P);
end

function [val] = d(theta, alpha, v_0, a, R_inv, X_w_p, M, P)
a_val = cell2mat(arrayfun(a,(1:M), theta*ones(1,M), alpha*ones(1,M), v_0*ones(1,M),'UniformOutput',false));
d = @(m, theta, alpha, v_0) (1/P * (a_val(:, m)' * X_w_p(:,m))' * (a_val(:, m)' * X_w_p(:,m)))/...
                                (a_val(:, m)' * a_val(:, m));
val = real(sum(arrayfun(d, (1:M), theta*ones(1,M), alpha*ones(1,M), v_0*ones(1,M))));
end