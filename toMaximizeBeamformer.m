%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: This function returns the function that we need to maximize for
% the beamformer estimator, feed it's output to MaximizeTheta.
% Inputs:
%   a: The function used to calculate the model, can be obtained from
%   model.m file.
%   R: The covariance matrix that is not used, but kept in the input for
%   compliance with functions of the same purpose.
%   X_w: The signal in the frequency domain.
%   M: The number of frequencies, in other files this is given by L.
%   P: The number of segments.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fun] = toMaximizeBeamformer(a, R, X_w, M, P)
    [~,K,~] = size(X_w);
    X_temp = squeeze(X_w(1,:,:));
    Rk = X_temp * X_temp';
    for p = 2:P
        X_temp = squeeze(X_w(p,:,:));
        Rk = 1/(p+1) * (p*Rk + X_temp * X_temp');
    end
    
    fun = @(theta, alpha, v_0) d(theta, alpha, v_0, a, Rk, M);
end

function [val] = d(theta, alpha, v_0, a, Rk, M)
a_val = cell2mat(arrayfun(a,(1:M), theta*ones(1,M), alpha*ones(1,M), v_0*ones(1,M),'UniformOutput',false));
d = @(m, theta, alpha, v_0) a_val(:, m)' * Rk * a_val(:, m);
val = real(sum(arrayfun(d, (1:M), theta*ones(1,M), alpha*ones(1,M), v_0*ones(1,M))));
end