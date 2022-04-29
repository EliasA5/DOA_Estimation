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