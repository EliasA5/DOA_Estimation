function [theta_estimate,theta_final] = Fisher_scoring(theta_0,Rv,v_0,alpha,K,X,iters,step_size,gamma,M,P,a_f,da_f,acc)
% This function estimates the DOA (theta) using Fisher's scoring

theta_estimate = zeros(1,iters);
theta_estimate(1) = theta_0;

pdf = zeros(1,iters);
% s = zeros(1,M);

Rv_inv = pinv(Rv);
X_w = squeeze(sum(X,1));

% a = @(m, theta) a_f(m,theta,alpha,v_0);
% da = @(m, theta) da_f(m,theta,alpha,v_0);

%--------------------------------------------------------------------------
for i = 1 : iters - 1

    theta = theta_estimate(i);
    A = cell2mat(arrayfun(a_f,1 : M,theta*ones(1,M),alpha*ones(1,M), v_0*ones(1,M),'uniformoutput',false));
    dA = cell2mat(arrayfun(da_f,1 : M,theta*ones(1,M),alpha*ones(1,M), v_0*ones(1,M),'uniformoutput',false));

    
    Up = diag(A' * Rv_inv * X_w);
    Down = diag(A' * Rv_inv * A);
    s = (Up ./ Down).' / P;
    Mue = A .* s;
    dMue = dA .* s;
    df_theta = sum(2 * real(diag((X_w - Mue)' * Rv_inv * dMue)),1);
    FIM = sum(2 * real(diag(dMue' * Rv_inv * dMue)));
    crb = FIM ^ (-1);
    ccrb = 2 - 2/(sqrt(crb + 1));

    pdf(i+1) = log_likelihood(X_w,Rv_inv,M,a_f,v_0,alpha,theta_estimate(i+1),P);

    theta_estimate(i+1) = wrapToPi(real(theta_estimate(i) + step_size * ccrb * df_theta));

    if (pdf(i+1) < pdf(i))
        step_size = step_size * gamma;
    end

    if (mod(abs(theta_estimate(i+1) - theta_estimate(i)),2*pi) < acc)
        theta_estimate(end) = theta_estimate(i+1);
        break;
    end
end
theta_final = theta_estimate(end);
end