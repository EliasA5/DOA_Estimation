function [theta_estimate] = Fisher_scoring(theta_0,Rv,f,v_0,alpha,K_3,K_1,X,iters,r_k,step_size,gamma,M,s)
% This function estimates the DOA (theta) using Fisher's scoring

theta_estimate = zeros(1,iters);
theta_estimate(1) = theta_0;

a = zeros(3*K_3 + K_1 , 1);
da = zeros(3*K_3 + K_1 , 1);
pdf = zeros(1,iters);

Rv_inv = pinv(Rv);
X_w = squeeze(sum(X,1));

sin_alpha = sin(alpha);
cos_alpha = cos(alpha);

for i = 1 : iters - 1
    
    FIM = 0;
    df_theta = 0;

    % calculate the a and a derivative 
    sin_theta = sin(theta_estimate(i));
    cos_theta = cos(theta_estimate(i));

    e_u = [sin_theta * sin_alpha ; cos_theta * sin_alpha ; cos_alpha];

    for j = 1 : 3 : 3 * K_3

        r_k_i = r_k(j,:);
        tau_3D = r_k_i * e_u/v_0;
        tau_3D_der = (1/v_0) * (r_k_i(1) * cos_theta * sin_alpha - r_k_i(2) * sin_theta * sin_alpha);
        e_i = exp(-1i * f * tau_3D);

        da(j) = sin_alpha * e_i * (cos_theta - 1i * f * sin_theta * tau_3D_der);
        da(j+1) = sin_alpha * e_i * (-sin_theta - 1i * f * cos_theta * tau_3D_der);
        da(j+2) = cos_alpha * e_i * (-1i * f * tau_3D_der);

        a(j) = sin_theta * sin_alpha * e_i;
        a(j+1) = cos_theta * sin_alpha * e_i;
        a(j+2) = cos_alpha * e_i;

    end

    for j = 3 * K_3 + 1 : 3 * K_3 + K_1

        r_k_i = r_k(j,:);
        tau_1D = r_k_i * e_u/v_0;
        tau_1D_der = (1/v_0) * (r_k_i(1) * cos_theta * sin_alpha - r_k_i(2) * sin_theta * sin_alpha);
        e_i = exp(-1i * f * tau_1D);

        da(j) = cos_alpha * e_i * (-1i * f * tau_1D_der);
        a(j) = cos_alpha * e_i;

    end

    for j = 1 : M

        mue = a * s(j);
        dmue = da * s(j);
        x_i = X_w(:,j);

        df_theta = df_theta + (x_i - mue)' * ...
        Rv_inv * dmue + (dmue' * Rv_inv * (x_i - mue)).';

        FIM = FIM + real(dmue' * Rv_inv * dmue);

    end

    FIM = FIM * 2;
    pdf(i+1) = log_likelihood(X_w,mue,Rv_inv,M);

    theta_estimate(i+1) = wrapToPi(real(theta_estimate(i) + step_size * FIM ^ (-1) * df_theta));

    if (pdf(i+1) < pdf(i))
        step_size = step_size * gamma;
    end

end
end