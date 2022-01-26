function [df_theta] = PDF_derivative(K_3,K_1,v_0,phi,omega,s,r_k,Rv,x,theta,da,M)
% This function calculates the derivative of the log of the PDF

a = zeros(3*K_3 + K_1 , 1);
R_inv = pinv(Rv);
df_theta = 0;
X_w_p = squeeze(sum(x,1));
sin_phi = sin(phi);
cos_phi = cos(phi);
sin_theta = sin(theta);
cos_theta = cos(theta);
e_u = [sin_theta*sin_phi ; cos_theta*sin_phi ; cos_phi];

for i = 1 : 3 : K_3
    r_k_i = r_k(i,:);
    tau_3D = r_k_i * e_u/v_0;
    a(i) = e_u(1)*exp(-1i*omega*tau_3D);
    a(i + 1) = e_u(2)*exp(-1i*omega*tau_3D);
    a(i + 2) = e_u(3)*exp(-1i*omega*tau_3D);
end

for i = K_3+1:K_3+K_1
    r_k_i = r_k(i,:);
    tau_1D = r_k_i * e_u/v_0;
    a(i) = cos_phi*exp(-1i*omega*tau_1D);
end

for i = 1:M
    mue = a * s(i);
    dmu = da * s(i);
    x_i = X_w_p(:,i);
    df_theta = df_theta + (x_i-mue)' * R_inv * dmu + (x_i-mue).' * R_inv.' * conj(dmu);

end


