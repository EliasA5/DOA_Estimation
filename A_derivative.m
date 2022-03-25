function [da_theta] = A_derivative(v_0,phi,r_k,K_3,K_1,theta,omega_m)
% This function calculates the Mean Mue's derivative
% theta and phi are in radians
da_theta = zeros(3*K_3 + K_1 ,1) ;
sin_phi = sin(phi);
cos_phi = cos(phi);
sin_theta = sin(theta);
cos_theta = cos(theta);
e_u = [sin_theta*sin_phi ; cos_theta*sin_phi ; cos_phi];

for i = 1: 3 : 3*K_3
    r_k_i = r_k(i,:);
    tau_3D = r_k_i * e_u/v_0;
    tau_3D_der = -1i*omega_m*sin_phi*(r_k_i(1)*cos_theta-r_k_i(2)*sin_theta)/v_0;
    e_i = exp(-1i*omega_m*tau_3D);
    da_theta(i) = sin_phi*e_i*(cos_theta+sin_theta*tau_3D_der);
    da_theta(i+1) = sin_phi*e_i*(-sin_theta+cos_theta*tau_3D_der);
    da_theta(i+2) = cos_phi*e_i*(tau_3D_der);
end

for i = 3*K_3+1 : 3*K_3+K_1
    r_k_i = r_k(i,:);
    tau_1D = r_k_i * e_u/v_0;
    tau_1D_der = -1i*omega_m*sin_phi*(r_k_i(1)*cos_theta-r_k_i(2)*sin_theta)/v_0;
    e_i = exp(-1i*omega_m*tau_1D);
    da_theta(i) = cos_phi*e_i*(tau_1D_der);
end
end
