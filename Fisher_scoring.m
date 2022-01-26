function [theta_estimate] = Fisher_scoring(theta,s,Rv,omega,v_0,phi,K_3,K_1,x,iter)
% This function estimates the DOA (theta) using Fisher's scoring
[K,M] = size(x);
theta_estimate = theta;
for i = 1:iter
    da =  A_derivative(v_0,phi,r_k,K_3,K_1,theta,omega_m);
    df_theta = PDF_derivative(K_3,K_1,v_0,phi,omega,s,r_k,Rv,x,theta,da,M);
    F = FIM(s,M,R_v,da); %scalar for only theta
    theta_estimate = theta_estimate + df_theta/F;
end

end