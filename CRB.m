function [crb] = CRB(v_0,phi,r_k,K_3,K_1,theta,omega_m,R_v,M,s)
%This function calculates the CRB
da = A_derivative(v_0,phi,r_k,K_3,K_1,theta,omega_m);
J = FIM(s,M,R_v,da);
crb = 1/J;
end