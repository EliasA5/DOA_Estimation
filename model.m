function [a] = model(r_m, K_1, K_3, w)
%gets matrix of locations r_m, each row is a 3 elements vector
%correcsponding to a sensor location
%K_1 the number of 1d sensors
%k_2 the number of 3d sensors
%w is a vector of frequencies
e_u = @(theta,alpha) [sin(theta)*sin(alpha); cos(theta)*cos(alpha); cos(alpha)];
e_z = @(theta, alpha) e_u(3);
u = @(theta, alpha, v_0) 1/v0 * e_u(theta,alpha);
tau = @(theta, alpha, v_0)  r_m * u(theta, alpha, v_0);
a = @(m, theta, alpha, v_0) [e_u(theta, alpha) .* exp(-1i*w(m) * tau(1:K_3)), e_z(theta,alpha) .* exp(-1i*w(m) * tau(K_3+1:end))];

end