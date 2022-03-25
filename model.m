function [a] = model(r_m, K_1, K_3, w)
%gets matrix of locations r_m, each row is a 3 elements vector
%correcsponding to a sensor location
%K_1 the number of 1d sensors
%k_2 the number of 3d sensors
%w is a vector of frequencies
e_u = @(theta,alpha) [sin(theta)*sin(alpha); cos(theta)*sin(alpha); cos(alpha)];
e_z = @(theta, alpha) cos(alpha);
u = @(theta, alpha, v_0) 1/v_0 * e_u(theta,alpha);
tau = @(r_ms, theta, alpha, v_0)  r_ms * u(theta, alpha, v_0);
if(K_3 ~= 0)
    a = @(m, theta, alpha, v_0) [repmat(e_u(theta, alpha), K_3, 1) .* exp(-1i*w(m) * tau(r_m(1:K_3*3, :), theta, alpha, v_0)); e_z(theta,alpha) .* exp(-1i*w(m) * tau(r_m(K_3*3+1:end, :), theta, alpha, v_0))];
else
    a = @(m, theta, alpha, v_0) [e_z(theta,alpha) .* exp(-1i*w(m) * tau(r_m, theta, alpha, v_0))];
end

end