function [a, da] = model(r_m, K_1, K_3, w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: calculates and returns function that model our problem.

% Inputs:
%   r_m: The distances matrix of size (3,K), where K is the number of
%   Sensors.
%   K_1: The number of 1d sensors.
%   k_2: The number of 3d sensors.
%   w: Is a vector of frequencies.

% Returns:
% a: a function that calculates the steering vector.
% da: a function that calculates the derivative of the steering vector.
% Inputs:
%   m: the index of the frequency in the w input vector.
%   theta: doa angle.
%   alpha: incidence angle.
%   v_0: the speed of the wave. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e_u = @(theta,alpha) [sin(theta)*sin(alpha); cos(theta)*sin(alpha); cos(alpha)];
e_z = @(theta, alpha) cos(alpha);
u = @(theta, alpha, v_0) 1/v_0 * e_u(theta,alpha);
tau = @(r_ms, theta, alpha, v_0)  r_ms * u(theta, alpha, v_0);
de_u = @(theta,alpha) [cos(theta)*sin(alpha); -sin(theta)*sin(alpha); 0];
dtau = @(r_ms, theta, alpha, v_0) 1/v_0 * r_ms * de_u(theta, alpha);
if(K_3 ~= 0)
    a = @(m, theta, alpha, v_0) [repmat(e_u(theta, alpha), K_3, 1) .* exp(-1i*w(m) * tau(r_m(1:K_3*3, :), theta, alpha, v_0));...
        e_z(theta,alpha) .* exp(-1i*w(m) * tau(r_m(K_3*3+1:end, :), theta, alpha, v_0))];
    da = @(m, theta, alpha, v_0) -1i*w(m)*dtau(r_m, theta, alpha, v_0) .* a(m, theta, alpha, v_0)...
        + [repmat(de_u(theta, alpha), K_3, 1) .* exp(-1i*w(m) * tau(r_m(1:K_3*3, :), theta, alpha, v_0)); zeros(K_1,1)];
else
    a = @(m, theta, alpha, v_0) [e_z(theta,alpha) .* exp(-1i*w(m) * tau(r_m, theta, alpha, v_0))];
    da = @(m, theta, alpha, v_0) -1i*w(m)*dtau(r_m, theta, alpha, v_0) .* a(m, theta, alpha, v_0);
end
end