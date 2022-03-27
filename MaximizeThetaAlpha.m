function [theta, alpha] = MaximizeThetaAlpha(fun, v_0, acc)
    a_vec = 0:acc:pi/2;
    t_vec = 0:acc*4:2*pi;
    V_0 = ones(size(t_vec)) * v_0;
    [Theta_vec, Alpha_vec, V0_vec] = meshgrid(t_vec, a_vec, V_0);
    res = arrayfun(fun, Theta_vec, Alpha_vec, V0_vec);
    [~, I] = max(res, [], 'all');
    [i, j, ~] = ind2sub(size(res), I);
    theta = theta_vec(i);
    alpha = alpha_vec(j);
end