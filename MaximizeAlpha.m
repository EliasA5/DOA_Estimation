function [val] = MaximizeAlpha(fun, theta, v_0, acc)
    a_vec = -pi/2:acc:pi/2;
    Thetas = ones(size(a_vec)) * theta;
    V_0 = ones(size(a_vec)) * v_0;
    max_vals = arrayfun(fun, Thetas, a_vec, V_0);
    [~, I] = max(max_vals);
    val = a_vec(I);

end