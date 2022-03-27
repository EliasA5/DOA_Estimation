function [val] = MaximizeTheta(fun, alpha, v_0, acc)
    t_vec = 0:acc:2*pi;
    Alphas = ones(size(t_vec)) * alpha;
    V_0 = ones(size(t_vec)) * v_0;
    max_vals = arrayfun(fun, t_vec, Alphas, V_0);
    [~, I] = max(max_vals);
    val = t_vec(I);
end