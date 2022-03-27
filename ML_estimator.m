
function [theta, alpha, v0] = ML_estimator(x, L, r_m, accuracy, theta_data, v0_data, alpha_data, to_maximize)

[M,N] = size(x); %M is number of sensors
K_3 = 3*floor((M - min([length(unique(rm(:,1))), length(unique(rm(:,2)))]))/3);
K_1 = M - K_3;
P = ceil(N/L);
X_w = zeros(P, M, L); %X_w(p,:,m) to address the pth segment and mth frequency
for i = 0:P-1
    X_w(i+1,:,:) = fft(x(:,i*L+1:min(N, (i+1)*L)), L, 2);
end
f = 40 * (-L/2:L/2-1)/L;

[a, ~] = model(r_m, K_1, K_3, f);
mcdRv = mcdcov(x.','cor', 1, 'plots', 0);
R = mcdRv.cov;
fu = toMaximizeMLE(a, R, X_w, M, P);
switch to_maximize
    case 'theta'
        theta = MaximizeTheta(fu, alpha_data, v0_data, accuracy);
        alpha = alpha_data;
        v0 = v0_data;
    case 'vel'
        theta = theta_data;
        alpha = alpha_data;
        v0 = MaximizeV(fun, theta_data, alpha_data, max(accuracy, 100));
    case 'alpha'
        theta = theta_data;
        alpha = MaximizeAlpha(fun, theta_data, v0_data, accuracy);
        v0 = v0_data;
    case 'all'
        [theta, alpha, v0] = Maximize3BySearch(func, accuracy, 0.1, 100);
    otherwise
        error('to_maximize can get the values: theta/vel/alpha/all')
        return
end
end



function [theta_max, alpha_max, v0_max] = Maximize3BySearch(func, theta_acc, alpha_acc, v0_acc)
    theta_vec = -pi:theta_acc:pi;
    alpha_vec = -pi:alpha_acc:0;
    v0_vec = 1000:v0_acc:10000;
    [Theta_vec, Alpha_vec, V0_vec] = meshgrid(theta_vec, alpha_vec, v0_vec);
    res = arrayfun(func, Theta_vec, Alpha_vec, V0_vec);
    [~, I] = max(res, [], 'all');
    [i, j, k] = ind2sub(size(res), I);
    theta_max = theta_vec(i);
    alpha_max = alpha_vec(j);
    v0_max = v0_vec(k);
end

function [val] = MaximizeV(fun, theta, alpha, acc)
    v_vec = 1000:acc:10000;
    Thetas = ones(size(t_vec)) * theta;
    Alphas = ones(size(t_vec)) * alpha;
    max_vals = arrayfun(fun, Thetas, Alphas, v_vec);
    [~, I] = max(max_vals);
    val = t_vec(I);
end

function [val] = MaximizeAlpha(fun, theta, v_0, acc)
    a_vec = -pi:acc:0;
    Thetas = ones(size(t_vec)) * theta;
    V_0 = ones(size(t_vec)) * v_0;
    max_vals = arrayfun(fun, Thetas, a_vec, V_0);
    [~, I] = max(max_vals);
    val = t_vec(I);

end



