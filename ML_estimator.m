close all
clc

verbose = false;
%oldpath = addpath('./LIBRA', '-end');
%https://wis.kuleuven.be/stat/robust/LIBRAfiles/LIBRA-home-orig
%system('conda activate obspy & python getData.py'); %uncomment to run the python data getter
data = load("data.mat");
x = data.data;
[M,N] = size(x); %M is number of sensors
L = 128;
P = ceil(N/L);
X_w = zeros(P, M, L); %X_w(p,:,m) to address the pth segment and mth frequency
for i = 0:P-1
    X_w(i+1,:,:) = fft(x(:,i*L+1:(i+1)*L), L, 2);
end
X_w_pr = fftshift(X_w, 3);
f = 40 * (-L/2:L/2-1)/L;

if(verbose)
    for p = (2:6)
        figure;
        for i = 1:M
            subplot(3,6,i);plot(f, abs(squeeze(X_w_pr(p,i,:))));title(['fft of channel ' num2str(i) newline 'at ' num2str(p) 'th segment']);
        end
    end
end

sum_X_w_p = @(X) squeeze(sum(X,1));
a = model(data.r_m, 16, 1, f);
mcdRv = mcdcov(x.','cor', 1, 'plots', 0);
R = mcdRv.cov;
fu = toMaximize(a, R, X_w, M, P);
%tic;[theta_max, alpha_max, v0_max] = Maximize3BySearch(fu, 0.1, 0.1, 1);toc;
thet = MaximizeTheta(fu, 1, 1, 0.1);


function [fun] = toMaximize(a, R, X_w, M, P)
    R_inv = inv(R);
    X_w_p = squeeze(sum(X_w,1));
    norm = @(m, theta, alpha, v_0) a(m, theta, alpha, v_0)' * R_inv * a(m, theta, alpha, v_0);
    upper = @(m, theta, alpha, v_0) a(m, theta, alpha, v_0)' * R_inv * X_w_p(:,m);
    f = @(m, theta, alpha, v_0) 1/P * upper(m, theta, alpha, v_0)' * upper(m, theta, alpha, v_0);
    fun = @(theta, alpha, v_0) 0;
    for m = 1:M
        fun = @(theta, alpha, v_0) fun(theta, alpha, v_0) + f(m, theta, alpha, v_0)/norm(m, theta, alpha, v_0);
    end
    fun = @(theta, alpha, v_0) real(fun(theta, alpha, v_0));
end

function [res] = applyVectorToFunction(func, theta_vec, alpha_vec, v0_vec)
    res = gpuArray(zeros(numel(theta_vec), numel(alpha_vec), numel(v0_vec)));
    for i = 1:numel(theta_vec)
        for j = 1:numel(alpha_vec)
            for k = 1:numel(v0_vec)
                res(i, j, k) = func(theta_vec(i), alpha_vec(j), v0_vec(k));
            end
        end
    end
    res = gather(res);

end


function [theta_max, alpha_max, v0_max] = Maximize3BySearch(func, theta_acc, alpha_acc, v0_acc)
    theta_vec = -pi:theta_acc:pi;
    alpha_vec = -pi:alpha_acc:0;
    v0_vec = 1:v0_acc:100;
    [Theta_vec, Alpha_vec, V0_vec] = meshgrid(theta_vec, alpha_vec, v0_vec);
    res = arrayfun(func, Theta_vec, Alpha_vec, V0_vec);
    [~, I] = max(res, [], 'all');
    [i, j, k] = ind2sub(size(res), I);
    theta_max = theta_vec(i);
    alpha_max = alpha_vec(j);
    v0_max = v0_vec(k);
end

function [val] = MaximizeTheta(fun, alpha, v_0, acc)
    t_vec = -pi:acc:pi;
    Alphas = ones(size(t_vec)) * alpha;
    V_0 = ones(size(t_vec)) * v_0;
    max_vals = arrayfun(fun, t_vec, Alphas, V_0);
    [~, I] = max(max_vals);
    val = t_vec(I);
end

function [val] = MaximizeV(fun, theta, alpha, acc)
    v_vec = 1:acc:100;
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


