close all
clear all
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
%mcdRv = mcdcov(x.','cor', 1, 'plots', 0);
%R = mcdRv.cov;
R_real = robustcov(x.');    % Covariance based on the assumption of correlation 
R_white = eye(M);   % no need for variance because it cancels out
theta_real = pi/4;
%------------------------------------------------------
fu_real = toMaximize(a, R_real, X_w, M, P);
fu_white = toMaximize(a, R_white, X_w, M, P);

theta_max_real = MaximizeTheta(fu_real, 1, 30, 0.001);
theta_max_white = MaximizeTheta(fu_white, 1, 30, 0.001);




%------------------------------------------------------
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

function [val] = MaximizeTheta(fun, alpha, v_0, acc)
    t_vec = -pi:acc:pi;
    val = -pi;
    curr_max = 0;
    for t = t_vec
        tmp = fun(t, alpha, v_0);
        if(tmp > curr_max)
            curr_max = tmp;
            val = t;
        end
    end
end

function [val] = MaximizeV(fun, theta, alpha, acc)
    v_vec = 1:acc:100;
    val = 1;
    curr_max = 0;
    for v = v_vec
        tmp = fun(theta, alpha, v);
        if(tmp > curr_max)
            curr_max = tmp;
            val = v;
        end
    end
end

function [val] = MaximizeAlpha(fun, theta, v_0, acc)
    a_vec = -pi:acc:0;
    val = -pi;
    curr_max = 0;
    for a = a_vec
        tmp = fun(theta, a, v_0);
        if(tmp > curr_max)
            curr_max = tmp;
            val = a;
        end
    end
end


