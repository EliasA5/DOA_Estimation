close all
clear all
clc

verbose = false;
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
f = 4 * (-L/2:L/2-1)/L;

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


