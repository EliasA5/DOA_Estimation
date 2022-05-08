%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: This function runs the estimator given by type on the input
% data, returns all 3 estimated parameters.
% Inputs:
%   x: The signal data in time domain, with (K,N) dimensions where K is the
%   number of sensors and N is the number of samples
%   L: The number of frequencies in each segment, this parameter is used to
%   perform DFT on x vector by dividing x into N/L segments and performing
%   DFT on each segment, preferably L is a power of 2.
%   r_m: The distances matrix of size (3,K).
%   accuracy: Used to define the accuracy of the estimator, for theta/alpha
%   this value can be of magnitude of 1e-3, for v0 estimation this value
%   needs to be on the order of 1e2.
%   theta_data: Starting value for theta, or real value when estimating the
%   other parameters.
%   v0_data: Same as above for v0.
%   alpha_data: Same as above for alpha.
%   to_maximize: This variable is used to decide what value to estimate,
%   can have the next values: ['theta', 'theta_alpha', 'vel', 'alpha',
%   'all'].
%   R: The covariance matrix, can either be supplied from real data by
%   estimating it using MSE or can be the identity matrix (assuming white
%   noise).
%   type: The type of estimator to use, can have the next values: ['MLE',
%   'MLE_WHITE', 'BEAMFORMER', 'FISHER_SCORING'], note that fisher scoring
%   only works with to_maximize=theta.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [theta, alpha, v0] = estimator(x, L, r_m, accuracy, theta_data, v0_data, alpha_data, to_maximize, R, type)

[M,N] = size(x); %M is number of sensors
K_3 = 0; %3*floor((M - min([length(unique(r_m(:,1))), length(unique(r_m(:,2)))]))/3);
K_1 = M - K_3;
P = ceil(N/L);
X_w = zeros(P, M, L); %X_w(p,:,m) to address the pth segment and mth frequency
for i = 0:P-1
    X_w(i+1,:,:) = fft(x(:,i*L+1:min(N, (i+1)*L)), L, 2);
end
f = 40 * (-L/2:L/2-1)/L;

iters = 1e3;
step_size = 1;
gamma = 0.95;
[a, da] = model(r_m, K_1, K_3, f);
% mcdRv = mcdcov(x.','cor', 1, 'plots', 0);
% R = mcdRv.cov;
switch type
    case 'MLE'
        fun = toMaximizeMLE(a, R, X_w, M, P);
    case 'MLE_WHITE'
        fun = toMaximizeMLEWhite(a, [], X_w, M, P);
    case 'BEAMFORMER'
        fun = toMaximizeBeamformer(a, [], X_w, M, P);
    case 'FISHER_SCORING'
        if(strcmp(to_maximize, 'theta'))
            [~, theta] = Fisher_scoring('',theta_data,R,v0_data,alpha_data,K_1+3*K_3,X_w,iters,step_size,gamma,L,P,a,da,accuracy);
            alpha = alpha_data;
            v0 = v0_data;
            return
        else
            error('fisher only works with to_maximize=theta')
        end
end

switch to_maximize
    case 'theta'
        theta = MaximizeTheta(fun, alpha_data, v0_data, accuracy);
        alpha = alpha_data;
        v0 = v0_data;
    case 'theta_alpha'
        [theta, alpha] = MaximizeThetaAlpha(fun, v0_data, accuracy/4);
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
    Thetas = ones(size(v_vec)) * theta;
    Alphas = ones(size(v_vec)) * alpha;
    max_vals = arrayfun(fun, Thetas, Alphas, v_vec);
    [~, I] = max(max_vals);
    val = v_vec(I);
end

function [val] = MaximizeAlpha(fun, theta, v_0, acc)
    a_vec = -pi:acc:0;
    Thetas = ones(size(a_vec)) * theta;
    V_0 = ones(size(a_vec)) * v_0;
    max_vals = arrayfun(fun, Thetas, a_vec, V_0);
    [~, I] = max(max_vals);
    val = a_vec(I);

end



