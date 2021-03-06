%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: runs the fisher scoring estimator on real data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
files = dir('./matFiles/*.mat');
estimated_theta = [];
real_thetas = [];
estimated_error_cyclic = [];
estimated_error_MSPE = [];
real_errors = [];
estimated_alphas = [];
L = 128;
alpha_accuracy = 0.01;
accuracy = 0.001;
j = 1;
limit = false;
%loop through all mat files
%f = waitbar(0,'Please wait...');
for file = files'
    load(fullfile(file.folder, file.name));
    %waitbar(j/length(files), f, append('working on: ', file.name, ', iter: ', string(j)));
    data_size = size(data);
    err_size = size(err);
    if(err_size(2) ~= 1)
        disp(append('something is wrong with: ', file.name))
        continue
    end
    %for each arrival
    parfor i=1:data_size(1)
        signal = squeeze(data(i,:,:));
        Rv = eye(height(signal));
        r_m = squeeze(distances(i,:,:));
        real_theta = doa(i)*pi/180;
        real_error = err(i)*pi/180;
        real_slowness = slow(i);
        alpha_data = 1;
        real_v0 = 1/real_slowness;
        [~, alpha_data, ~] = estimator(signal, L, r_m, alpha_accuracy, real_theta, real_v0, alpha_data, 'alpha', Rv, 'MLE_WHITE');
        estimated_alphas = [estimated_alphas, alpha_data];
        [theta_est, alpha_est, v0_est] = estimator(signal, L, r_m, accuracy, real_theta, real_v0, alpha_data, 'theta', Rv, 'FISHER_SCORING');
        estimated_theta = [estimated_theta, theta_est];
        real_thetas = [real_thetas, real_theta];
        estimated_error_cyclic = [estimated_error_cyclic, MSPE(real_theta, theta_est, 'cyclic')];
        estimated_error_MSPE = [estimated_error_MSPE, MSPE(real_theta, theta_est, 'MSPE')];
        real_errors = [real_errors, real_error];
    end
    j = j+1;
    if(limit && j == 3), break; end
end
%close(f)
res = dir('./res/Fisher_scoring_simulation_real_white_*.mat');
save(append('./res/Fisher_scoring_simulation_real_white_', string(length(res)+1)), 'estimated_theta','real_thetas','estimated_error_cyclic','estimated_error_MSPE','real_errors','estimated_alphas');
delete(gcp);
clear all;
