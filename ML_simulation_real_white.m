%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: runs the ML estimator on real data assuming white noise.
% Note: Since this simulation assumes white noise we can take Rv to be the
% identity matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_dir = ["./res/"];
simulation_property = ["filter_1_response_16secs", "filter_2_response_16secs", "filter_3_response_16secs", "response_16secs"];

for sim_prop = simulation_property
files = dir(append('./matFiles_' , sim_prop, '/*.mat'));
estimated_theta = [];
real_thetas = [];
estimated_error_cyclic = [];
estimated_error_MSPE = [];
real_errors = [];
estimated_alphas = [];
snrs = [];
L = 128;
alpha_accuracy = 0.01;
accuracy = 0.001;
j = 0;
limit = false;
types = ["MLE_WHITE", "BEAMFORMER", "FISHER_SCORING", "FISHER_SCORING_ORIG"];
%loop through all mat files
for type = types
tic;
for file = files'
    load(fullfile(file.folder, file.name));
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
        snrs = [snrs, snr(i)];
        alpha_data = pi/3;
        real_v0 = 1/real_slowness;
        [~, alpha_data, ~] = estimator(signal, L, r_m, alpha_accuracy, real_theta, real_v0, alpha_data, 'alpha', Rv, 'MLE_WHITE');
        estimated_alphas = [estimated_alphas, alpha_data];
        [theta_est, alpha_est, v0_est] = estimator(signal, L, r_m, accuracy, real_theta, real_v0, alpha_data, 'theta', Rv, type);
        estimated_theta = [estimated_theta, theta_est];
        real_thetas = [real_thetas, real_theta];
        estimated_error_cyclic = [estimated_error_cyclic, MSPE(real_theta, theta_est, 'cyclic')];
        estimated_error_MSPE = [estimated_error_MSPE, MSPE(real_theta, theta_est, 'MSPE')];
        real_errors = [real_errors, real_error];
        j = j+1;
    end

    if(limit && j == 3), break; end
end
%close(f)
toc;
fprintf("finished %d samples using %s, sample type: %s\n", length(real_thetas), type, sim_prop);
res = dir(append(result_dir, type, '_simulation_real_white_results_', sim_prop, '_*.mat'));
save(append(result_dir, type, '_simulation_real_white_results_', sim_prop, '_', string(length(res)+1)), 'estimated_theta','real_thetas','estimated_error_cyclic','estimated_error_MSPE','real_errors','estimated_alphas', 'snrs');
estimated_theta = [];
real_thetas = [];
estimated_error_cyclic = [];
estimated_error_MSPE = [];
real_errors = [];
estimated_alphas = [];
snrs = [];
end

end
delete(gcp);
clear all;
