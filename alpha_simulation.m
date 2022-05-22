close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: runs the ML estimator on real data assuming white noise.
% Note: Since this simulation assumes white noise we can take Rv to be the
% identity matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_dir = ["./res/"]; %#ok<NBRAK> 
simulation_property = ["40secs"]; %#ok<NBRAK> 

for sim_prop = simulation_property
files = dir(append('./matFiles_' , sim_prop, '/*.mat'));
files = files';
estimated_theta = [];
real_thetas = [];
estimated_error_cyclic = [];
estimated_error_MSPE = [];
real_errors = [];
estimated_alphas = [];
snrs = [];

estimated_theta_event = [];
real_thetas_event = [];
estimated_error_cyclic_event = [];
estimated_error_MSPE_event = [];
real_errors_event = [];
estimated_alphas_event = [];
snrs_event = [];


L = 128;
alpha_accuracy = 0.01;
accuracy = 0.001;
j = 0;
limit = false;
types = ["MLE_WHITE", "BEAMFORMER", "FISHER_SCORING", "FISHER_SCORING_ORIG"];
%loop through all mat files
for type = types
%fprintf(2, "on type: %s\n", type);
tic;
for file = files(1:5)
    %fprintf(2, "on file: %s\n", file.name);
    load(fullfile(file.folder, file.name));
    data_size = size(data);
    err_size = size(err);
    if(err_size(2) ~= 1)
        disp(append('something is wrong with: ', file.name))
        continue
    end
    %for each arrival
    for i=1:data_size(1)
        signal = squeeze(data(i,:,:));
        Rv = eye(height(signal));
        r_m = squeeze(distances(i,:,:));
        real_theta = doa(i)*pi/180;
        real_error = err(i)*pi/180;
        real_slowness = slow(i);
        real_v0 = 1/real_slowness;
        alpha_data = pi/3;
        [~, alpha_est, ~] = estimator(signal, L, r_m, alpha_accuracy, real_theta, real_v0, alpha_data, 'alpha', Rv, 'MLE_WHITE');
        arr = [linspace(-pi/2+0.1, pi/2-0.1, 11), alpha_est];
        parfor alpha_data_index=1:length(arr)
            estimated_alphas = [estimated_alphas, arr(alpha_data_index)];
            [theta_est, ~, v0_est] = estimator(signal, L, r_m, accuracy, real_theta, real_v0, arr(alpha_data_index), 'theta', Rv, type);
            estimated_theta = [estimated_theta, theta_est];
            real_thetas = [real_thetas, real_theta];
            estimated_error_cyclic = [estimated_error_cyclic, MSPE(real_theta, theta_est, 'cyclic')];
            estimated_error_MSPE = [estimated_error_MSPE, MSPE(real_theta, theta_est, 'MSPE')];
            real_errors = [real_errors, real_error];
            snrs = [snrs, snr(i)];
        end
        estimated_theta_event = [estimated_theta_event; estimated_theta];
        real_thetas_event = [real_thetas_event; real_thetas];
        estimated_error_cyclic_event = [estimated_error_cyclic_event; estimated_error_cyclic];
        estimated_error_MSPE_event = [estimated_error_MSPE_event; estimated_error_MSPE];
        real_errors_event = [real_errors_event; real_errors];
        estimated_alphas_event = [estimated_alphas_event; estimated_alphas];
        snrs_event = [snrs_event; snrs];

        estimated_theta = [];
        real_thetas = [];
        estimated_error_cyclic = [];
        estimated_error_MSPE = [];
        real_errors = [];
        estimated_alphas = [];
        snrs = [];

        j = j+1;
    end

    if(limit && j == 3), break; end
end
%close(f)
toc;

estimated_theta = estimated_theta_event;
real_thetas = real_thetas_event;
estimated_error_cyclic = estimated_error_cyclic_event;
estimated_error_MSPE = estimated_error_MSPE_event;
real_errors = real_errors_event;
estimated_alphas = estimated_alphas_event;
snrs = snrs_event;

sim_prop = "alpha_effect";
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
