close all
clc

oldpath = addpath('./LIBRA', '-end');
files = dir('./matFiles/*.mat');
estimated_theta = [];
real_thetas = [];
estimated_error = [];
real_errors = [];
L = 128;
accuracy = 0.001;
alpha_data = 1;
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
    end
    %for each arrival
    parfor i=1:data_size(1)
        signal = squeeze(data(i,:,:));
        r_m = squeeze(distances(i,:,:));
        real_theta = doa(i)*pi/180;
        real_error = err(i);
        real_slowness = slow(i);
        [theta_est, alpha_est, v0_est] = ML_estimator(signal, L, r_m, accuracy, real_theta, 1/real_slowness, alpha_data, 'theta');
        estimated_theta = [estimated_theta, theta_est];
        real_thetas = [real_thetas, real_theta];
        estimated_error = [estimated_error, mod(sqrt(abs(real_theta^2-theta_est^2)), 2*pi)];
        real_errors = [real_errors, real_error];
    end
    j = j+1;
    if(limit && j == 10), break; end
end
%close(f)
save('results', 'estimated_theta','real_thetas','estimated_error','real_errors');
delete(gcp);
clear all;