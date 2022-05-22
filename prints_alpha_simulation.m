%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Prints the results of the files named:
% [method_name]_simulation_real[_white]?.m
% Instructions: load the output file of simulation then run
% this function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_pics = false;
samples_for_mean = 30;
result_dir = ["./res/"];
simulation_properties = ["alpha_effect"];
getFullDir = @(file) append(file.folder, '\', file.name);
getInfo = @(error) sprintf("\nmean: %f, std: %f", mean(error), std(error));

%fig1=figure;title("SNR histogram");
%histogram(snrs);

for simulation_property = simulation_properties

fig_beamformer = figure;sgtitle(append("alpha effect ", "beamformer"));
ylabel("err");xlabel("alpha");
fig_mle_whites = figure;sgtitle(append("alpha effect ", "MLE white"));
ylabel("err");xlabel("alpha");
fig_periodic_fisher_scorings = figure;sgtitle(append("alpha effect ", "Periodic fisher scoring"));
ylabel("err");xlabel("alpha");
fig_original_fisher_scorings = figure;sgtitle(append("alpha effect ", "regular fisher scoring"));
ylabel("err");xlabel("alpha");

beamformers = dir(append(result_dir, 'BEAMFORMER_simulation_real_white_results_', simulation_property, '*'));
mle_whites = dir(append(result_dir, 'MLE_WHITE_simulation_real_white_results_', simulation_property, '*'));
periodic_fisher_scorings = dir(append(result_dir, 'FISHER_SCORING_simulation_real_white_results_', simulation_property, '*'));
original_fisher_scorings = dir(append(result_dir, 'FISHER_SCORING_ORIG_simulation_real_white_results_', simulation_property, '*'));

est_alphas_mean = [];
est_errors_mean = [];
est_alpha = [];
est_error = [];

for i=1:length(beamformers)
load(getFullDir(beamformers(i)));
est_alphas_mean = [est_alphas_mean; mean(estimated_alphas(1:max(1, end-1), :), 1)];
est_errors_mean = [est_errors_mean; mean(estimated_error_cyclic(1:max(1, end-1), :), 1)];
est_alpha = [est_alpha, estimated_alphas(end)];
est_error = [est_error, estimated_error_cyclic(end)];
end
figure(fig_beamformer); hold on;
stem(mean(est_alphas_mean,1), mean(est_errors_mean,1), 'o');
stem(est_alpha, est_error, 'x');
legend("constant alphas", "estimated alphas")

est_alphas_mean = [];
est_errors_mean = [];
est_alpha = [];
est_error = [];

for i=1:length(mle_whites)
load(getFullDir(mle_whites(i)));
est_alphas_mean = [est_alphas_mean; mean(estimated_alphas(1:max(1, end-1), :), 1)];
est_errors_mean = [est_errors_mean; mean(estimated_error_cyclic(1:max(1, end-1), :), 1)];
est_alpha = [est_alpha, estimated_alphas(end)];
est_error = [est_error, estimated_error_cyclic(end)];
end
figure(fig_mle_whites); hold on;
stem(mean(est_alphas_mean,1), mean(est_errors_mean,1), 'o');
stem(est_alpha, est_error, 'x');
legend("constant alphas", "estimated alphas")

est_alphas_mean = [];
est_errors_mean = [];
est_alpha = [];
est_error = [];

for i=1:length(periodic_fisher_scorings)
load(getFullDir(periodic_fisher_scorings(i)));
est_alphas_mean = [est_alphas_mean; mean(estimated_alphas(1:max(1, end-1), :), 1)];
est_errors_mean = [est_errors_mean; mean(estimated_error_cyclic(1:max(1, end-1), :), 1)];
est_alpha = [est_alpha, estimated_alphas(end)];
est_error = [est_error, estimated_error_cyclic(end)];
end
figure(fig_periodic_fisher_scorings); hold on;
stem(mean(est_alphas_mean,1), mean(est_errors_mean,1), 'o');
stem(est_alpha, est_error, 'x');
legend("constant alphas", "estimated alphas")

est_alphas_mean = [];
est_errors_mean = [];
est_alpha = [];
est_error = [];

for i=1:length(original_fisher_scorings)
load(getFullDir(original_fisher_scorings(i)));
est_alphas_mean = [est_alphas_mean; mean(estimated_alphas(1:max(1, end-1), :), 1)];
est_errors_mean = [est_errors_mean; mean(estimated_error_cyclic(1:max(1, end-1), :), 1)];
est_alpha = [est_alpha, estimated_alphas(end)];
est_error = [est_error, estimated_error_cyclic(end)];
end
figure(fig_original_fisher_scorings); hold on;
stem(mean(est_alphas_mean,1), mean(est_errors_mean,1), 'o');
stem(est_alpha, est_error, 'x');
legend("constant alphas", "estimated alphas")


% fig.Position(3:4) = [1600, 480];
% movegui(fig, 'north');

end


fig_beamformer.WindowState = "maximize";
fig_mle_whites.WindowState = "maximize";
fig_periodic_fisher_scorings.WindowState = "maximize";
fig_original_fisher_scorings.WindowState = "maximize";

if (save_pics)
for i=1:3
    name = append("simulation_real_white_1min_", string(i));
    saveas(figure(i),append('./pics/',name),'png')
end
end
