%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Prints the results of the files named:
% [method_name]_simulation_real[_white]?.m
% Instructions: load the output file of simulation then run
% this function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_pics = false;
samples_for_mean = 30;
result_dir = ["./res/"];
simulation_properties = ["all_filters"];%["filter_1_response_16secs", "filter_2_response_16secs", "filter_3_response_16secs", "response_16secs"];%["filter_1_16secs", "filter_2_16secs", "filter_3_16secs"];
getFullDir = @(file) append(file.folder, '\', file.name);
getInfo = @(error) sprintf("\nmean: %f, std: %f", rad2deg(mean(error)), rad2deg(std(error)));



for simulation_property = simulation_properties
fig2 = figure;sgtitle(append("error histrograms ", simulation_property));

beamformers = dir(append(result_dir, 'BEAMFORMER_simulation_real_white_results_', simulation_property, '*'));
mle_whites = dir(append(result_dir, 'MLE_WHITE_simulation_real_white_results_', simulation_property, '*'));
periodic_fisher_scorings = dir(append(result_dir, 'FISHER_SCORING_simulation_real_white_results_', simulation_property, '*'));
original_fisher_scorings = dir(append(result_dir, 'FISHER_SCORING_ORIG_simulation_real_white_results_', simulation_property, '*'));

if(~isempty(beamformers))
load(getFullDir(beamformers(end))); %uncomment and replace i with requested results
figure(fig2);
subplot(4,2,1);histogram(estimated_error_cyclic, 0:0.1:4.2);title("cyclic Beamformer" + getInfo(estimated_error_cyclic));
xlabel("bin");ylabel("error value")
subplot(4,2,2);histogram(estimated_error_MSPE, 0:0.2:10);title("MSPE Beamformer" + getInfo(estimated_error_MSPE));
xlabel("bin");ylabel("error value")
end

if(~isempty(mle_whites))
load(getFullDir(mle_whites(end)));
figure(fig2);
subplot(4,2,3);histogram(estimated_error_cyclic, 0:0.1:4.2);title("cyclic ML" + getInfo(estimated_error_cyclic));
xlabel("bin");ylabel("error value")
subplot(4,2,4);histogram(estimated_error_MSPE, 0:0.2:10);title("MSPE ML" + getInfo(estimated_error_MSPE));
xlabel("bin");ylabel("error value")
end

if(~isempty(periodic_fisher_scorings))
load(getFullDir(periodic_fisher_scorings(end)));
figure(fig2);
subplot(4,2,5);histogram(estimated_error_cyclic, 0:0.1:4.2);title("cyclic Periodic Fisher Scoring" + getInfo(estimated_error_cyclic));
xlabel("bin");ylabel("error value")
subplot(4,2,6);histogram(estimated_error_MSPE, 0:0.2:10);title("MSPE Periodic Fisher Scoring" + getInfo(estimated_error_MSPE));
xlabel("bin");ylabel("error value")
end

if(~isempty(original_fisher_scorings))
load(getFullDir(original_fisher_scorings(end)));
figure(fig2);
subplot(4,2,7);histogram(estimated_error_cyclic, 0:0.1:4.2);title("cyclic Fisher Scoring" + getInfo(estimated_error_cyclic));
xlabel("bin");ylabel("error value")
subplot(4,2,8);histogram(estimated_error_MSPE, 0:0.2:10);title("MSPE Fisher Scoring" + getInfo(estimated_error_MSPE));
xlabel("bin");ylabel("error value")
end
% fig.Position(3:4) = [1600, 480];
% movegui(fig, 'north');
fig2.WindowState = "maximize";

end

snrs(snrs>80) = [];
fig1=figure;sgtitle("SNR histograms");
subplot(1,2,1);histogram(snrs);xlabel("bin-edges");ylabel("count");title("in abs value");
subplot(1,2,2);histogram(db(snrs, "power"));xlabel("bin-edges_{db}");ylabel("count");title("in db");
fig1.WindowState = "maximize";

if (save_pics)
for i=1:7
    name = append("simulation_real_white_timings_", string(i));
    saveas(figure(i),append('./pics/',name),'png')
end
end
