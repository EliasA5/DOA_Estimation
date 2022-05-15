

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Prints the results of the files named:
% [method_name]_simulation_real[_white]?.m
% Instructions: load the output file of simulation then run
% this function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('./res/BEAMFORMER_simulation_real_white_results_3.mat'); %uncomment and replace i with requested results
save_pics = false;
samples_for_mean = 30;
fig2 = figure;sgtitle("error histrograms");
subplot(4,2,1);histogram(estimated_error_cyclic);title("cyclic beamformer");
xlabel("bin");ylabel("error value")
subplot(4,2,2);histogram(estimated_error_MSPE);title("MSPE beamformer");
xlabel("bin");ylabel("error value")

fig = figure;sgtitle("errors")
%% figure to errors
subplot(1,4,1);title('beamformer');
cyclic = mean(reshape(estimated_error_cyclic(1:length(estimated_error_cyclic)-mod(length(estimated_error_cyclic),samples_for_mean)), samples_for_mean, []), 1);
MSPE = mean(reshape(estimated_error_MSPE(1:length(estimated_error_cyclic)-mod(length(estimated_error_cyclic),samples_for_mean)), samples_for_mean, []), 1);
real_err = mean(reshape(real_errors(1:length(real_errors)-mod(length(estimated_error_cyclic),samples_for_mean)), samples_for_mean, []), 1);
hold on;plot(cyclic);plot(MSPE);plot(real_err);
xlabel('# of sample');ylabel('mean error');
legend('cyclic', 'MSPE', 'real');

load('./res/MLE_WHITE_simulation_real_white_results_2.mat');
figure(fig2);
subplot(4,2,3);histogram(estimated_error_cyclic);title("cyclic ML");
xlabel("bin");ylabel("error value")
subplot(4,2,4);histogram(estimated_error_MSPE);title("MSPE ML");
xlabel("bin");ylabel("error value")

figure(fig);
subplot(1,4,2);title('ML real data assuming white noise');
cyclic = mean(reshape(estimated_error_cyclic(1:length(estimated_error_cyclic)-mod(length(estimated_error_cyclic),samples_for_mean)), samples_for_mean, []), 1);
MSPE = mean(reshape(estimated_error_MSPE(1:length(estimated_error_cyclic)-mod(length(estimated_error_cyclic),samples_for_mean)), samples_for_mean, []), 1);
real_err = mean(reshape(real_errors(1:length(real_errors)-mod(length(estimated_error_cyclic),samples_for_mean)), samples_for_mean, []), 1);
hold on;plot(cyclic);plot(MSPE);plot(real_err);
xlabel('# of sample');ylabel('mean error');
legend('cyclic', 'MSPE', 'real');

load('./res/FISHER_SCORING_simulation_real_white_results_2.mat');
figure(fig2);
subplot(4,2,5);histogram(estimated_error_cyclic);title("cyclic fisher scoring");
xlabel("bin");ylabel("error value")
subplot(4,2,6);histogram(estimated_error_MSPE);title("MSPE fisher scoring");
xlabel("bin");ylabel("error value")

figure(fig);
subplot(1,4,3);title('Fisher scoring real data assuming white noise');
cyclic = mean(reshape(estimated_error_cyclic(1:length(estimated_error_cyclic)-mod(length(estimated_error_cyclic),samples_for_mean)), samples_for_mean, []), 1);
MSPE = mean(reshape(estimated_error_MSPE(1:length(estimated_error_cyclic)-mod(length(estimated_error_cyclic),samples_for_mean)), samples_for_mean, []), 1);
real_err = mean(reshape(real_errors(1:length(real_errors)-mod(length(estimated_error_cyclic),samples_for_mean)), samples_for_mean, []), 1);
hold on;plot(cyclic);plot(MSPE);plot(real_err);
xlabel('# of sample');ylabel('mean error');
legend('cyclic', 'MSPE', 'real');

load('./res/FISHER_SCORING_ORIG_simulation_real_white_results_1.mat');
figure(fig2);
subplot(4,2,7);histogram(estimated_error_cyclic);title("cyclic orig fisher scoring");
xlabel("bin");ylabel("error value")
subplot(4,2,8);histogram(estimated_error_MSPE);title("MSPE orig fisher scoring");
xlabel("bin");ylabel("error value")

figure(fig);
subplot(1,4,4);title('Orig Fisher scoring real data assuming white noise');
cyclic = mean(reshape(estimated_error_cyclic(1:length(estimated_error_cyclic)-mod(length(estimated_error_cyclic),samples_for_mean)), samples_for_mean, []), 1);
MSPE = mean(reshape(estimated_error_MSPE(1:length(estimated_error_cyclic)-mod(length(estimated_error_cyclic),samples_for_mean)), samples_for_mean, []), 1);
real_err = mean(reshape(real_errors(1:length(real_errors)-mod(length(estimated_error_cyclic),samples_for_mean)), samples_for_mean, []), 1);
hold on;plot(cyclic);plot(MSPE);plot(real_err);
xlabel('# of sample');ylabel('mean error');
legend('cyclic', 'MSPE', 'real');

fig.Position(3:4) = [1600, 480];
movegui(fig, 'north');
if (save_pics)
for i=1:8
    name = num2str(i);
    saveas(figure(i),append('./pics/',name),'png')
end
end