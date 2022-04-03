
close all;
%load('./res/ML_simulation_results_i.mat'); %uncomment and replace i with requested results
save_pics = false;
%% figure to show the convergance of the methods W-W, C-C
figure;
hold on; grid on;
plot(SNR,theta_og*ones(size(SNR)),'LineWidth', 2)
plot(SNR,ThetaEst_MLE_white_white); plot(SNR,ThetaEst_MLE_colored_colored);
plot(SNR,ThetaEst_fisher_white_white); plot(SNR,ThetaEst_fisher_colored_colored);
legend('Original', 'MLE W-W', 'MLE C-C', 'Fisher W-W', 'Fisher C-C')
xlabel('SNR'); ylabel('\theta'); title('Comparison of the Estimations')
hold off
set(gca,'Xscale','log')
%% figure to show the convergance of the methods W-C
figure;
hold on; grid on;
plot(SNR,theta_og*ones(size(SNR)),'LineWidth', 2)
plot(SNR,ThetaEst_MLE_white_colored); plot(SNR,ThetaEst_MLE_colored_white);
plot(SNR,ThetaEst_fisher_white_colored); plot(SNR,ThetaEst_fisher_colored_white);
legend('Original', 'MLE W-C', 'MLE C-W', 'Fisher W-C', 'Fisher C-W')
xlabel('SNR'); ylabel('\theta'); title('Comparison of the Estimations')
hold off
set(gca,'Xscale','log')
%% figure to compare the MSPE to all kinds of CRB W-W, C-C
figure;
subplot(3,1,1)
hold on; grid on;
plot(SNR,CRB_white_reg,'LineWidth', 2);plot(SNR,CRB_colored_reg,'LineWidth', 2);
plot(SNR,RMSPE_MLE_white_white); plot(SNR,RMSPE_MLE_colored_colored);
plot(SNR,RMSPE_fisher_white_white); plot(SNR,RMSPE_fisher_colored_colored);
legend('CRB W', 'CRB C', 'MLE W-W', 'MLE C-C', 'Fisher W-W', 'Fisher C-C')
xlabel('SNR'); ylabel('Regular MSPE'); title('regular CRB vs. MSPE of the Estimations')
hold off
set(gca,'Xscale','log')
set(gca,'Yscale','log')



subplot(3,1,2)
hold on; grid on;
plot(SNR,CRB_white_cyc1,'LineWidth', 2);plot(SNR,CRB_colored_cyc2,'LineWidth', 2);
plot(SNR,RMSPE_MLE_white_white); plot(SNR,RMSPE_MLE_colored_colored);
plot(SNR,RMSPE_fisher_white_white); plot(SNR,RMSPE_fisher_colored_colored);
legend('CCRB_1 W', 'CCRB_1 C', 'MLE W-W', 'MLE C-C', 'Fisher W-W', 'Fisher C-C')
xlabel('SNR'); ylabel('Regular MSPE'); title('cyclic 1 CRB vs. MSPE of the Estimations')
hold off
set(gca,'Xscale','log')
set(gca,'Yscale','log')


subplot(3,1,3)
hold on; grid on;
plot(SNR,CRB_white_cyc2,'LineWidth', 2);plot(SNR,CRB_colored_cyc2,'LineWidth', 2);
plot(SNR,RMSPE_MLE_white_white); plot(SNR,RMSPE_MLE_colored_colored);
plot(SNR,RMSPE_fisher_white_white); plot(SNR,RMSPE_fisher_colored_colored);
legend('CCRB_2 W', 'CCRB_2 C', 'MLE W-W', 'MLE C-C', 'Fisher W-W', 'Fisher C-C')
xlabel('SNR'); ylabel('Regular MSPE'); title('cyclic 2 CRB vs. MSPE of the Estimations')
hold off
set(gca,'Xscale','log')
set(gca,'Yscale','log')
%% figure to compare the MSPE to all kinds of CRB W-C
figure;
subplot(3,1,1)
hold on; grid on;
plot(SNR,CRB_white_reg,'LineWidth', 2);plot(SNR,CRB_colored_reg,'LineWidth', 2);
plot(SNR,RMSPE_MLE_white_colored); plot(SNR,RMSPE_MLE_colored_white);
plot(SNR,RMSPE_fisher_white_colored); plot(SNR,RMSPE_fisher_colored_white);
legend('CRB W', 'CRB C', 'MLE W-C', 'MLE C-W', 'Fisher W-C', 'Fisher C-W')
xlabel('SNR'); ylabel('Regular MSPE'); title('regular CRB vs. MSPE of the Estimations')
hold off
set(gca,'Xscale','log')
set(gca,'Yscale','log')



subplot(3,1,2)
hold on; grid on;
plot(SNR,CRB_white_cyc1,'LineWidth', 2);plot(SNR,CRB_colored_cyc2,'LineWidth', 2);
plot(SNR,RMSPE_MLE_white_colored); plot(SNR,RMSPE_MLE_colored_white);
plot(SNR,RMSPE_fisher_white_colored); plot(SNR,RMSPE_fisher_colored_white);
legend('CCRB_1 W', 'CCRB_1 C', 'MLE W-C', 'MLE C-W', 'Fisher W-C', 'Fisher C-W')
xlabel('SNR'); ylabel('Regular MSPE'); title('cyclic 1 CRB vs. MSPE of the Estimations')
hold off
set(gca,'Xscale','log')
set(gca,'Yscale','log')


subplot(3,1,3)
hold on; grid on;
plot(SNR,CRB_white_cyc2,'LineWidth', 2);plot(SNR,CRB_colored_cyc2,'LineWidth', 2);
plot(SNR,RMSPE_MLE_white_colored); plot(SNR,RMSPE_MLE_colored_white);
plot(SNR,RMSPE_fisher_white_colored); plot(SNR,RMSPE_fisher_colored_white);
legend('CCRB_2 W', 'CCRB_2 C', 'MLE W-C', 'MLE C-W', 'Fisher W-C', 'Fisher C-W')
xlabel('SNR'); ylabel('Regular MSPE'); title('cyclic 2 CRB vs. MSPE of the Estimations')
hold off
set(gca,'Xscale','log')
set(gca,'Yscale','log')

%% figure to compare the cyclic MSPE to all kinds of CRB W-W, C-C
figure;
subplot(3,1,1)
hold on; grid on;
plot(SNR,CRB_white_reg,'LineWidth', 2);plot(SNR,CRB_colored_reg,'LineWidth', 2);
plot(SNR,CyclicErr_MLE_white_white); plot(SNR,CyclicErr_MLE_colored_colored);
plot(SNR,CyclicErr_fisher_white_white); plot(SNR,CyclicErr_fisher_colored_colored);
legend('CRB W', 'CRB C', 'MLE W-W', 'MLE C-C', 'Fisher W-W', 'Fisher C-C')
xlabel('SNR'); ylabel('Cyclic Error'); title('Regular CRB vs. CMSPE of the Estimations')
hold off;
set(gca,'Xscale','log');set(gca,'Yscale','log')

subplot(3,1,2)
hold on; grid on;
plot(SNR,CRB_white_cyc1,'LineWidth', 2);plot(SNR,CRB_colored_cyc2,'LineWidth', 2);
plot(SNR,CyclicErr_MLE_white_white); plot(SNR,CyclicErr_MLE_colored_colored);
plot(SNR,CyclicErr_fisher_white_white); plot(SNR,CyclicErr_fisher_colored_colored);
legend('CRBC_1 W', 'CRBC_1 C', 'MLE W-W', 'MLE C-C', 'Fisher W-W', 'Fisher C-C')
xlabel('SNR'); ylabel('Cyclic Error'); title('cyclic 1 CRB vs. CMSPE of the Estimations')
hold off
set(gca,'Xscale','log');set(gca,'Yscale','log')

subplot(3,1,3)
hold on; grid on;
plot(SNR,CRB_white_cyc2,'LineWidth', 2);plot(SNR,CRB_colored_cyc2,'LineWidth', 2);
plot(SNR,CyclicErr_MLE_white_white); plot(SNR,CyclicErr_MLE_colored_colored);
plot(SNR,CyclicErr_fisher_white_white); plot(SNR,CyclicErr_fisher_colored_colored);
legend('CRBC_2 W', 'CRBC_2 C', 'MLE W-W', 'MLE C-C', 'Fisher W-W', 'Fisher C-C')
xlabel('SNR'); ylabel('Cyclic Error'); title('cyclic 2 CRB vs. CMSPE of the Estimations')
hold off
set(gca,'Xscale','log');set(gca,'Yscale','log')

%% figure to compare the cyclic MSPE to all kinds of CRB W-C
figure;
subplot(3,1,1)
hold on; grid on;
plot(SNR,CRB_white_reg,'LineWidth', 2);plot(SNR,CRB_colored_reg,'LineWidth', 2);
plot(SNR,CyclicErr_MLE_white_colored); plot(SNR,CyclicErr_MLE_colored_white);
plot(SNR,CyclicErr_fisher_white_colored); plot(SNR,CyclicErr_fisher_colored_white);
legend('CRB W', 'CRB C', 'MLE W-C', 'MLE C-W', 'Fisher W-C', 'Fisher C-W')
xlabel('SNR'); ylabel('Cyclic Error'); title('Regular CRB vs. CMSPE of the Estimations')
hold off;
set(gca,'Xscale','log');set(gca,'Yscale','log')

subplot(3,1,2)
hold on; grid on;
plot(SNR,CRB_white_cyc1,'LineWidth', 2);plot(SNR,CRB_colored_cyc2,'LineWidth', 2);
plot(SNR,CyclicErr_MLE_white_colored); plot(SNR,CyclicErr_MLE_colored_white);
plot(SNR,CyclicErr_fisher_white_colored); plot(SNR,CyclicErr_fisher_colored_white);
legend('CRBC_1 W', 'CRBC_1 C', 'MLE W-C', 'MLE C-W', 'Fisher W-C', 'Fisher C-W')
xlabel('SNR'); ylabel('Cyclic Error'); title('cyclic 1 CRB vs. CMSPE of the Estimations')
hold off
set(gca,'Xscale','log');set(gca,'Yscale','log')

subplot(3,1,3)
hold on; grid on;
plot(SNR,CRB_white_cyc2,'LineWidth', 2);plot(SNR,CRB_colored_cyc2,'LineWidth', 2);
plot(SNR,CyclicErr_MLE_white_colored); plot(SNR,CyclicErr_MLE_colored_white);
plot(SNR,CyclicErr_fisher_white_colored); plot(SNR,CyclicErr_fisher_colored_white);
legend('CRBC_2 W', 'CRBC_2 C', 'MLE W-C', 'MLE C-W', 'Fisher W-C', 'Fisher C-W')
xlabel('SNR'); ylabel('Cyclic Error'); title('cyclic 2 CRB vs. CMSPE of the Estimations')
hold off
set(gca,'Xscale','log');set(gca,'Yscale','log')

if (save_pics)
for i=1:6
    name = num2str(i);
    saveas(figure(i),append('./pics/',name),'png')
end
end