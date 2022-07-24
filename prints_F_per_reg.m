
close all;
load('./res/Fvers_10.mat'); %uncomment and replace i with requested results
save_pics = false;
%% figure to show the nember of iterations for each methods W-W, C-C
% figure;
% hold on; grid on;
% plot(SNR,Iter_F_1_white_white); 
% plot(SNR,Iter_F_2_white_white); 
% plot(SNR,Iter_F_3_white_white); 
% legend( 'Reg W-W',  'Per1 W-W','Per2 W-W')
% xlabel('SNR'); ylabel('Iterations'); title('Comparison of the Estimations')
% hold off
% 
% figure
% hold on
% plot(SNR,Iter_F_1_colored_colored);
% plot(SNR,Iter_F_2_colored_colored);
% plot(SNR,Iter_F_3_colored_colored);
% legend('Reg C-C', 'Per1 C-C', 'Per2 C-C')
% xlabel('\theta'); ylabel('Iterations'); title('Comparison of the Estimations')
% hold off
% % set(gca,'Xscale','log')
%% figure to show the convergance of the methods W-W, C-C
figure;
hold on; grid on;
plot(SNR,theta_og*ones(size(SNR)),'LineWidth', 2)
plot(SNR,ThetaEst_F_1_white_white); %plot(SNR,ThetaEst_F_1_colored_colored);
plot(SNR,ThetaEst_F_2_white_white); %plot(SNR,ThetaEst_F_2_colored_colored);
plot(SNR,ThetaEst_F_3_white_white); %plot(SNR,ThetaEst_F_3_colored_colored);
legend('Original', 'Reg W-W',  'Per1 W-W','Per2 W-W')
xlabel('SNR'); ylabel('\theta'); title('Comparison of the Estimations')
hold off
set(gca,'Xscale','log')

figure;
hold on; grid on;
plot(SNR,theta_og*ones(size(SNR)),'LineWidth', 2)
plot(SNR,ThetaEst_F_1_colored_colored);
plot(SNR,ThetaEst_F_2_colored_colored);
plot(SNR,ThetaEst_F_3_colored_colored);
legend('Original', 'Reg C-C', 'Per1 C-C', 'Per2 C-C')
xlabel('SNR'); ylabel('\theta'); title('Comparison of the Estimations')
hold off
set(gca,'Xscale','log')
%% figure to show the convergance of the methods W-C
figure;
hold on; grid on;
plot(SNR,theta_og*ones(size(SNR)),'LineWidth', 2)
plot(SNR,ThetaEst_F_1_white_colored); %plot(SNR,ThetaEst_F_1_colored_white);
plot(SNR,ThetaEst_F_2_white_colored); %plot(SNR,ThetaEst_F_2_colored_white);
plot(SNR,ThetaEst_F_3_white_colored); %plot(SNR,ThetaEst_F_3_colored_white);
legend('Original', 'Reg W-C','Per1 W-C','Per2 W-C')
xlabel('SNR'); ylabel('\theta'); title('Comparison of the Estimations')
hold off
set(gca,'Xscale','log')

figure;
hold on; grid on;
plot(SNR,theta_og*ones(size(SNR)),'LineWidth', 2)
plot(SNR,ThetaEst_F_1_colored_white);
plot(SNR,ThetaEst_F_2_colored_white);
plot(SNR,ThetaEst_F_3_colored_white);
legend('Original','Reg C-W','Per1 C-W','Per2 C-W')
xlabel('SNR'); ylabel('\theta'); title('Comparison of the Estimations')
hold off
set(gca,'Xscale','log')
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% figure to compare the MSPE to all kinds of CRB W-W, C-C
figure;
hold on; grid on;
plot(SNR,CRB_white_reg,'LineWidth', 2);plot(SNR,CRB_white_cyc1,'LineWidth', 2);plot(SNR,CRB_white_cyc2,'LineWidth', 2);
plot(SNR,RMSPE_F_1_white_white); %plot(SNR,RMSPE_F_1_colored_colored);
plot(SNR,RMSPE_F_2_white_white); %plot(SNR,RMSPE_F_2_colored_colored);
plot(SNR,RMSPE_F_3_white_white); %plot(SNR,RMSPE_F_3_colored_colored);
legend('CRB W', 'CCRB1 W','CCRB2 W','Reg W-W','Per1 W-W', 'Per2 W-W')
xlabel('SNR'); ylabel('Regular MSPE'); title('MSPE of the Estimations & CRB (Cyclic and Regular')
hold off
set(gca,'Xscale','log')
set(gca,'Yscale','log')

figure;
hold on; grid on;
plot(SNR,CRB_colored_reg,'LineWidth', 2);plot(SNR,CRB_colored_cyc1,'LineWidth', 2);plot(SNR,CRB_colored_cyc2,'LineWidth', 2);
plot(SNR,RMSPE_F_1_colored_colored);
plot(SNR,RMSPE_F_2_colored_colored);
plot(SNR,RMSPE_F_3_colored_colored);
legend('CRB C','CCRB1 C','CCRB2 C', 'Reg C-C', 'Per1 C-C', 'Per2 C-C')
xlabel('SNR'); ylabel('Regular MSPE'); title('MSPE of the Estimations & CRB (Cyclic and Regular')
hold off
set(gca,'Xscale','log')
set(gca,'Yscale','log')
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% figure to compare the MSPE to all kinds of CRB W-C
figure;
hold on; grid on;
plot(SNR,CRB_colored_reg,'LineWidth', 2);plot(SNR,CRB_colored_cyc1,'LineWidth', 2);plot(SNR,CRB_colored_cyc2,'LineWidth', 2);
plot(SNR,RMSPE_F_1_white_colored); % White_colored means that the model is colored but the signal is white
plot(SNR,RMSPE_F_2_white_colored);% plot(SNR,RMSPE_F_2_colored_white);
plot(SNR,RMSPE_F_3_white_colored); %plot(SNR,RMSPE_F_3_colored_white);
legend('CRB C','CCRB_1 C','CCRB_2 C','Reg W-C', 'Per1 W-C','Per2 W-C')
xlabel('SNR'); ylabel('Regular MSPE'); title('MSPE of the Estimations & CRB (Cyclic and Regular')
hold off
set(gca,'Xscale','log')
set(gca,'Yscale','log')

figure;
hold on; grid on;
plot(SNR,CRB_white_reg,'LineWidth', 2);plot(SNR,CRB_white_cyc1,'LineWidth', 2);plot(SNR,CRB_white_cyc2,'LineWidth', 2);
plot(SNR,RMSPE_F_1_colored_white); % colored_white means that the model is white but the signal is colored
plot(SNR,RMSPE_F_2_colored_white);
plot(SNR,RMSPE_F_3_colored_white);
legend('CRB W','CCRB_1 W','CCRB_2 W','Reg C-W','Per1 C-W','Per2 C-W')
xlabel('SNR'); ylabel('Regular MSPE'); title('MSPE of the Estimations & CRB (Cyclic and Regular')
hold off
set(gca,'Xscale','log')
set(gca,'Yscale','log')

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% figure to compare the cyclic MSE to all kinds of CRB W-W, C-C
figure;
hold on; grid on;
plot(SNR,CRB_white_reg,'LineWidth', 2);plot(SNR,CRB_white_cyc1,'LineWidth', 2);plot(SNR,CRB_white_cyc2,'LineWidth', 2);
plot(SNR,CyclicErr_F_1_white_white); %plot(SNR,CyclicErr_F_1_colored_colored);
plot(SNR,CyclicErr_F_2_white_white); %plot(SNR,CyclicErr_F_2_colored_colored);
plot(SNR,CyclicErr_F_3_white_white); %plot(SNR,CyclicErr_F_3_colored_colored);
legend('CRB W','CRBC_1 W','CRBC_2 W','Reg W-W','Per1 W-W', 'Per2 W-W')
xlabel('SNR'); ylabel('Cyclic Error'); title('CMSPE of the Estimations & CRB (Cyclic and Regular')
hold off;
set(gca,'Xscale','log');
set(gca,'Yscale','log')

figure;
hold on; grid on;
plot(SNR,CRB_colored_reg,'LineWidth', 2);plot(SNR,CRB_colored_cyc1,'LineWidth', 2);plot(SNR,CRB_colored_cyc2,'LineWidth', 2);
plot(SNR,CyclicErr_F_1_colored_colored);
plot(SNR,CyclicErr_F_2_colored_colored);
plot(SNR,CyclicErr_F_3_colored_colored);
legend('CRB C', 'CRBC_1 C','CRBC_2 C','Reg C-C', 'Per1 C-C', 'Per2 C-C')
xlabel('SNR'); ylabel('Cyclic Error'); title('CMSPE of the Estimations & CRB (Cyclic and Regular')
hold off;
set(gca,'Xscale','log');
set(gca,'Yscale','log')

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% figure to compare the cyclic MSE to all kinds of CRB W-C
figure;
hold on; grid on;
plot(SNR,CRB_colored_reg,'LineWidth', 2);plot(SNR,CRB_colored_cyc1,'LineWidth', 2);plot(SNR,CRB_colored_cyc2,'LineWidth', 2);
plot(SNR,CyclicErr_F_1_white_colored); % White_colored means that the model is colored but the signal is white
plot(SNR,CyclicErr_F_2_white_colored); %plot(SNR,CyclicErr_F_2_colored_white);
plot(SNR,CyclicErr_F_3_white_colored); %plot(SNR,CyclicErr_F_3_colored_white);
legend('CRB C','CRBC_1 C','CRBC_2 C','Reg W-C','Per1 W-C', 'Per2 W-C')
xlabel('SNR'); ylabel('Cyclic Error'); title('CMSPE of the Estimations & CRB (Cyclic and Regular')
hold off;
set(gca,'Xscale','log');set(gca,'Yscale','log')

figure;
hold on; grid on;
plot(SNR,CRB_white_reg,'LineWidth', 2);plot(SNR,CRB_white_cyc1,'LineWidth', 2);plot(SNR,CRB_white_cyc2,'LineWidth', 2);
plot(SNR,CyclicErr_F_1_colored_white); % colored_white means that the model is white but the signal is colored
plot(SNR,CyclicErr_F_2_colored_white);
plot(SNR,CyclicErr_F_3_colored_white);
legend('CRB W', 'CRBC_1 W','CRBC_2 W','Reg C-W', 'Per1 C-W', 'Per2 C-W')
xlabel('SNR'); ylabel('Cyclic Error'); title('CMSPE of the Estimations & CRB (Cyclic and Regular')
hold off
set(gca,'Xscale','log');set(gca,'Yscale','log')

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% figure to compare the MCE to all kinds of CRB W-W, C-C
figure;
hold on; grid on;
plot(SNR,CRB_white_reg,'LineWidth', 2);plot(SNR,CRB_white_cyc1,'LineWidth', 2);plot(SNR,CRB_white_cyc2,'LineWidth', 2);
plot(SNR,MCE_F_1_white_white); %plot(SNR,MCE_F_1_colored_colored);
plot(SNR,MCE_F_2_white_white); %plot(SNR,MCE_F_2_colored_colored);
plot(SNR,MCE_F_3_white_white); %plot(SNR,MCE_F_2_colored_colored);
legend('CRB W','CRBC_1 W','CRBC_2 W','Reg W-W', 'Per1 W-W','Per2 W-W')
xlabel('SNR'); ylabel('Cyclic Error'); title('MCE of the Estimations & CRB (Cyclic and Regular')
hold off;
set(gca,'Xscale','log');
set(gca,'Yscale','log')


figure;
hold on; grid on;
plot(SNR,CRB_colored_reg,'LineWidth', 2);plot(SNR,CRB_colored_cyc1,'LineWidth', 2);plot(SNR,CRB_colored_cyc2,'LineWidth', 2);
plot(SNR,MCE_F_1_colored_colored);
plot(SNR,MCE_F_2_colored_colored);
plot(SNR,MCE_F_2_colored_colored);
legend('CRB C', 'CRBC_1 C','CRBC_2 C','Reg C-C', 'Per1 C-C', 'Per2 C-C')
xlabel('SNR'); ylabel('Cyclic Error'); title('MCE of the Estimations & CRB (Cyclic and Regular')
hold off;
set(gca,'Xscale','log');
set(gca,'Yscale','log')

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% figure to compare the cyclic MSPE to all kinds of CRB W-C
figure;
hold on; grid on;
plot(SNR,CRB_colored_reg,'LineWidth', 2);plot(SNR,CRB_colored_cyc1,'LineWidth', 2);plot(SNR,CRB_colored_cyc2,'LineWidth', 2);
plot(SNR,MCE_F_1_white_colored); % White_colored means that the model is colored but the signal is white
plot(SNR,MCE_F_2_white_colored); %plot(SNR,MCE_F_2_colored_white);
plot(SNR,MCE_F_3_white_colored); %plot(SNR,MCE_F_3_colored_white);
legend('CRB C','CRBC_1 C','CRBC_2 C','Reg W-C','Per1 W-C', 'Per2 W-C')
xlabel('SNR'); ylabel('Cyclic Error'); title('MCE of the Estimations & CRB (Cyclic and Regular')
hold off;
set(gca,'Xscale','log');
set(gca,'Yscale','log')

figure;
hold on; grid on;
plot(SNR,CRB_white_reg,'LineWidth', 2);plot(SNR,CRB_white_cyc1,'LineWidth', 2);plot(SNR,CRB_white_cyc1,'LineWidth', 2);
plot(SNR,MCE_F_1_colored_white); % colored_white means that the model is white but the signal is colored
plot(SNR,MCE_F_2_colored_white);
plot(SNR,MCE_F_3_colored_white);
legend('CRB W', 'CRBC_1 W','CRBC_2 W','Reg C-W', 'Per1 C-W', 'Per2 C-W')
xlabel('SNR'); ylabel('Cyclic Error'); title('MCE of the Estimations & CRB (Cyclic and Regular')
hold off;
set(gca,'Xscale','log');
set(gca,'Yscale','log')

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%% figure to compare the Bias of the estimators
figure;
hold on; grid on;
plot(SNR,abs(ThetaEst_F_1_colored_colored_mean));plot(SNR,abs(ThetaEst_F_1_colored_white_mean));
plot(SNR,abs(ThetaEst_F_1_white_white_mean)); plot(SNR,abs(ThetaEst_F_1_white_colored_mean));
legend('C-C', 'C-W', 'W-W', 'W-W')
xlabel('SNR'); ylabel('Bias'); title('Bias of Regular Fisher-scoring')
hold off;
set(gca,'Xscale','log');set(gca,'Yscale','log')

figure; 
hold on; grid on;
plot(SNR,abs(ThetaEst_F_2_colored_colored_mean));plot(SNR,abs(ThetaEst_F_2_colored_white_mean));
plot(SNR,abs(ThetaEst_F_2_white_white_mean)); plot(SNR,abs(ThetaEst_F_2_white_colored_mean));
legend('C-C', 'C-W', 'W-W', 'W-W')
xlabel('SNR'); ylabel('Bias'); title('Bias of Periodic 1 Fisher-scoring')
hold off;
set(gca,'Xscale','log');set(gca,'Yscale','log')

figure; 
hold on; grid on;
plot(SNR,abs(ThetaEst_F_3_colored_colored_mean));plot(SNR,abs(ThetaEst_F_3_colored_white_mean));
plot(SNR,abs(ThetaEst_F_3_white_white_mean)); plot(SNR,abs(ThetaEst_F_3_white_colored_mean));
legend('C-C', 'C-W', 'W-W', 'W-W')
xlabel('SNR'); ylabel('Bias'); title('Bias of Periodic 2 Fisher-scoring')
hold off;
set(gca,'Xscale','log');set(gca,'Yscale','log')

%%
figure; 
hold on; grid on;
plot(SNR,abs(ThetaEst_F_1_white_white_mean));
plot(SNR,abs(ThetaEst_F_2_white_white_mean));
plot(SNR,abs(ThetaEst_F_3_white_white_mean));
legend('reg', 'per1', 'per2')
xlabel('SNR'); ylabel('Bias'); title('Bias of verions for Fisher-scoring')
hold off;
set(gca,'Xscale','log');set(gca,'Yscale','log')
%%
if (save_pics)
for i=1:6
    name = num2str(i);
    saveas(figure(i),append('./pics/',name),'png')
end
end