close all
clear all
clc
%-------------------------------------------------------
data = load("data.mat","data");
x = data.data;
[K,N] = size(x); %k is number of sensors
c = 0.5;
hcn = dsp.ColoredNoise('InverseFrequencyPower',1,'SamplesPerFrame',N,'NumChannels',K); %colored noise

%-------------------------------------------------------
 psi_whiteN = zeros(1,150); % Test for white noise
 psi_coloredN = zeros(1,150); % Test for colored noise
false_neg = zeros(1,10);
false_pos = zeros(1,10);
caclc_Rv = 0;

% for i  = 1:10
%     for j = 1:150
%         x_white = randn(size(x));
%         x_color = hcn().';
%     end 
% end
% 
% 

for j =1:10
    for i=1:150
        x_ = randn(size(x));              %x is a white noise
        [Rv,psi_white] = NoiseTest1(x_,K,N);
         psi_whiteN(i) = psi_white;
        if (psi_white < c)
            false_pos(j) = false_pos(j) + 1;
        end
        x_ = hcn().';                     %x is a colored noise
        [Rv,psi_color] = NoiseTest1(x_,K,N);
        psi_coloredN(i) = psi_color;
        if (psi_color > c)
            false_neg(j) = false_neg(j) + 1;
        end
    end
end
% 
% false_neg = false_neg/150;
% false_pos = false_pos/150;
% 
% %Building the ROC
% figure;
% plot(false_pos , (1-false_neg));
% xlabel("False Positive");
% ylabel("True Positive");
% title("ROC")


% x_white = randn(size(x));
% [Rv,psi_white] = NoiseTest1(x_white,K,N);
% figure;
% heatmap(Rv);
% 
% x_color = hcn().';
% [Rv,psi_color] = NoiseTest1(x_color,K,N);
% figure;
% heatmap(Rv);
% 
% [Rv_data,psi_data] = NoiseTest1(x,K,N);
% figure;
% heatmap(Rv_data);

function [Rv,psi_] = NoiseTest1(x,K,N) 
 % X = K * N
    Rv = zeros(K,K);
    for i =1:N
       Rv = Rv + x(:,i)*x(:,i).';
    end
    Rv = Rv/N;
    psi_ = (det(Rv)/((trace(Rv)/K)^K));
    psi_ = -db(psi_);
end

