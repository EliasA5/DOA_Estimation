close all
clear all
clc
%-------------------------------------------------------
data = load("data.mat","data");
x = data.data;
[K,N] = size(x); %k is number of sensors
c_white = 0.978; %the threshholds for the ML estimator for estimating if the noise is colored
c_colored = 0.4;
hcn = dsp.ColoredNoise('InverseFrequencyPower',1,'SamplesPerFrame',N,'NumChannels',K); %colored noise

%-------------------------------------------------------
 psi_whiteN = zeros(1,150); % Test for white noise
 psi_coloredN = zeros(1,150); % Test for colored noise
true_pos = zeros(1,20);
false_pos = zeros(1,20);
caclc_Rv = 0;


for j =1:20
    for i=1:150
        x_ = randn(size(x));              %x is a white noise
        [Rv,psi_white] = NoiseTest1(x_,K,N);
         psi_whiteN(i) = psi_white;
        if (psi_white < c_white)
            false_pos(j) = false_pos(j) + 1;
        end
        x_ = hcn().';                     %x is a colored noise
        [Rv,psi_color] = NoiseTest1(x_,K,N);
        psi_coloredN(i) = psi_color;
        if (psi_color < c_colored)
            true_pos(j) = true_pos(j) + 1;
        end
    end
  
end
% 
true_pos = true_pos/150;
false_pos = false_pos/150;


%Building the ROC
figure;
hold on;
stem(false_pos , true_pos,'X');
xlabel("False Positive");
ylabel("True Positive");
title("ROC")
xlim([0 1])
ylim([0 1])
plot([0:0.1:1],[0:0.1:1]) %comparing our model to the random classifier
hold off;

function [Rv,psi_] = NoiseTest1(x,K,N) 
 % X = K * N
    Rv = zeros(K,K);
    for i =1:N
       Rv = Rv + x(:,i)*x(:,i).';
    end
    Rv = Rv/N;
    psi_ = (det(Rv)/((trace(Rv)/K)^K));
end

