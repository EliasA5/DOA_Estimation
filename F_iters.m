close all
clear all
clc
%--------------------------------------------------------------------------
% This scripts checks the effect of the number of iterations of the Fisher's scoring ethod
%%
K_1 = 15;
K_3 = 1;
alpha = pi/4;
v_0 = 1;  % P waves velocities can be from 6 Km/sec to 11 Km/sec
SNR = 10;
sigma_source = 10;
sigma_noise = 1;  % for white and colored noise
data = load("data.mat");
rm = data.r_m;   % r_m is loaded form the data matrix we extrcted from the getData python script
M = 100;   % for these simlutaions we decided to use one segement
P = 1;
f = 3;      % the omega_m is constant at this point for all the frquencies

%%
Tests = 10;
iters = [100 , 1000 , 2000 , 3000 , 4000];    % the number of iterations for the Fisher's scoring (second method)

RMSPE_colored_100 = zeros(Tests,1);
RMSPE_white_100 = zeros(Tests,1);

RMSPE_colored_1000 = zeros(Tests,1);
RMSPE_white_1000 = zeros(Tests,1);

RMSPE_colored_2000 = zeros(Tests,1);
RMSPE_white_2000 = zeros(Tests,1);

RMSPE_colored_3000 = zeros(Tests,1);
RMSPE_white_3000 = zeros(Tests,1);

RMSPE_colored_4000 = zeros(Tests,1);
RMSPE_white_4000 = zeros(Tests,1);

theta_og = -pi +2*pi/Tests : 2*pi/Tests : pi;
%theta_og = zeros(Tests,1);

for i = 1:Tests

    theta = theta_og(i); % theta is in [-pi,pi]

    [X_colored,~,Rv_colored,~] = synData(rm, theta, alpha, v_0, sigma_source, sigma_noise, M, 'colored', f, K_1, K_3, P);
    [X_white,s,Rv_white,a] = synData(rm, theta, alpha, v_0, sigma_source, sigma_noise, M, 'white', f, K_1, K_3, P);

    % Fisher's scoring for colored noise
    theta_estimate = 0;
    
    for j = 1:4000
        da =  A_derivative(v_0,alpha,rm,K_3,K_1,theta,f);
        df_theta = PDF_derivative(K_3,K_1,v_0,alpha,f,s,rm,Rv_colored,X_colored,theta,da,M);
        F = FIM(s,M,Rv_colored,da); %scalar for only theta
        theta_estimate = theta_estimate + df_theta/F;

        if (j==100)
            theta_colored = real(theta_estimate);
%             if abs(theta_colored) > pi
%                   theta_colored = mod(theta_colored,pi);
%             end
            RMSPE_colored_100(i) =  sqrt(mean((theta-theta_colored).^2));
        end

        if (j==1000)
            theta_colored = real(theta_estimate);
%             if abs(theta_colored) > pi
%                   theta_colored = mod(theta_colored,pi);
%             end
            RMSPE_colored_1000(i) =  sqrt(mean((theta-theta_colored).^2));
        end

        if (j==2000)
            theta_colored = real(theta_estimate);
%             if abs(theta_colored) > pi
%                   theta_colored = mod(theta_colored,pi);
%             end
            RMSPE_colored_2000(i) =  sqrt(mean((theta-theta_colored).^2));
        end

        if (j==3000)
            theta_colored = real(theta_estimate);
%             if abs(theta_colored) > pi
%                   theta_colored = mod(theta_colored,pi);
%             end
            RMSPE_colored_3000(i) =  sqrt(mean((theta-theta_colored).^2));
        end
    end
    theta_colored = real(theta_estimate);
%     if abs(theta_colored) > pi
%           theta_colored = mod(theta_colored,pi);
%     end
    RMSPE_colored_4000(i) =  sqrt(mean((theta-theta_colored).^2));
    

    % Fisher's scoring for white noise
    theta_estimate = 0;
    
    for j = 1:4000
        da =  A_derivative(v_0,alpha,rm,K_3,K_1,theta,f);
        df_theta = PDF_derivative(K_3,K_1,v_0,alpha,f,s,rm,Rv_white,X_white,theta,da,M);
        F = FIM(s,M,Rv_white,da); %scalar for only theta
        theta_estimate = theta_estimate + df_theta/F;

        if (j==100)
            theta_white = real(theta_estimate);
%             if abs(theta_white) > pi
%                   theta_white = mod(theta_white,pi);
%             end
            RMSPE_white_100(i) =  sqrt(mean((theta-theta_white).^2));
        end

        if (j==1000)
            theta_white = real(theta_estimate);
%             if abs(theta_white) > pi
%                   theta_white = mod(theta_white,pi);
%             end
            RMSPE_white_1000(i) =  sqrt(mean((theta-theta_white).^2));
        end

        if (j==2000)
            theta_white = real(theta_estimate);
%             if abs(theta_white) > pi
%                   theta_white = mod(theta_white,pi);
%             end
            RMSPE_white_2000(i) =  sqrt(mean((theta-theta_white).^2));
        end

        if (j==3000)
            theta_white = real(theta_estimate);
%             if abs(theta_white) > pi
%                   theta_white = mod(theta_white,pi);
%             end
            RMSPE_white_3000(i) =  sqrt(mean((theta-theta_white).^2));
        end

    end
    theta_white = real(theta_estimate);
%     if abs(theta_white) > pi
%           theta_white = mod(theta_white,pi);
%     end
    RMSPE_white_4000(i) =  sqrt(mean((theta-theta_white).^2));

end

figure;
hold on;
stem(theta_og,RMSPE_colored_100)
stem(theta_og,RMSPE_colored_1000)
stem(theta_og,RMSPE_colored_2000)
stem(theta_og,RMSPE_colored_3000)
stem(theta_og,RMSPE_colored_4000)
hold off;
grid on;
legend('100' , '1000' , '2000' , '3000' , '4000');
title('Colored')

figure;
hold on;
stem(theta_og,RMSPE_white_100)
stem(theta_og,RMSPE_white_1000)
stem(theta_og,RMSPE_white_2000)
stem(theta_og,RMSPE_white_3000)
stem(theta_og,RMSPE_white_4000)
hold off;
grid on;
legend('100' , '1000' , '2000' , '3000' , '4000');
title('white')




%--------------------------------------------------------------------------
