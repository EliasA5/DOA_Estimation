close all
clear all
clc
%--------------------------------------------------------------------------
K_1 = 15;
K_3 = 1;
alpha = pi/4;
v_0 = 9600;  % P waves velocities can be from 6 Km/sec to 11 Km/sec
SNR = 10;
sigma_source = 10;
sigma_noise = 1;  % for white and colored noise
data = load("data.mat");
r_k = data.r_m;   % r_m is loaded form the data matrix we extrcted from the getData python script
M = 100;   % for these simlutaions we decided to use one segement
P = 1;
f = 3;      % the omega_m is constant at this point for all the frquencies
%--------------------------------------------------------------------------
Tests = 10;

theta_og = -pi +2*pi/Tests : 2*pi/Tests : pi;
theta_0 = pi/4; % starting estimate at 45 deg

iters = 500;

theta_white_f = zeros(Tests,1);
theta_colored_f = zeros(Tests,1);

step_size_white = 1;
step_size_colored = 1;
gamma = 0.9;

for t = 1 : Tests
    
    [X_colored,s_colored,Rv_colored,~] = synData(r_k, theta_og(t), alpha, v_0, sigma_source, sigma_noise, M, 'colored', f, K_1, K_3, P);
    [X_white,s_white,Rv_white,~] = synData(r_k, theta_og(t), alpha, v_0, sigma_source, sigma_noise, M, 'white', f, K_1, K_3, P);

    theta_white = Fisher_scoring(theta_0,Rv_white,f,v_0,alpha,K_3,K_1,X_white,iters,r_k,step_size_white,gamma,M,s_white);
    theta_colored = Fisher_scoring(theta_0,Rv_colored,f,v_0,alpha,K_3,K_1,X_colored,iters,r_k,step_size_colored,gamma,M,s_colored);

    theta_white_f(t) = (theta_white(end));
    theta_colored_f(t) = (theta_colored(end));

    MSP_white(t) = MSPE(theta_white_f(t), theta_og(t));
    MSP_colored(t) = MSPE(theta_colored_f(t), theta_og(t));

%     figure;
%     hold on
%     plot(theta_colored)
%     plot(theta_white)
%     hold off
%     legend('colored' , 'white')

end

bias_white = wrapToPi((theta_white_f - theta_og.'));
bias_colored = wrapToPi((theta_colored_f - theta_og.'));

% MSP_white = wrapToPi((theta_white_f - theta_og.') .^ 2);
% MSP_colored = wrapToPi((theta_colored_f - theta_og.') .^ 2);

% MSP_white = MSPE(theta_white_f, theta_og.');
% MSP_colored = MSPE(theta_colored_f, theta_og.');

figure;
hold on;
plot(theta_og,'*')
plot(theta_white_f,'o')
stem(theta_colored_f)
legend('Original', 'White Estimate', 'Colored Estimate')
xlabel('Test')
ylabel('\theta')
hold off

figure;
hold on;
stem(bias_white)
stem(bias_colored)
legend('White Estimate','Colored Estimate')
xlabel('Test')
ylabel('Bias')
hold off

figure;
hold on;
stem(MSP_white)
stem(MSP_colored)
legend('White Estimate','Colored Estimate')
xlabel('Test')
ylabel('MSE')
hold off

%%
angles = linspace(0, 2*pi, 500);
radius = 20;
CenterX = 50;
CenterY = 40;
x = radius * cos(angles) + CenterX;
y = radius * sin(angles) + CenterY;

x_og = radius * cos(theta_og) + CenterX;
x_white = radius * cos(theta_white_f) + CenterX;
x_colored = radius * cos(theta_colored_f) + CenterX;

y_og = radius * sin(theta_og) + CenterY;
y_white = radius * sin(theta_white_f) + CenterY;
y_colored = radius * sin(theta_colored_f) + CenterY;

figure
hold on
grid on
axis equal;
plot(x, y, 'k-', 'LineWidth', 1);
plot(x_og, y_og, 'rx','LineWidth', 2)
plot(x_white, y_white, 'bo','LineWidth', 1.5)
plot(x_colored, y_colored, 'go','LineWidth', 1.5)
xlabel('X', 'FontSize', 10);
ylabel('Y', 'FontSize', 10);
legend('','Original','White','Colored')



