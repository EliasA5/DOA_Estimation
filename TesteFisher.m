close all
clear all
clc
%--------------------------------------------------------------------------
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
iters = 5000;
step = 1;
gamma = 0.8;

theta_colored = zeros(Tests,iters);
theta_white= zeros(Tests,iters);


theta_og = -pi +2*pi/Tests : 2*pi/Tests : pi;
%theta_og = zeros(Tests,1);
%theta_og = rand*2*pi - pi;
epsilon = 0.1;

for i = 1:Tests
    
    theta_0 = 0;
    
    
    [X_colored,Rv_colored,~] = synData(rm, theta_og(i), alpha, v_0, sigma_source, sigma_noise, M, 'colored', f, K_1, K_3, P);
    [X_white,Rv_white,~] = synData(rm, theta_og(i), alpha, v_0, sigma_source, sigma_noise, M, 'white', f, K_1, K_3, P);

%     Rv_colored = EstimateRv(X_colored,3*K_3+K_1,M);
%     Rv_white = EstimateRv(X_white,3*K_3+K_1,M);


    
    theta_colored(i,:) = wrapToPi(real(Fisher_scoring(theta_0,Rv_colored,f,v_0,alpha,K_3,K_1,X_colored,iters,rm,step,gamma)));
    theta_white(i,:) = wrapToPi(real(Fisher_scoring(theta_0,Rv_white,f,v_0,alpha,K_3,K_1,X_white,iters,rm,step,gamma)));

    figure
    hold on
    plot(theta_white(i,:))
    plot(theta_colored(i,:))
    legend('white','colored')
    hold off
end

figure;
hold on;
plot(theta_og,'*')
plot(theta_white(:,iters),'o')
stem(theta_colored(:,iters))
legend('Original','White Estimate','Colored Estimate')
xlabel('Test')
ylabel('\theta')
hold off


%%
function [Rv] = EstimateRv(x,K,N) 

  f = @(k) k.' * k;
  g = @(k) k * k.';
  sigma_squared_array = zeros(1,N);
  Rv_array = zeros(K,K,N);

  for i = 1:N
    k = x(:,i);
    sigma_squared_array(i) = f(k);
    Rv_array(:,:,i) = g(k);

  end

  %sigma_squared = 1/(K*N) * sum(sigma_squared_array);
  Rv = 1/N * sum(Rv_array,3);

end
