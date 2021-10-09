close all
clear all
clc

%for every time step (1,2,3,...,9600) we calculate the estimation of sigma and R by using
%all the 17 sensors 


N1 = 1; % start measurement, end measurement
N2 = 9600;
N = N2-N1+1 ; % number of past measurements
data = load("data.mat","data");
x = data.data;
[K,N] = size(x); %k is number of sensors

f = @(k) k.' * k;
g = @(k) k * k.';
sigma_squared_array = zeros(1,N);
Rv_array = zeros(K,K,N);

for i = N1:N2
  k = x(:,i);
  sigma_squared_array(i) = f(k);
  Rv_array(:,:,i) = g(k);

end

sigma_squared = 1/(K*N) * sum(sigma_squared_array);
Rv = 1/N * sum(Rv_array,3);

psi_ = ( det(Rv) / ((trace(Rv)/K)^K) );
test = -N*log(psi_);

c = 10^(-10);

A = sigma_squared/test




