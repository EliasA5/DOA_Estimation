close all
clear all

%for every time step (1,2,3,...,9600) we calculate the estimation of sigma and R by using
%all the 17 sensors 

K = 17; %k is number of sensors
N1 = 1; % start measurement, end measurement
N2 = 9600;
N = N2-N1+1 ; % number of past measurements
data = load("data.mat","data");
x = data.data;

f = @(k) k.' * k;
g = @(k) k * k.';
sigma_squared = 0;
Rv = 0;

for i = N1:N2
  k = x(:,i);

  sigma_squared = sigma_squared + f(k);
  Rv = Rv + g(k);

end

sigma_squared = 1/(K*N) * sigma_squared;
Rv = 1/N * Rv;

psi_ = ( det(Rv) / ((trace(Rv)/K)^K) );
test = -N*log(psi_);

c = 10^(-10)

A = sigma_squared/test




