close all
clear all


K = 2; %k is number of sensors
N1 = 1; % start measurement, end measurement
N2 = 9600;
N = N2-N1+1 ; % number of past measurements

for i = 1:10
  x = randn(K,N);

  f = @(k) k.' * k;
  sigma_squared = 0;
  for k = x
    sigma_squared = sigma_squared + f(k);
   end
  sigma_squared = 1/(K*N) * sigma_squared;

  g = @(k) k * k.';
  Rv = 0;
  for k = x
    Rv = Rv + g(k);
  end
  Rv = 1/N * Rv;

  psi_ = ( det(Rv) / (trace(Rv)/K)^K ) ^ (-N/2);
  test = 2*log(psi_)
end



