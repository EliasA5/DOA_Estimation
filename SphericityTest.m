close all
clear all
clc

%for every time step (1,2,3,...,9600) we calculate the estimation of sigma and R by using
%all the 19 sensors 

%system('conda activate obspy & python getData.py'); %uncomment to run the python data getter
data = load("data.mat","data");
x = data.data;
[K,N] = size(x); %k is number of sensors
c = 10^-10;
hcn = dsp.ColoredNoise('InverseFrequencyPower',1,'SamplesPerFrame',N,'NumChannels',K);
%x = hcn().';
psi_whiteN = zeros(1,150);
psi_coloredN = zeros(1,150);
for i=1:150
    x_ = randn(size(x));
    [Rv,psi_] = NoiseTest(x_,K,N,10^(-10));
    psi_whiteN(i) = psi_;
    x_ = hcn().';
    [Rv,psi_] = NoiseTest(x_,K,N,10^(-10));
    psi_coloredN(i) = psi_;
end

[Rv,psi_] = NoiseTest(x,K,N,10^(-10));
figure;
heatmap(db(Rv), 'Colormap', bone);

%oldpath = addpath('C:\Users\Elias\Desktop\project\LIBRA', '-end');
mcdRv = mcdcov(x.','cor', 1, 'plots', 0);
mcdCov = mcdRv.cov;
figure;
heatmap(db(mcdCov), 'Colormap', bone);


function [Rv,psi_] = NoiseTest(x,K,N,c) 

  f = @(k) k.' * k;
  g = @(k) k * k.';
  sigma_squared_array = zeros(1,N);
  Rv_array = zeros(K,K,N);

  for i = 1:N
    k = x(:,i);
    sigma_squared_array(i) = f(k);
    Rv_array(:,:,i) = g(k);

  end

  sigma_squared = 1/(K*N) * sum(sigma_squared_array);
  Rv = 1/N * sum(Rv_array,3);

  psi_ = ( det(Rv) / ((trace(Rv)/K)^K) )^(N/2);
  %figure;
  %imagesc(db(co));
end

%useful code
%heatmap(db(Rv), 'Colormap', bone); %, 'Colorlimits', [-350, -270]
%hcn = dsp.ColoredNoise('InverseFrequencyPower',1,'SamplesPerFrame',9600,'NumChannels',19);



