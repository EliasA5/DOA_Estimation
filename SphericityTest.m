close all
clear all
clc

%for every time step (1,2,3,...,9600) we calculate the estimation of sigma and R by using
%all the 19 sensors 
for i=1:1
    %system('conda activate obspy & python getData.py'); %uncomment to run the python data getter
    data = load("6.mat","data");
    x = data.data;
    %x(6,:) = [];
    [K,N] = size(x); %k is number of sensors
    c_white = 0.978; %the threshholds for the ML estimator for estimating if the noise is colored
    c_colored = 0.4;
    [Rv,psi_color] = NoiseTest1(x,K,N);
    if (psi_color < c_colored)
        test_noise = 1;
    elseif (psi_color > c_white)
        tes_noise = -1;
    end

end


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

  psi_ = ( det(Rv) / ((trace(Rv)/K)^K) );
  %figure;
  %imagesc(db(co));
end

function [Rv,psi_] = NoiseTest1(x,K,N) 
 % X = K * N
    Rv = zeros(K,K);
    for i =1:N
       Rv = Rv + x(:,i)*x(:,i).';
    end
    Rv = Rv/N;
    psi_ = (det(Rv)/((trace(Rv)/K)^K));
end

%useful code
%heatmap(db(Rv), 'Colormap', bone); %, 'Colorlimits', [-350, -270]
%hcn = dsp.ColoredNoise('InverseFrequencyPower',1,'SamplesPerFrame',9600,'NumChannels',19);



