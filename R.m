%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Calculate the covariance matrix from real data saved in
% signal_*.mat files.
% Problems: the dataset doesn't always contain all the sensors, so for
% signals with less than 25 sensors this script will not run correctly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
files = dir('./Geres/matFiles/signal_*.mat');
r = zeros(25,25);
i = 0;
for file = files'
load(fullfile(file.folder, file.name));
r = 1/(i+1) * (i*r + R_fromData(data));
i = i+1;
end
