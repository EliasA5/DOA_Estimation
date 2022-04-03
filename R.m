
files = dir('./Geres/matFiles/signal_*.mat');
r = zeros(25,25);
i = 0;
for file = files'
load(fullfile(file.folder, file.name));
r = 1/(i+1) * (i*r + R_fromData(data));
i = i+1;
end
