
files = dir('./matFiles/signal_*.mat');

function [Rv] = R_fromDataa(x)
x = double(x);
col = x(:,1);
Rv =  col * col.';
for i =2:length(x)
    col = x(:,1);
    Rv = 1/(i+1) * (i*Rv + col * col.');
end

end