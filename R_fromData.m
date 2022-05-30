%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: this fucntion calculates the covariance matrix using MSE
% estimator from the data supplied by x, preferably x has large amounts of
% data for better estimation of the matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Rv] = R_fromData(x)
x = double(x);
col = x(:,1);
Rv =  col * col.';
for i =2:length(x)
    col = x(:,1);
    Rv = 1/(i+1) * (i*Rv + col * col.');
end

end