function [crb] = CRB(type,v_0,alpha,theta,Rv,M,X,a_f,da_f,P)
%This function calculates the CRB
% depending on the type of the bound ('cyclic 1', 'cyclic 2','regular')
Rv_inv = pinv(Rv);
X_w = squeeze(sum(X,1));
s = zeros(1,M);
FIM = 0;
for j = 1 : M
    x_i = X_w(:,j);
    % calculate the a and a derivative 
    a = a_f(j, theta, alpha, v_0);
    da = da_f(j, theta, alpha, v_0);

    % estimating s(m) using the MLE of the deterministic signal
    down = a' * Rv_inv * a;   
    s_j = a' * Rv_inv * x_i / P;   % per the formula only x is a variable of p
    s(j) = s_j / down;

    dmue = da * s(j);
    
    FIM = FIM + real(dmue' * Rv_inv * dmue);
end

FIM = FIM * 2;

switch(type)
    
    case 'cyclic 1'
        crb =  FIM ^ (-1);
        crb = 2 - 2/(sqrt(crb + 1));

    case 'cyclic 2'
        crb =  FIM ^ (-1);
        crb = 2 + 2 * (1 - sqrt(1 + 2 * crb))/crb;

    case 'regular'
        crb = FIM ^ (-1);

end


