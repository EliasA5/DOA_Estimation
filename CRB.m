function [crb] = CRB(type,v_0,alpha,theta,Rv,M,X,a_f,da_f,P)
%This function calculates the CRB
% depending on the type of the bound ('cyclic 1', 'cyclic 2','regular')
Rv_inv = pinv(Rv);      
X_w = squeeze(sum(X,1));
A = cell2mat(arrayfun(a_f,1 : M,theta*ones(1,M), alpha*ones(1,M), v_0*ones(1,M),'uniformoutput',false));
dA = cell2mat(arrayfun(da_f,1 : M,theta*ones(1,M), alpha*ones(1,M), v_0*ones(1,M),'uniformoutput',false));
Up = diag(A' * Rv_inv * X_w);
Down = diag(A' * Rv_inv * A);
s = (Up ./ Down).' / P;
dMue = dA .* s;
FIM = sum(2 * real(diag(dMue' * Rv_inv * dMue)));

switch(type)
    
    case 'cyclic 1'
        crb =  FIM ^ (-1);
        crb = 2 - 2/(sqrt(crb + 1));

    case 'cyclic 2'
        crb =  FIM ^ (-1);
        crb = 2  + 2 / crb - (2 * sqrt(1 + 2 * crb)) / crb;

    case 'regular'
        crb = FIM ^ (-1);

end


