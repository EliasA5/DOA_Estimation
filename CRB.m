function [crb] = CRB(type,v_0,alpha,theta,Rv,M,X,a_f,da_f,P)
%This function calculates the CRB
% depending on the type of the bound ('cyclic 1', 'cyclic 2','regular')
% v_0 is the speed of the wave
% alpha is the incidence angle
% theta is the doa angle
% Rv is the correlation matrix, note that we take the pseudo inverse of it,
% since our correlation matrix is non-singular
% M the number of frequencies in the signal
% X is the signal in frequency, which has P segments and each segment
% contains M samples
% a_f is the first function returned from model.m (the steering vector)
% da_f is the second function returned from model.m (derivative of the
% steering vector)
% P the number of segments in the signal, unused

%returns the value of the CRB bound
Rv_inv = pinv(Rv);
X_w = squeeze(sum(X,1));
A = cell2mat(arrayfun(a_f,1 : M,theta*ones(1,M), alpha*ones(1,M), v_0*ones(1,M),'uniformoutput',false));
dA = cell2mat(arrayfun(da_f,1 : M,theta*ones(1,M), alpha*ones(1,M), v_0*ones(1,M),'uniformoutput',false));
    
Up = diag(A' * Rv_inv * X_w);
Down = diag(A' * Rv_inv * A);
s = (Up ./ Down).';
dMue = dA .* s;
FIM = sum(2 * real(diag(dMue' * Rv_inv * dMue)));

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


