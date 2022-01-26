function [F] = FIM(s,M,R_v,da)
%This function calculates the  FIM 
R_v_inv = pinv(R_v);
F = 0;
for i = 1:M
    dm = da*s(i);
    F = F + real(dm' * R_v_inv * dm);
end
F = F * 2;
end