function [Xt Abar sigma2_e] = fgn_pred(Xts,H,sigma2)
% [Xt Abar sigma2_e] = fgn_pred(Xts,H,sigma2)
% 
% Estimate Xt from Xt-1, ..., Xt-k for an fGn series with hurst exponent H 
% and variance sigma2, where input Xts is a vector of form [Xt-1, Xt-2, ..., Xt-k]
%
% Xt = sum(Abar .* Xts)


k = length(Xts);

n = 0:k ;
r = sigma2/2*(abs(n-1).^(2*H) + abs(n+1).^(2*H) - 2*abs(n).^(2*H)) ;
Rbar = r(2:end)' ; 
temp = r(1:end-1) ; 
Rbarbar = toeplitz(temp') ; 
Abar = inv(Rbarbar)*Rbar  ;

sigma2_e = r(1) - sum(Abar'.*Rbar');

Abar = Abar';

Xt = sum(Abar .* Xts);