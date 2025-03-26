function w = akaikeWeights(AIC,scaling)
% w = akaikeWeights(AIC,scaling)
%
% Given an array of AIC values, computes Akaike weights for each AIC.
%
% If scaling is set to 0, then output is just L(models|data) without the
% weighted scaling s.t. they sum to 1.

if ~exist('scaling','var') || isempty(scaling), scaling = 1; end

minAIC = min(AIC);

if scaling
    w = exp(-.5*(AIC-minAIC)) / sum(exp(-.5*(AIC-minAIC)));
else
    w = exp(-.5*(AIC-minAIC));
end