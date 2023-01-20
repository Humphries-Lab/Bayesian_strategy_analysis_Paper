function [alpha,beta] = set_Beta_prior(strPrior)

% SET_BETA_PRIOR set requested prior parameters for Beta distribution
% [A,B] = SET_BETA_PRIOR(P) sets the Beta distribution to the specified prior 
%   in P: 'Uniform'; 'Jeffreys'
%
% Initial version: 22/1/2019
% Mark Humphries

switch strPrior
    case 'Uniform'
        alpha = 1; beta = 1;
    case 'Jeffreys'
        alpha = 0.5; beta = 0.5;
    otherwise
        error('Unknown prior distribution for the Beta distribution')
end