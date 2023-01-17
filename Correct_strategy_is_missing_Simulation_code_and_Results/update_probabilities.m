function [MAP, variance] = update_probabilities(trial_number,gamma,x,alpha0,beta0)

% update the MAP probability and variance for a strategy for this trial
% given: (trial number, gammma, x, alpha0, beta0),
% where x is the binary sequences of success (1) or failure (0) to
% execute the strategy up to trial_number

trial_sequence = [1:trial_number]';
 
% compute s* and f* (Eqs 8 and 9 in the paper)
s_star = sum(gamma.^(trial_number - trial_sequence) .* x);
f_star = sum(gamma.^(trial_number - trial_sequence) .* (1-x));

% update alpha, beta (text after Eq 10 in the paper)
alpha = alpha0 + s_star;
beta = beta0 + f_star;


% compute MAP probabilities & precision
MAP = (alpha - 1) / (alpha + beta - 2);
variance = (alpha*beta) / ((alpha+beta).^2 .* (alpha+beta+1));
% precision = 1/variance;

