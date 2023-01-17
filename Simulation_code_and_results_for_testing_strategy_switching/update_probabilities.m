function [p_current, p_new] = update_probabilities(trial_number,gamma,x_1,x_2,alpha0,beta0)

% update the probabilities of each strategy for this trial
% given: (trial number, gammma, x_1, x_2),
% where x_1, x_2 are the binary sequences of success (1) or failure (2) to
% execute strategies 1 or 2 up to trial_number
%
% 3/10/2022: original version
% 7/11/2022: updated to use nansum to handle nans when trial type is null
% Mark Humphries

trial_sequence = [1:trial_number]';
 
% compute s* and f* (Eqs 8 and 9 in the paper)
s_star_current = nansum(gamma.^(trial_number - trial_sequence) .* x_1);
s_star_new = nansum(gamma.^(trial_number - trial_sequence) .* x_2);
f_star_current = nansum(gamma.^(trial_number - trial_sequence) .* (1-x_1));
f_star_new = nansum(gamma.^(trial_number - trial_sequence) .* (1 - x_2));

% update alpha, beta (text after Eq 10 in the paper)
alpha_current = alpha0 + s_star_current;
alpha_new = alpha0 + s_star_new;
beta_current = beta0 + f_star_current;
beta_new = beta0 + f_star_new;

% compute MAP probabilities
p_current = (alpha_current - 1) / (alpha_current + beta_current - 2);
p_new = (alpha_new - 1) / (alpha_new + beta_new - 2);