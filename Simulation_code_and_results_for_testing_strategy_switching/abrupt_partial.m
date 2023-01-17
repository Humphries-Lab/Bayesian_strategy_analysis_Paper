function [next_x_current,next_x_new] = abrupt_partial(trial_number)

% abrupt partial: immediately executes the new strategy on every trial 
% after trial t, but in N=2 choice experiment, choice is consistent with 
% previous strategy on half of trials
% Modelled as an alternating sequence of 0s, 1s

next_x_new = 1;

% if an even-numbered trial, then choice is consistent with both strategies
if rem(trial_number,2) == 0
    next_x_current = 1;
else
    next_x_current = 0;
end