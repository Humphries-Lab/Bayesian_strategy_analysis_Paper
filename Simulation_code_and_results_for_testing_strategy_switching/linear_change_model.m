function [p1,p2] = linear_change_model(t_since_switch,probability_pair,slope)

% compute linear change in probability of strategies

if t_since_switch < (probability_pair(1)-probability_pair(2))/slope
    % if not yet reached the target probability, then
    % update
    p1 = probability_pair(1) - slope*t_since_switch; % decrease larger probability
    p2 = probability_pair(2) + slope*t_since_switch;  % increase larger probability
else
    % assign target probabilities: switched from original
    p1 = probability_pair(2);
    p2 = probability_pair(1);
end
