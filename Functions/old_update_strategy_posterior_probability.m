function [alpha_now,beta_now] = update_strategy_posterior_probability(trial_type,decay_rate,alpha_previous,beta_previous)

% UPDATE_STRATEGY_POSTERIOR_PROBABILITY updates the Beta distribution parameters
% [ALPHA_NEW,BETA_NEW] = UPDATE_STRATEGY_POSTERIOR_PROBABILITY(TYPE,DECAY_RATE,ALPHA,BETA)
% updates the prior estimates of the Beta distribution parameters (ALPHA,
% BETA), according to the trial type of the current trial: a failure,
% success, or null attempt to execute the strategy. Prior evidence is
% decayed at DECAY_RATE (parameter "gamma")
%
% Mark Humphries 31/3/2022

switch trial_type
    case "success"
        alpha_now = decay_rate * alpha_previous + 1;
        beta_now = decay_rate * beta_previous;
    case "failure"
        alpha_now = decay_rate * alpha_previous;
        beta_now = decay_rate * beta_previous + 1;        
    case "null"
        alpha_now = decay_rate * alpha_previous;
        beta_now = decay_rate * beta_previous;           
end