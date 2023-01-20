function [alpha,beta,success_total,failure_total] = update_strategy_posterior_probability(trial_type,decay_rate,success_total,failure_total,alpha_zero,beta_zero,varargin)

% UPDATE_STRATEGY_POSTERIOR_PROBABILITY updates the Beta distribution parameters
% [ALPHA,BETA,S_TOTAL,F_TOTAL] = UPDATE_STRATEGY_POSTERIOR_PROBABILITY(TYPE,DECAY_RATE,S_TOTAL,F_TOTAL,ALPHA_ZERO,BETA_ZERO)
% updates the Beta distribution parameters (ALPHA, BETA), according to:
% TYPE: a string giving the trial type of the current trial: a "failure",
% "success", or "null" attempt to execute the strategy, as returned by the
% strategy model functions.
% DECAY_RATE: scalar, the decay weight of previous trials' events (parameter "gamma")
% S_TOTAL, F_TOTAL: scalars, respectively the running total of previous successes
% and failures, decayed. Set S_TOTAL = F_TOTAL = 0 for the first trial.
% ALPHA_ZERO, BETA_ZERO: scalars, the priors for alpha and beta  
%
% If trial is NULL, then ALPHA=BETA=NaN by default
% 
% ... = UPDATE_STRATEGY_POSTERIOR_PROBABILITY(...,'DecayNull') will decay
% the success and failure totals during null trials too, and return the
% updated ALPHA and BETA.
%
% Mark Humphries 2/4/2022

blnDecayNull = 0;   % default is to not decay during Null trials

if nargin >= 7 && strcmp(varargin{1},'DecayNull') == 1
    blnDecayNull = 1;
end

switch trial_type
    case "success"
        % update success total, and decay both prior sets of evidence
        success_total = decay_rate * success_total + 1;
        failure_total = decay_rate * failure_total;
    case "failure"
        success_total = decay_rate * success_total;
        failure_total = decay_rate * failure_total + 1;   
    case "null"
        if blnDecayNull
            success_total = decay_rate * success_total;
            failure_total = decay_rate * failure_total;
        end
end

% update Beta distribution parameter estimates
if strcmp(trial_type,"null") && ~blnDecayNull
    % return NaNs for null trials by default
    alpha = NaN;
    beta = NaN;
else
    alpha = alpha_zero + success_total; 
    beta = beta_zero + failure_total;
end