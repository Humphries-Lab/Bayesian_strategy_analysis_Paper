function [alpha_per_trial,beta_per_trial,MAP_per_trial,precision_per_trial] = interpolate_null_trials(alpha_timeseries,beta_timeseries,alpha_zero,beta_zero)

% INTERPOLATE_NULL_TRIALS replaces NaN trials with last evidence update
% [ALPHA_PER_TRIAL,BETA_PER_TRIAL,MAP_PER_TRIAL,PRECISION_PER_TRIAL] = 
%       INTERPOLATE_NULL_TRIALS(ALPHA_TIMESERIES,BETA_TIMESERIES,ALPHA_ZERO,BETA_ZERO)
% give:
%   ALPHA_TIMESERIES: N-element vector, the time-series of the alpha parameter per 
%   trial that contains NaNs for Null trials (i.e. those where the 
%   strategy's probability could not be updated)
%   BETA_TIMESERIES: N-element vector, the time-series of the beta parameter per 
%   trial that contains NaNs for Null trials 
%   ALPHA_ZERO, BETA_ZERO: scalars, the priors for alpha and beta 
%
% Returns:
%   ALPHA_PER_TRIAL,BETA_PER_TRIAL = N-element vectors, which are interpolated 
%   by replacing NaNs in the original time-series with the values from the 
%   last evidence update. For the first set of trials with no update, the prior 
%   probabilities are used.
%   MAP_PER_TRIAL,PRECISION_PER_TRIAL: N-element vectors of MAP probability
%   and precision per trial, derived from the interpolated (alpha,beta) values
%
% Note:
%   this is constant interpolation, assigning a constant to all trials with
%   NaNs; alternative approaches could consider linear interpolation etc
% 
% 4/5/2022: initial version
% Mark Humphries

%% set up storage
N = numel(alpha_timeseries);

if N ~= numel(beta_timeseries)
    error('alpha and beta time-series are different lengths');
end

alpha_per_trial = zeros(N,1);
beta_per_trial = zeros(N,1);
MAP_per_trial = zeros(N,1);
precision_per_trial = zeros(N,1);

%% interpolate alpha,beta
alpha_entries = find(~isnan(alpha_timeseries)); % indices of non-NaN entries

alpha_per_trial(1:alpha_entries(1)-1) = alpha_zero;  % trials before first update are set to the prior 
for index_entry = 1:numel(alpha_entries)-1
    % iterate through queue of non-NaN entries to update
    % for all indices between this non-NaN entry and the next one, update
    % to the current non-NaN value
    alpha_per_trial(alpha_entries(index_entry):alpha_entries(index_entry+1)-1) ... 
        = alpha_timeseries(alpha_entries(index_entry)); 
end
alpha_per_trial(alpha_entries(end):end) = alpha_timeseries(alpha_entries(end));  % trials after last update are set to the value of last update

% beta
beta_entries = find(~isnan(beta_timeseries)); % indices of non-NaN entries

beta_per_trial(1:beta_entries(1)-1) = beta_zero;  % trials before first update are set to the prior 
for index_entry = 1:numel(beta_entries)-1
    % iterate through queue of non-NaN entries to update
    % for all indices between this non-NaN entry and the next one, update
    % to the current non-NaN value
    beta_per_trial(beta_entries(index_entry):beta_entries(index_entry+1)-1) ... 
        = beta_timeseries(beta_entries(index_entry)); 
end
beta_per_trial(beta_entries(end):end) = beta_timeseries(beta_entries(end));


%% compute MAP, precision using interpolated values
for index_trials = 1:N
    MAP_per_trial(index_trials) = Summaries_of_Beta_distribution(alpha_per_trial(index_trials),beta_per_trial(index_trials),'MAP');
    precision_per_trial(index_trials) = Summaries_of_Beta_distribution(alpha_per_trial(index_trials),beta_per_trial(index_trials),'Precision');
end

