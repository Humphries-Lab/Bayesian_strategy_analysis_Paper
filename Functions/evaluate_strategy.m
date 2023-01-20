function trial_type = evaluate_strategy(strStrategy,trial_data)

% EVALUATE_STRATEGY helper function for update of specified strategy
% TYPE = EVALUATE_STRATEGY(STRATEGY,TRIAL_DATA) executes the strategy model 
% in string STRATEGY given TRIAL_DATA, the Table of behavioural data up to 
% the current trial.
%
% Returns current trial's TYPE for the specified strategy:
% "success", "failure", or "null.
%
% The string STRATEGY must match the name of a strategy model function in
% the Strategy_models/ folder.
% 
% Mark Humphries 5/4/2022

fcn_handle = str2func(char(strStrategy)); % turn string into function
trial_type = fcn_handle(trial_data); % execute strategy function and return trial type