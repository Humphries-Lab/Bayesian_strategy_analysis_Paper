function trial_type = win_stay_spatial(trial_data)

% WIN_STAY_SPATIAL checks if subject chose same spatial option after a reward
% TYPE = WIN_STAY_SPATIAL(TRIAL_DATA) takes the Table of data TRIAL_DATA up to the current trial, and
% returns the TYPE ('success','failure','null')
%
% Mark Humphries 31/3/2022

% default is that trial did not meet criteria for win-stay 
trial_type = "null";

number_trials = size(trial_data,1);

% if there are more than one trial, and was rewarded on the previous trial then candidate for win-stay
if number_trials > 1 && trial_data.Reward(end-1) == "yes"
    
    if trial_data.Choice(end) == trial_data.Choice(end-1)
        % chose same spatial option, so is a success
        trial_type = "success";
    else
        % got reward but chose a different option
        trial_type = "failure"; 
    end
end
