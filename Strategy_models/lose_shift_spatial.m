function trial_type = lose_shift_spatial(trial_data)

% LOSE_SHIFT_SPATIAL checks if subject chose different spatial option after no reward
% TYPE = LOSE_SHIFT_SPATIAL(TRIAL_DATA) takes the Table of data TRIAL_DATA up to the current trial, and
% returns the TYPE ('success','failure','null')
%
% Mark Humphries 31/3/2022

% default is that trial did not meet criteria for lose-shift 
trial_type = "null";
number_trials = size(trial_data,1);

% if more than one trial, and was not rewarded on the previous trial then candidate for lose-shift
if number_trials > 1 && trial_data.Reward(end-1) == "no" 
    
    if trial_data.Choice(end) ~= trial_data.Choice(end-1)
        % chose different spatial option, so is a success
        trial_type = "success";
    else
        % got reward but chose a different option
        trial_type = "failure"; 
    end
end
