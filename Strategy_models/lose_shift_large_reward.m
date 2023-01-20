function trial_type = lose_shift_large_reward(trial_data)

% LOSE_SHIFT_LARGE_REWARD checks if subject shifted choice based on reward
% value
% TYPE = LOSE_SHIFT_LARGE_REWARD(TRIAL_DATA) takes the Table of data TRIAL_DATA up to the current trial, and
% returns the TYPE ('success','failure','null')
%
% 31/3/2022 Initial version
% 20/4/2022 Add complete set of conditions
% Mark Humphries 

% default is that trial did not meet criteria for win-stay 
trial_type = "null";
number_trials = size(trial_data,1);

% if more than 1 trial AND did not get large reward on the previous trial then is
% candidate for lose-shift
if number_trials > 1 && trial_data.MaxReward(end-1) == "no"
    if trial_data.Choice(end) ~= trial_data.Choice(end-1)
        % if swtiched choice of option on the current trial
        trial_type = "success";        
    else
        % same choice
        trial_type = "failure"; 
    end
end