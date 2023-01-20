function trial_type = win_stay_large_reward(trial_data)

% WIN_STAY_LARGE_REWARD checks if subject chose option again after a large reward 
% TYPE = WIN_STAY_LARGE_REWARD(TRIAL_DATA) takes the Table of data TRIAL_DATA up to the current trial, and
% returns the TYPE ('success','failure','null')
%
% 31/3/2022 Initial version
% 20/4/2022 Add complete set of conditions
% Mark Humphries 

% default is that trial did not meet criteria for win-stay 
trial_type = "null";
number_trials = size(trial_data,1);

% if more than one trial AND got the large reward on the previous trial then candidate for win-stay
if number_trials > 1 && trial_data.MaxReward(end-1) == "yes"
    if trial_data.Choice(end) == trial_data.Choice(end-1)
        % if made the same choice of option on the curren trial
        trial_type = "success";        
    else
        % switched choice
        trial_type = "failure"; 
    end
end

