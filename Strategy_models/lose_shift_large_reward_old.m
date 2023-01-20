function trial_type = lose_shift_large_reward(trial_data)

% LOSE_SHIFT_CUED checks if subject shifted choice based on the cue, after no
% reward
% TYPE = LOSE_SHIFT_CUED(TRIAL_DATA) takes the Table of data TRIAL_DATA up to the current trial, and
% returns the TYPE ('success','failure','null')
%
% NOTE: this is lose-shift based on the cue, so is either:
% (a) chose cued option when unrewarded, then chose uncued option on next trial
% (b) chose uncued option when unrewarded, then chose cued option on next
% trial
% 
% 31/3/2022 Initial version
% 20/4/2022 Add complete set of conditions
% Mark Humphries 

% default is that trial did not meet criteria for win-stay 
trial_type = "null";
number_trials = size(trial_data,1);

% if more than 1 trial AND was NOT rewarded on the previous trial then is
% candidate for lose-shift
if number_trials > 1 && trial_data.MaxReward(end-1) == "no"
    
    if trial_data.Choice(end-1) == trial_data.OptReward(end-1) && ...
            trial_data.Choice(end) ~= trial_data.OptReward(end)       
        % if chose cued option on previous trial AND chose an 
        % uncued option on this trial
        trial_type = "success";
    elseif trial_data.Choice(end-1) ~= trial_data.OptReward(end-1) && ...
            trial_data.Choice(end) == trial_data.OptReward(end)       
        % if chose an uncued option on the previous trial AND did choose the
        % cued option on this trial
        trial_type = "success";
    else
        % chose the same cue-based option on the current trial
        trial_type = "failure"; 
    end
end
