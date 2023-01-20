function trial_type = win_stay_cued(trial_data)

% WIN_STAY_CUED checks if subject chose cued option again after a reward on
% cued option
% TYPE = WIN_STAY_CUED(TRIAL_DATA) takes the Table of data TRIAL_DATA up to the current trial, and
% returns the TYPE ('success','failure','null')
%
% 31/3/2022 Initial version
% 20/4/2022 Add complete set of conditions
% Mark Humphries 

% default is that trial did not meet criteria for win-stay 
trial_type = "null";
number_trials = size(trial_data,1);

% if more than one trial AND was rewarded on the previous trial then candidate for win-stay
if number_trials > 1 && trial_data.Reward(end-1) == "yes"
    if trial_data.Choice(end-1) == trial_data.CuePosition(end-1) && trial_data.Choice(end) == trial_data.CuePosition(end)       
        % if chose the cued option on previous trial AND chose an
        % cued option on this trial too
        trial_type = "success";
    elseif trial_data.Choice(end-1) ~= trial_data.CuePosition(end-1) && ...
            trial_data.Choice(end) ~= trial_data.CuePosition(end)       
        % if did NOT choose the cued option on the previous trial AND chose an
        % uncued option on this trial too
        trial_type = "success";    
    else
        % switched option based on cue
        trial_type = "failure"; 
    end
end
