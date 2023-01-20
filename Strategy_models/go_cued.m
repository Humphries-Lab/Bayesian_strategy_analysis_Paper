function trial_type = go_cued(trial_data)

% GO_CUED checks if subject chose the cued option on this trial
% TYPE = GO_CUED(TRIAL_DATA) takes the Table of data TRIAL_DATA up to the current trial, and
% returns the TYPE ('success','failure')
%
% Mark Humphries 31/3/2022

% check only the current trial
% if the chosen option matched the cued option
if trial_data.Choice(end) == trial_data.CuePosition(end)
    trial_type = "success";
else 
    trial_type = "failure";
end