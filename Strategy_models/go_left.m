function trial_type = go_left(trial_data)

% GO_LEFT checks if subject chose the left option on this trial
% TYPE = GO_LEFT(TRIAL_DATA) takes the Table of data TRIAL_DATA up to the current trial, and
% returns the TYPE ('success','failure')
%
% Mark Humphries 31/3/2022

% check only the current trial
if trial_data.Choice(end) == "left"
    trial_type = "success";
else 
    trial_type = "failure";
end