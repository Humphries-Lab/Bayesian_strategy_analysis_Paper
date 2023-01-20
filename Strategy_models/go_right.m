function trial_type = go_right(trial_data)

% GO_RIGHT checks if subject chose the right-hand option on this trial
% TYPE = GO_RIGHT(TRIAL_DATA) takes the Table of data TRIAL_DATA up to the current trial, and
% returns the TYPE ('success','failure')
%
% Mark Humphries 31/3/2022

% check only the current trial
if trial_data.Choice(end) == "right"
    trial_type = "success";
else 
    trial_type = "failure";
end