function trial_type = sticky(trial_data)

% STICKY checks if subject chose same spatial option on consecutive trials
% TYPE = STICKY(TRIAL_DATA) takes the Table of data TRIAL_DATA up to the current trial, and
% returns the TYPE ('success','failure')
%
% Mark Humphries 31/3/2022

number_trials = size(trial_data,1);
    
if number_trials > 1 && trial_data.Choice(end) == trial_data.Choice(end-1)
    % chose same spatial option on consecutive trials, so is a success
    trial_type = "success";
else
    % chose a different option on this trial
    trial_type = "failure"; 
end
