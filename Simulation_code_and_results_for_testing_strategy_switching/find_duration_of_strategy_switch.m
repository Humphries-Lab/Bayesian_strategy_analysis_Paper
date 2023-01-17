%% script to test duration of strategy switching
% basic model: the current strategy is executed for t trials, from trial 1
% a new strategy is then adopted.
% For different models of how the new strategy is adopted (abrupt, gradual
% etc), we ask:
% How many trials does it take for p(new) > p(current)?
% [where p(strategy) is MAP estimate from Beta distribution]
% 
% Mark Humphries 27/9/2022

clearvars;

%% key parameters
% priors
alpha0 = 1; 
beta0 = 1;

% evidence decay - range to test
gamma = 0.5:0.1:1;

% duration of current strategy
switching_trial = [5:5:50]; % range of possible durations

%% for each gamma, duration and model, find stopping trial
for iGamma = 1:numel(gamma)
    
    for iDurations = 1:numel(switching_trial)
        % models for initial sequence of trial outcomes (x)
        % worst-case scenario: only executed current strategy for all t trials
        x_current = ones(switching_trial(iDurations),1);
        x_new = zeros(switching_trial(iDurations),1);
    
        % for each trial after the switch, check for p(new) > p(current)
        p_new = 0; p_current = inf; % force true check
        trial_number = switching_trial(iDurations);
        while p_current >= p_new
            trial_number = trial_number + 1;
            % update trial outcomes according to model
            % [next_x_current,next_x_new] = abrupt_switch;
            [next_x_current,next_x_new] = abrupt_partial(trial_number);
            x_current = [x_current; next_x_current];
            x_new = [x_new; next_x_new];

            trial_sequence = [1:trial_number]';
            % compute s* and f* (Eqs 8 and 9 in the paper)
            s_star_current = sum(gamma(iGamma).^(trial_number - trial_sequence) .* x_current);
            s_star_new = sum(gamma(iGamma).^(trial_number - trial_sequence) .* x_new);
            f_star_current = sum(gamma(iGamma).^(trial_number - trial_sequence) .* (1-x_current));
            f_star_new = sum(gamma(iGamma).^(trial_number - trial_sequence) .* (1 - x_new));

            % update alpha, beta (text after Eq 10 in the paper)
            alpha_current = alpha0 + s_star_current;
            alpha_new = alpha0 + s_star_new;
            beta_current = beta0 + f_star_current;
            beta_new = beta0 + f_star_new;

            % compute MAP probabilities
            p_current = (alpha_current - 1) / (alpha_current + beta_current - 2);
            p_new = (alpha_new - 1) / (alpha_new + beta_new - 2);
        end
        % store result
        number_of_trials(iGamma,iDurations) = trial_number - switching_trial(iDurations);
    end
end

plot(switching_trial,number_of_trials)
xlabel('switching trial')
ylabel('number of trials to win')
exportPPTfig(gcf,'AbruptSwitch',[pwd '\'])
