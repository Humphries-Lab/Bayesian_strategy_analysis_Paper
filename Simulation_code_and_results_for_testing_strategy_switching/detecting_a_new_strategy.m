%% script to look at time until new strategy is detected
% i.e. how many trials does it take to detect strategy_new
% starting from no competing strategy
% 
% Define detection as 1-1/n
%
% This is general script, but we're looking at it notably for new
% conditional strategies
%
% Mark Humphries 7/11/2022

clearvars; close all

%% key parameters
% priors
alpha0 = 1; 
beta0 = 1;

% how many options in task?
N=2;

% evidence decay - range to test
gamma = 0.5:0.1:1;

% how many trials after reaching criterion in which to check stability
stability_window = 50;

% probability of meeting condition for strategy
p_meet_condition = [0.1:0.1:1];  % 1=deterministic 

% slope of probability change (how quickly reach p1=1, from p1=0)
m = [0.05 0.1 0.2 1];  % m=1 is abrupt switch to the new strategy

% repeats of each pair
number_of_repeats = 50; % for each pair of probabilities, how many repeats

%% run simulations
for iCondition = 1:numel(p_meet_condition)
    % for each condition probability
        
    for iGamma = 1:numel(gamma) 
        % for each gamma value
        for iSlope = 1:numel(m)
            % for each slope value
          
            for iRepeats = 1:number_of_repeats
                p_new = 0; p_current = inf; % force true check
                trial_number = 0; x_1 = [];
                while p_new <= 1-(1/N)
                    trial_number = trial_number+1;
                    % update probability of strategy
                    p1 = m(iSlope)* trial_number;
                    if p1 > 1, p1 = 1; end
                    
                    % draw new success/fail/null
                    if rand > p_meet_condition(iCondition)
                        x_1 = [x_1; nan];
                        p_new = p_new;  % just to emphasise that p_new is unchanged
                    else
                        x_1 = [x_1; rand < p1];
                        % update MAP probabilities for new strategy (enter 0 for
                        % second strategy here)
                        [p_new, ~] = update_probabilities(trial_number,gamma(iGamma),x_1,0,alpha0,beta0);
  
                    end
                end
                Results.trials_to_detect(iCondition,iGamma,iSlope,iRepeats) = trial_number;
                
                % do stability check
                p_new_window = p_new;
                for iTrial = 1:stability_window
                   trial_number = trial_number + 1;

                   % update probability of strategy
                    p1 = m(iSlope)* trial_number;
                    if p1 > 1, p1 = 1; end

                    % draw new success/fail/null
                    if rand > p_meet_condition(iCondition)
                        x_1 = [x_1; nan];
                        if iTrial == 1
                            p_new_window(iTrial) = p_new;
                        else
                        p_new_window(iTrial) = p_new_window(iTrial-1); 
                        end
                    else
                        x_1 = [x_1; rand < p1];
                        % update MAP probabilities for new strategy (enter 0 for
                        % second strategy here)
                        [p_new_window(iTrial), ~] = update_probabilities(trial_number,gamma(iGamma),x_1,0,alpha0,beta0);
  
                    end
                end
                % proportion of trials in which condition is met
                Results.stability(iCondition,iGamma,iSlope,iRepeats) = sum(p_new_window > 1-(1/N)) / stability_window;
            end
        end
    end
end
% store means and STD over repeats
Results.mean_Trials_to_detect = mean(Results.trials_to_detect,4);
Results.std_Trials_to_detect = std(Results.trials_to_detect,0,4);
Results.mean_stability = mean(Results.stability,4);
Results.std_stability = std(Results.stability,0,4);

save('Results_for_detecting_new_strategy','Results','p_meet_condition','N','gamma','stability_window','m')
    
%% plot results

close all

%% abrupt switches...
% ixAbrupt = find(m == 1);
% figure
% plot(p_meet_condition, squeeze(Results.mean_Trials_to_detect(:,:,ixAbrupt)));
% set(gca,'yLim',[0 11])
% xlabel('P(meet condition)')
% ylabel('Trials to detect strategy')

%% loop over m....
legend_text = arrayfun(@num2str, gamma, 'UniformOutput', 0);
for iSlope = 1:numel(m)
    figure
    plot(p_meet_condition, squeeze(Results.mean_Trials_to_detect(:,:,iSlope)));
        set(gca,'yLim',[0 inf])
    xlabel('P(meet condition)')
    ylabel('Trials to detect strategy')
    title(['Slope = ' num2str(m(iSlope))])
    if iSlope==numel(m)
        legend(legend_text);
    end
    exportPPTfig(gcf,['Conditional_trials_to_detect_m_' num2str(m(iSlope)) '.png'],[pwd '\'],[10 15 10 10])
    
    figure
    plot(p_meet_condition, squeeze(Results.mean_stability(:,:,iSlope)));
        set(gca,'yLim',[0 inf])
    xlabel('P(meet condition)')
    ylabel('Stability after detecting strategy')
    title(['Slope = ' num2str(m(iSlope))])
    if iSlope==numel(m)
        legend(legend_text);
    end  
    exportPPTfig(gcf,['Conditional_stability_m_' num2str(m(iSlope)) '.png'],[pwd '\'],[10 15 10 10])
    
end

