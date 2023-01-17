%% detection of strategy switching, when strategies are stochastic and change gradually (linear model)
% basic model: two strategies are adopted with probability p1 > p2 for t
% trials; then switched to p2 > p1.
% The switch is linear after trial t: p1 decreases at rate m; p2 increases
% at the same rate.
%
% Here: we explore effects of different rates m, for the plateau case (i.e.
% t is some huge number)
%
% For different pairs (p1,p2) before and after switch, we ask:
% How many trials does it take for p(strategy2) > p(strategy)?
% [where p(strategy) is MAP estimate from Beta distribution]
% 
% Mark Humphries 29/9/2022

clearvars; close all

%% key parameters
% priors
alpha0 = 1; 
beta0 = 1;

% evidence decay - range to test
gamma = 0.5:0.1:1;

% plateau effects: duration of strategy p1 > p2
switching_trial = 500; % range of possible durations

% how many trials after reaching equality in which to check stability
stability_window = 50;

% range of probabilities
p_higher = 0.6:0.1:1;
p_lower = 0:0.1:0.5;

pairs = combvec(p_higher,p_lower)'; % combvec from Neural Networks toolbox?

% slope of probability change
m = [0.01:0.01:0.05];

% repeats of each pair
number_of_repeats = 50; % for each pair of probabilities, how many repeats

%% for each gamma, slope, and probabilities, find stopping trial

for iPair = 1:size(pairs,1)
    % for each pair of probabilities...
        
    for iGamma = 1:numel(gamma) 
        % for each gamma value
        for iSlope = 1:numel(m)
            % for each slope value
            expected_crossover(iPair,iSlope) = (pairs(iPair,1)-pairs(iPair,2))/(2*m(iSlope)); % trial since t at which p1 < p2
          
            for iRepeats = 1:number_of_repeats
                % generate initial sequence of trial outcomes (x) according to
                % initial probability pair
                x_1 = rand(switching_trial,1) < pairs(iPair,1);  % higher probability in column 1
                x_2 = rand(switching_trial,1) < pairs(iPair,2);

                % for each trial after the switch, check for p(new) > p(current)
                p_new = 0; p_current = inf; % force true check
                trial_number = switching_trial;
                while p_current >= p_new
                    trial_number = trial_number + 1;
                    t_since_switch = trial_number - switching_trial;  
                    
                    % update trial outcomes according to model
                    % 1. linear change in probability
                    [p1,p2] = linear_change_model(t_since_switch,pairs(iPair,:),m(iSlope));

                    % find new trial outcome using those probabilities for
                    % each strategy
                    x_1 = [x_1; rand < p1];
                    x_2 = [x_2; rand < p2]; 
                    
                    % update MAP probabilities for each strategy
                    [p_current, p_new] = update_probabilities(trial_number,gamma(iGamma),x_1,x_2,alpha0,beta0);
    
                end
                % store result
                Results(iPair).trials_since_t(iGamma,iSlope,iRepeats) = trial_number - switching_trial;
                Results(iPair).trials_since_crossover(iGamma,iSlope,iRepeats) = trial_number - switching_trial - expected_crossover(iPair,iSlope);
                
                % do stability check
                for iTrial = 1:stability_window
                    trial_number = trial_number + 1;
                    t_since_switch = trial_number - switching_trial;  

                    % update trial outcomes according to model
                    [p1,p2] = linear_change_model(t_since_switch,pairs(iPair,:),m(iSlope));

                    % find new trial outcome using those probabilities for each strategy
                    x_1 = [x_1; rand < p1];
                    x_2 = [x_2; rand < p2]; 

                    % update MAP probabilities for each strategy
                    [p_current_window(iTrial), p_new_window(iTrial)] = update_probabilities(trial_number,gamma(iGamma),x_1,x_2,alpha0,beta0);
                end
                % proportion of trials in which condition is met
                Results(iPair).stability(iGamma,iSlope,iRepeats) = sum(p_new_window > p_current_window) / stability_window;
            end
        end
    end
    % store means and STD over repeats
    Results(iPair).mean_Trial_since_t = mean(Results(iPair).trials_since_t,3);
    Results(iPair).std_Trial_since_t = std(Results(iPair).trials_since_t,0,3);
    Results(iPair).mean_Trial_since_crossover = mean(Results(iPair).trials_since_crossover,3);
    Results(iPair).std_Trial_since_crossover = std(Results(iPair).trials_since_crossover,0,3);
    Results(iPair).mean_stability = mean(Results(iPair).stability,3);
    Results(iPair).std_stability = std(Results(iPair).stability,0,3);

end

save('Results_for_slope_effects','Results','pairs','p_higher','p_lower','gamma','switching_trial','m')


%% plot results
close all;

%% effects of slope versus number of trials at plateau

% deterministic strategies
pair_index = find(pairs(:,1) == 1 & pairs(:,2) == 0);
plot_slope_vs_trials(m,Results,pair_index,pairs)

figure
plot(m,Results(pair_index).mean_stability')
xlabel('slope')
ylabel('stability')
exportPPTfig(gcf,['Independent_stability_for_' num2str(pairs(pair_index,1)) '_' num2str(pairs(pair_index,2)) '_crossover.png'],[pwd '\'],[10 15 10 10])


% best-case here
pair_index = find(pairs(:,1) == 0.9 & pairs(:,2) == 0.1);
plot_slope_vs_trials(m,Results,pair_index,pairs)

figure
plot(m,Results(pair_index).mean_stability')
xlabel('slope')
ylabel('stability')
exportPPTfig(gcf,['Independent_stability_for_' num2str(pairs(pair_index,1)) '_' num2str(pairs(pair_index,2)) '_crossover.png'],[pwd '\'],[10 15 10 10])


% worst-case here
pair_index = find(pairs(:,1) == 0.6 & pairs(:,2) == 0.5);
plot_slope_vs_trials(m,Results,pair_index,pairs)

figure
plot(m,Results(pair_index).mean_stability')
xlabel('slope')
ylabel('stability')
exportPPTfig(gcf,['Independent_stability_for_' num2str(pairs(pair_index,1)) '_' num2str(pairs(pair_index,2)) '_crossover.png'],[pwd '\'],[10 15 10 10])


%% functon to do plotting
function plot_slope_vs_trials(slope,Results,pair_index,pairs)
 figure
    subplot(1,2,1)
    plot(slope,Results(pair_index).mean_Trial_since_crossover(1:end-1,:)')  % leave off gamma = 1; transpose to plot for slope (rows of matrix)
    title(['Pair ' num2str(pairs(pair_index,:)) ' crossover'])
    xlabel('slope')
    ylabel('Mean trials to win since crossover')

    subplot(1,2,2)
    plot(slope,Results(pair_index).std_Trial_since_crossover(1:end-1,:)')
    xlabel('slope')
    ylabel('STD trials to win since crossover')
    exportPPTfig(gcf,['Plateau_slope_for_' num2str(pairs(pair_index,1)) '_' num2str(pairs(pair_index,2)) '_crossover.png'],[pwd '\'],[10 15 20 10])

    figure
    subplot(1,2,1)
    plot(slope,Results(pair_index).mean_Trial_since_t(1:end-1,:)')  % leave off gamma = 1; transpose to plot for slope (rows of matrix)
    title(['Pair ' num2str(pairs(pair_index,:)) ' since t'])

    xlabel('slope')
    ylabel('Mean trials to win since t')

    subplot(1,2,2)
    plot(slope,Results(pair_index).std_Trial_since_t(1:end-1,:)')
    xlabel('slope')
    ylabel('STD trials to win since t')
    exportPPTfig(gcf,['Plateau_slope_for_' num2str(pairs(pair_index,1)) '_' num2str(pairs(pair_index,2)) 'since_t.png'],[pwd '\'],[10 15 20 10])

end