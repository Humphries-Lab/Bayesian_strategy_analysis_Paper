%% script to test duration of strategy switching, when strategies are stochastic and switch abruptly
% basic model: two strategies are adopted with probability p1 > p2 for t
% trials; then switched to p2 > p1.
% For different pairs (p1,p2) before and after switch, we ask:
% How many trials does it take for p(strategy2) > p(strategy)?
% [where p(strategy) is MAP estimate from Beta distribution]
% 
% Mark Humphries 27/9/2022

clearvars; close all

%% key parameters
% priors
alpha0 = 1; 
beta0 = 1;

% evidence decay - range to test
gamma = 0.5:0.1:1;

% duration of strategy p1 > p2
switching_trial = [5:5:50 500]; % range of possible durations

% range of probabilities
p_higher = 0.6:0.1:1;
p_lower = 0:0.1:0.5;

pairs = combvec(p_higher,p_lower)'; % combvec from Neural Networks toolbox?

number_of_repeats = 50; % for each pair of probabilities, how many repeats

stability_window = 50;

%% for each gamma, duration and probabilities, find stopping trial
Pairs_Mean = zeros(numel(p_higher),numel(p_lower),numel(gamma));
Pairs_STD = Pairs_Mean;

for iPair = 1:size(pairs,1)
    % for each pair of probabilities
    for iGamma = 1:numel(gamma) 
        % for each gamma value
        for iDurations = 1:numel(switching_trial)
            % for each duration until the switch
            
            for iRepeats = 1:number_of_repeats
                % generate initial sequence of trial outcomes (x) according to
                % initial probability pair
                x_1 = rand(switching_trial(iDurations),1) < pairs(iPair,1);  % higher probability in column 1
                x_2 = rand(switching_trial(iDurations),1) < pairs(iPair,2);

                % for each trial after the switch, check for p(new) > p(current)
                p_new = 0; p_current = inf; % force true check
                trial_number = switching_trial(iDurations);
                while p_current >= p_new
                    trial_number = trial_number + 1;
                    % update trial outcomes according to model
                    x_1 = [x_1; rand < pairs(iPair,2)];  % NOTE: switch in probabilities
                    x_2 = [x_2; rand < pairs(iPair,1)]; 

                     % update MAP probabilities for detecting strategies
                    [p_current, p_new] = update_probabilities(trial_number,gamma(iGamma),x_1,x_2,alpha0,beta0);

                end
                % store result
                Results(iPair).number_of_trials(iGamma,iDurations,iRepeats) = trial_number - switching_trial(iDurations);
                
                
                % do stability check
                p_current_window = zeros(stability_window,1); p_new_window = p_current_window;
                for iTrial = 1:stability_window
                    trial_number = trial_number + 1;

                    % update trial outcomes according to model
                    x_1 = [x_1; rand < pairs(iPair,2)];  % NOTE: switch in probabilities
                    x_2 = [x_2; rand < pairs(iPair,1)]; 

                    % update MAP probabilities for each strategy
                    [p_current_window(iTrial), p_new_window(iTrial)] = update_probabilities(trial_number,gamma(iGamma),x_1,x_2,alpha0,beta0);
                end
                % proportion of trials in which condition is met
                Results(iPair).stability(iGamma,iDurations,iRepeats) = sum(p_new_window > p_current_window) / stability_window;
                
                
            end
        end
    end
    mean_Trial = mean(Results(iPair).number_of_trials,3);
    std_Trial = std(Results(iPair).number_of_trials,0,3);
    mean_stability = mean(Results(iPair).stability,3);
    std_stability = std(Results(iPair).stability,0,3);

    
    if pairs(iPair,1) == 0.8 && pairs(iPair,2) == 0.4
        % for this pair, plot mean over repeats
        figure
        subplot(1,2,1)
        title(['Mean for pair ' num2str(pairs(iPair,:))])
        plot(switching_trial(1:end-1),mean_Trial(:,1:end-1))
        xlabel('switching trial')
        ylabel('mean number of trials to win')

        subplot(1,2,2)
        title(['STD for pair ' num2str(pairs(iPair,:))])
        plot(switching_trial(1:end-1),std_Trial(:,1:end-1))
        xlabel('switching trial')
        ylabel('STD of number of trials to win')
        exportPPTfig(gcf,['ExampleStochasticSwitch_' num2str(pairs(iPair,1)) '_' num2str(pairs(iPair,2)) '.png'],[pwd '\'],[10 15 20 10])
    end
    
    % store summary statistics per pair as longest-duration (i.e. ~plateau)
    ixHigher = find(p_higher == pairs(iPair,1));
    ixLower = find(p_lower == pairs(iPair,2));
    
    Summary.Pairs_Mean(ixHigher,ixLower,:) = mean_Trial(:,end);
    Summary.Pairs_STD(ixHigher,ixLower,:) = std_Trial(:,end);
    Summary.Pairs_Mean_stability(ixHigher,ixLower,:) = mean_stability(:,end);
    Summary.Pairs_STD_stability(ixHigher,ixLower,:) = std_stability(:,end);
    
end

save('Results_AbruptStochastic_variable_t','Results','Summary','pairs','p_higher','p_lower','gamma','switching_trial')

%% plot plateau results per gamma
close all;

lower_bounds = [-0.1, 0.5];

% get colormap
cmap = brewermap(15,'BuPu');  % colormap for all gamma except gamma=1
cmap_nodecay = brewermap(15,'Greys');

% find maxima and minima
set_of_means = Summary.Pairs_Mean(:,:,1:end-1); % omit gamma = 1
set_of_STD = Summary.Pairs_STD(:,:,1:end-1);
[min_mean,max_mean] = bounds(set_of_means(:));
[min_std,max_std] = bounds(set_of_STD(:));

for iGamma = 1:numel(gamma)-1 % omit no-decay: plot separately
    figure
    subplot(121),
    colormap(cmap)
    plotMatrix(p_lower,p_higher,squeeze(Summary.Pairs_Mean(:,:,iGamma)),lower_bounds)
    set(gca,'CLim',[min_mean, max_mean]);  % scale colormap the same in all figures
    title(['Mean for \gamma = ' num2str(gamma(iGamma))])
    axis square
    colorbar
    
    subplot(122),
    colormap(cmap)
    plotMatrix(p_lower,p_higher,squeeze(Summary.Pairs_STD(:,:,iGamma)),lower_bounds)
    set(gca,'CLim',[min_std, max_std]);
    title(['STD for \gamma = ' num2str(gamma(iGamma))])
    axis square
    colorbar
    
    exportPPTfig(gcf,['MeanSTD_trials_StochasticSwitch_' num2str(gamma(iGamma)) '.png'],[pwd '\'],[10 15 12 6])
end

% plot gamma = 1 result
figure
subplot(121),
colormap(cmap_nodecay)
plotMatrix(p_lower,p_higher,squeeze(Summary.Pairs_Mean(:,:,end)),lower_bounds)
title(['Mean for \gamma = ' num2str(gamma(end))])
axis square
colorbar

subplot(122),
colormap(cmap_nodecay)
plotMatrix(p_lower,p_higher,squeeze(Summary.Pairs_STD(:,:,end)),lower_bounds)
title(['STD for \gamma = ' num2str(gamma(end))])
axis square
colorbar
exportPPTfig(gcf,['MeanSTD_trials_StochasticSwitch_Gamma1'],[pwd '\'],[10 15 15 7])

%% plots of stability
% find maxima and minima
set_of_means = Summary.Pairs_Mean_stability(:);
set_of_STD = Summary.Pairs_STD_stability(:);
[min_mean,max_mean] = bounds(set_of_means(:));
[min_std,max_std] = bounds(set_of_STD(:));

for iGamma = 1:numel(gamma)
    figure
    subplot(121),
    colormap(cmap)
    plotMatrix(p_lower,p_higher,squeeze(Summary.Pairs_Mean_stability(:,:,iGamma)),lower_bounds)
    set(gca,'CLim',[min_mean, max_mean]);  % scale colormap the same in all figures
    title(['Mean for \gamma = ' num2str(gamma(iGamma))])
    axis square
    colorbar
    
    subplot(122),
    colormap(cmap)
    plotMatrix(p_lower,p_higher,squeeze(Summary.Pairs_STD_stability(:,:,iGamma)),lower_bounds)
    set(gca,'CLim',[min_std, max_std]);
    title(['STD for \gamma = ' num2str(gamma(iGamma))])
    axis square
    colorbar
    
    exportPPTfig(gcf,['MeanSTD_stability_StochasticSwitch_' num2str(gamma(iGamma)) '.png'],[pwd '\'],[10 15 12 6])
end