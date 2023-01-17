%% detection of strategy switching, when strategies are stochastic and change gradually (linear model)
% basic model: two strategies are adopted with probability p1 > p2 for t
% trials; then switched to p2 > p1.
% The switch is linear after trial t: p1 decreases at rate m; p2 increases
% at the same rate.
%
% For different pairs (p1,p2) before and after switch, we ask:
% How many trials does it take for p(strategy2) > p(strategy)?
% [where p(strategy) is MAP estimate from Beta distribution]
% 
% This script looks at one slope for the linear model;
%
% Mark Humphries 29/9/2022

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

% slope of probability change
m = 0.03;

% repeats of each pair
number_of_repeats = 50; % for each pair of probabilities, how many repeats

% length of window to compute stability
stability_window = 50;

%% for each gamma, duration and probabilities, find stopping trial
Pairs_Mean_since_t = zeros(numel(p_higher),numel(p_lower),numel(gamma));
Pairs_STD_since_t = Pairs_Mean_since_t;

for iPair = 1:size(pairs,1)
    % for each pair of probabilities...
    
    expected_crossover(iPair) = (pairs(iPair,1)-pairs(iPair,2))/(2*m); % trial since t at which p1 < p2
    
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
                    t_since_switch = trial_number - switching_trial(iDurations);

                    % update probabilities of executing each strategy
                    [p1,p2] = linear_change_model(t_since_switch,pairs(iPair,:),m);

                    % find new trial outcome using those probabilities for
                    % each strategy
                    x_1 = [x_1; rand < p1];
                    x_2 = [x_2; rand < p2]; 
                    
                    % update MAP probabilities for detecting strategies
                    [p_current, p_new] = update_probabilities(trial_number,gamma(iGamma),x_1,x_2,alpha0,beta0);

                end
                % store result
                Results(iPair).trials_since_t(iGamma,iDurations,iRepeats) = trial_number - switching_trial(iDurations);
                Results(iPair).trials_since_crossover(iGamma,iDurations,iRepeats) = trial_number - switching_trial(iDurations) - expected_crossover(iPair);
                
                % do stability check
                p_current_window = zeros(stability_window,1); p_new_window = p_current_window;
                for iTrial = 1:stability_window
                    trial_number = trial_number + 1;
                    t_since_switch = trial_number - switching_trial(iDurations);  

                    % update trial outcomes according to model
                    [p1,p2] = linear_change_model(t_since_switch,pairs(iPair,:),m);

                    % find new trial outcome using those probabilities for each strategy
                    x_1 = [x_1; rand < p1];
                    x_2 = [x_2; rand < p2]; 

                    % update MAP probabilities for each strategy
                    [p_current_window(iTrial), p_new_window(iTrial)] = update_probabilities(trial_number,gamma(iGamma),x_1,x_2,alpha0,beta0);
                end
                % proportion of trials in which condition is met
                Results(iPair).stability(iGamma,iDurations,iRepeats) = sum(p_new_window > p_current_window) / stability_window;
            end
        end
    end
    
    mean_Trial_since_t = mean(Results(iPair).trials_since_t,3);
    std_Trial_since_t = std(Results(iPair).trials_since_t,0,3);
    mean_Trial_since_crossover = mean(Results(iPair).trials_since_crossover,3);
    std_Trial_since_crossover = std(Results(iPair).trials_since_crossover,0,3);
    mean_stability = mean(Results(iPair).stability,3);
    std_stability = std(Results(iPair).stability,0,3);

    if pairs(iPair,1) == 0.8 && pairs(iPair,2) == 0.4
        % for this pair, plot mean over repeats
        figure
        subplot(1,2,1)
        title(['Mean for pair ' num2str(pairs(iPair,:))])
        plot(switching_trial(1:end-1),mean_Trial_since_t(:,1:end-1))
        xlabel('switching trial')
        ylabel('Mean trials to win since t')

        subplot(1,2,2)
        title(['STD for pair ' num2str(pairs(iPair,:))])
        plot(switching_trial(1:end-1),std_Trial_since_t(:,1:end-1))
        xlabel('switching trial')
        ylabel('STD trials to win since t')
        exportPPTfig(gcf,['ExampleStochasticSwitch_' num2str(pairs(iPair,1)) '_' num2str(pairs(iPair,2)) '.png'],[pwd '\'],[10 15 20 10])
    end
    
    % store summary statistics per pair as longest-duration (i.e. ~plateau)
    ixHigher = find(p_higher == pairs(iPair,1));
    ixLower = find(p_lower == pairs(iPair,2));
    
    Summary.Pairs_Mean_since_t(ixHigher,ixLower,:) = mean_Trial_since_t(:,end);
    Summary.Pairs_STD_since_t(ixHigher,ixLower,:) = std_Trial_since_t(:,end);
    Summary.Pairs_Mean_since_crossover(ixHigher,ixLower,:) = mean_Trial_since_crossover(:,end);
    Summary.Pairs_STD_since_crossover(ixHigher,ixLower,:) = std_Trial_since_crossover(:,end);
    Summary.Pairs_Mean_stability(ixHigher,ixLower,:) = mean_stability(:,end);
    Summary.Pairs_STD_stability(ixHigher,ixLower,:) = std_stability(:,end);
end

fname = ['Results_Gradual_variable_t_with_slope_' num2str(m) '.mat'];
save(fname,'Results','Summary','pairs','p_higher','p_lower','gamma','switching_trial','m')

pause

%% plot plateau results per gamma, for trials since t
close all

lower_bounds = [-0.1, 0.5];

% get colormap
cmap = brewermap(15,'BuPu');  % colormap for all gamma except gamma=1
cmap_nodecay = brewermap(15,'Greys');

% find maxima and minima
set_of_means = Summary.Pairs_Mean_since_t(:,:,1:end-1); % omit gamma = 1
set_of_STD = Summary.Pairs_STD_since_t(:,:,1:end-1);
[min_mean,max_mean] = bounds(set_of_means(:));
[min_std,max_std] = bounds(set_of_STD(:));

for iGamma = 1:numel(gamma)-1 % omit no-decay: plot separately
    figure
    subplot(121),
    colormap(cmap)
    plotMatrix(p_lower,p_higher,squeeze(Summary.Pairs_Mean_since_t(:,:,iGamma)),lower_bounds);
    set(gca,'CLim',[min_mean, max_mean]);  % scale colormap the same in all figures
    title(['Mean for \gamma = ' num2str(gamma(iGamma))])
    axis square
    colorbar
    
    subplot(122),
    colormap(cmap)
    plotMatrix(p_lower,p_higher,squeeze(Summary.Pairs_STD_since_t(:,:,iGamma)),lower_bounds);
    set(gca,'CLim',[min_std, max_std]);
    title(['STD for \gamma = ' num2str(gamma(iGamma))])
    axis square
    colorbar
    
    %exportPPTfig(gcf,['PlateauGradualStochasticSwitch_' num2str(gamma(iGamma)) '_since_t.png'],[pwd '\'],[10 15 12 6])
end

% plot gamma = 1 result
figure
subplot(121),
colormap(cmap_nodecay)
    plotMatrix(p_lower,p_higher,squeeze(Summary.Pairs_Mean_since_t(:,:,end)),lower_bounds);
title(['Mean for \gamma = ' num2str(gamma(end))])
axis square
colorbar

subplot(122),
colormap(cmap_nodecay)
plotMatrix(p_lower,p_higher,squeeze(Summary.Pairs_STD_since_t(:,:,end)),lower_bounds);
title(['STD for \gamma = ' num2str(gamma(end))])
axis square
colorbar
%exportPPTfig(gcf,['PlateauGradualStochasticSwitch_Gamma1_since_t'],[pwd '\'],[10 15 15 7])

%% plot plateau results per gamma, for trials since crossover

% get colormap
cmap = brewermap(15,'BuPu');  % colormap for all gamma except gamma=1
cmap_nodecay = brewermap(15,'Greys');

% find maxima and minima
set_of_means = Summary.Pairs_Mean_since_crossover(:,:,1:end-1); % omit gamma = 1
set_of_STD = Summary.Pairs_STD_since_crossover(:,:,1:end-1);
[min_mean,max_mean] = bounds(set_of_means(:));
[min_std,max_std] = bounds(set_of_STD(:));

for iGamma = 1:numel(gamma)-1 % omit no-decay: plot separately
    figure
    subplot(121),
    colormap(cmap)
    plotMatrix(p_lower,p_higher,squeeze(Summary.Pairs_Mean_since_crossover(:,:,iGamma)),lower_bounds);
    set(gca,'CLim',[min_mean, max_mean]);  % scale colormap the same in all figures
    title(['Mean for \gamma = ' num2str(gamma(iGamma))])
    axis square
    colorbar
    
    subplot(122),
    colormap(cmap)
    plotMatrix(p_lower,p_higher,squeeze(Summary.Pairs_STD_since_crossover(:,:,iGamma)),lower_bounds);
    set(gca,'CLim',[min_std, max_std]);
    title(['STD for \gamma = ' num2str(gamma(iGamma))])
    axis square
    colorbar
    
    exportPPTfig(gcf,['PlateauGradualStochasticSwitch_' num2str(gamma(iGamma)) '_since_crossover.png'],[pwd '\'],[10 15 12 6])
end

% plot gamma = 1 result
figure
subplot(121),
colormap(cmap_nodecay)
plotMatrix(p_lower,p_higher,squeeze(Summary.Pairs_Mean_since_t(:,:,end)),lower_bounds);
title(['Mean for \gamma = ' num2str(gamma(end))])
axis square
colorbar

subplot(122),
colormap(cmap_nodecay)
plotMatrix(p_lower,p_higher,squeeze(Summary.Pairs_STD_since_t(:,:,end)),lower_bounds);
title(['STD for \gamma = ' num2str(gamma(end))])
axis square
colorbar
exportPPTfig(gcf,['PlateauGradualSwitch_Gamma1_since_crossover'],[pwd '\'],[10 15 15 7])

%% stability results...
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
    
    exportPPTfig(gcf,['MeanSTD_stability_GradualSwitch_' num2str(gamma(iGamma)) '.png'],[pwd '\'],[10 15 12 6])
end

