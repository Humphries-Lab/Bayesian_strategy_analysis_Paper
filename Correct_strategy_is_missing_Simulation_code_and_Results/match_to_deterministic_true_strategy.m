%% script to check convergence of P(strategy) on P(match)
% for deterministic true strategy
%
% 14/10/2022 Mark Humphries

clearvars; close all;

%% parameters

% decay rate
gamma = 0.5:0.1:1;

% prior
alpha0 = 1;
beta0 = 1;

% duration
t = 50;

% P(match) to test
p_match = [0.25, 0.5, 0.75];

% how many repetitions
number_of_repeats = 50;

%% compute MAP and precision for each tested strategy


% for each P(match)
for iMatch = 1:numel(p_match)
    % for each decay rate
    for iGamma = 1:numel(gamma)
        % for each repetition
        for iRepeat = 1:number_of_repeats
            %   draw sequence of x(t)
            x_Match = rand(t,1) < p_match(iMatch);

            %   compute MAP and precision
            for iTrial = 1:t
                [Results(iMatch,iGamma).MAP_for_match(iRepeat,iTrial), Results(iMatch,iGamma).variance(iRepeat,iTrial)] = update_probabilities(iTrial,gamma(iGamma),x_Match(1:iTrial),alpha0,beta0);
            end

        end
    end
end

%% compute MAP and precision for true strategy...
x_True = ones(t,1);
for iGamma = 1:numel(gamma)
    for iTrial = 1:t
        [Results(numel(p_match)+1,iGamma).MAP_for_match(iTrial), Results(numel(p_match)+1,iGamma).variance(iTrial)] = update_probabilities(iTrial,gamma(iGamma),x_True(1:iTrial),alpha0,beta0);
    end    
end

%% plot results

match_colours = brewermap(3,'Dark2');
true_colour = [0.7 0.7 0.7];

% plot an example run of each
ixRepeat = 1;

for iGamma = 1:numel(gamma)
    
    figure
    for iMatch = 1:numel(p_match)
        % plot target
        line([0 t],[p_match(iMatch) p_match(iMatch)],'Color',match_colours(iMatch,:),'LineStyle','--')

        % plot match
        shadedErrorBar(1:t,Results(iMatch,iGamma).MAP_for_match(ixRepeat,:),Results(iMatch,iGamma).variance(ixRepeat,:),'lineprops',{'Color',match_colours(iMatch,:)});        
    end
    % plot True
    shadedErrorBar(1:t,Results(numel(p_match)+1,iGamma).MAP_for_match,Results(numel(p_match)+1,iGamma).variance,'lineprops',{'Color',true_colour});        
   
    xlabel('Trials')
    ylabel('Probability of match')
    title(['Example time-series for gamma = ' num2str(gamma(iGamma))])
    % exportPPTfig(gcf,['example_matching_of_p_strategy_' num2str(gamma(iGamma)) '.png'],[pwd '\'],[10 15 10 10])

end

% MAP averaged over all repeats
for iGamma = 1:numel(gamma)
    figure
    for iMatch = 1:numel(p_match)
        % plot target
        line([0 t],[p_match(iMatch) p_match(iMatch)],'Color',match_colours(iMatch,:),'LineStyle','--')

        % plot match
        mean_MAP = mean(Results(iMatch,iGamma).MAP_for_match);
        std_MAP = std(Results(iMatch,iGamma).MAP_for_match);
        shadedErrorBar(1:t,mean_MAP,std_MAP,'lineprops',{'Color',match_colours(iMatch,:)});
    end
    % plot True
    shadedErrorBar(1:t,Results(numel(p_match)+1,iGamma).MAP_for_match,Results(numel(p_match)+1,iGamma).variance,'lineprops',{'Color',true_colour});        
    
    xlabel('Trials')
    ylabel('Probability of match')
    title(['Average + STD time-series for gamma = ' num2str(gamma(iGamma))])
    
    % exportPPTfig(gcf,['Mean_matching_of_p_strategy_' num2str(gamma(iGamma)) '.png'],[pwd '\'],[10 15 10 10])

end


%% scaling of MAP variation with gamma
figure
line([gamma(1) gamma(end)],[0.25 0.25],'Color',[0.7 0.7 0.7],'LineStyle','--')
line([gamma(1) gamma(end)],[0.25 0.25]./2,'Color',[0.7 0.7 0.7],'LineStyle','--')
hold on
for iMatch = 1:numel(p_match)
    for iGamma = 1:numel(gamma)
        % compute STD at end trial
        std_MAP_per_gamma(iGamma) = std(Results(iMatch,iGamma).MAP_for_match(:,end));
    end
    plot(gamma,std_MAP_per_gamma,'.-','Color',match_colours(iMatch,:)); 
end
xlabel('Gamma')
ylabel('STD of MAP')

% exportPPTfig(gcf,['Scalingof_MAP_std_with_Gamma.png'],[pwd '\'],[10 15 10 10])


save MAP_is_P_match Results gamma p_match t