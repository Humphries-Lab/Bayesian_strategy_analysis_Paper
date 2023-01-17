%%% script to plot SI Figure of MAP estimation of p(match)
clearvars; close all

run figure_properties.m

if ispc
    simOutputPath = '';
else
    simOutputPath = '/Users/mqbssmhg/Dropbox/Projects/Bayesian strategy analysis/Correct strategy is missing/';
end

% load data
load([simOutputPath 'MAP_is_P_match']);


%% panel: examples of MAP estimation at two values of gamma

gamma_plot = [0.5, 0.9];

% MAP averaged over all repeats
for iGamma = 1:numel(gamma)
    if any(gamma(iGamma) == gamma_plot)
        figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.square]); 

        for iMatch = 1:numel(p_match)
            % plot target
            line([0 t],[p_match(iMatch) p_match(iMatch)],'Color',colourmaps.match_colours(iMatch,:),'LineStyle','--')

            % plot match
            mean_MAP = mean(Results(iMatch,iGamma).MAP_for_match);
            std_MAP = std(Results(iMatch,iGamma).MAP_for_match);
            shadedErrorBar(1:t,mean_MAP,std_MAP,'lineprops',{'Color',colourmaps.match_colours(iMatch,:)});
        end
        % plot True
        shadedErrorBar(1:t,Results(numel(p_match)+1,iGamma).MAP_for_match,Results(numel(p_match)+1,iGamma).variance,'lineprops',{'Color',colourmaps.true_colour});        

        xlabel('Trials')
        ylabel('Probability of match')
        axis tight
        
        FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);

        print([exportpath 'MAP_estimation_for_Gamma_' num2str(gamma(iGamma)) '.png'],'-dpng')
        print([exportpath 'MAP_estimation_for_Gamma_' num2str(gamma(iGamma)) '.svg'],'-dsvg')
    end
end

%% plot scaling of MAP variation with gamma

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.square]); 
line([gamma(1) gamma(end)],[0.25 0.25],'Color',[0.7 0.7 0.7],'LineStyle','--')
line([gamma(1) gamma(end)],[0.25 0.25]./2,'Color',[0.7 0.7 0.7],'LineStyle','--')
hold on
for iMatch = 1:numel(p_match)
    for iGamma = 1:numel(gamma)
        % compute STD at end trial
        std_MAP_per_gamma(iGamma) = std(Results(iMatch,iGamma).MAP_for_match(:,end));
    end
    plot(gamma,std_MAP_per_gamma,'.-','Color',colourmaps.match_colours(iMatch,:)); 
end
xlabel('Evidence decay \gamma')
ylabel('S.D. of MAP')
FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);

print([exportpath 'MAP_variation_scaling_with_gamma'],'-dpng')
print([exportpath 'MAP_variation_scaling_with_gamma'],'-dsvg')
