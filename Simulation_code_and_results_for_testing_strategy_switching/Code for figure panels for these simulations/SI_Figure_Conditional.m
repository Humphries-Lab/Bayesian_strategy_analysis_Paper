%%% script to plot results for detecting the conditional strategy - SI Figure
clearvars; close all

run figure_properties.m

if ispc
    simOutputPath = 'C:\Users\lpzmdh\Dropbox\Projects\Bayesian strategy analysis\Analytical_changes_of_distribution\';
else
    simOutputPath = '/Users/mqbssmhg/Dropbox/Projects/Bayesian strategy analysis/Analytical_changes_of_distribution/';
end

load([simOutputPath 'Results_for_detecting_new_strategy']);

%% dependence of detection on p(meet condition)
for iSlope = 1:numel(m)
    figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.square]); 
    set(gca,'ColorOrder',colourmaps.gamma_lines,'NextPlot', 'replacechildren');
    plot(p_meet_condition, squeeze(Results.mean_Trials_to_detect(:,:,iSlope)));
            set(gca,'yLim',[0 inf])
    xlabel('P(meet condition)')
    ylabel('Trials')
    FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);

    print([exportpath 'trials_to_detect_conditional' num2str(m(iSlope)) '.png'],'-dpng')
    print([exportpath 'trials_to_detect_conditional' num2str(m(iSlope)) '.svg'],'-dsvg')

end

