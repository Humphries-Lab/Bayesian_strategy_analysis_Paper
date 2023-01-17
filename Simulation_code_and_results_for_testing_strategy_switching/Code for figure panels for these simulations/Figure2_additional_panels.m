%%% script to plot results from
clearvars; close all

run figure_properties.m

if ispc
    simOutputPath = '';
else
    simOutputPath = '/Users/mqbssmhg/Dropbox/Projects/Bayesian strategy analysis/Analytical_changes_of_distribution/';
end

%% panel: effect of trials-before-switch on detection speed [abrupt switch, deterministic]

load([simOutputPath 'Results_AbruptStochastic_variable_t']);

ixPair = find(pairs(:,1) == 1 & pairs(:,2) == 0);
Trials_to_detection = mean(Results(ixPair).number_of_trials,3);

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.square]); 
set(gca,'ColorOrder',colourmaps.gamma_lines,'NextPlot', 'replacechildren');
plot(switching_trial(1:end-1),Trials_to_detection(:,1:end-1),'LineWidth',widths.plot);
xlabel('Trials before switch')
ylabel('Switch detected')
FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);

print([exportpath 'abrupt_switch_affect_of_pre-trials'],'-dpng')
print([exportpath 'abrupt_switch_affect_of_pre-trials'],'-dsvg')


%% panel: effect of slope of change on detection speed [gradual switch, deterministic]
% error-shading for SD?

load([simOutputPath 'Results_for_slope_effects']);

ixPair = find(pairs(:,1) == 1 & pairs(:,2) == 0);
mean_Trials_to_detection = Results(ixPair).mean_Trial_since_crossover(1:end-1,:)'; % leave off gamma = 1; transpose to plot for slope (rows of matrix)
std_Trials_to_detection = Results(ixPair).std_Trial_since_crossover(1:end-1,:)';
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.square]); 
% plot(m,mean_Trials_to_detection)  
for iGamma = 1:numel(gamma)-1
    shadedErrorBar(m,mean_Trials_to_detection(iGamma,:),std_Trials_to_detection(iGamma,:),'lineprops',{'Color',colourmaps.gamma_lines(iGamma,:)})  
end
axis tight
xlabel('Probability change')
ylabel({'Switch detected'})
FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);

print([exportpath 'abrupt_switch_gradual_detection_speed'],'-dpng')
print([exportpath 'abrupt_switch_gradual_detection_speed'],'-dsvg')


%% panel: stability of switch
% error-shading for SD?

mean_stability = 100*Results(ixPair).mean_stability(1:end-1,:)'; % leave off gamma = 1; transpose to plot for slope (rows of matrix)
std_stability = 100*Results(ixPair).std_stability(1:end-1,:)';
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.square]); 
% plot(m,mean_Trials_to_detection)  
for iGamma = 1:numel(gamma)-1
    shadedErrorBar(m,mean_stability(iGamma,:),std_stability(iGamma,:),'lineprops',{'Color',colourmaps.gamma_lines(iGamma,:)})  
end
axis tight
xlabel('Probability change')
ylabel({'Stability (%)'})
FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);

print([exportpath 'abrupt_switch_gradual_detection_stability'],'-dpng')
print([exportpath 'abrupt_switch_gradual_detection_stability'],'-dsvg')
