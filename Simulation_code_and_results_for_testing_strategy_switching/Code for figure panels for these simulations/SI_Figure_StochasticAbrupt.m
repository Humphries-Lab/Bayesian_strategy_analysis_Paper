%%% script to plot results for the stochastic, abrupt model in SI Figure
clearvars; close all

run figure_properties.m

if ispc
    simOutputPath = 'C:\Users\lpzmdh\Dropbox\Projects\Bayesian strategy analysis\Analytical_changes_of_distribution\';
else
    simOutputPath = '/Users/mqbssmhg/Dropbox/Projects/Bayesian strategy analysis/Analytical_changes_of_distribution/';
end

%% panel: effect of trials-before-switch on detection speed [abrupt switch, deterministic]

load([simOutputPath 'Results_AbruptStochastic_variable_t']);

ixPair = find(pairs(:,1) == 0.9 & pairs(:,2) == 0.1);
Trials_to_detection = mean(Results(ixPair).number_of_trials,3);

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.square]); 
set(gca,'ColorOrder',colourmaps.gamma_lines,'NextPlot', 'replacechildren');
plot(switching_trial(1:end-1),Trials_to_detection(:,1:end-1),'LineWidth',widths.plot);
xlabel('Trials before switch')
ylabel('Switch detected')
FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);

print([exportpath 'abrupt_switch_affect_of_pre-trials_pt9_pt1'],'-dpng')
print([exportpath 'abrupt_switch_affect_of_pre-trials_pt9_pt1'],'-dsvg')


ixPair = find(pairs(:,1) == 0.6 & pairs(:,2) == 0.5);
Trials_to_detection = mean(Results(ixPair).number_of_trials,3);

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.square]); 
set(gca,'ColorOrder',colourmaps.gamma_lines,'NextPlot', 'replacechildren');
plot(switching_trial(1:end-1),Trials_to_detection(:,1:end-1),'LineWidth',widths.plot);
xlabel('Trials before switch')
ylabel('Switch detected')
FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);

print([exportpath 'abrupt_switch_affect_of_pre-trials_pt6_pt5'],'-dpng')
print([exportpath 'abrupt_switch_affect_of_pre-trials_pt6_pt5'],'-dsvg')

%% panels: mean & STD times per gamma
lower_bounds = [-0.1, 0.5];  % to set limits on x and y axis for pseudo-colour plots

% find maxima and minima to scale all colormaps to the same range
set_of_means = Summary.Pairs_Mean(:,:,1:end-1); % omit gamma = 1
set_of_STD = Summary.Pairs_STD(:,:,1:end-1);
[min_mean,max_mean] = bounds(set_of_means(:));
[min_std,max_std] = bounds(set_of_STD(:));

for iGamma = 1:numel(gamma)-1 % omit no-decay: plot separately
    % MEANS
    figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.smallmatrices]); 
    colormap(colourmaps.trials)
    plotMatrix(p_lower,p_higher,squeeze(Summary.Pairs_Mean(:,:,iGamma)),lower_bounds);
    set(gca,'CLim',[min_mean, max_mean]);  % scale colormap the same in all figures
    axis square
    %xlabel('Probability 2')
    %ylabel('Probability 1')
    FormatFig_For_Export(gcf,fontsize-2,fontname,widths.axis);

    print([exportpath 'abrupt_switch_mean_per_probability_gamma_' num2str(gamma(iGamma)) '.png'],'-dpng')
    print([exportpath 'abrupt_switch_mean_per_probability_gamma_' num2str(gamma(iGamma)) '.svg'],'-dsvg')

    % STD
    figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.smallmatrices ]); 
    colormap(colourmaps.trials)
    plotMatrix(p_lower,p_higher,squeeze(Summary.Pairs_STD(:,:,iGamma)),lower_bounds);
    set(gca,'CLim',[min_std, max_std]);  % scale colormap the same in all figures
    axis square
    %xlabel('Probability 2')
    %ylabel('Probability 1')
    FormatFig_For_Export(gcf,fontsize-2,fontname,widths.axis);

    print([exportpath 'abrupt_switch_STD_per_probability_gamma_' num2str(gamma(iGamma)) '.png'],'-dpng')
    print([exportpath 'abrupt_switch_STD_per_probability_gamma_' num2str(gamma(iGamma)) '.svg'],'-dsvg')
   
end
% separate colorbars
% MEAN
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.smallmatrices ]);
colormap(colourmaps.trials)
set(gca,'CLim',[min_mean, max_mean]);  % scale colormap the same in all figures
hc = colorbar('location','north');
set(gca,'Visible','off')
pos = hc.Position;
hc.Position = [pos(1) pos(2) pos(3) figsize.colourbar_height];
FormatFig_For_Export(gcf,fontsize,fontname,widths.axis); 
print([exportpath 'colorbar_abrupt_switch_mean_per_probability'],'-dpng')
print([exportpath 'colorbar_abrupt_switch_mean_per_probability'],'-dsvg')

% separate colorbars
% STD
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.smallmatrices ]);
colormap(colourmaps.trials)
set(gca,'CLim',[min_std, max_std]);  % scale colormap the same in all figures
hc = colorbar('location','north');
set(gca,'Visible','off')
pos = hc.Position;
hc.Position = [pos(1) pos(2) pos(3) figsize.colourbar_height];
FormatFig_For_Export(gcf,fontsize,fontname,widths.axis); 
print([exportpath 'colorbar_abrupt_switch_STD_per_probability'],'-dpng')
print([exportpath 'colorbar_abrupt_switch_STD_per_probability'],'-dsvg')


%% panels: mean and STD stability per gamma

set_of_means = Summary.Pairs_Mean_stability(:);
set_of_STD = Summary.Pairs_STD_stability(:);
[min_mean,max_mean] = bounds(set_of_means(:));
[min_std,max_std] = bounds(set_of_STD(:));

for iGamma = 1:numel(gamma)
    % Means
    figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.smallmatrices ]); 
    colormap(colourmaps.stability)
    plotMatrix(p_lower,p_higher,squeeze(Summary.Pairs_Mean_stability(:,:,iGamma)),lower_bounds);
    set(gca,'CLim',[min_mean, max_mean]);  % scale colormap the same in all figures
    axis square
   % xlabel('Initial probability of strategy 2')
    %ylabel('Initial probability of strategy 1')
    FormatFig_For_Export(gcf,fontsize-2,fontname,widths.axis);

    print([exportpath 'abrupt_switch_mean_Stability_per_probability_gamma_' num2str(gamma(iGamma)) '.png'],'-dpng')
    print([exportpath 'abrupt_switch_mean_Stability_per_probability_gamma_' num2str(gamma(iGamma)) '.svg'],'-dsvg')
    
    % STD
    figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.smallmatrices ]); 
    colormap(colourmaps.stability)
    plotMatrix(p_lower,p_higher,squeeze(Summary.Pairs_STD_stability(:,:,iGamma)),lower_bounds);
    set(gca,'CLim',[min_std, max_std]);  % scale colormap the same in all figures
    axis square
    %xlabel('Initial probability of strategy 2')
    %ylabel('Initial probability of strategy 1')
    FormatFig_For_Export(gcf,fontsize-2,fontname,widths.axis);

    print([exportpath 'abrupt_switch_STD_Stability_per_probability_gamma_' num2str(gamma(iGamma)) '.png'],'-dpng')
    print([exportpath 'abrupt_switch_STD_Stability_per_probability_gamma_' num2str(gamma(iGamma)) '.svg'],'-dsvg')
end

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.smallmatrices ]);
colormap(colourmaps.stability)
set(gca,'CLim',[min_mean, max_mean]);  % scale colormap the same in all figures
hc = colorbar('location','north');
set(gca,'Visible','off')
pos = hc.Position;
hc.Position = [pos(1) pos(2) pos(3) figsize.colourbar_height];
FormatFig_For_Export(gcf,fontsize,fontname,widths.axis); 
print([exportpath 'colorbar_abrupt_switch_mean_Stability_per_probability'],'-dpng')
print([exportpath 'colorbar_abrupt_switch_mean_Stability_per_probability'],'-dsvg')

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.smallmatrices ]);
colormap(colourmaps.stability)
set(gca,'CLim',[min_std, max_std]);  % scale colormap the same in all figures
hc = colorbar('location','north');
set(gca,'Visible','off')
pos = hc.Position;
hc.Position = [pos(1) pos(2) pos(3) figsize.colourbar_height];
FormatFig_For_Export(gcf,fontsize,fontname,widths.axis); 
print([exportpath 'colorbar_abrupt_switch_STD_Stability_per_probability'],'-dpng')
print([exportpath 'colorbar_abrupt_switch_STD_Stability_per_probability'],'-dsvg')

