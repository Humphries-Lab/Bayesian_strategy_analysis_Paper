clearvars; close all;

% strategies to assess
strategies = ["go_right", "go_cued", "go_left", "go_uncued",...
                "win_stay_spatial","lose_shift_spatial","win_stay_cued","lose_shift_cued",...
                   "alternate","sticky"];
nstrategy = numel(strategies);
strategies_label = {'Go Right','Go Cued','Go Left','Go Uncued','Win-Stay-Spatial',...
    'Lose-Shift-Spatial','Win-Stay-Cued','Lose-Shift-Cued', 'Alternate','Sticky'};

fontsize = 7;
axlinewidth = 0.5;
figpath = 'Figures\';

% load strategy profiles data for synthetic data
load('Processed_data/SyntheticDataStrategies_seed_VaryingDecay.mat')
load('Processed_data/SyntheticData_seed.mat') % Data contain 3 vector [light choice reward];

fields = fieldnames(Output{1});
strategy = {'Go Right','Alternate','Lose-Shift-Cue','Go Cued','Lose-Shift-Choice'};

% number of trials for each strategy
ntrial = 500;
nstr = ntrial/length(strategy);

barsize = [5 5 20 3.5*4];

indgm = find(decay_rate == .9);

%% plot the rule strategies for every gamma
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize)
for st = 1:nstrategy  
    heatmap = [];
    for dec = 1 : length(decay_rate)
        if sum(contains(fieldnames(Output{indgm}.(fields{st})),'MAPprob_interpolated'))
            heatmap = [heatmap; Output{dec}.(fields{st}).MAPprob_interpolated'];
        else
            heatmap = [heatmap; Output{dec}.(fields{st}).MAPprobability'];
        end
    end
    subplot(5,2,st)
    imagesc(heatmap); hold on
    colorbar()
    plot([0:nstr:ntrial; 0:nstr:ntrial],[0 length(decay_rate)],'k--','LineWidth', 1.2); hold on
    text(ntrial+2,indgm,'\leftarrow')
    set(gca,'FontName','Helvetica','FontSize',fontsize);
    set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth); 
    set(gca,'YTick',1:10:length(decay_rate))
    set(gca,'YTickLabel',decay_rate(1:10:end))
    title([strategies_label{st}])
%     title('Synthetic data')
%     if st == 1;
    ylabel('Decay rate (\gamma)')
%     end
    if (st == nstrategy-1) || (st == nstrategy); xlabel('Trial'); end
end
print([figpath 'figS3_MAP_Heatmap'],'-depsc')
saveas(gcf,[figpath, 'figS3_MAP_Heatmap.png'])

