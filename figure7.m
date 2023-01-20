clearvars;
close all

gamma_rule = 0.9; % forgetting rate for this analysis

load('Processed_data/Trial_Struct_Ukey_2018-06-05.mat')
load('Processed_data\MonkeyDataStrategies.mat')
figpath = 'Figures\';

% rename strategies for this task
strategies_label = {'Go Right','Go High-Prob','Go Left','Go Low-Prob','Win-Stay-Spatial', ...
    'Lose-Shift-Spatial','Repeat-High-Prob','Shift-To-High-Prob', 'Alternate','Sticky',...
    'Repeat-Large-Reward','Shift-After-Small-Reward'};
fields = fieldnames(Output);
nstrategy = numel(fields);

% visualization parameters
cmapstrategies = [brewermap(4,'Set2'); brewermap(8,'Set1')];
fontsize = 7;
axlinewidth = 0.5;

% define the side with higher probability
rewPright = Trial_Struct.rewProbRight;
rewPleft = Trial_Struct.rewProbLeft;
[maxP, maxSide]=max([rewPright'; rewPleft']);
maxSide(maxSide==2)=0;

% param for alignining around probability change
nbefore = 18;
nafter = 10;
x = -nbefore:nafter;
% find index of change in probability
idx_change = find(diff(Trial_Struct.rewProbRight)~=0)+1;
% align Rule strategies to probability change point
rulestrategy_probchange = cell(nstrategy,1);
for ic = 1:length(idx_change)-1
    for st = 1:nstrategy
        if sum(contains(fieldnames(Output.(fields{st})),'MAPprob_interpolated'))
            rulestrategy_probchange{st} = [rulestrategy_probchange{st} Output.(fields{st}).MAPprob_interpolated(idx_change(ic)-nbefore:idx_change(ic)+nafter)];
        else
            rulestrategy_probchange{st} = [rulestrategy_probchange{st} Output.(fields{st}).MAPprobability(idx_change(ic)-nbefore:idx_change(ic)+nafter)];
        end
    end
end


%% plot Rule Strategies
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[5 5 10 3.5])
for str = 1:3
    plot(Output.(fields{str}).MAPprobability,'Color',cmapstrategies(str,:)); hold on
%     text(length(Sbj.Rule.MAPts(:,str))+1,Sbj.Rule.MAPts(end,str),...
%         nameRstrategy{str},'Color',cmapRule(str,:),'FontSize',fontsize)
end
new_maxSide=maxSide;
new_maxSide(maxSide==0)=.49;
% plot right vs left probability with different color reflecting right left
% strategy
idx_right = find(maxSide==1);
idx_left = find(maxSide==0);
plot(idx_right,maxSide(idx_right),'.','Color',cmapstrategies(1,:)); hold on
plot(idx_left,maxSide(idx_left),'.','Color',cmapstrategies(3,:)); hold on
plot([idx_change idx_change],[0 1],'--','Color',[.6 .6 .6]); hold on
% plot(1:Trial_Struct.nTrials,maxSide,'k.'); hold on
plot([0 length(Output.(fields{str}).MAPprobability)],[0.5 0.5],'--','Color',[.6 .6 .6]); hold on
xlim([0 Trial_Struct.nTrials])
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
ylabel('Probability')
xlabel('Trial')
% ylim([0.45 1])

set(gcf,'PaperPositionMode','auto')
print([figpath 'fig7_monkey_rulestrategy'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'fig7_monkey_rulestrategy.png'])


%% plot Rule strategies around probability change trial
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[5 5 4 3.5])
for st = 1:3 
    % add shaded area for 1 SD
    y1 = mean(rulestrategy_probchange{st},2,'omitnan') + std(rulestrategy_probchange{st},1,2,'omitnan')./sqrt(size(rulestrategy_probchange{st},2));
    y2 = mean(rulestrategy_probchange{st},2,'omitnan') - std(rulestrategy_probchange{st},1,2,'omitnan')./sqrt(size(rulestrategy_probchange{st},2));
    h = patch([x fliplr(x)], [y1' fliplr(y2')],cmapstrategies(st,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x,mean(rulestrategy_probchange{st},2,'omitnan'),'Color',cmapstrategies(st,:)); hold on
    text(-9,0.3-0.1*(st-1),strategies_label{st},'Color',cmapstrategies(st,:),...
        'FontSize',fontsize); hold on
end
plot([-nbefore nafter],[0.5 0.5],'--','Color',[.6 .6 .6]); hold on
plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
ylabel('Probability')
xlabel('Trial')

set(gcf,'PaperPositionMode','auto')
print([figpath 'fig7_monkey_rulestrategy_rulechange'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'fig7_monkey_rulestrategy_rulechange.png'])


%% plot the proportion of dominant trials for win-stay strategies
ws = [7 11]; % win-stay strategy to test
ls = [8 12]; % lose-shift strategy to test
% find max(MAP) and max(Precision) to identify dominant strategy
matrixMAP = [Output.(fields{ws(1)}).MAPprob_interpolated Output.(fields{ws(2)}).MAPprob_interpolated];
maxval = max(matrixMAP,[],2);
[map_row,map_col] = find([Output.(fields{ws(1)}).MAPprob_interpolated Output.(fields{ws(2)}).MAPprob_interpolated] == maxval);

matrixPrecision = [Output.(fields{ws(1)}).precision_interpolated Output.(fields{ws(2)}).precision_interpolated];% sbj(2).Explore.Precisionts];
maxvalue = max(matrixPrecision,[],2);
[pre_row,pre_col] = find(matrixPrecision == maxvalue);

%% Build a matrix where max(MAP) and max(Precision) match. If so 
% strategy is recorded as 1, otherwise zeros

[~,idx_ws] = max(matrixMAP,[],2);

uni_val = unique(idx_ws);
prop_dominant_str = zeros(length(uni_val),1);
% estimate binomial confidence interval
N = 100;
for iu = 1:length(uni_val)
    prop_dominant_str(iu) = sum(idx_ws==uni_val(iu))./Trial_Struct.nTrials; 
    [~,pci(iu,:)] = binofit(prop_dominant_str(iu),N);

end
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[5 5 4 3.5])
b = bar(prop_dominant_str,'EdgeColor','none'); hold on
b.FaceColor = 'flat';
b.CData(1,:) = cmapstrategies(ws(1),:);
b.CData(2,:) = cmapstrategies(ws(2),:);
plot([1 1],pci(1,:)+prop_dominant_str(1),'-','Color',cmapstrategies(ws(1),:),'LineWidth',1.5); hold on
plot([2 2],pci(2,:)+prop_dominant_str(2),'-','Color',cmapstrategies(ws(2),:),'LineWidth',1.5); hold on
text(1,max(pci(1,:)+prop_dominant_str(1))+.05,strategies_label{ws(1)},'Color',...
    cmapstrategies(ws(1),:),'FontSize',fontsize); hold on
text(1,max(pci(2,:)+prop_dominant_str(2))+.05,strategies_label{ws(2)},'Color',...
    cmapstrategies(ws(2),:),'FontSize',fontsize); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
ylabel('Prop of trials')
set(gca,'xticklabel',[])
title('Dominant Probability')
ylim([0 0.8])

set(gcf,'PaperPositionMode','auto')
print([figpath 'fig7_monkey_propTrialWS'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'fig7_monkey_propTrialWS.png'])


%% plot the proportion of dominant trials for lose-shift strategies
matrixMAP_ls = [Output.(fields{ls(1)}).MAPprob_interpolated Output.(fields{ls(2)}).MAPprob_interpolated];
[~,idx_ls] = max(matrixMAP_ls,[],2);
uni_val = unique(idx_ls);
prop_dominant_str = zeros(length(uni_val),1);
% estimate binomial confidence interval
N = 100;
for iu = 1:length(uni_val)
    prop_dominant_str(iu) = sum(idx_ls==uni_val(iu))./Trial_Struct.nTrials; 
    [~,pci(iu,:)] = binofit(prop_dominant_str(iu),N);

end
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[5 5 4 3.5])
b = bar(prop_dominant_str,'EdgeColor','none'); hold on
b.FaceColor = 'flat';
b.CData(1,:) = cmapstrategies(ls(1),:);
b.CData(2,:) = cmapstrategies(ls(2),:);
plot([1 1],pci(1,:)+prop_dominant_str(1),'-','Color',cmapstrategies(ls(1),:),'LineWidth',1.5); hold on
plot([2 2],pci(2,:)+prop_dominant_str(2),'-','Color',cmapstrategies(ls(2),:),'LineWidth',1.5); hold on
text(1,max(pci(1,:)+prop_dominant_str(1))+.05,strategies_label{ls(1)},'Color',...
    cmapstrategies(ls(1),:),'FontSize',fontsize); hold on
text(1,max(pci(2,:)+prop_dominant_str(2))+.05,strategies_label{ls(2)},'Color',...
    cmapstrategies(ls(2),:),'FontSize',fontsize); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
ylabel('Prop of trials')
title('Dominant Probability')
set(gca,'xticklabel',[])

set(gcf,'PaperPositionMode','auto')
print([figpath 'fig7_monkey_propTrialLS'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'fig7_monkey_propTrialLS.png'])

%% plot win-stay and lose-shift High-prob vs High-Reward
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[5 5 10 3.5*2])
subplot(2,1,1)
for str = ws 
    plot(Output.(fields{str}).MAPprob_interpolated,'Color',cmapstrategies(str,:)); hold on
    text(length(Output.(fields{str}).MAPprob_interpolated)+1,Output.(fields{str}).MAPprob_interpolated(end),...
        strategies_label{str},'Color',cmapstrategies(str,:),'FontSize',fontsize)
end
plot(idx_right,maxSide(idx_right),'.','Color',cmapstrategies(1,:)); hold on
plot(idx_left,maxSide(idx_left),'.','Color',cmapstrategies(3,:)); hold on
plot([idx_change idx_change],[0 1],'--','Color',[.6 .6 .6]); hold on
% plot(1:Trial_Struct.nTrials,maxSide,'k.'); hold on
plot([0 Trial_Struct.nTrials],[0.5 0.5],'--','Color',[.6 .6 .6]); hold on
xlim([0 Trial_Struct.nTrials])
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
ylabel('Probability')

subplot(2,1,2)
for str = ls
    plot(Output.(fields{str}).MAPprob_interpolated,'Color',cmapstrategies(str,:)); hold on
    text(Trial_Struct.nTrials+1,Output.(fields{str}).MAPprob_interpolated(end),...
        strategies_label{str},'Color',cmapstrategies(str,:),'FontSize',fontsize)
end   
plot(idx_right,maxSide(idx_right),'.','Color',cmapstrategies(1,:)); hold on
plot(idx_left,maxSide(idx_left),'.','Color',cmapstrategies(3,:)); hold on
plot([idx_change idx_change],[0 1],'--','Color',[.6 .6 .6]); hold on
% plot(1:Trial_Struct.nTrials,maxSide,'k.'); hold on
plot([0 Trial_Struct.nTrials],[0.5 0.5],'--','Color',[.6 .6 .6]); hold on
xlim([0 Trial_Struct.nTrials])
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
ylabel('Probability')
xlabel('Trial')
% ylim([0.45 1])

set(gcf,'PaperPositionMode','auto')
print([figpath 'fig7_monkey_WS_LS'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'fig7_monkey_WS_LS.png'])
