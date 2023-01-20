clearvars; 
close all
%% load Adrien data

% load strategy profiles data for 4 rats from Peyarache et al., 2009
% dataset
load('Processed_data/PeyracheDataStrategies.mat')
figpath = 'Figures\';

% Header = {'Animal', 'Session', 'Direction', 'Reward', 'Rule', 'Light', 'SesName', 'Learning'};
load('Processed_data\SummaryDataTable_AllSessions.mat');
Data = double(Data);
Data(Data(:,4)==-1,4)=0;

strategies_label = {'Go Right','Go Cued','Go Left','Go Uncued','Win-Stay-Spatial', ...
    'Lose-Shift-Spatial','Win-Stay-Cued','Lose-Shift-Cued', 'Alternate','Sticky'};

nRats = length(Output);  % number of rats
fields = fieldnames(Output{1});
nstrategy = numel(fields);

% visualise choices
cmapStrategy = [brewermap(4,'Set2'); brewermap(6,'Dark2')];
cmapRchange = brewermap(2,'Set1');
fontsize = 7;
axlinewidth = 0.5;
barsize = [5 5 10 3.5];

nbefore = 22;
nafter = 30;

% collect MAPs values around rule switch trial. Both Ego-->Cue and
% Cue-->Ego
ego_cue = cell(nstrategy,1); % first cell column for shift right-->cue
cue_ego = cell(nstrategy,1);
record_ruleChange = [];
 
for iR = 1:nRats
%     rules = unique(Rat(iR).rules);
    idx_rat = find(Data(:,1)==iR);
    rchange = find(diff(Data(idx_rat,5))~=0)+1;
    if iR == 1; rchange(1:2)=[]; end
    rulesRC = Data(idx_rat(rchange),5);
    record_ruleChange = [record_ruleChange; rulesRC];

    if intersect([2 4],rulesRC)
        nrc = rchange((rulesRC==2 | rulesRC==4));
        for str = 1: nstrategy
            for r = 1:length(nrc)
                if sum(contains(fieldnames(Output{iR}.(fields{str})),'MAPprob_interpolated'))
                    ego_cue{str} = [ego_cue{str} Output{iR}.(fields{str}).MAPprob_interpolated(nrc(r)-nbefore:nrc(r)+nafter)];                
                else
                    ego_cue{str} = [ego_cue{str} Output{iR}.(fields{str}).MAPprobability(nrc(r)-nbefore:nrc(r)+nafter)];                
                end
            end
        end
    end
    
    if intersect([1 3],rulesRC)
        nrc = rchange((rulesRC==1 | rulesRC==3));
        for str = 1: nstrategy
            for r = 1:length(nrc)
                if sum(contains(fieldnames(Output{iR}.(fields{str})),'MAPprob_interpolated'))
                    cue_ego{str} = [cue_ego{str} Output{iR}.(fields{str}).MAPprob_interpolated(nrc(r)-nbefore:nrc(r)+nafter)];
                else
                    cue_ego{str} = [cue_ego{str} Output{iR}.(fields{str}).MAPprobability(nrc(r)-nbefore:nrc(r)+nafter)];
                end
            end
        end
    end
end


%% plot strategy profiles around rule change point for Adrien data. 
% First subplot show the succesfull strategy for the ego rule in the shift from Ego --> Cue
% and the shift from Cue --> Ego
str_r = find(contains(strategies_label,'Go Right'));
str_l = find(contains(strategies_label,'Go Left'));
% Second subplot show the succesful strategy for cue rule in the shift from Cue --> Ego
% and Ego-->Cue
str_c = find(contains(strategies_label,'Go Cued'));
% str_c = find(contains(strRuleStrategy,'Go Cue'));

cuer = record_ruleChange(record_ruleChange==2 | record_ruleChange==4);
right_rat = find(cuer==2);
left_rat = find(cuer==4);

egor = record_ruleChange(record_ruleChange==1 | record_ruleChange==3);
cueright_rat = find(egor==1);
cueleft_rat = find(egor==3);

successE_egocue = [ego_cue{str_r}(:,right_rat) ego_cue{str_l}(:,left_rat)];
successE_cueego = [cue_ego{str_l}(:,cueleft_rat)];

successC_egocue = [ego_cue{str_c}(:,right_rat) ego_cue{str_c}(:,left_rat)];
successC_cueego = [cue_ego{str_c}(:,cueleft_rat)];


x = -nbefore+1:nafter+1;
 
%%****************************************
fig1 = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize);
subplot(121)
curve1 = nanmean(successE_egocue,2)-nanstd(successE_egocue,1,2)./sqrt(nRats);
curve2 = nanmean(successE_egocue,2)+nanstd(successE_egocue,1,2)./sqrt(nRats);
h = fill([x, fliplr(x)], [curve1', fliplr(curve2')], cmapRchange(1,:)); hold on
set(h,'facealpha',.3,'LineStyle','none')
plot(x,nanmean(successE_egocue,2),'-','Color',cmapRchange(1,:)); hold on
text(x(end)-10, curve2(end)+.1,'go spatial','Color',cmapRchange(1,:),'FontSize',fontsize)

curve1 = nanmean(successC_egocue,2)-nanstd(successC_egocue,1,2)./sqrt(nRats);
curve2 = nanmean(successC_egocue,2)+nanstd(successC_egocue,1,2)./sqrt(nRats);
h = fill([x, fliplr(x)], [curve1', fliplr(curve2')], cmapRchange(2,:)); hold on
set(h,'facealpha',.3,'LineStyle','none')
plot(x,nanmean(successC_egocue,2),'-','Color',cmapRchange(2,:)); hold on
text(x(end)-10, curve1(end)-.1,'go cued','Color',cmapRchange(2,:),'FontSize',fontsize)

plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([x(1) x(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
title('Spatial \rightarrow Cued');
xlabel('Trial')
ylabel('Probability')
ylim([0.3 1])

subplot(122)
curve1 = nanmean(successE_cueego,2)-nanstd(successE_cueego,1,2)./sqrt(nRats);
curve2 = nanmean(successE_cueego,2)+nanstd(successE_cueego,1,2)./sqrt(nRats);
h = fill([x, fliplr(x)], [curve1', fliplr(curve2')], cmapRchange(1,:)); hold on
set(h,'facealpha',.3,'LineStyle','none')
plot(x,nanmean(successE_cueego,2),'-','Color',cmapRchange(1,:)); hold on
% text(x(end)-10, curve2(end)+.1,'go ego','Color',cmapRchange(1,:),'FontSize',fontsize)

curve1 = nanmean(successC_cueego,2)-nanstd(successC_cueego,1,2)./sqrt(nRats);
curve2 = nanmean(successC_cueego,2)+nanstd(successC_cueego,1,2)./sqrt(nRats);
h = fill([x, fliplr(x)], [curve1', fliplr(curve2')], cmapRchange(2,:)); hold on
set(h,'facealpha',.3,'LineStyle','none')
plot(x,nanmean(successC_cueego,2),'-','Color',cmapRchange(2,:)); hold on

plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([x(1) x(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
title('Cued \rightarrow Spatial')
xlabel('Trial')
ylabel('Probability')
ylim([0.3 1])

set(gcf,'PaperPositionMode','auto')
print([figpath 'fig4_YmazeRuleStr_rulechange'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'fig4_YmazeRuleStr_rulechange.png'])


%% visualize individual variability
ymaze_E_egocue = successE_egocue(nbefore:end,:);
successE_egocue = [];
successE_egocue = ymaze_E_egocue;

ymaze_C_egocue = successC_egocue(nbefore:end,:);
successC_egocue = [];
successC_egocue = ymaze_C_egocue;

ymaze_E_cueego = successE_cueego(nbefore:end,:);
successE_cueego = [];
successE_cueego = ymaze_E_cueego;

ymaze_C_cueego = successC_cueego(nbefore:end,:);
successC_cueego = [];
successC_cueego = ymaze_C_cueego;

len_egocue = size(successE_egocue,2);
jit = randn(len_egocue,1)*.1; jit2 = randn(size(successE_cueego,2),1)*.1;
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[5 5 7 3.5]);
subplot(1,2,1)
plot(1+jit,sum(successE_egocue>.5)./length(successE_egocue),'wo','MarkerSize',4,...
    'MarkerFaceColor',cmapRchange(1,:)); hold on
plot(2+jit,sum(successC_egocue>.5)./length(successC_egocue),'wo','MarkerSize',4,...
    'MarkerFaceColor',cmapRchange(2,:)); hold on
plot([1+jit'; 2+jit'],[(sum(successE_egocue>.5)./length(successE_egocue)); ...
    (sum(successC_egocue>.5)./length(successC_egocue))],'-','Color',[.7 .7 .7]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
plot([.5 2.5], [.5 .5],'--','Color',[.6 .6 .6]); hold on
xlim([.5 2.5])
xticks(1:2)
ylim([0 1])
xticklabels({'go spatial','go cued'})
xtickangle(45)
title('Spatial \rightarrow Cued')
ylabel('Prop of trial >0.5')

subplot(1,2,2)
plot(2+jit2,sum(successE_cueego>.5)./length(successE_cueego),'wo','MarkerSize',4,...
    'MarkerFaceColor',cmapRchange(1,:)); hold on
plot(1+jit2,sum(successC_cueego>.5)./length(successC_cueego),'wo','MarkerSize',4,...
    'MarkerFaceColor',cmapRchange(2,:)); hold on
plot([2+jit2'; 1+jit2'],[(sum(successE_cueego>.5)./length(successE_cueego)); ...
    (sum(successC_cueego>.5)./length(successC_cueego))],'-','Color',[.7 .7 .7]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
plot([.5 2.5], [.5 .5],'--','Color',[.6 .6 .6]); hold on
xlim([.5 2.5])
xticks(1:2)
ylim([0 1])
xticklabels({'go cued','go spatial'})
xtickangle(45)
title('Cued \rightarrow Spatial')
ylabel('Prop of trial >0.5')

set(gcf,'PaperPositionMode','auto')
print([figpath 'fig4_Ymaze_ProbTrial'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'fig4_Ymaze_ProbTrial.png'])

%% plot explore strategies around rule change trial
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize);
subplot(121)
for str = 5:nstrategy
    curve1 = nanmean(ego_cue{str},2)-nanstd(ego_cue{str},1,2)./sqrt(nRats);
    curve2 = nanmean(ego_cue{str},2)+nanstd(ego_cue{str},1,2)./sqrt(nRats);
    h = fill([x, fliplr(x)], [curve1', fliplr(curve2')],cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x,nanmean(ego_cue{str},2),'-','Color',cmapStrategy(str,:)); hold on
    plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
    plot([x(1) x(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
    set(gca,'FontName','Helvetica','FontSize',fontsize);
    set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
    y = nanmean(ego_cue{str},2);
    text(x(end)+1,y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
xlabel('Trial')
ylabel('Probability')
title('Spatial-->Cued')
ylim([0.3 1])
xlim([x(1) x(end)+20])

subplot(122)
for str = 5:nstrategy
    curve1 = nanmean(cue_ego{str},2)-nanstd(cue_ego{str},1,2)./sqrt(nRats);
    curve2 = nanmean(cue_ego{str},2)+nanstd(cue_ego{str},1,2)./sqrt(nRats);
    h = fill([x, fliplr(x)], [curve1', fliplr(curve2')], cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x,nanmean(cue_ego{str},2),'-','Color',cmapStrategy(str,:)); hold on
    plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
    plot([x(1) x(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
    set(gca,'FontName','Helvetica','FontSize',fontsize);
    set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
    y = nanmean(cue_ego{str},2);
    text(x(end)+1,y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
xlabel('Trial')
title('Cued-->Spatial')
ylim([0.3 1])
xlim([x(1) x(end)+20])

set(gcf,'PaperPositionMode','auto')
print([figpath 'fig4_SI_YmazeExploreStr_rulechange'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'fig4_SI_YmazeExploreStr_rulechange.png'])

%% plot explore strategies around rule change trial
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize);
subplot(121)
for str = 5:6
    curve1 = nanmean(ego_cue{str},2)-nanstd(ego_cue{str},1,2)./sqrt(nRats);
    curve2 = nanmean(ego_cue{str},2)+nanstd(ego_cue{str},1,2)./sqrt(nRats);
    h = fill([x, fliplr(x)], [curve1', fliplr(curve2')],cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x,nanmean(ego_cue{str},2),'-','Color',cmapStrategy(str,:)); hold on
    plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
    plot([x(1) x(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
    set(gca,'FontName','Helvetica','FontSize',fontsize);
    set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
    y = nanmean(ego_cue{str},2);
    text(x(end)+1,y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
xlabel('Trial')
ylabel('Probability')
title('Spatial-->Cued')
ylim([0.2 1])
xlim([x(1) x(end)+20])

subplot(122)
for str = 5:6 %nExploreStrategy
    curve1 = nanmean(cue_ego{str},2)-nanstd(cue_ego{str},1,2)./sqrt(nRats);
    curve2 = nanmean(cue_ego{str},2)+nanstd(cue_ego{str},1,2)./sqrt(nRats);
    h = fill([x, fliplr(x)], [curve1', fliplr(curve2')], cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x,nanmean(cue_ego{str},2),'-','Color',cmapStrategy(str,:)); hold on
    plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
    plot([x(1) x(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
    set(gca,'FontName','Helvetica','FontSize',fontsize);
    set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
    y = nanmean(cue_ego{str},2);
    text(x(end)+1,y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
xlabel('Trial')
title('Cued-->Spatial')
ylim([0.2 1])
xlim([x(1) x(end)+20])

set(gcf,'PaperPositionMode','auto')
print([figpath 'fig5_YmazeWSLSchoice_rulechange'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'fig5_YmazeWSLSchoice_rulechange.png'])

%% plot explore strategies around rule change trial
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize);
subplot(121)
for str = 7:8
    curve1 = nanmean(ego_cue{str},2)-nanstd(ego_cue{str},1,2)./sqrt(nRats);
    curve2 = nanmean(ego_cue{str},2)+nanstd(ego_cue{str},1,2)./sqrt(nRats);
    h = fill([x, fliplr(x)], [curve1', fliplr(curve2')],cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x,nanmean(ego_cue{str},2),'-','Color',cmapStrategy(str,:)); hold on
    plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
    plot([x(1) x(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
    set(gca,'FontName','Helvetica','FontSize',fontsize);
    set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
    y = nanmean(ego_cue{str},2);
    text(x(end)+1,y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
xlabel('Trial')
ylabel('Probability')
title('Spatial-->Cued')
ylim([0.3 1])
xlim([x(1) x(end)+20])

subplot(122)
for str = 7:8 %nExploreStrategy
    curve1 = nanmean(cue_ego{str},2)-nanstd(cue_ego{str},1,2)./sqrt(nRats);
    curve2 = nanmean(cue_ego{str},2)+nanstd(cue_ego{str},1,2)./sqrt(nRats);
    h = fill([x, fliplr(x)], [curve1', fliplr(curve2')], cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x,nanmean(cue_ego{str},2),'-','Color',cmapStrategy(str,:)); hold on
    plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
    plot([x(1) x(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
    set(gca,'FontName','Helvetica','FontSize',fontsize);
    set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
    y = nanmean(cue_ego{str},2);
    text(x(end)+1,y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
xlabel('Trial')
title('Cued-->Spatial')
ylim([0.3 1])
xlim([x(1) x(end)+20])

set(gcf,'PaperPositionMode','auto')
print([figpath 'fig5_YmazeWSLScue_rulechange'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'Fig5_YmazeWSLScue_rulechange.png'])

%% load Tobias data (only control animal)
load('Processed_data/LeverPressDataStrategies.mat')
fields = fieldnames(Output{1});
 
% this dataset contain as header = {'Animal ID', 'Ses ID', 'Dir', 'Reward', ...
% 'Rule', 'Light', 'Date', 'Learn tr', 'Group', 'Reaction Time'};
load('Processed_data/SummaryDataTable_AllSessions_LeverPress.mat')

% for this task the rules are numbered as follow:
% right = 1 
% left = 2
% cue = 3
nSbj = length(Output);

nbefore = 30;
nafter = 70;
x = -nbefore+1:nafter+1;
rc = NaN(nSbj,3);
%% record the rule change trial and the rules before and after change for 
% plot of MAP values around rule change 
sbj_name = unique(Datas(:,1));
for iS = 1:nSbj
    idx_sbj = find(Datas(:,1)==sbj_name(iS));
    rchange = find(diff(Datas(idx_sbj,5))~=0)+1;
    rulesRC = Datas(idx_sbj(rchange),5);
    % record rule change trial and the rule before and after shift for each subject
    if ~isempty(rchange)
        rc(iS,:) = [rchange Datas(idx_sbj(rchange-1),5) rulesRC];
    end
    if length(rulesRC)>1; disp('Multiple rule changes'); return; end
end

rulechange_str = cell(nstrategy,1);
for str = 1: nstrategy
    rulechange_str{str} = NaN(nafter+nbefore+1,nSbj);
    for iS = 1:nSbj
        if ~isnan(rc(iS,1))
            if sum(contains(fieldnames(Output{iS}.(fields{str})),'MAPprob_interpolated'))
                rulechange_str{str}(:,iS) = Output{iS}.(fields{str}).MAPprob_interpolated(rc(iS,1)-nbefore:rc(iS,1)+nafter);
            else
                rulechange_str{str}(:,iS) = Output{iS}.(fields{str}).MAPprobability(rc(iS,1)-nbefore:rc(iS,1)+nafter);
            end
        end
    end
end


%% plot strategy profiles around rule change point for Tobias data. 
% First subplot show the succesfull strategy for the ego rule in the shift from Ego --> Cue
% and the shift from Cue --> Ego

% identify subject with shift from right or left to cue
[right_sbj, ~] = find(ismember(rc(:,2:3),[1 3],'rows'));
[left_sbj,~] = find(ismember(rc(:,2:3),[2 3],'rows'));

% identify subject with shift from cue to ego right or left
[cueright_sbj, ~] = find(ismember(rc(:,2:3),[3 1],'rows'));
[cueleft_sbj,~] = find(ismember(rc(:,2:3),[3 2],'rows'));

% collect the MAPs around rule change for succesful ego strategy for each subject
successE_egocue = [rulechange_str{str_r}(:,right_sbj) rulechange_str{str_l}(:,left_sbj)]; 
successE_cueego = [rulechange_str{str_r}(:,cueright_sbj) rulechange_str{str_l}(:,cueleft_sbj)]; 

% collect the MAPs around rule change for succesful cue strategy for each subject
successC_egocue = [rulechange_str{str_c}(:,right_sbj) rulechange_str{str_c}(:,left_sbj)]; 
successC_cueego = [rulechange_str{str_c}(:,cueright_sbj) rulechange_str{str_c}(:,cueleft_sbj)]; 

% x = -nbefore:nafter;
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize);
subplot(121)
curve1 = nanmean(successE_egocue,2)-nanstd(successE_egocue,1,2)./sqrt(nRats);
curve2 = nanmean(successE_egocue,2)+nanstd(successE_egocue,1,2)./sqrt(nRats);
h = fill([x, fliplr(x)], [curve1', fliplr(curve2')], cmapRchange(1,:)); hold on
set(h,'facealpha',.3,'LineStyle','none')
plot(x,nanmean(successE_egocue,2),'-','Color',cmapRchange(1,:)); hold on
text(x(end)-20, curve2(end)+.1,'go spatial','Color',cmapRchange(1,:),'FontSize',fontsize)

curve1 = nanmean(successC_egocue,2)-nanstd(successC_egocue,1,2)./sqrt(nRats);
curve2 = nanmean(successC_egocue,2)+nanstd(successC_egocue,1,2)./sqrt(nRats);
h = fill([x, fliplr(x)], [curve1', fliplr(curve2')], cmapRchange(2,:)); hold on
set(h,'facealpha',.3,'LineStyle','none')
plot(x,nanmean(successC_egocue,2),'-','Color',cmapRchange(2,:)); hold on
text(x(end)-20, curve1(end)-.1,'go cued','Color',cmapRchange(2,:),'FontSize',fontsize)

plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([x(1) x(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
title('Spatial \rightarrow Cued');
xlabel('Trial')
ylabel('Probability')
ylim([0.3 1.1])
xlim([-nbefore nafter])

subplot(122)
curve1 = nanmean(successE_cueego,2)-nanstd(successE_cueego,1,2)./sqrt(nRats);
curve2 = nanmean(successE_cueego,2)+nanstd(successE_cueego,1,2)./sqrt(nRats);
h = fill([x, fliplr(x)], [curve1', fliplr(curve2')], cmapRchange(1,:)); hold on
set(h,'facealpha',.3,'LineStyle','none')
plot(x,nanmean(successE_cueego,2),'-','Color',cmapRchange(1,:)); hold on

curve1 = nanmean(successC_cueego,2)-nanstd(successC_cueego,1,2)./sqrt(nRats);
curve2 = nanmean(successC_cueego,2)+nanstd(successC_cueego,1,2)./sqrt(nRats);
h = fill([x, fliplr(x)], [curve1', fliplr(curve2')], cmapRchange(2,:)); hold on
set(h,'facealpha',.3,'LineStyle','none')
plot(x,nanmean(successC_cueego,2),'-','Color',cmapRchange(2,:)); hold on

plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([x(1) x(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
title('Cued \rightarrow Spatial')
xlabel('Trial')
ylabel('Probability')
ylim([0.3 1.1])
xlim([-nbefore nafter])

set(gcf,'PaperPositionMode','auto')
print([figpath 'fig4_LeverpressRuleStr_rulechange'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'fig4_LeverpressRuleStr_rulechange.png'])


%% visualize individual variability
ymaze_E_egocue = successE_egocue(nbefore:end,:);
successE_egocue = [];
successE_egocue = ymaze_E_egocue;

ymaze_C_egocue = successC_egocue(nbefore:end,:);
successC_egocue = [];
successC_egocue = ymaze_C_egocue;

ymaze_E_cueego = successE_cueego(nbefore:end,:);
successE_cueego = [];
successE_cueego = ymaze_E_cueego;

ymaze_C_cueego = successC_cueego(nbefore:end,:);
successC_cueego = [];
successC_cueego = ymaze_C_cueego;

len_egocue = size(successE_egocue,2);
jit = randn(len_egocue,1)*.1; jit2 = randn(size(successE_cueego,2),1)*.1;
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[5 5 7 3.5]);
subplot(1,2,1)
plot(1+jit,sum(successE_egocue>.5)./length(successE_egocue),'wo','MarkerSize',4,...
    'MarkerFaceColor',cmapRchange(1,:)); hold on
plot(2+jit,sum(successC_egocue>.5)./length(successC_egocue),'wo','MarkerSize',4,...
    'MarkerFaceColor',cmapRchange(2,:)); hold on
plot([1+jit'; 2+jit'],[(sum(successE_egocue>.5)./length(successE_egocue)); ...
    (sum(successC_egocue>.5)./length(successC_egocue))],'-','Color',[.7 .7 .7]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
plot([.5 2.5], [.5 .5],'--','Color',[.6 .6 .6]); hold on
xlim([.5 2.5])
xticks(1:2)
ylim([0.4 1])
xticklabels({'go spatial','go cued'})
xtickangle(45)
title('Spatial \rightarrow Cued')
ylabel('Prop of trial >0.5')

subplot(1,2,2)
plot(2+jit2,sum(successE_cueego>.5)./length(successE_cueego),'wo','MarkerSize',4,...
    'MarkerFaceColor',cmapRchange(1,:)); hold on
plot(1+jit2,sum(successC_cueego>.5)./length(successC_cueego),'wo','MarkerSize',4,...
    'MarkerFaceColor',cmapRchange(2,:)); hold on
plot([2+jit2'; 1+jit2'],[(sum(successE_cueego>.5)./length(successE_cueego)); ...
    (sum(successC_cueego>.5)./length(successC_cueego))],'-','Color',[.7 .7 .7]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
plot([.5 2.5], [.5 .5],'--','Color',[.6 .6 .6]); hold on
xlim([.5 2.5])
ylim([0.4 1])
xticks(1:2)
xticklabels({'go cued','go spatial'})
xtickangle(45)
title('Cued \rightarrow Spatial')
ylabel('Prop of trial >0.5')

set(gcf,'PaperPositionMode','auto')
print([figpath 'fig4_Leverpress_ProbTrial'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'fig4_Leverpress_ProbTrial.png'])

%% plot explore strategies around rule change trials. Both from Ego to Cue 
% and viceversa
nsb = length([left_sbj;right_sbj]);
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize);
subplot(121)
for str = 5:nstrategy
    curve1 = nanmean(rulechange_str{str}(:,[left_sbj;right_sbj]),2)-nanstd(rulechange_str{str}(:,[left_sbj;right_sbj]),1,2)./sqrt(nsb);
    curve2 = nanmean(rulechange_str{str}(:,[left_sbj;right_sbj]),2)+nanstd(rulechange_str{str}(:,[left_sbj;right_sbj]),1,2)./sqrt(nsb);
    h = fill([x, fliplr(x)], [curve1', fliplr(curve2')], cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x,nanmean(rulechange_str{str}(:,[left_sbj;right_sbj]),2),'-','Color',cmapStrategy(str,:)); hold on
    y = nanmean(rulechange_str{str}(:,[left_sbj;right_sbj]),2);
    text(x(end),y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([x(1) x(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
xlabel('Trial')
ylabel('Probability')
ylim([0.3 1])
title('Spatial-->Cued')
xlim([x(1) x(end)+20])

nsb = length([cueleft_sbj;cueright_sbj]);
subplot(122)
for str = 5:nstrategy
    curve1 = nanmean(rulechange_str{str}(:,[cueright_sbj; cueleft_sbj]),2)-nanstd(rulechange_str{str}(:,[cueright_sbj; cueleft_sbj]),1,2)./sqrt(nsb);
    curve2 = nanmean(rulechange_str{str}(:,[cueright_sbj; cueleft_sbj]),2)+nanstd(rulechange_str{str}(:,[cueright_sbj; cueleft_sbj]),1,2)./sqrt(nsb);
    h = fill([x, fliplr(x)], [curve1', fliplr(curve2')], cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x,nanmean(rulechange_str{str}(:,[cueright_sbj; cueleft_sbj]),2),'-','Color',cmapStrategy(str,:)); hold on
    y = nanmean(rulechange_str{str}(:,[cueright_sbj; cueleft_sbj]),2);
    text(x(end),y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([x(1) x(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
xlabel('Trial')
% ylabel('P(go cue)')
ylim([0.3 1])
title('Cued-->Spatial')
xlim([x(1) x(end)+20])

set(gcf,'PaperPositionMode','auto')
print([figpath 'fig4_SI_LeverpressExploreStr_rulechange'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'fig4_SI_LeverpressExploreStr_rulechange.png'])

%% plot explore strategies around rule change trials. Both from Ego to Cue 
% and viceversa
nsb = length([left_sbj;right_sbj]);
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize);
subplot(121)
for str = 5:6 %nExploreStrategy
    curve1 = nanmean(rulechange_str{str}(:,[left_sbj;right_sbj]),2)-nanstd(rulechange_str{str}(:,[left_sbj;right_sbj]),1,2)./sqrt(nsb);
    curve2 = nanmean(rulechange_str{str}(:,[left_sbj;right_sbj]),2)+nanstd(rulechange_str{str}(:,[left_sbj;right_sbj]),1,2)./sqrt(nsb);
    h = fill([x, fliplr(x)], [curve1', fliplr(curve2')], cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x,nanmean(rulechange_str{str}(:,[left_sbj;right_sbj]),2),'-','Color',cmapStrategy(str,:)); hold on
    y = nanmean(rulechange_str{str}(:,[left_sbj;right_sbj]),2);
    text(x(end),y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([x(1) x(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
xlabel('Trial')
ylabel('Probability')
ylim([0.2 1])
title('Spatial-->Cued')
xlim([x(1) x(end)+20])

nsb = length([cueleft_sbj;cueright_sbj]);
subplot(122)
for str = 5:6 %nExploreStrategy
    curve1 = nanmean(rulechange_str{str}(:,[cueright_sbj; cueleft_sbj]),2)-nanstd(rulechange_str{str}(:,[cueright_sbj; cueleft_sbj]),1,2)./sqrt(nsb);
    curve2 = nanmean(rulechange_str{str}(:,[cueright_sbj; cueleft_sbj]),2)+nanstd(rulechange_str{str}(:,[cueright_sbj; cueleft_sbj]),1,2)./sqrt(nsb);
    h = fill([x, fliplr(x)], [curve1', fliplr(curve2')], cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x,nanmean(rulechange_str{str}(:,[cueright_sbj; cueleft_sbj]),2),'-','Color',cmapStrategy(str,:)); hold on
    y = nanmean(rulechange_str{str}(:,[cueright_sbj; cueleft_sbj]),2);
    text(x(end),y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([x(1) x(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
xlabel('Trial')
% ylabel('P(go cue)')
ylim([0.2 1])
title('Cued-->Spatial')
xlim([x(1) x(end)+20])

set(gcf,'PaperPositionMode','auto')
print([figpath 'fig5_LeverpressWSLSchoice_rulechange'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'fig5_LeverpressWSLSchoice_rulechange.png'])

%% plot explore strategies around rule change trials. Both from Ego to Cue 
% and viceversa
nsb = length([left_sbj;right_sbj]);
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize);
subplot(121)
for str = 7:8
    curve1 = nanmean(rulechange_str{str}(:,[left_sbj;right_sbj]),2)-nanstd(rulechange_str{str}(:,[left_sbj;right_sbj]),1,2)./sqrt(nsb);
    curve2 = nanmean(rulechange_str{str}(:,[left_sbj;right_sbj]),2)+nanstd(rulechange_str{str}(:,[left_sbj;right_sbj]),1,2)./sqrt(nsb);
    h = fill([x, fliplr(x)], [curve1', fliplr(curve2')], cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x,nanmean(rulechange_str{str}(:,[left_sbj;right_sbj]),2),'-','Color',cmapStrategy(str,:)); hold on
    y = nanmean(rulechange_str{str}(:,[left_sbj;right_sbj]),2);
    text(x(end),y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([x(1) x(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
xlabel('Trial')
ylabel('Probability')
ylim([0.2 1])
title('Spatial-->Cued')
xlim([x(1) x(end)+20])

nsb = length([cueleft_sbj;cueright_sbj]);
subplot(122)
for str = 7:8
    curve1 = nanmean(rulechange_str{str}(:,[cueright_sbj; cueleft_sbj]),2)-nanstd(rulechange_str{str}(:,[cueright_sbj; cueleft_sbj]),1,2)./sqrt(nsb);
    curve2 = nanmean(rulechange_str{str}(:,[cueright_sbj; cueleft_sbj]),2)+nanstd(rulechange_str{str}(:,[cueright_sbj; cueleft_sbj]),1,2)./sqrt(nsb);
    h = fill([x, fliplr(x)], [curve1', fliplr(curve2')], cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x,nanmean(rulechange_str{str}(:,[cueright_sbj; cueleft_sbj]),2),'-','Color',cmapStrategy(str,:)); hold on
    y = nanmean(rulechange_str{str}(:,[cueright_sbj; cueleft_sbj]),2);
    text(x(end),y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([x(1) x(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
xlabel('Trial')
% ylabel('P(go cue)')
ylim([0.2 1])
title('Cued-->Spatial')
xlim([x(1) x(end)+20])

set(gcf,'PaperPositionMode','auto')
print([figpath 'fig5_LeverpressWSLScue_rulechange'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'fig5_LeverpressWSLScue_rulechange.png'])
