clearvars; clc
close all
%% load Adrien data
% load strategy profiles data for 4 rats from Peyarache et al., 2009
% dataset
load('Processed_data/PeyracheDataStrategies.mat')
figpath = 'Figures\';

strategies_label = {'Go Right','Go Cued','Go Left','Go Uncued','Win-Stay-Spatial', ...
    'Lose-Shift-Spatial','Win-Stay-Cued','Lose-Shift-Cued', 'Alternate','Sticky'};

nRats = length(Output);  % number of rats
fields = fieldnames(Output{1});
nstrategy = numel(fields);

% visualise choices
cmapStrategy = [brewermap(4,'Set2'); brewermap(6,'Dark2')];
grey = [.6 .6 .6];
orange = [1 .45 0];
fontsize = 7;
axlinewidth = 0.5;
barsize = [5 5 10 3.5];
siz = 4;

%% plot strategy profiles around learning trial
% Header = {'Animal', 'Session', 'Direction', 'Reward', 'Rule', 'Light', 'SesName', 'Learning'};
load('Processed_data\SummaryDataTable_AllSessions.mat');
Data = double(Data);
Data(Data(:,4)==-1,4)=0;
learn = load('Processed_data\LearningTrials.txt');

learndata = NaN(size(learn,1),size(Data,2));
trialn = NaN(size(learn,1),1);
learning_from_start = NaN(size(learn,1),1);
% identify learning trials based on rule learned
for ses = 1:size(learn,1)
    sesidx = find(Data(:,7)==learn(ses,1)); % learning session. TrialsID in the learning session 
    trialn(ses) = sesidx(learn(ses,2)); % trial number in which the learning happen (from Data matrix). 
    ratid = unique(Data(sesidx,1)); % find which rat it is
    if numel(ratid)>1; disp('error'); break; end
    firstTrial = find(Data(:,1)==ratid,1,'first'); % find first trial for the current rat
    learning_from_start(ses) = trialn(ses) - firstTrial - 1;% trial number in which the learning happen from the beginning of the current animal training 
    learndata(ses,:) = Data(sesidx(learn(ses,2)),:); % all the data relative to learning trial
    if learndata(ses,end); continue; else; disp('Error'); break; end       
end
% add to learndata a column with the number of learning trial in
% progressive value
learndata = [learndata trialn learn(:,2) learning_from_start];
% add to learndata a column with the learning from last rulechange
newcolumn = NaN(size(learndata,1),1);
for ln = 1: size(learndata,1)
    iR = learndata(ln,1);
    idx = find(Data(:,1)==iR);
    rulechange = find(diff(Data(idx,5)) ~=0);
    rctrial = find(rulechange<learndata(ln, 11),1,'last');
       if rctrial
           newcolumn(ln) = learndata(ln,11)-rulechange(rctrial);
       else
           newcolumn(ln) = learndata(ln,11);
       end
end
learndata = [learndata newcolumn];

nbefore_learn = 9;
nafter_learn = 9;
egostr = find((contains(fields,'go_right') | contains(fields,'go_left')));
cuestr = find((contains(fields,'go_cued') | contains(fields,'go_uncued')));
trials_len = NaN(nRats,1);
learntrial = NaN(size(learndata,1),1);
ego_strategy = cell(nstrategy,1); cue_strategy = cell(nstrategy,1);
for iR = 1: nRats
    trials_len(iR) = size(Data(Data(:,1)==iR,:),1);
    rat_idx = find(learndata(:,1)==iR);
    learntrial(rat_idx) = trialn(rat_idx);
    if iR >1
        learntrial(rat_idx) = trialn(rat_idx)-sum(trials_len(1:iR-1));
    end

    % select MAP values around learning ego rule 
    egol = find(ismember(learndata(rat_idx,5),egostr));
        for k = 1:length(egol)
            for s = 1:nstrategy
                if sum(contains(fieldnames(Output{iR}.(fields{s})),'MAPprob_interpolated'))
                    ego_strategy{s} = [ego_strategy{s} Output{iR}.(fields{s}).MAPprob_interpolated(learntrial(rat_idx(egol(k)))-nbefore_learn:learntrial(rat_idx(egol(k)))+nafter_learn)];
                else
                    ego_strategy{s} = [ego_strategy{s} Output{iR}.(fields{s}).MAPprobability(learntrial(rat_idx(egol(k)))-nbefore_learn:learntrial(rat_idx(egol(k)))+nafter_learn)];
                end
            end
        end
    
%     % select MAP values around learning cue rule 
    cuel = find(ismember(learndata(rat_idx,5),cuestr));
        for k = 1:length(cuel)
            for s = 1:nstrategy
                if sum(contains(fieldnames(Output{iR}.(fields{s})),'MAPprob_interpolated'))
                    cue_strategy{s} = [cue_strategy{s} Output{iR}.(fields{s}).MAPprob_interpolated(learntrial(rat_idx(cuel(k)))-nbefore_learn:learntrial(rat_idx(cuel(k)))+nafter_learn)];
                else
                    cue_strategy{s} = [cue_strategy{s} Output{iR}.(fields{s}).MAPprobability(learntrial(rat_idx(cuel(k)))-nbefore_learn:learntrial(rat_idx(cuel(k)))+nafter_learn)];
                end
            end
        end
    
end

%% plot explore strategies around learning of cue and ego rule for Ymaze data (learning criterion)
x_learn = -nbefore_learn+1:nafter_learn+1;

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize);
subplot(121)
for str = 5:nstrategy 
    curve1 = mean(ego_strategy{str},2,'omitnan')-std(ego_strategy{str},1,2,'omitnan')./sqrt(nRats);
    curve2 = mean(ego_strategy{str},2,'omitnan')+std(ego_strategy{str},1,2,'omitnan')./sqrt(nRats);
    h = fill([x_learn, fliplr(x_learn)], [curve1', fliplr(curve2')],cmapStrategy(str,:)); hold on
    plot(x_learn,mean(ego_strategy{str},2,'omitnan'),'-','Color',cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
    plot([x_learn(1) x_learn(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
    set(gca,'FontName','Helvetica','FontSize',fontsize);
    set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
    y = mean(ego_strategy{str},2,'omitnan');
%     text(x_learn(end)+1,y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
xlabel('Trial')
ylabel('Probability')
title('Learning Spatial (from criterion)')
ylim([0.3 .85])
xlim([x_learn(1) x_learn(end)+5])

subplot(122)
for str = 5:nstrategy 
    curve1 = mean(cue_strategy{str},2,'omitnan')-std(cue_strategy{str},1,2,'omitnan')./sqrt(nRats);
    curve2 = mean(cue_strategy{str},2,'omitnan')+std(cue_strategy{str},1,2,'omitnan')./sqrt(nRats);
    h = fill([x_learn, fliplr(x_learn)], [curve1', fliplr(curve2')], cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x_learn,mean(cue_strategy{str},2,'omitnan'),'-','Color',cmapStrategy(str,:)); hold on
    plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
    plot([x_learn(1) x_learn(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
    set(gca,'FontName','Helvetica','FontSize',fontsize);
    set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
    y = mean(cue_strategy{str},2,'omitnan');
    text(x_learn(end)+1,y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
xlabel('Trial')
title('Learning Cued')
ylim([0.3 .85])
xlim([x_learn(1) x_learn(end)+5])

set(gcf,'PaperPositionMode','auto')
print([figpath 'Fig3_SI_YmazeExploreStr_criterion'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'Fig3_SI_YmazeExploreStr_criterion.png'])

%% plot 
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize);
subplot(121)
for str = 5:6
    curve1 = mean(ego_strategy{str},2,'omitnan')-std(ego_strategy{str},1,2,'omitnan')./sqrt(nRats);
    curve2 = mean(ego_strategy{str},2,'omitnan')+std(ego_strategy{str},1,2,'omitnan')./sqrt(nRats);
    h = fill([x_learn, fliplr(x_learn)], [curve1', fliplr(curve2')],cmapStrategy(str,:)); hold on
    plot(x_learn,mean(ego_strategy{str},2,'omitnan'),'-','Color',cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
    plot([x_learn(1) x_learn(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
    set(gca,'FontName','Helvetica','FontSize',fontsize);
    set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
    y = mean(ego_strategy{str},2,'omitnan');
    text(x_learn(end)+1,y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
xlabel('Trial')
ylabel('Probability')
title('Learning Spatial (from criterion)')
ylim([0.2 1])
xlim([x_learn(1) x_learn(end)+5])

subplot(122)
for str = 7:8
    curve1 = mean(cue_strategy{str},2,'omitnan')-std(cue_strategy{str},1,2,'omitnan')./sqrt(nRats);
    curve2 = mean(cue_strategy{str},2,'omitnan')+std(cue_strategy{str},1,2,'omitnan')./sqrt(nRats);
    h = fill([x_learn, fliplr(x_learn)], [curve1', fliplr(curve2')], cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x_learn,mean(cue_strategy{str},2,'omitnan'),'-','Color',cmapStrategy(str,:)); hold on
    plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
    plot([x_learn(1) x_learn(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
    set(gca,'FontName','Helvetica','FontSize',fontsize);
    set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
    y = mean(cue_strategy{str},2,'omitnan');
    text(x_learn(end)+1,y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
xlabel('Trial')
title('Learning Cued')
ylim([0.2 1])
xlim([x_learn(1) x_learn(end)+5])

set(gcf,'PaperPositionMode','auto')
print([figpath 'fig5_Ymaze_WSLS_criterion'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'fig5_Ymaze_WSLS_criterion.png'])


%% identify learning trial based on strategy profile (i.e. last trial after 
% which the successful strategy is above chance till the end of the session
lbefore = nbefore_learn-2;
lafter = nafter_learn;
x=-lbefore:lafter;
init = 9;
fine = 7;
learningstrategy = cell(nRats,1);
lstrategy = NaN(nRats, 4, lbefore+lafter+1);
for iR = 1 :nRats
   ratidx = find(Data(:,1)==iR);
   rulechange = find(diff(Data(ratidx,5)) ~=0);
   nses = length(unique(Data(ratidx,2)));
   nrule = length(unique(Data(ratidx,5)));
   for r = 1:nrule
       % find sessions that have this rule
       ridx = find(Data(ratidx,5)==r);
       sesidx = unique(Data(ratidx(ridx),2));
       for s = 1:length(sesidx)
           sestr = find(Data(ratidx,2)==sesidx(s));
           ltr = find(Output{iR}.(fields{r}).MAPprobability(sestr(init:end-fine))<=0.5,1,'last');
           if ltr & (Output{iR}.(fields{r}).MAPprobability(sestr(end-fine:end))>0.5) %>init & ltr<fine
%                ltrial = sestr(ltr);
               ltrial = sestr(ltr)+init;
           elseif sum(Output{iR}.(fields{r}).MAPprobability(sestr)>0.5)==length(sestr)
               [~, ltr] = min(Output{iR}.(fields{r}).MAPprobability(sestr(1:end-fine)));
               ltrial = sestr(ltr);
           else
               ltrial = [];
           end
           
           if ltrial
               rctrial = find(rulechange<ltrial,1,'last');
               if rctrial
                   learningstrategy{iR}(r,:) = [iR r sesidx(s) ltrial ltrial-rulechange(rctrial)];
               else
                   learningstrategy{iR}(r,:) = [iR r sesidx(s) ltrial ltrial];
               end
                   lstrategy(iR,r,:) = Output{iR}.(fields{r}).MAPprobability(ltrial-lbefore:ltrial+lafter);
               break
           end
       end
   end
end

lstrategydata = [];
for iR = 1:nRats
    lstrategydata = [lstrategydata; learningstrategy{iR}];
end

% store the explore strategies around learning from strategy
explorelstrategy = NaN(nrule,size(lstrategydata,1),length(x));
for ln = 1:size(lstrategydata,1)
    for st = 1:nstrategy
        if sum(contains(fieldnames(Output{lstrategydata(ln,1)}.(fields{st})),'MAPprob_interpolated'))
            explorelstrategy(ln,st,:) = Output{lstrategydata(ln,1)}.(fields{st}).MAPprob_interpolated(lstrategydata(ln,4)-lbefore:lstrategydata(ln,4)+lafter);
        else
            explorelstrategy(ln,st,:) = Output{lstrategydata(ln,1)}.(fields{st}).MAPprobability(lstrategydata(ln,4)-lbefore:lstrategydata(ln,4)+lafter);
        end
    end  
end

%% plot Explore strategies around learning from strategy
% store ego and cue rule in new dataset with learning defined by strategy
ego = find((lstrategydata(:,2)==1) | (lstrategydata(:,2)==3));
cue = find((lstrategydata(:,2)==2) | (lstrategydata(:,2)==4));
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize);
subplot(121)
for str = 5:nstrategy 
    curve1 = mean(squeeze(explorelstrategy(ego,str,:)),'omitnan')-std(squeeze(explorelstrategy(ego,str,:)),1,'omitnan')./sqrt(nRats);
    curve2 = mean(squeeze(explorelstrategy(ego,str,:)),'omitnan')+std(squeeze(explorelstrategy(ego,str,:)),1,'omitnan')./sqrt(nRats);
    h = fill([x, fliplr(x)], [curve1, fliplr(curve2)],cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x,mean(squeeze(explorelstrategy(ego,str,:)),'omitnan'),'-','Color',cmapStrategy(str,:)); hold on
    plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
    plot([x_learn(1) x_learn(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
    set(gca,'FontName','Helvetica','FontSize',fontsize);
    set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
    y = mean(squeeze(explorelstrategy(ego,str,:)),'omitnan');
    text(x_learn(end)+1,y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
xlabel('Trial')
ylabel('Probability')
title('Learning Spatial (from strategy)')
ylim([0.2 1])
xlim([x_learn(1) x_learn(end)+5])

subplot(122)
for str = 5:nstrategy
    curve1 = mean(squeeze(explorelstrategy(cue,str,:)),'omitnan')-std(squeeze(explorelstrategy(cue,str,:)),1,'omitnan')./sqrt(nRats);
    curve2 = mean(squeeze(explorelstrategy(cue,str,:)),'omitnan')+std(squeeze(explorelstrategy(cue,str,:)),1,'omitnan')./sqrt(nRats);
    h = fill([x, fliplr(x)], [curve1, fliplr(curve2)], cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x,mean(squeeze(explorelstrategy(cue,str,:)),'omitnan'),'-','Color',cmapStrategy(str,:)); hold on
    plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
    plot([x_learn(1) x_learn(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
    set(gca,'FontName','Helvetica','FontSize',fontsize);
    set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
    y = mean(squeeze(explorelstrategy(cue,str,:)),'omitnan');
    text(x_learn(end)+1,y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
xlabel('Trial')
title('Learning Cued')
ylim([0.2 1])
xlim([x_learn(1) x_learn(end)+5])

set(gcf,'PaperPositionMode','auto')
print([figpath 'Fig3_SI_YmazeExploreStr_strategy'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'Fig3_SI_YmazeExploreStr_strategy.png'])

%% plot WSLS strategies around learning (from strategy)
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize);
subplot(121)
for str = 5:6
    curve1 = mean(squeeze(explorelstrategy(ego,str,:)),'omitnan')-std(squeeze(explorelstrategy(ego,str,:)),1,'omitnan')./sqrt(nRats);
    curve2 = mean(squeeze(explorelstrategy(ego,str,:)),'omitnan')+std(squeeze(explorelstrategy(ego,str,:)),1,'omitnan')./sqrt(nRats);
    h = fill([x, fliplr(x)], [curve1, fliplr(curve2)],cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x,mean(squeeze(explorelstrategy(ego,str,:)),'omitnan'),'-','Color',cmapStrategy(str,:)); hold on
    plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
    plot([x_learn(1) x_learn(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
    set(gca,'FontName','Helvetica','FontSize',fontsize);
    set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
    y = mean(squeeze(explorelstrategy(ego,str,:)),'omitnan');
    text(x_learn(end)+1,y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
xlabel('Trial')
ylabel('Probability')
title('Learning Spatial (from strategy)')
ylim([0.2 1])
xlim([x_learn(1) x_learn(end)+5])

subplot(122)
for str = 7:8
    curve1 = mean(squeeze(explorelstrategy(cue,str,:)),'omitnan')-std(squeeze(explorelstrategy(cue,str,:)),1,'omitnan')./sqrt(nRats);
    curve2 = mean(squeeze(explorelstrategy(cue,str,:)),'omitnan')+std(squeeze(explorelstrategy(cue,str,:)),1,'omitnan')./sqrt(nRats);
    h = fill([x, fliplr(x)], [curve1, fliplr(curve2)], cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x,mean(squeeze(explorelstrategy(cue,str,:)),'omitnan'),'-','Color',cmapStrategy(str,:)); hold on
    plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
    plot([x_learn(1) x_learn(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
    set(gca,'FontName','Helvetica','FontSize',fontsize);
    set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
    y = mean(squeeze(explorelstrategy(cue,str,:)),'omitnan');
    text(x_learn(end)+1,y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
xlabel('Trial')
title('Learning Cued')
ylim([0.2 1])
xlim([x_learn(1) x_learn(end)+5])

set(gcf,'PaperPositionMode','auto')
print([figpath 'fig5_Ymaze_WSLS_strategy'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'fig5_Ymaze_WSLS_strategy.png'])

%% plot learning trial (from criterion and strategy) for each rule from the beginning of the training
jit = randn(size(learndata,1),1)*0.1-.15;
jit_str = randn(size(lstrategydata,1),1)*0.1+.15;
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[5 5 4 4]);
plot(learndata(:,5)+jit,learndata(:,11),'o','Color',grey,'MarkerSize', siz,...
    'MarkerFaceColor',grey,'MarkerEdgeColor','w'); hold on
plot(lstrategydata(:,2)+jit_str,lstrategydata(:,4),'o','Color',orange,...
    'MarkerSize', siz, 'MarkerFaceColor',orange,'MarkerEdgeColor','w'); hold on
% alpha .6
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
xticks(1:4)
xticklabels({'Go right','Go cued','Go left','Go un-cue'})
xtickangle(45)
text(.7,350,'Criterion','Color',grey,'FontSize',fontsize);
text(.7,300,'Strategy','Color',orange,'FontSize',fontsize);
xlim([0.5 4.5])
xlabel('Rule')
ylabel('Learning trial')


%% plot learning criterion compared to learning strategy from rule change points
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[5 5 4 4]);
plot(learndata(:,5)+jit,learndata(:,12),'o','Color',grey,'MarkerSize', siz,...
    'MarkerFaceColor',grey,'MarkerEdgeColor','w'); hold on
plot(lstrategydata(:,2)+jit_str,lstrategydata(:,5),'o','Color',orange,...
    'MarkerSize', siz, 'MarkerFaceColor',orange,'MarkerEdgeColor','w'); hold on
% alpha .6
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
xticks(1:4)
xticklabels({'Go right','Go cued','Go left','Go uncued'})
xtickangle(45)
text(2.3,170,'Criterion','Color',grey,'FontSize',fontsize);
text(2.3,140,'Strategy','Color',orange,'FontSize',fontsize);
xlim([0.5 4.5])
xlabel('Rule')
ylabel('Learning trial')

%% Stats
[~,p] = kstest2(learndata(:,12),lstrategydata(:,5));

print([figpath 'Fig3_YmazeLearningTrial_criterion_strategy'],'-depsc')
saveas(gcf,[figpath, 'Fig3_YmazeLearningTrial_criterion_strategy.png'])


%% plot successful rule strategy around learning from criterion and strategy
idx_ego = find((learndata(:,5)==1) | (learndata(:,5)==3));
idx_cue = find(learndata(:,5)==2); 
ego_learnstrategy = []; cue_learnstrategy = [];
for t = 1:length(idx_ego)
    ego_learnstrategy = [ego_learnstrategy ego_strategy{learndata(idx_ego(t),5)}(:,t)];
end
for t = 1:length(idx_cue)
    cue_learnstrategy = [cue_learnstrategy cue_strategy{learndata(idx_cue(t),5)}(:,t)];
end

fig2 = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize);
subplot(121)
curve1 = mean(ego_learnstrategy,2,'omitnan')-std(ego_learnstrategy,1,2,'omitnan')./sqrt(size(ego_learnstrategy,2));
curve2 = mean(ego_learnstrategy,2,'omitnan')+std(ego_learnstrategy,1,2,'omitnan')./sqrt(size(ego_learnstrategy,2));
h = fill([x_learn, fliplr(x_learn)], [curve1', fliplr(curve2')],grey); hold on
set(h,'facealpha',.3,'LineStyle','none')
plot(x_learn,mean(ego_learnstrategy,2,'omitnan'),'-','Color',grey); hold on

egolearn = squeeze([lstrategy(:,1,:); lstrategy(:,3,:)]);
egolearn(any(isnan(egolearn), 2), :) = [];
curve1 = mean(egolearn,'omitnan')-std(egolearn,1,'omitnan')./sqrt(size(egolearn,1));
curve2 = mean(egolearn,'omitnan')+std(egolearn,1,'omitnan')./sqrt(size(egolearn,1));
h = fill([x, fliplr(x)], [curve1, fliplr(curve2)], orange); hold on
set(h,'facealpha',.3,'LineStyle','none')
plot(x,mean(egolearn,'omitnan'),'-','Color',orange); hold on

plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([x(1) x(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
text(2,max(mean(ego_learnstrategy,2,'omitnan'))+.06,'Criterion','Fontsize',fontsize,'Color',grey);
text(2,max(mean(ego_learnstrategy,2,'omitnan'))+.14,'Strategy','Fontsize',fontsize,'Color',orange);
xlabel('Trial')
ylabel({'Probability of','correct strategy'})
ylim([.35 1])
xlim([x_learn(1) x_learn(end)])
title('Learning Spatial')

subplot(122)
curve1 = mean(cue_learnstrategy,2,'omitnan')-std(cue_learnstrategy,1,2,'omitnan')./sqrt(size(cue_learnstrategy,2));
curve2 = mean(cue_learnstrategy,2,'omitnan')+std(cue_learnstrategy,1,2,'omitnan')./sqrt(size(cue_learnstrategy,2));
h = fill([x_learn, fliplr(x_learn)], [curve1', fliplr(curve2')], grey); hold on
set(h,'facealpha',.3,'LineStyle','none')
plot(x_learn,mean(cue_learnstrategy,2,'omitnan'),'-','Color',grey); hold on

cuelearn = squeeze([lstrategy(:,2,:); lstrategy(:,4,:)]);
cuelearn(any(isnan(cuelearn), 2), :) = [];
curve1 = mean(cuelearn,'omitnan')-std(cuelearn,1,'omitnan')./sqrt(size(cuelearn,1));
curve2 = mean(cuelearn,'omitnan')+std(cuelearn,1,'omitnan')./sqrt(size(cuelearn,1));
h = fill([x, fliplr(x)], [curve1, fliplr(curve2)], orange); hold on
set(h,'facealpha',.3,'LineStyle','none')
plot(x,mean(cuelearn,'omitnan'),'-','Color',orange); hold on

plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([x(1) x(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
xlabel('Trial')
% ylabel('P(successful strategy)')
ylim([.35 1])
xlim([x_learn(1) x_learn(end)])
title('Learning Cued')

set(gcf,'PaperPositionMode','auto')
print([figpath 'Fig3_YmazeLearning_criterion_strategy'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'Fig3_YmazeLearning_criterion_strategy.png'])


%% ****************************************************************************
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

nbefore = nbefore_learn;
nafter = nafter_learn;
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

%% plot strategy profiles around learning trial for Tobias data
% first find the learning trials as the first of 10 consecutive correct
% trials. Omit trials (Sbj.choice=0) will not count for the identification 
% of the learning trial

learn = NaN(nSbj,3); %first column contain learning first rule, second column
% learning second rule from beginning of training, third col learning second 
% rule from rule change trial
for iS = 1:nSbj
    %% find learning of first rule
    if isnan(rc(iS,1)); rc(iS,1) = length(Output{iS}.(fields{1}).alpha); end   
    
    idx_sbj = find(Datas(:,1)==sbj_name(iS));
    omit = find(Datas(idx_sbj(1:rc(iS,1)-1),3)==0);
    % remove omission trial to reward array and rename
    reward = Datas(idx_sbj(1:rc(iS,1)-1),4);
    reward(omit) = [];
    % find all the index with 10 consecutive correct (reward = 1)
    ind10cor = strfind(reward',ones(1,10));
    
    if sum(ind10cor>10)
        possible_learn = ind10cor(ind10cor>10);
        idx_omit = find(omit<=possible_learn(1));
        learn(iS,1) = possible_learn(1)+length(idx_omit);
    end
    
    %% find learning second rule
    if rc(iS,1)<length(Datas(idx_sbj,3))
        clear omit reward
        omit = find(Datas(idx_sbj(rc(iS,1):end),3)==0);
        % remove omission trial to reward array and rename
        reward = Datas(idx_sbj(rc(iS,1):end),4);
        reward(omit) = [];
        ind10cor = strfind(reward',ones(1,10));
        if ind10cor
            idx_omit = find(omit<=ind10cor(1));
            learn(iS,2) = ind10cor(1)+length(idx_omit)+rc(iS,1)-1;
            learn(iS,3) = ind10cor(1)+length(idx_omit);
        end
    end
end

%% identify learning based on Strategy Analysis
learn_leverpress = NaN(nSbj,3);
%% match the lever-press rules (1=right,2=left, 3=cue) with the succesfull 
% strategy (1=go right, 2=go cue, 3=go left) for each subject
rules = rc(:,2:3); % first colum is first rule, second is second rule
rules(rc(:,2:3)==2)=3;
rules(rc(:,2:3)==3)=2;
for iS = 1 :nSbj
    sbj_idx = find(Datas(:,1)==sbj_name(iS)); % identify in Datas the current subject    
    if ~isnan(rc(iS,2)) 
        %% ******************************
        %% find learning of first rule
        trial = find(Output{iS}.(fields{rules(iS,1)}).MAPprobability(1:rc(iS,1))<=.5,1,'last');
        if all(Output{iS}.(fields{rules(iS,1)}).MAPprobability(1:rc(iS,1))>.5,'all')
            [~,min_idx] = min(Output{iS}.(fields{rules(iS,1)}).MAPprobability(1:rc(iS,1)));
                learn_leverpress(iS,1) = min_idx;
        elseif trial %>init & trial < rc(iS,1)-fine
            learn_leverpress(iS,1) = trial;
        end
    end
    %% ******************************
    %% find learning of second rule
    if rc(iS,1)<length(Datas(sbj_idx,3)) & ~isnan(rc(iS,3)) 
        postrial = find(Output{iS}.(fields{rules(iS,2)}).MAPprobability(rc(iS,1):end)<=.5,1,'last');
        % length training secodn rule
        len_train_secondrule = length(Output{iS}.(fields{rules(iS,2)}).MAPprobability(rc(iS,1):end)); 
       if postrial  %< len_train_secondrule-fine
           learn_leverpress(iS,2) =  postrial+rc(iS,1);
           learn_leverpress(iS,3) =  postrial;
       elseif all(Output{iS}.(fields{rules(iS,2)}).MAPprobability(rc(iS,1):end)>.5,'all')
           [~,min_idx] = min(Output{iS}.(fields{rules(iS,2)}).MAPprobability(rc(iS,1):end));
%            if min_idx>init & min_idx<len_train_secondrule-fine
               learn_leverpress(iS,2) = min_idx + rc(iS,1); 
               learn_leverpress(iS,3) = min_idx; 
%            end
       end
    end
            
    learn_leverpress(learn_leverpress<8)=8;        
            
end

%% Plot learning first and second rule (learning from criterion and strategy)
% Learning from beginning of training
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[5 5 8 3.5]);
subplot(1,2,1)
jit = randn(size(learn,1),1)*0.07-.15;
jit_str = randn(size(learn_leverpress,1),1)*0.07+.15;
% plot([ones(size(learn,1),1).*rules(:,1)+jit ones(size(learn_leverpress,1),1).*rules(:,1)+jit_str]',...
%     [learn(:,1) learn_leverpress(:,1)]','-','Color',[.6 .6 .6]); hold on
plot(rules(:,1)+jit,learn(:,1),'o','MarkerSize', siz, 'MarkerFaceColor',grey,'MarkerEdgeColor','w'); hold on
plot(rules(:,1)+jit_str,learn_leverpress(:,1),'o','MarkerSize', siz, 'MarkerFaceColor',orange,'MarkerEdgeColor','w'); hold on
% add stats *******************
p1 = signrank(learn(:,1), learn_leverpress(:,1));
p2 = signrank(learn(:,3), learn_leverpress(:,3));
text(1.3,max(learn(:,1))+50,['p=',num2str(round(p1,2))],'FontSize',fontsize);
% *****************************
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
text(.7,300,'Criterion','FontSize',fontsize,'Color',grey);
text(.7,250,'Strategy','FontSize',fontsize,'Color',orange);
xticks(1:3)
xticklabels({'Go right','Go cued','Go left'})
xtickangle(45)
xlim([.5 3.5])
ylabel('Learning trial')
title('First rule')
subplot(1,2,2)
plot(rules(:,2)+jit,learn(:,3),'o','MarkerSize', siz, 'MarkerFaceColor',grey,'MarkerEdgeColor','w'); hold on
plot(rules(:,2)+jit_str,learn_leverpress(:,3),'o','MarkerSize', siz, 'MarkerFaceColor',orange,'MarkerEdgeColor','w'); hold on
text(1.3,650,['p=',num2str(round(p2,2))],'FontSize',fontsize);
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
xticks(1:3)
xticklabels({'Go right','Go cued','Go left'})
xtickangle(45)
xlim([.5 3.5])
title('Second rule')

print([figpath 'Fig3_LeverpressLearningTrial_criterion_strategy'],'-depsc')
saveas(gcf,[figpath, 'Fig3_LeverpressLearningTrial_criterion_strategy.png'])

%% record rule strategies around learning (criterion) of first and second rule
% nafter_learn = 9;
learning_str = cell(nstrategy,2);
for str = 1: nstrategy
    learning_str{str,1} = NaN(nafter_learn+nbefore_learn+1,nSbj);
    learning_str{str,2} = NaN(nafter_learn+nbefore_learn+1,nSbj);
    for iS = 1:nSbj
        if ~isnan(learn(iS,1))
            if sum(contains(fieldnames(Output{iS}.(fields{str})),'MAPprob_interpolated'))
                learning_str{str,1}(:,iS) = Output{iS}.(fields{str}).MAPprob_interpolated(learn(iS,1)-nbefore_learn:learn(iS,1)+nafter_learn);
            else
                learning_str{str,1}(:,iS) = Output{iS}.(fields{str}).MAPprobability(learn(iS,1)-nbefore_learn:learn(iS,1)+nafter_learn);
            end
        end
        
        if ~isnan(learn(iS,2))
            if sum(contains(fieldnames(Output{iS}.(fields{str})),'MAPprob_interpolated'))
                learning_str{str,2}(:,iS) = Output{iS}.(fields{str}).MAPprob_interpolated(learn(iS,2)-nbefore_learn:learn(iS,2)+nafter_learn);
            else
                learning_str{str,2}(:,iS) = Output{iS}.(fields{str}).MAPprobability(learn(iS,2)-nbefore_learn:learn(iS,2)+nafter_learn);
            end
        end        
    end
end 

%% record rule strategies around learning (from strategy profile) of first and second rule
% nafter_learn = 9;
learning_strategy = cell(nstrategy,2);
for str = 1: nstrategy
    learning_strategy{str,1} = NaN(nafter_learn+nbefore_learn+1,nSbj);
    learning_strategy{str,2} = NaN(nafter_learn+nbefore_learn+1,nSbj);
    for iS = 1:nSbj
        if ~isnan(learn_leverpress(iS,1)) && learn_leverpress(iS,1)-nbefore_learn>0
            if sum(contains(fieldnames(Output{iS}.(fields{str})),'MAPprob_interpolated'))
                learning_strategy{str,1}(:,iS) = Output{iS}.(fields{str}).MAPprob_interpolated(learn_leverpress(iS,1)-nbefore_learn:learn_leverpress(iS,1)+nafter_learn);
            else
                learning_strategy{str,1}(:,iS) = Output{iS}.(fields{str}).MAPprobability(learn_leverpress(iS,1)-nbefore_learn:learn_leverpress(iS,1)+nafter_learn);
            end
        end
        
        if ~isnan(learn_leverpress(iS,2)) && learn_leverpress(iS,2)+nafter_learn<=length(Output{iS}.(fields{str}).alpha)
            if sum(contains(fieldnames(Output{iS}.(fields{str})),'MAPprob_interpolated'))
                learning_strategy{str,2}(:,iS) = Output{iS}.(fields{str}).MAPprob_interpolated(learn_leverpress(iS,2)-nbefore_learn:learn_leverpress(iS,2)+nafter_learn);
            else
                learning_strategy{str,2}(:,iS) = Output{iS}.(fields{str}).MAPprobability(learn_leverpress(iS,2)-nbefore_learn:learn_leverpress(iS,2)+nafter_learn);
            end
        end        
    end
end 


% identify subject with shift from right or left to cue
[right_sbj, ~] = find(ismember(rc(:,2:3),[1 3],'rows'));
[left_sbj,~] = find(ismember(rc(:,2:3),[2 3],'rows'));

% identify subject with shift from cue to ego right or left
[cueright_sbj, ~] = find(ismember(rc(:,2:3),[3 1],'rows'));
[cueleft_sbj,~] = find(ismember(rc(:,2:3),[3 2],'rows'));

str_r = find(contains(strategies_label,'Go Right'));
str_l = find(contains(strategies_label,'Go Left'));
% Second subplot show the succesful strategy for cue rule in the shift from Cue --> Ego
% and Ego-->Cue
str_c = find(contains(strategies_label,'Go Cued'));
% str_c = find(contains(strRuleStrategy,'Go Cue'));

% collect the MAPs around learning the ego rule for each
learnE_egocue = [learning_str{str_r,1}(:,right_sbj) learning_str{str_l,1}(:,left_sbj)]; 
learnE_cueego = [learning_str{str_r,2}(:,cueright_sbj) learning_str{str_l,2}(:,cueleft_sbj)]; 
learnEgo = [learnE_egocue learnE_cueego];

% collect the MAPs around learning the cue rule for each subject
learnC_egocue = [learning_str{str_c,2}(:,right_sbj) learning_str{str_c,2}(:,left_sbj)]; 
learnC_cueego = [learning_str{str_c,1}(:,cueright_sbj) learning_str{str_c,1}(:,cueleft_sbj)]; 
learnCue = [learnC_egocue learnC_cueego];

% collect the MAPs around learning (from strategy) the ego rule for each
learnstrategyE_egocue = [learning_strategy{str_r,1}(:,right_sbj) learning_strategy{str_l,1}(:,left_sbj)]; 
learnstrategyE_cueego = [learning_strategy{str_r,2}(:,cueright_sbj) learning_strategy{str_l,2}(:,cueleft_sbj)]; 
learnstrategyEgo = [learnstrategyE_egocue learnstrategyE_cueego];

% collect the MAPs around learning (from strategy) the cue rule for each subject
learnstrategyC_egocue = [learning_strategy{str_c,2}(:,right_sbj) learning_strategy{str_c,2}(:,left_sbj)]; 
learnstrategyC_cueego = [learning_strategy{str_c,1}(:,cueright_sbj) learning_strategy{str_c,1}(:,cueleft_sbj)]; 
learnstrategyCue = [learnstrategyC_egocue learnstrategyC_cueego];

%% plot explore strategies around learning trials for Tobias data
explorelearning_ego = cell(nstrategy,1);
explorelearning_cue = cell(nstrategy,1);
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize);
subplot(121)
for str = 5:nstrategy
    explorelearning_ego{str} = [learning_str{str,1}(:,[right_sbj;left_sbj]) learning_str{str,2}(:,[cueright_sbj;cueleft_sbj])];   
    curve1 = mean(explorelearning_ego{str},2,'omitnan')-std(explorelearning_ego{str},1,2,'omitnan')./sqrt(nSbj);
    curve2 = mean(explorelearning_ego{str},2,'omitnan')+std(explorelearning_ego{str},1,2,'omitnan')./sqrt(nSbj);
    h = fill([x_learn, fliplr(x_learn)], [curve1', fliplr(curve2')], cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x_learn,mean(explorelearning_ego{str},2,'omitnan'),'-','Color',cmapStrategy(str,:)); hold on
    y = mean(explorelearning_ego{str},2,'omitnan');
%     text(x_learn(end),y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([x_learn(1) x_learn(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
xlabel('Trial')
ylabel('Probability')
title('Learning Spatial (from criterion)')
ylim([.2 1])
xlim([x_learn(1) x_learn(end)+5])

subplot(122)
for str = 5:nstrategy
    explorelearning_cue{str} = [learning_str{str,2}(:,[right_sbj;left_sbj]) learning_str{str,1}(:,[cueright_sbj;cueleft_sbj])]; 
    
    curve1 = mean(explorelearning_cue{str},2,'omitnan')-std(explorelearning_cue{str},1,2,'omitnan')./sqrt(nSbj);
    curve2 = mean(explorelearning_cue{str},2,'omitnan')+std(explorelearning_cue{str},1,2,'omitnan')./sqrt(nSbj);
    h = fill([x_learn, fliplr(x_learn)], [curve1', fliplr(curve2')], cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x_learn,mean(explorelearning_cue{str},2,'omitnan'),'-','Color',cmapStrategy(str,:)); hold on
    y = mean(explorelearning_cue{str},2,'omitnan');
    text(x_learn(end),y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([x_learn(1) x_learn(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
xlabel('Trial')
% ylabel('P(explore strategy)')
ylim([.2 1])
xlim([x_learn(1) x_learn(end)+5])
title('Learning Cued')

set(gcf,'PaperPositionMode','auto')
print([figpath 'Fig3_SI_LeverpressExploreStr_criterion'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'Fig3_SI_LeverpressExploreStr_criterion.png'])

%% plot WSLS around learning (from criterion)
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize);
subplot(121)
for str = 5:6
    explorelearning_ego{str} = [learning_str{str,1}(:,[right_sbj;left_sbj]) learning_str{str,2}(:,[cueright_sbj;cueleft_sbj])];   
    curve1 = mean(explorelearning_ego{str},2,'omitnan')-std(explorelearning_ego{str},1,2,'omitnan')./sqrt(nSbj);
    curve2 = mean(explorelearning_ego{str},2,'omitnan')+std(explorelearning_ego{str},1,2,'omitnan')./sqrt(nSbj);
    h = fill([x_learn, fliplr(x_learn)], [curve1', fliplr(curve2')], cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x_learn,mean(explorelearning_ego{str},2,'omitnan'),'-','Color',cmapStrategy(str,:)); hold on
    y = mean(explorelearning_ego{str},2,'omitnan');
%     text(x_learn(end),y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([x_learn(1) x_learn(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
xlabel('Trial')
ylabel('Probability')
title('Learning Spatial (from criterion)')
ylim([.2 1])
xlim([x_learn(1) x_learn(end)+5])

subplot(122)
for str = 7:8
    explorelearning_cue{str} = [learning_str{str,2}(:,[right_sbj;left_sbj]) learning_str{str,1}(:,[cueright_sbj;cueleft_sbj])]; 
    
    curve1 = mean(explorelearning_cue{str},2,'omitnan')-std(explorelearning_cue{str},1,2,'omitnan')./sqrt(nSbj);
    curve2 = mean(explorelearning_cue{str},2,'omitnan')+std(explorelearning_cue{str},1,2,'omitnan')./sqrt(nSbj);
    h = fill([x_learn, fliplr(x_learn)], [curve1', fliplr(curve2')], cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x_learn,mean(explorelearning_cue{str},2,'omitnan'),'-','Color',cmapStrategy(str,:)); hold on
    y = mean(explorelearning_cue{str},2,'omitnan');
    text(x_learn(end),y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([x_learn(1) x_learn(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
xlabel('Trial')
% ylabel('P(explore strategy)')
ylim([.2 1])
xlim([x_learn(1) x_learn(end)+5])
title('Learning Cued')

set(gcf,'PaperPositionMode','auto')
print([figpath 'fig5_Leverpress_WSLS_criterion'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'fig5_Leverpress_WSLS_criterion.png'])

%% plot explore strategies around learning trials (from strategy analysis) for Tobias data
explorelearning_ego = cell(nstrategy,1);
explorelearning_cue = cell(nstrategy,1);
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize);
subplot(121)
for str = 5:nstrategy
    explorelearning_ego{str} = [learning_strategy{str,1}(:,[right_sbj;left_sbj]) learning_strategy{str,2}(:,[cueright_sbj;cueleft_sbj])];   
    curve1 = mean(explorelearning_ego{str},2,'omitnan')-std(explorelearning_ego{str},1,2,'omitnan')./sqrt(nSbj);
    curve2 = mean(explorelearning_ego{str},2,'omitnan')+std(explorelearning_ego{str},1,2,'omitnan')./sqrt(nSbj);
    h = fill([x_learn, fliplr(x_learn)], [curve1', fliplr(curve2')], cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x_learn,mean(explorelearning_ego{str},2,'omitnan'),'-','Color',cmapStrategy(str,:)); hold on
    y = mean(explorelearning_ego{str},2,'omitnan');
%     text(x_learn(end),y(end),strExploreStrategy{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([x_learn(1) x_learn(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
xlabel('Trial')
ylabel('Probability')
title('Learning Spatial (from strategy)')
ylim([.2 1])
xlim([x_learn(1) x_learn(end)+5])

subplot(122)
for str = 5:nstrategy
    explorelearning_cue{str} = [learning_strategy{str,2}(:,[right_sbj;left_sbj]) learning_strategy{str,1}(:,[cueright_sbj;cueleft_sbj])]; 
    
    curve1 = mean(explorelearning_cue{str},2,'omitnan')-std(explorelearning_cue{str},1,2,'omitnan')./sqrt(nSbj);
    curve2 = mean(explorelearning_cue{str},2,'omitnan')+std(explorelearning_cue{str},1,2,'omitnan')./sqrt(nSbj);
    h = fill([x_learn, fliplr(x_learn)], [curve1', fliplr(curve2')], cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x_learn,mean(explorelearning_cue{str},2,'omitnan'),'-','Color',cmapStrategy(str,:)); hold on
    y = mean(explorelearning_cue{str},2,'omitnan');
    text(x_learn(end),y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([x_learn(1) x_learn(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
xlabel('Trial')
% ylabel('P(explore strategy)')
ylim([.2 1])
xlim([x_learn(1) x_learn(end)+5])
title('Learning Cued')

set(gcf,'PaperPositionMode','auto')
print([figpath 'Fig3_SI_LeverpressExploreStr_strategy'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'Fig3_SI_LeverpressExploreStr_strategy.png'])

%% plot WSLS around the learning (from strategy)
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize);
subplot(121)
for str = 5:6
    explorelearning_ego{str} = [learning_strategy{str,1}(:,[right_sbj;left_sbj]) learning_strategy{str,2}(:,[cueright_sbj;cueleft_sbj])];   
    curve1 = mean(explorelearning_ego{str},2,'omitnan')-std(explorelearning_ego{str},1,2,'omitnan')./sqrt(nSbj);
    curve2 = mean(explorelearning_ego{str},2,'omitnan')+std(explorelearning_ego{str},1,2,'omitnan')./sqrt(nSbj);
    h = fill([x_learn, fliplr(x_learn)], [curve1', fliplr(curve2')], cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x_learn,mean(explorelearning_ego{str},2,'omitnan'),'-','Color',cmapStrategy(str,:)); hold on
    y = mean(explorelearning_ego{str},2,'omitnan');
%     text(x_learn(end),y(end),strExploreStrategy{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([x_learn(1) x_learn(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
xlabel('Trial')
ylabel('Probability')
title('Learning Spatial (from strategy)')
ylim([.2 1])
xlim([x_learn(1) x_learn(end)+5])

subplot(122)
for str = 7:8
    explorelearning_cue{str} = [learning_strategy{str,2}(:,[right_sbj;left_sbj]) learning_strategy{str,1}(:,[cueright_sbj;cueleft_sbj])]; 
    
    curve1 = mean(explorelearning_cue{str},2,'omitnan')-std(explorelearning_cue{str},1,2,'omitnan')./sqrt(nSbj);
    curve2 = mean(explorelearning_cue{str},2,'omitnan')+std(explorelearning_cue{str},1,2,'omitnan')./sqrt(nSbj);
    h = fill([x_learn, fliplr(x_learn)], [curve1', fliplr(curve2')], cmapStrategy(str,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x_learn,mean(explorelearning_cue{str},2,'omitnan'),'-','Color',cmapStrategy(str,:)); hold on
    y = mean(explorelearning_cue{str},2,'omitnan');
    text(x_learn(end),y(end),strategies_label{str},'Fontsize',fontsize,'Color',cmapStrategy(str,:));
end
plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([x_learn(1) x_learn(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
xlabel('Trial')
% ylabel('P(explore strategy)')
ylim([.2 1])
xlim([x_learn(1) x_learn(end)+5])
title('Learning Cued')

set(gcf,'PaperPositionMode','auto')
print([figpath 'fig5_Leverpress_WSLS_strategy'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'fig5_Leverpress_WSLS_strategy.png'])



x_learn = -nbefore_learn+1:nafter_learn+1;

%% learning Tobias compare criterion and strategy learning
fig2 = figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize);
subplot(121)
curve1 = mean(learnEgo,2,'omitnan')-std(learnEgo,1,2,'omitnan')./sqrt(size(learnEgo,2));
curve2 = mean(learnEgo,2,'omitnan')+std(learnEgo,1,2,'omitnan')./sqrt(size(learnEgo,2));
h = fill([x_learn, fliplr(x_learn)], [curve1', fliplr(curve2')], grey); hold on
set(h,'facealpha',.3,'LineStyle','none')
plot(x_learn,mean(learnEgo,2,'omitnan'),'-','Color',grey); hold on

curve1 = mean(learnstrategyEgo,2,'omitnan')-std(learnstrategyEgo,1,2,'omitnan')./sqrt(size(learnstrategyEgo,2));
curve2 = mean(learnstrategyEgo,2,'omitnan')+std(learnstrategyEgo,1,2,'omitnan')./sqrt(size(learnstrategyEgo,2));
h = fill([x_learn, fliplr(x_learn)], [curve1', fliplr(curve2')], orange); hold on
set(h,'facealpha',.3,'LineStyle','none')
plot(x_learn,mean(learnstrategyEgo,2,'omitnan'),'-','Color',orange); hold on

plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([x(1) x(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
text(-7,max(mean(learnEgo,2,'omitnan')),'Criterion','Fontsize',fontsize,'Color',grey);
text(-7,max(mean(learnCue,2,'omitnan'))-.1,'Strategy','Fontsize',fontsize,'Color',orange);
xlabel('Trial')
ylabel({'Probability of','correct strategy'})
ylim([.35 1])
xlim([x_learn(1) x_learn(end)])
title('Learning Ego')

subplot(122)
curve1 = mean(learnCue,2,'omitnan')-std(learnCue,1,2,'omitnan')./sqrt(size(learnCue,2));
curve2 = mean(learnCue,2,'omitnan')+std(learnCue,1,2,'omitnan')./sqrt(size(learnCue,2));
h = fill([x_learn, fliplr(x_learn)], [curve1', fliplr(curve2')], grey); hold on
set(h,'facealpha',.3,'LineStyle','none')
plot(x_learn,mean(learnCue,2,'omitnan'),'-','Color',grey); hold on

curve1 = mean(learnstrategyCue,2,'omitnan')-std(learnstrategyCue,1,2,'omitnan')./sqrt(size(learnstrategyCue,2));
curve2 = mean(learnstrategyCue,2,'omitnan')+std(learnstrategyCue,1,2,'omitnan')./sqrt(size(learnstrategyCue,2));
h = fill([x_learn, fliplr(x_learn)], [curve1', fliplr(curve2')],orange); hold on
set(h,'facealpha',.3,'LineStyle','none')
plot(x_learn,mean(learnstrategyCue,2,'omitnan'),'-','Color',orange); hold on

plot([0 0],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([x(1) x(end)],[.5 .5],'--','Color',[.6 .6 .6]); hold on

set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
xlabel('Trial')
% ylabel('P(successful strategy)')
ylim([.35 1])
xlim([x_learn(1) x_learn(end)])
title('Learning Cued')

set(gcf,'PaperPositionMode','auto')
print([figpath 'Fig3_LeverpressLearning_Criterion_strategy'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'Fig3_LeverpressLearning_criterion_strategy.png'])
