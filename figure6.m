clearvars; close all;

% common parameters
fontsize = 7;
axlinewidth = .5;

% colors = distinguishable_colors(ngroups);
colors = [0 0 1; 1 0 0; 0 1 0; 0 0 0];
figpath = 'Figures\';
nTrials = 30;


% strategy translation for Human Instrumental Learning (Pessiglione et al.,
% 2006)
humanStr_label = ["Go Bottom","Go Gain/Avoid Loss","Go Top","Avoid Gain/Go Loss",...
    "Win-Stay-Location","Lose-Shift-Location","Win-Stay-Stimulus",...
    "Lose-Shift-Stimulus","Alternate","Sticky"];

% visualise choices
cmapPrior = brewermap(2,'PuOr');
cmapUnif = brewermap(2,'Oranges');
cmapJeff = brewermap(2,'Purples');
cmaptask = brewermap(3,'Set1');

barsize = [5 5 10 3.5];

%% rules
% rule(1) = 'Gain';
% rule(2) = 'Loss';
% rule(3) = 'Look';
trial_type = ["go_gain", "avoid_loss", "look"];

strategy = find(humanStr_label=='Go Gain/Avoid Loss'); % look for 'go gain/avoid loss' strategy
bottom = find(humanStr_label=='Go Bottom'); % look for 'go gain/avoid loss' strategy
top = find(humanStr_label=='Go Top'); % look for 'go gain/avoid loss' strategy
trtype = find((trial_type == 'go_gain') | (trial_type == 'avoid_loss'));

%% for an example subject show Rule startegies profile and stability 
% for changes in prior
uni = load('Processed_data/HumanDataStrategies.mat');
jef = load('Processed_data/HumanDataStrategies_Jeffreys.mat');

x = 1:nTrials;
nsbj = length(uni.Output);
learn = zeros(nsbj,2); 
bias = NaN(nsbj,1);

fields = fieldnames(uni.Output{1}.go_gain);

for sbj = 1:nsbj
    lg = find(uni.Output{sbj}.go_gain.(fields{strategy}).MAPprobability<=.5,1,'last');
    ll = find(uni.Output{sbj}.avoid_loss.(fields{strategy}).MAPprobability<=.5,1,'last');
    if ~isempty(lg)
        learn(sbj,1) = lg;
    end
    if ~isempty(ll)
        learn(sbj,2) = ll;
    end
end
learn(learn==nTrials)=NaN;

%% plot MAP profile for succesful strategy for 'Gain' and 'Loss' trials.
% Compare the profiles for different prior. Highlight the learning point as
% the last trial after which the strategy profile is above chance
sbj = 19;
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize)
subplot(121)
plot(uni.Output{sbj}.go_gain.(fields{strategy}).MAPprobability,'-','Color',cmapUnif(1,:)); hold on
plot(jef.Output{sbj}.go_gain.(fields{strategy}).MAPprobability,'.','Color',cmapUnif(2,:)); hold on
plot([x(1) x(end)], [0.5 0.5], '--', 'Color', [.6 .6 .6]); hold on
text(12,0.15,'Uniform','Color',cmapUnif(1,:),'FontSize',fontsize);
h=text(learn(sbj,1), .55,'\leftarrow');
set(h,'Rotation',90);
text(12,0.35,'Jeffreys','Color',cmapUnif(2,:),'FontSize',fontsize);
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
ylabel('Probability')
title('Go Gain')
xlabel('Trial')
ylim([0 1])

subplot(122)
plot(uni.Output{sbj}.avoid_loss.(fields{strategy}).MAPprobability,'-','Color',cmapJeff(1,:)); hold on
plot(jef.Output{sbj}.avoid_loss.(fields{strategy}).MAPprobability,'.','Color',cmapJeff(2,:)); hold on
plot([x(1) x(end)], [0.5 0.5], '--', 'Color', [.6 .6 .6]); hold on
h=text(learn(sbj,2), .45,'\leftarrow');
set(h,'Rotation',-90);
text(12,0.6,'Uniform','Color',cmapJeff(1,:),'FontSize',fontsize);
text(12,0.75,'Jeffreys','Color',cmapJeff(2,:),'FontSize',fontsize);
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
ylim([0 1])
title('Avoid Loss')
xlabel('Trial')

set(gcf,'PaperPositionMode','auto')
print([figpath 'fig6_exampleSbj_UniJef'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'fig6_exampleSbj_UniJef.png'])

%% plot difference between Uniform and Jeffreys prior for example subject
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize)
subplot(121)
plot(uni.Output{sbj}.go_gain.(fields{strategy}).MAPprobability-...
    jef.Output{sbj}.go_gain.(fields{strategy}).MAPprobability,'-','Color',cmapUnif(1,:)); hold on
text(10,0.03,'Uniform-Jeffreys','Color',cmapUnif(1,:),'FontSize',fontsize);
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
ylabel({'Difference','Uniform-Jeffreys'})
title('Go Gain')
xlabel('Trial')

subplot(122)
plot(uni.Output{sbj}.avoid_loss.(fields{strategy}).MAPprobability-...
    jef.Output{sbj}.avoid_loss.(fields{strategy}).MAPprobability,'-','Color',cmapJeff(1,:)); hold on
text(10,0.07,'Uniform-Jeffreys','Color',cmapJeff(1,:),'FontSize',fontsize);
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
title('Avoid Loss')
xlabel('Trial')

set(gcf,'PaperPositionMode','auto')
print([figpath 'figS6_exampleSbj_UniJef'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'figS6_exampleSbj_UniJef.png'])

%% Get a list of all files in the folder with the desired file name pattern.
% filePattern = fullfile('Processed Data/Musa data/', '*.mat'); % Change to whatever pattern you need.
% theFiles = dir(filePattern);

gogain = []; %cell(ngroups,1);
avoidloss = []; %cell(ngroups,1);

jef_gogain = [];
jef_avoidloss = [];

explore_gain = cell(6,1);
explore_loss = cell(6,1);

for sbj = 1 :nsbj
    
        %% record 'go gain' strategy for gain trials
        gogain = [gogain; uni.Output{sbj}.go_gain.(fields{strategy}).MAPprobability']; % sbj(1).Rule.MAPts(:,strategy)'];   
        jef_gogain = [jef_gogain; jef.Output{sbj}.go_gain.(fields{strategy}).MAPprobability']; % sbj(1).Rule.MAPts(:,strategy)'];   
        %% record 'avoid loss' strategy for loss trial 
        avoidloss = [avoidloss; uni.Output{sbj}.avoid_loss.(fields{strategy}).MAPprobability']; %sbj(2).Rule.MAPts(:,strategy)'];
        jef_avoidloss = [jef_avoidloss; jef.Output{sbj}.avoid_loss.(fields{strategy}).MAPprobability']; %sbj(2).Rule.MAPts(:,strategy)'];
        
        %% record proportion of bias trials
        bias(sbj) = sum(uni.Output{sbj}.look.(fields{bottom}).MAPprobability>0.5)./nTrials;
        bias(sbj) = max([bias(sbj) 1-bias(sbj)]); 
        
        %% record explore strategies for subjects that learned (gain and loss trials)
        for str = 5:length(humanStr_label)
            if ~isnan(learn(sbj,1))
%                 explore_gain{ng} = [sbj(1).Explore.MAPts];
                if sum(contains(fieldnames(uni.Output{sbj}.go_gain.(fields{str})),'MAPprob_interpolated'))
                    explore_gain{str-4} = [explore_gain{str-4}; uni.Output{sbj}.go_gain.(fields{str}).MAPprob_interpolated'];
                else
                    explore_gain{str-4} = [explore_gain{str-4}; uni.Output{sbj}.go_gain.(fields{str}).MAPprobability'];
                end
            end

            if ~isnan(learn(sbj,2))
                if sum(contains(fieldnames(uni.Output{sbj}.avoid_loss.(fields{str})),'MAPprob_interpolated'))
                    explore_loss{str-4} = [explore_loss{str-4}; uni.Output{sbj}.avoid_loss.(fields{str}).MAPprob_interpolated'];
                else
                    explore_loss{str-4} = [explore_loss{str-4}; uni.Output{sbj}.avoid_loss.(fields{str}).MAPprobability'];
                end
            end
        end
        
end

%% Figure showing average +- SEM of 'Go Gain' strategy for 'Gain trials'
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize)
subplot(131)
y1 = mean(gogain,'omitnan') + std(gogain,1,'omitnan')./sqrt(size(gogain,1));
y2 = mean(gogain,'omitnan') - std(gogain,1,'omitnan')./sqrt(size(gogain,1));
h = patch([x fliplr(x)]', [y1 fliplr(y2)]',cmapPrior(1,:)); hold on
set(h,'facealpha',.3,'LineStyle','none')
plot(x,mean(gogain,'omitnan'),'-','Color',cmapPrior(1,:)); hold on
y1 = mean(avoidloss,'omitnan') + std(avoidloss,1,'omitnan')./sqrt(size(avoidloss,1));
y2 = mean(avoidloss,'omitnan') - std(avoidloss,1,'omitnan')./sqrt(size(avoidloss,1));
h = patch([x fliplr(x)], [y1 fliplr(y2)],cmapPrior(2,:)); hold on
set(h,'facealpha',.3,'LineStyle','none')
plot(x,mean(avoidloss,'omitnan'),'-','Color',cmapPrior(2,:)); hold on

plot([0 nTrials],[0.5 0.5],'k--'); hold on    
text(10,0.55,'Gain','Color',cmapPrior(1,:),'FontSize',fontsize);
text(10,0.35,'Loss','Color',cmapPrior(2,:),'FontSize',fontsize);
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
ylim([.2 .8])
title('Go Gain/Avoid Loss')
xlabel('Trial')
ylabel('Probability')

% Figure showing learning trial for gain and loss trials
jit = randn(size(learn,1),1)*.05;
siz = 3;
subplot(133)
boxplot(learn,'Colors',[.6 .6 .6]); hold on
plot(jit+1,learn(:,1),'o','Markersize',siz,'Color',cmapPrior(1,:)); hold on
plot(jit+2,learn(:,2),'o','Markersize',siz,'Color',cmapPrior(2,:)); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
ylabel('Learning trial')
xticks([1 2])
xticklabels({'Gain','Loss'})

subplot(132)
percnolearn(1) = sum(~isnan(learn(:,1)))./nsbj;
percnolearn(2) = sum(~isnan(learn(:,2)))/nsbj;

% estimate binomial confidence interval
N = 100;
[~,pci] = binofit(percnolearn(1),N);
[~,pci_loss] = binofit(percnolearn(2),N);

% plot the proportion of subject that didn't learn for each group
b = bar(percnolearn,'EdgeColor','none'); hold on
b.FaceColor = 'flat';
b.CData(1,:) = cmapPrior(1,:);
b.CData(2,:) = cmapPrior(2,:);
plot([1 1],pci+percnolearn(1),'-','Color',cmapPrior(1,:),'LineWidth',1.5); hold on
plot([2 2],pci_loss+percnolearn(2),'-','Color',cmapPrior(2,:),'LineWidth',1.5); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
% ylim([0 1])
ylabel('Prop of learning sbjs')
xticks([1 2])
xticklabels({'Gain','Loss'})

set(gcf,'PaperPositionMode','auto')
print([figpath 'fig6_group_learning'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'fig6_group_learning.png'])

%% plot probability of go gain/avoid loss for jeffreys prior and difference 
% between uniform and jeffreys

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize)
subplot(131)
y1 = mean(jef_gogain,'omitnan') + std(jef_gogain,1,'omitnan')./sqrt(size(gogain,1));
y2 = mean(jef_gogain,'omitnan') - std(jef_gogain,1,'omitnan')./sqrt(size(gogain,1));
h = patch([x fliplr(x)]', [y1 fliplr(y2)]',cmapPrior(1,:)); hold on
set(h,'facealpha',.3,'LineStyle','none')
plot(x,mean(jef_gogain,'omitnan'),'-','Color',cmapPrior(1,:)); hold on
y1 = mean(jef_avoidloss,'omitnan') + std(jef_avoidloss,1,'omitnan')./sqrt(size(avoidloss,1));
y2 = mean(jef_avoidloss,'omitnan') - std(jef_avoidloss,1,'omitnan')./sqrt(size(avoidloss,1));
h = patch([x fliplr(x)], [y1 fliplr(y2)],cmapPrior(2,:)); hold on
set(h,'facealpha',.3,'LineStyle','none')
plot(x,mean(jef_avoidloss,'omitnan'),'-','Color',cmapPrior(2,:)); hold on

plot([0 nTrials],[0.5 0.5],'k--'); hold on    
text(10,0.55,'Gain','Color',cmapPrior(1,:),'FontSize',fontsize);
text(10,0.35,'Loss','Color',cmapPrior(2,:),'FontSize',fontsize);
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
ylim([.2 .8])
title('Go Gain/Avoid Loss')
xlabel('Trial')
ylabel('Probability')

subplot(132)
for sbj = 1:nsbj
    plot(uni.Output{sbj}.go_gain.(fields{strategy}).MAPprobability-...
        jef.Output{sbj}.go_gain.(fields{strategy}).MAPprobability,'-','Color',cmapUnif(1,:)); hold on
end
text(2,0.1,'Uniform-Jeffreys','Color',cmapUnif(1,:),'FontSize',fontsize);
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
% ylabel({'Difference','Uniform-Jeffreys'})
title('Go Gain')
xlabel('Trial')

subplot(133)
for sbj = 1:nsbj
    plot(uni.Output{sbj}.avoid_loss.(fields{strategy}).MAPprobability-...
        jef.Output{sbj}.avoid_loss.(fields{strategy}).MAPprobability,'-','Color',cmapJeff(1,:)); hold on
end
text(2,0.1,'Uniform-Jeffreys','Color',cmapJeff(1,:),'FontSize',fontsize);
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
title('Avoid Loss')
xlabel('Trial')

set(gcf,'PaperPositionMode','auto')
print([figpath 'figS6_gain_loss_DeltaUniJeff'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'figS6_gain_loss_DeltaUniJeff.png'])

%% Plot explore startegies for the gain and loss trials for subjects that did learn
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[5 5 10 6.5])
for iS = 1:6
    subplot(2,3,iS)
    y1 = mean(explore_gain{iS},1,'omitnan') + std(explore_gain{iS},1,1,'omitnan')./sqrt(size(explore_gain{iS},1));
    y2 = mean(explore_gain{iS},1,'omitnan') - std(explore_gain{iS},1,1,'omitnan')./sqrt(size(explore_gain{iS},1));
    h = patch([x fliplr(x)], [y1 fliplr(y2)],cmapPrior(1,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x,mean(explore_gain{iS},'omitnan'),'-','Color',cmapPrior(1,:)); hold on
    y1 = mean(explore_loss{iS},'omitnan') + std(explore_loss{iS},1,1,'omitnan')./sqrt(size(explore_loss{iS},1));
    y2 = mean(explore_loss{iS},'omitnan') - std(explore_loss{iS},1,1,'omitnan')./sqrt(size(explore_loss{iS},1));
    h = patch([x fliplr(x)], [y1 fliplr(y2)],cmapPrior(2,:)); hold on
    set(h,'facealpha',.3,'LineStyle','none')
    plot(x,mean(explore_loss{iS},'omitnan'),'-','Color',cmapPrior(2,:)); hold on

    plot([0 nTrials],[0.5 0.5],'--','Color',[.6 .6 .6]); hold on
    set(gca,'FontName','Helvetica','FontSize',fontsize);
    set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
    if iS == 1
        text(15,0.32,'Gain','Color',cmapPrior(1,:),'FontSize',fontsize);
        text(15,0.2,'Loss','Color',cmapPrior(2,:),'FontSize',fontsize);
        ylabel('Probability')
    end
    if iS == 4
        ylabel('Probability')
    end

    title(humanStr_label{iS+4})
    if iS >3; xlabel('Trial'); end
    ylim([0 1])
end

set(gcf,'PaperPositionMode','auto')
print([figpath 'fig6_group_explore'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'fig6_group_explore.png'])

winningStr = NaN(nsbj,nTrials); 
winningStr_gain = NaN(nsbj,nTrials); 

% matrixMAPs = []; matrixPrecision = [];
for k = 1 :nsbj
        %% collect all MAPs and Precision for every subject during 'loss' trial
        matrixMAPs = []; matrixPrecision = [];
        for str = 1:3
            matrixMAPs = [matrixMAPs uni.Output{k}.avoid_loss.(fields{str}).MAPprobability];
            matrixPrecision = [matrixPrecision uni.Output{k}.avoid_loss.(fields{str}).precision];
        end
            
        maxval = max(matrixMAPs,[],2);
        [row,col] = find(matrixMAPs == maxval);

        maxvalue = max(matrixPrecision,[],2);   
        [prow,pcol] = find(matrixPrecision == maxvalue);

        %% Build a matrix where for every gamma and every trial the winning 
        % strategy is recorded as 1, otherwise zeros
        A = [[1:30]' ones(30,1)*2];
        B = [prow, pcol]; % coordinates for max precision strategies
        C = [row, col]; % coordinates for max MAPs strategies
        %% find where max precision and max MAP match the expected succesful strategy
        for t = 1:size(A,1)
            if ismember(A(t,:),B,'row') & ismember(A(t,:),C,'row')
                winningStr(k,A(t,1)) = 1;
            end 
        end

        %% collect all MAPs and Precision for every subject during 'gain' trial
        matrixMAPs = []; matrixPrecision = [];
        for str = 1:3
            matrixMAPs = [matrixMAPs uni.Output{k}.go_gain.(fields{str}).MAPprobability];
            matrixPrecision = [matrixPrecision uni.Output{k}.go_gain.(fields{str}).precision];
        end
        
        maxval = max(matrixMAPs,[],2);
        [row,col] = find(matrixMAPs == maxval);

        maxvalue = max(matrixPrecision,[],2);   
        [prow,pcol] = find(matrixPrecision == maxvalue);

        %% Build a matrix where for every gamma and every trial the winning 
        % strategy is recorded as 1, otherwise zeros
        A = [[1:30]' ones(30,1)*2];
        B = [prow, pcol]; % coordinates for max precision strategies
        C = [row, col]; % coordinates for max MAPs strategies
        %% find where max precision and max MAP match the expected succesful strategy
        for t = 1:size(A,1)
            if ismember(A(t,:),B,'row') & ismember(A(t,:),C,'row')
                winningStr_gain(k,A(t,1)) = 1;
            end 
        end
end

percWin(:,1) = sum(~isnan(winningStr_gain),2)./30;
percWin(:,2) = sum(~isnan(winningStr),2)./30;

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize)
subplot(131)
violinplot(bias,ones(length(bias),1),'EdgeColor',cmaptask(3,:),...
        'ViolinColor',cmaptask(3,:)); hold on
plot([.5 1.5], [.5 .5],'--','Color', [.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
% ylim([0 50])
ylabel('Prop of biased trial')
xticks(1)
xticklabels({''})

bias(isnan(bias))=[];

subplot(132)
brob = robustfit(bias,percWin(:,1)); 
rsq = corr(percWin(:,1),brob(1)+brob(2)*bias);
plot(bias,percWin(:,1),'o','MarkerFaceColor',cmapPrior(1,:),'MarkerEdgeColor','w'); hold on
plot(bias,brob(1)+brob(2)*bias,'-','Color',cmapPrior(1,:)); hold on
text(.6,.2,['r = ',num2str(round(rsq,2))],'Color',cmapPrior(1,:),'FontSize',fontsize);
plot([0 1],[.5 .5],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
ylabel('Prop winning strategy')
xlabel('Bias')
ylim([0 1])
xlim([0.5 1])

subplot(133)
brob = robustfit(bias,percWin(:,2)); 
rsq = corr(percWin(:,2),brob(1)+brob(2)*bias);
plot(bias,percWin(:,2),'o','MarkerFaceColor',cmapPrior(2,:),'MarkerEdgeColor','w'); hold on
plot(bias,brob(1)+brob(2)*bias,'-','Color',cmapPrior(2,:)); hold on
text(.6,.2,['r = ',num2str(round(rsq,2))],'Color',cmapPrior(2,:),'FontSize',fontsize);
plot([0 1],[.5 .5],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
% xlabel('Prop winning strategy')
xlabel('Bias')
ylim([0 1])
xlim([0.5 1])

set(gcf,'PaperPositionMode','auto')
print([figpath 'figS6_group_bias_winningStr'],'-depsc','-r0','-painters')
saveas(gcf,[figpath, 'figS6_group_bias_winningStr.png'])

