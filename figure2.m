clearvars; close all;

% strategies to assess
strategies = ["go_right", "go_cued", "go_left", "go_uncued",...
                "win_stay_spatial","lose_shift_spatial","win_stay_cued","lose_shift_cued",...
                   "alternate","sticky"];

nstrategy = numel(strategies);

strategies_label = {'Go Right','Go Cued','Go Left','Go Uncued','Win-Stay-Spatial',...
    'Lose-Shift-Spatial','Win-Stay-Cued','Lose-Shift-Cued', 'Alternate','Sticky'};

% decay rate
% decay_rate = .9;

% visualise choices
cmapStrategy = [brewermap(4,'Paired'); brewermap(6,'Dark2')];
fontsize = 7;
axlinewidth = 0.5;
figpath = 'Figures\';

% load strategy profiles data for synthetic data
load('Processed_data/SyntheticDataStrategies_seed_VaryingDecay.mat')
load('Processed_data/SyntheticData_seed.mat') % Data contain 3 vector [light choice reward];

fields = fieldnames(Output{1});

% rule strategies assessed
% strRuleStrategy = {'Go Right','Go Cued','Go Left'};
rule = 'Cue'; % It can be direction rule or cue rule
ntrial = 500;
strategy = ["Go Right","Alternate","Lose-Shift-Cued","Go Cued","Lose-Shift-Spatial"];
% number of trials for each strategy
nstr = ntrial/length(strategy);

choice = Data(:,2);
choice(choice==0)=-1;
cue = Data(:,1)+1;
cue(cue==2)=-1;
Data(Data(:,3)==0,3)=-1;
ccue = [.6 .2 0]; % color for cue position
cumreward = cumsum(Data(:,3));
cumchoice = cumsum(choice);
cumcue = cumsum(cue);

% strategy fragmentation
% [M,I_rule] = max([Rat.Rule.MAPts],[],2);

barsize = [5 5 11 3.5];
barsize3sb = [5 5 11 3.5*3];
barsize6sb = [5 5 11*2 3.5*3];

indgm = find(decay_rate == .9);
%% plot the raw behavioral data (synthetic data) 
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize)
plot(cumreward,'-','Color',[.6 .6 .6]); hold on
plot(cumsum(choice),'k'); hold on
plot(cumsum(cue),'-','Color',ccue); hold on
plot([0:nstr:ntrial; 0:nstr:ntrial],[min(cumsum(cue)) max(cumreward)],'k--'); hold on
for str = 1:length(strategy)
    if strcmp(strategy{str},'Go Cue')
        text(30+nstr*(str-1), max(cumreward)+20, strategy{str},'FontSize',fontsize) 
    else
        text(5+nstr*(str-1), max(cumreward)+20, strategy{str},'FontSize',fontsize) 
    end
end
text(length(choice)+3,cumchoice(end)+8,'choice','FontSize',fontsize)
text(length(choice)+3,cumreward(end)-8,'reward','Color',[.6 .6 .6],'FontSize',fontsize)
text(length(choice)+3,cumcue(end),'cue','Color',ccue,'FontSize',fontsize)
text(length(choice)+50,max([cumchoice; cumreward])*1/4,'\downarrow','FontSize',20); hold on
text(length(choice)+50,max([cumchoice; cumreward])*3/4,'\uparrow','FontSize',20); hold on
text(length(choice)+75,max([cumchoice; cumreward])*1/4,'left','FontSize',fontsize); hold on
text(length(choice)+75,max([cumchoice; cumreward])*3/4,'right','FontSize',fontsize); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
ylabel('Cumulative distribution')
title('')
xlabel('Trials')
xlim([0 length(choice)+70])

print([figpath 'fig2_rawBehav'],'-depsc')
saveas(gcf,[figpath, 'fig2_rawBehav.png'])

%% plot the 5 strategies implemented
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize)
% sim_strategy = [Output{indgm}.Rule.MAPts(:,:,indgm) Output{indgm}.Explore.MAPts(:,:,indgm)];
% tested_strategies = [strRuleStrategy strExploreStrategy];
% cmap = brewermap(length(tested_strategies)-2,'Set1');
for str = 1:nstrategy
   if  ismember(strategies_label{str},strategy)
       if sum(contains(fieldnames(Output{indgm}.(fields{str})),'MAPprob_interpolated' ) ) %exist(Output{indgm}.(fields{str}).MAPprob_intepolated)==1
           plot(Output{indgm}.(fields{str}).MAPprob_interpolated,'-','Color',cmapStrategy(str,:)); hold on
           text(ntrial,Output{indgm}.(fields{str}).MAPprob_interpolated(end)+randn*0.07,strategies_label(str),...
                'Color',cmapStrategy(str,:),'FontSize',fontsize); hold on
       else
           plot(Output{indgm}.(fields{str}).MAPprobability,'-','Color',cmapStrategy(str,:)); hold on
           text(ntrial,Output{indgm}.(fields{str}).MAPprobability(end)+randn*0.07,strategies_label(str),...
                'Color',cmapStrategy(str,:),'FontSize',fontsize); hold on
       end
       
   end
end
plot([0:nstr:ntrial; 0:nstr:ntrial],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([0 ntrial],[0.5 0.5],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
ylim([-.01 1.01]) 
xlim([0 length(choice)+70])
ylabel('Probability')
xlabel('Trials')
print([figpath 'fig2_MAPs'],'-depsc')
saveas(gcf,[figpath, 'fig2_MAPs.png'])

%% plot precision for the strategies applyed by the agent
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize)

for str = 1:nstrategy
    precision = Output{indgm}.(fields{str}).precision;
   if  ismember(strategies_label{str},strategy)
       if sum(contains(fieldnames(Output{indgm}.(fields{str})),'precision_interpolated' ) ) 
           semilogy(1:ntrial,Output{indgm}.(fields{str}).precision_interpolated,'-',...
               'Color',cmapStrategy(str,:)); hold on
           text(ntrial,log(Output{indgm}.(fields{str}).precision_interpolated(end)),strategies_label{str},...
                'Color',cmapStrategy(str,:),'FontSize',fontsize); hold on 
       else
          semilogy(1:ntrial,precision,'-',...
               'Color',cmapStrategy(str,:)); hold on
          text(ntrial,log(precision(end)),strategies_label{str},...
               'Color',cmapStrategy(str,:),'FontSize',fontsize); hold on 
       end
   end
end
semilogy([0:nstr:ntrial; 0:nstr:ntrial],[min(precision(:)) 10^2.5],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
ylabel('log_{10}(Precision)')
xlabel('Trials')
xlim([0 length(choice)+70])
print([figpath 'fig2_Precision'],'-depsc')
saveas(gcf,[figpath, 'fig2_Precision.png'])

%% build a matrix with the implemented strategy
block = ntrial/length(strategy);
implementedStr = zeros(length(strategies_label),ntrial);
for st = 1:length(strategy)
    idx = find(strcmp(strategy(st), strategies_label));
    implementedStr(idx,1+(st-1)*block:block*st) = 1;
    
end


%% identify winning strategy combining between dominant strategy and 
% strategy with higher Precision for each gamma
% [precision,idx_precision] = max([Rat.Rule.Precisionts Rat.Explore.Precisionts],[],2);
[irow,icol] = find(implementedStr == 1);
winningStr = zeros(length(decay_rate),ntrial); 
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[5 5 10 3.5])
for g = 1:length(decay_rate)
    % g = indgm;
    matrixMAPs = []; matrixPrecision = [];
    for str = 1:nstrategy
        if sum(contains(fieldnames(Output{indgm}.(fields{str})),'MAPprob_interpolated' ) ) 
            matrixMAPs = [matrixMAPs; Output{g}.(fields{str}).MAPprob_interpolated'];
            matrixPrecision = [matrixPrecision; Output{g}.(fields{str}).precision_interpolated'];
        else
            matrixMAPs = [matrixMAPs; Output{g}.(fields{str}).MAPprobability'];
            matrixPrecision = [matrixPrecision; Output{g}.(fields{str}).precision'];
        end
    end
    maxval = max(matrixMAPs);
    [row,col] = find(matrixMAPs == maxval);

    % find indices with MAP<0.5 and change precision for those indices to
    % zero
    map = find(matrixMAPs<=0.5);
    matrixPrecision(map) = 0;
    maxvalue = max(matrixPrecision);   
    [prow,pcol] = find(matrixPrecision == maxvalue);

    %% Build a matrix where for every gamma and every trial the winning 
    % strategy is recorded as 1, otherwise zeros
    clear A B C
    A = [irow, icol]; % coordinates for implemented strategies
    B = [prow, pcol]; % coordinates for max precision strategies
    C = [row, col]; % coordinates for max MAPs strategies

    %% trying with for loop
    for t = 1:ntrial
        if ismember(A(t,:),B,'row') & ismember(A(t,:),C,'row')
%             winningStr(g,A(t,1)) = 1;
            winningStr(g,t) = 1;
        end       
    end

    if g == indgm
        plot(icol,irow,'+','Color',[1 .6 0]); hold on
        plot(pcol,prow,'o','Color',[.6 .6 .6],'MarkerSize',3); hold on %,'MarkerSize',3); hold on
        plot(col,row,'k.','MarkerSize',3); hold on
        plot([100:100:ntrial; 100:100:ntrial],[0 9],'--','Color',[.7 .7 .7]); hold on
        set(gca,'FontName','Helvetica','FontSize',fontsize);
        set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
        yticks(1:size(matrixMAPs,2))
        yticklabels(strategies_label)
        xlabel('Trial')
    end
    
    for str = 1:length(strategy)
        if strcmp(strategy{st},'Go Cue')
            text(30+nstr*(str-1), nstrategy+0.5, strategy{str},'FontSize',fontsize) 
        else
            text(5+nstr*(str-1), nstrategy+.5, strategy{str},'FontSize',fontsize) 
        end
    end
end
    
text(420,3,'Agent strategy','FontSize',fontsize,'Color',[1 .6 0]); hold on
text(420,2,'Max(Precision)','FontSize',fontsize,'Color',[.6 .6 .6]); hold on
text(420,1,'Max(Probability)','FontSize',fontsize); hold on

print([figpath 'fig2_WinningStrat'],'-depsc')
saveas(gcf,[figpath, 'fig2_WinningStrat.png'])
%%

successrate = sum(winningStr,2)./ntrial;
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[5 5 11 7])
subplot(1,6,1:5)
imagesc(winningStr); hold on
yticks(1:10:length(decay_rate))
yticklabels(decay_rate(1:10:end))
colormap(1-gray); %colorbar()
text(-50,indgm,'\rightarrow')
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
ylabel('Decay rate (\gamma)')
xlabel('Trial')
% print([figpath 'fig2_WinningStrat_Gamma'],'-depsc')
% saveas(gcf,[figpath, 'fig2_WinningStrat_Gamma.png'])

NtrialTosucces = (ntrial-sum(winningStr,2))/length(strategy);

idx = find(successrate == max(successrate),1,'first');
thresh = successrate(idx)*.95;
signval = find(successrate>=thresh);
x = [decay_rate(signval(1)) decay_rate(signval(end))];
% figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[5 5 1 7])
subplot(1,6,6)
curve1 = [0 0];
curve2 = [max(successrate) max(successrate)];
h = fill([curve1, fliplr(curve2)], [x, fliplr(x)],  [.6 .6 .6]); hold on
plot(successrate,decay_rate,'k'); hold on
% set(h,'facealpha',.3,'LineStyle','none')
hold on
% set(gca, 'XDir','reverse')
set(gca, 'YDir','reverse')
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
xlabel('Accuracy')
% ylabel('gamma')
% xlim([.2 .9])

print([figpath 'fig2_Accuracy_Gamma'],'-depsc')
saveas(gcf,[figpath, 'fig2_Accuracy_Gamma.png'])

