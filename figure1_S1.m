clearvars; close all;

% load summary behavioural data for rats in Peyrache et al., 2009
% Nat Neurosci paper
% Header = {'Animal', 'Session', 'Direction', 'Reward', 'Rule', 'Light', 'SesName', 'Learning'};
load('Processed_data\SummaryDataTable_AllSessions.mat');
Data = double(Data);
Data(Data(:,4)==-1,4)=0;

nRats = length(unique(Data(:,1)));  % number of rats
subjects = unique(Data(:,1));

% strategies to assess
strategies = ["go_right", "go_cued", "go_left", "go_uncued",...
                "win_stay_spatial","lose_shift_spatial","win_stay_cued","lose_shift_cued",...
                   "alternate","sticky"];

nstrategy = numel(strategies);

strategies_label = {'Go Right','Go Cued','Go Left','Go Uncued','Win-Stay-Spatial', ...
    'Lose-Shift-Spatial','Win-Stay-Cued','Lose-Shift-Cued', 'Alternate','Sticky'};

% visualise choices
cmapStrategy = [brewermap(4,'Set2'); brewermap(6,'Dark2')];
Clearn = [.9 .9 0]; % color for learning sessions
fontsize = 7;
axlinewidth = 0.5;
figpath = 'Figures\';

% load strategy profiles data for 4 rats from Peyarache et al., 2009
% dataset
load('Processed_data/PeyracheDataStrategies.mat')

session = cell(nRats,1);
% Loop to identify each animal
for iR = 1:nRats
    % find index of each animal
    index = find(Data(:,1)==subjects(iR));
    % find index of session separation
    sessionsID = Data(index,2);
    session{iR} = find(diff(sessionsID) ~=0)+1;  % assumes rat ID is a string the XLS sheet, and so will be in cell array once loaded
end


%% plot beta distributions of rule strategies for example trial of Rat iR  
trial = 247; % example trial for posterior distributions
iR = 2; % identify the example subject
fields = fieldnames(Output{iR});
x = 0:0.001:1;

barsize0 = [5 5 3.2 11];
barsize = [5 5 11 3.2];

%% example plot of the three main rule strategies
nstr = 3;
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize0)
subplot(311)
for str = 1:nstr      
    alpha = Output{iR}.(fields{str}).alpha(trial-2);
    beta = Output{iR}.(fields{str}).beta(trial-2);
    y = betapdf(x,alpha,beta);
    [~,iy] = sort(y,'descend'); % sort values, then max will be first row of indices
    estimator = x(iy(1));  % maximum a posterioi estimate
    plot(x,y,'Color',cmapStrategy(str,:)); hold on
    plot([estimator estimator],[0 y(iy(1))],'--','Color',cmapStrategy(str,:)); hold on
    set(gca,'FontName','Helvetica','FontSize',fontsize);
    set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
%     text(estimator,-.15,'MAP','Color',cmapStrategy(str,:),'FontSize',fontsize); 
    text(estimator-.1,y(iy(1))+.15,strategies_label{str},'Color',cmapStrategy(str,:),'FontSize',fontsize)
    ylabel('PDF')
    title('Trial t-2')
end

subplot(312)
for str = 1:nstr
    alpha = Output{iR}.(fields{str}).alpha(trial);
    beta = Output{iR}.(fields{str}).beta(trial);
    y = betapdf(x,alpha,beta);
    [~,iy] = sort(y,'descend'); % sort values, then max will be first row of indices
    estimator = x(iy(1));  % maximum a posterioi estimate
    plot(x,y,'Color',cmapStrategy(str,:)); hold on
    plot([estimator estimator],[0 y(iy(1))],'--','Color',cmapStrategy(str,:)); hold on
    set(gca,'FontName','Helvetica','FontSize',fontsize);
    set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
    ylabel('PDF')
    title('Trial t')
end

subplot(313)
for str = 1:nstr
    alpha = Output{iR}.(fields{str}).alpha(trial+2);
    beta = Output{iR}.(fields{str}).beta(trial+2);
    y = betapdf(x,alpha,beta);
    [~,iy] = sort(y,'descend'); % sort values, then max will be first row of indices
    estimator = x(iy(1));  % maximum a posterioi estimate
    plot(x,y,'Color',cmapStrategy(str,:)); hold on
    plot([estimator estimator],[0 y(iy(1))],'--','Color',cmapStrategy(str,:)); hold on
    set(gca,'FontName','Helvetica','FontSize',fontsize);
    set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
    ylabel('PDF')
    xlabel('P(strategy)')
    title('Trial t+2')
end

print([figpath 'fig1_examplesPDF'],'-depsc')
saveas(gcf,[figpath, 'fig1_examplesPDF.png'])
% exportfig(gcf,[figpath 'Fig1_examplesPDF'],'Color',color,'Format',format,'Resolution',dpi)

%% plot raw behavioural data
% compute cumulative distribution for reward (-1 for no reward, +1 for
% reward)
index = find(Data(:,1)==iR);
reward = Data(index,4);
reward(reward==0) = -1;
cumreward = cumsum(reward);
% compute cumulative distribution for choices (-1 for left, +1 for
% right)
choice = Data(index,3)+1; % choice is 0=right, 1=left
choice(choice==2) = -1; % I want right=1, left=-1
cumchoice = cumsum(choice);
cue = Data(index,6)+1; % light cue (0=right, 1=left)
cue(cue==2)=-1; % I want right=1, left=-1
cumcue = cumsum(cue);
ccue = [.6 .2 0]; % color for cue position

%% plot raw behavioural data
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize)
plot(cumchoice,'k'); hold on
plot(cumreward,'Color',[.6 .6 .6]); hold on
plot(cumsum(cue),'-','Color',ccue); hold on
text(length(reward)+1,cumchoice(end),'choice','FontSize',fontsize)
text(length(reward)+1,cumreward(end),'reward','Color',[.6 .6 .6],'FontSize',fontsize)
text(length(choice)+3,cumcue(end),'cue','Color',ccue,'FontSize',fontsize)
for str = 1:nstrategy 
    
    tr_id = find(Data(index,5)==str);
    text(tr_id,ones(size(tr_id))*(max([cumchoice; cumreward])+20),'\bullet' ,'Color',cmapStrategy(str,:)); hold on
end
text(trial,max([cumchoice; cumreward])-10,'\downarrow'); hold on
text(length(reward)+35,max([cumchoice; cumreward])*1/3,'\downarrow','FontSize',20); hold on
text(length(reward)+35,max([cumchoice; cumreward])*2/3,'\uparrow','FontSize',20); hold on
text(length(reward)+50,max([cumchoice; cumreward])*1/3,'left','FontSize',fontsize); hold on
text(length(reward)+50,max([cumchoice; cumreward])*2/3,'right','FontSize',fontsize); hold on
plot([session{iR} session{iR}],[min([cumchoice; cumreward]) max([cumchoice; cumreward])],'--','Color',[.6 .6 .6]); hold on
text(length(reward)+15,max([cumchoice; cumreward])+20,'rule','FontSize',fontsize); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
ylabel('Cumulative distribution') 
xlim([0 length(reward)+50])

print([figpath 'fig1_rawBehav'],'-depsc')
saveas(gcf,[figpath, 'fig1_rawBehav.png'])

%% plot MAP for rule strategy with gamma = 0.9
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize)
for str = 1:nstr
    plot(Output{iR}.(fields{str}).MAPprobability,'Color',cmapStrategy(str,:)); hold on
    text(length(Output{iR}.(fields{str}).MAPprobability)+1,Output{iR}.(fields{str}).MAPprobability(end),...
        strategies_label{str},'Color',cmapStrategy(str,:),'FontSize',fontsize)
end
text(trial,max(Output{iR}.(fields{str}).MAPprobability(trial))+0.2,'\downarrow'); hold on
plot([0 length(Output{iR}.(fields{str}).MAPprobability)],[0.5 0.5],'--','Color',[.6 .6 .6]); hold on
plot([session{iR} session{iR}],[0 1],'--','Color',[.6 .6 .6]); hold on
xlim([0 length(reward)+50])
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
ylabel('Probability')
xlabel('Trial')

print([figpath 'fig1_ruleMAP'],'-depsc')
saveas(gcf,[figpath, 'fig1_ruleMAP.png'])

%% plot confidence in estimation (1/var)
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize)
nTrials = length(Output{iR}.(fields{str}).precision);
for str = 1:nstr+1
    lin(str) = semilogy(1:nTrials,Output{iR}.(fields{str}).precision,'.-',...
        'Color',cmapStrategy(str,:)); hold on
    text(length(Output{iR}.(fields{str}).precision)+1,Output{iR}.(fields{str}).precision(end),...
        strategies_label{str},'Color',cmapStrategy(str,:),'FontSize',fontsize)
end
text(trial,max(Output{iR}.(fields{str}).precision(trial))+70,'\downarrow'); hold on
semilogy([session{iR} session{iR}],[min(Output{iR}.(fields{str}).precision) ...
    max(Output{iR}.(fields{str}).precision)],'--','Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
legend(lin(3:4),'Go Left & Go Right','Go Cued & Go Uncued','Location','southeast')
xlabel('Trial')
ylabel('log_{10}(Precision)')
xlim([0 length(reward)+50])

print([figpath 'fig1_precision'],'-depsc')
saveas(gcf,[figpath, 'fig1_precision.png'])

%% plot explore strategies
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize)

for str = [5 6 8] %1:nExploreStrategy 
    plot(Output{iR}.(fields{str}).MAPprob_interpolated,'Color',cmapStrategy(str,:)); hold on
    text(length(Output{iR}.(fields{str}).MAPprob_interpolated)+1,Output{iR}.(fields{str}).MAPprob_interpolated(end)+randn*0.1,...
        strategies_label{str},'Color',cmapStrategy(str,:),'FontSize',fontsize)
end
text(trial,max(Output{iR}.(fields{str}).MAPprob_interpolated(trial))+0.6,'\downarrow'); hold on
plot([session{iR} session{iR}],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([0 length(reward)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
% text(length(reward)+15,1,'rule','FontSize',fontsize); hold on
xlim([0 length(reward)+50])
% ylim([0.4 1])
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
ylabel('Probability')
xlabel('Trial')

print([figpath 'fig1_exploreMAP'],'-depsc')
saveas(gcf,[figpath, 'fig1_exploreMAP.png'])
%% #########################################

% load strategy profiles data for 4 rats from Peyarache et al., 2009
% dataset
load('Processed_data/PeyracheDataStrategies_NoDecay.mat')

% plot MAP for Strategy Analysis without no decay
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize)

for str = 1:nstr
    plot(Output{iR}.(fields{str}).MAPprobability,'Color',cmapStrategy(str,:)); hold on
    text(length(Output{iR}.(fields{str}).MAPprobability)+1,Output{iR}.(fields{str}).MAPprobability(end),...
        strategies_label{str},'Color',cmapStrategy(str,:),'FontSize',fontsize)
end
text(trial,max(Output{iR}.(fields{str}).MAPprobability(trial))+0.5,'\downarrow'); hold on
plot([session{iR} session{iR}],[0 1],'--','Color',[.6 .6 .6]); hold on
plot([0 length(reward)],[.5 .5],'--','Color',[.6 .6 .6]); hold on
xlim([0 length(reward)+50])
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
ylabel('Probability')
xlabel('Trial')

print([figpath 'figS1_ruleMAP_nodecay'],'-depsc')
saveas(gcf,[figpath, 'figS1_ruleMAP_nodecay.png'])

