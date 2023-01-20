clearvars; close all;

% load behavioural data for rats in Peyrache te al., 2009
% Nat Neurosci paper
% Header = {'Animal', 'Session', 'Direction', 'Reward', 'Rule', 'Light', 'SesName', 'Learning'};
load('Processed_data\SummaryDataTable_AllSessions.mat');
Data = double(Data);
Data(Data(:,4)==-1,4)=0;

subjects = unique(Data(:,1));
nRats = length(subjects);  % number of rats

% strategies to assess
strategies = ["go_right", "go_cued", "go_left", "go_uncued",...
                "win_stay_spatial","lose_shift_spatial","win_stay_cued","lose_shift_cued",...
                   "alternate","sticky"];

nstrategy = numel(strategies);

strategies_label = {'Go Right','Go Cued','Go Left','Go Uncued','Win-Stay-Spatial', ...
    'Lose-Shift-Spatial','Win-Stay-Cued','Lose-Shift-Cued', 'Alternate','Sticky'};

% visualise choices
cmapStrategy = [brewermap(4,'Set2'); brewermap(6,'Dark2')];
fontsize = 7;
axlinewidth = 0.5;
figpath = 'Figures\';

% load strategy profiles data for 4 rats from Peyarache et al., 2009
% dataset
uni = load('Processed_data/PeyracheDataStrategies.mat');
jef = load('Processed_data/PeyracheDataStrategies_Jeffreys.mat');

session = cell(nRats,1);
% Loop to identify each animal
for iR = 1:nRats
    % find index of each animal
    index = find(Data(:,1)==subjects(iR));
    % find index of session separation
    sessionsID = Data(index,2);
    session{iR} = find(diff(sessionsID) ~=0)+1;  % assumes rat ID is a string the XLS sheet, and so will be in cell array once loaded
end

barsize = [5 5 4 3.5];

%% plot beta distributions of rule strategies for example trial of Rat iR  
trial = 167; % example trial for posterior distributions
iR = 2;
x = 0:0.001:1;
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize)
%% plot example of prior distribution
% case 'Uniform'
alpha = 1; beta = 1;
y_uni = betapdf(x,alpha,beta);
% case 'Jeffreys'
alpha = 0.5; beta = 0.5;
y_jef = betapdf(x,alpha,beta);

% subplot(122)
plot(x,y_uni,'k'); hold on
plot(x,y_jef,'Color',[.6 .6 .6]); hold on
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
text(.1, 2,'Uniform prior','FontSize',fontsize)
text(.1, 7,'Jeffreys prior','Color', [.6 .6 .6],'FontSize',fontsize)
ylabel('PDF')
xlabel('P(strategy)')

print([figpath 'figS2_priors'],'-depsc')
saveas(gcf,[figpath, 'figS2_priors.png'])

%% plot MAPs for uniform and jeffreys prior
barsize = [5 5 11 3.5];
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

nstr = 3;
fields = fieldnames(uni.Output{iR});

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize)
for str = 1:nstr
    plot(jef.Output{iR}.(fields{str}).MAPprobability,'Color',cmapStrategy(str,:)); hold on
    text(length(jef.Output{iR}.(fields{str}).MAPprobability)+1,jef.Output{iR}.(fields{str}).MAPprobability(end),...
        strategies_label{str},'Color',cmapStrategy(str,:),'FontSize',fontsize)
end
for str = 1:nstr+1
    tr_id = find(Data(index,5)==str);
    text(tr_id,ones(size(tr_id))+0.1,'\bullet' ,'Color',cmapStrategy(str,:)); hold on
end
text(trial,max(jef.Output{iR}.(fields{str}).MAPprobability(trial))+0.2,'\downarrow'); hold on
plot([0 length(jef.Output{iR}.(fields{str}).MAPprobability)],[0.5 0.5],'--','Color',[.6 .6 .6]); hold on
plot([session{iR} session{iR}],[0 1],'--','Color',[.6 .6 .6]); hold on
xlim([0 length(reward)+50])
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
ylabel('Probability')
xlabel('Trial')

print([figpath 'figS2_MAP_Jefpriors'],'-depsc')
saveas(gcf,[figpath, 'figS2_MAP_Jefpriors.png'])

%% plot difference between Uniform and Jeffrey's prior
figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',barsize)
for str = 1:nstr
    plot(uni.Output{iR}.(fields{str}).MAPprobability-jef.Output{iR}.(fields{str}).MAPprobability,'Color',cmapStrategy(str,:)); hold on
    text(length(jef.Output{iR}.(fields{str}).MAPprobability)+1,(uni.Output{iR}.(fields{str}).MAPprobability(end)-jef.Output{iR}.(fields{str}).MAPprobability(end))+randn*0.02,...
        strategies_label{str},'Color',cmapStrategy(str,:),'FontSize',fontsize)
end
text(trial,max((uni.Output{iR}.(fields{str}).MAPprobability(trial)-jef.Output{iR}.(fields{str}).MAPprobability(trial)))+0.2,'\downarrow'); hold on
plot([session{iR} session{iR}],[-.1 .1],'--','Color',[.6 .6 .6]); hold on
xlim([0 length(reward)+50])
set(gca,'FontName','Helvetica','FontSize',fontsize);
set(gca,'Box','off','TickDir','out','LineWidth',axlinewidth);
ylabel({'Difference', 'Uniform - Jeffreys'})
xlabel('Trial')

print([figpath 'figS2_DeltaMAP_UniJef'],'-depsc')
saveas(gcf,[figpath, 'Fig1_SI_DeltaMAP_UniJef.png'])

