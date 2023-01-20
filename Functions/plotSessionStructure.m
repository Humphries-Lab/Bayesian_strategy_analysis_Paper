function plotSessionStructure(axes_handle,number_of_trials,new_session_trials,rule_change_trials,strRules,varargin)

% PLOTSESSIONSTRUCTURE plots structure of task to axes in figure
% PLOTSESSIONSTRUCTURE(AH,N_TRIALS,N_SESSIONS,RULE_CHANGE,RULES) will add
% the rule and session structure of the task to the set of axes specified
% by the handle AH, given:
% N_TRIALS: total number of trials for this subject to be plotted
% N_SESSIONS: the trial numbers for the start of each session
% RULE_CHANGE: the trial numbers for the rule switches (in order)
% RULES: a string array of the sequence of experienced rules 
% 
%
% PLOTSESSIONSTRUCTURE(...,'patch') will colour code sessions by coloured
% patches rather than vertical lines
%
% Mark Humphries

axes(axes_handle);  % make sure current axis is focus
blnPatch = 0;

if nargin >= 6 && strcmp(varargin{1},'patch')
    blnPatch = 1;
end

%% define ranges to plot things in
ymax = axes_handle.YLim(2);
yrange = axes_handle.YLim(2) - axes_handle.YLim(1);
yabove = ymax + 0.1*yrange;
ytext = ymax + (0.1*yrange/2);

%% sessions: plot as faces

face = 1:4;
fcolor = [0.95 0.95 1];

sessions = [1; new_session_trials; number_of_trials];
for iSh = 1:numel(sessions)-1
    if blnPatch
        % coloured patches on alternating sessions
        if rem(1+iSh,2)  % if an odd number
            vert = [sessions(iSh) 0; sessions(iSh+1) 0; sessions(iSh+1) ymax; sessions(iSh) ymax]; 
            patch('Faces',face,'Vertices',vert,'FaceColor',fcolor,'EdgeColor','none');
        end
    else
        % or just draw vertical line
        line([sessions(iSh) sessions(iSh)],[0 ymax],'Color',[0.7 0.7 0.7])
    end
end

%% rules: plot as coloured bars

changes = [1; rule_change_trials; number_of_trials]; % add first and last trials to create ends of patches
rcolors = repmat([0 0 0; 0.6 0.6 0.6],ceil(numel(changes)/2),1);

% % add error is Rat 1 is being plotted
% if Rat(iR).ID == 1 % mistake; 2nd "change" is within-session accident
%     line([Rat(iR).ruleChange(2) Rat(iR).ruleChange(2)],[0 1],'Color',[0.8 0.8 0.8]);
%     changes(3) = []; % delete from plotting array
% end

for iC = 1:numel(changes)-1
    % line([changes(iC) changes(iC+1)],[ymax ymax],'Color',rcolors(iC,:),'LineWidth',15,'Clipping','off')
    vert = [changes(iC) ymax; changes(iC+1) ymax; changes(iC+1) yabove; changes(iC) yabove];
    patch('Faces',face,'Vertices',vert,'FaceColor',rcolors(iC,:),'EdgeColor','none','Clipping','off');   
    
    
    % add label 
    % add "char" casting for old MATLAB versions where text() function
    % can't handle strings
    text(changes(iC)+5,ytext,char(strRules(iC)),'Color',[1 1 1],'Fontsize',9);
end

