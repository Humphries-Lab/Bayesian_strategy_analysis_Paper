%% script to plot SI Figure of what p(match) is expected to be
%%
clearvars; close all

run figure_properties.m

if ispc
    simOutputPath = 'C:\Users\lpzmdh\Dropbox\Projects\Bayesian strategy analysis\Correct strategy is missing\';
else
    simOutputPath = '/Users/mqbssmhg/Dropbox/Projects/Bayesian strategy analysis/Correct strategy is missing/';
end

% load([simOutputPath 'Results_for_Missing_Strategies_True_is_Deterministic_N_2']);
load([simOutputPath 'Results_for_Missing_Strategies_True_is_Deterministic_N_3']);

%% panel A: distribution of p(match) when True & Test strategy are deterministic

plot_cdfs(P.true_deterministic_test_deterministic(:),[],n);
print([exportpath 'Test_deterministic_True_deterministic_P_match_N_' num2str(n)],'-dpng')
print([exportpath 'Test_deterministic_True_deterministic_P_match_N_' num2str(n)'],'-dsvg')

%% panel B: distribution of p(match) when True is conditional and Test strategy is deterministic

plot_cdfs(P.true_cond_deterministic_test_deterministic(:),P.true_cond_deterministic_test_deterministic_no_equivalence(:),n);
print([exportpath 'Test_deterministic_True_conditional_P_match_N_' num2str(n)],'-dpng')
print([exportpath 'Test_deterministic_True_conditional_P_match_N_' num2str(n)'],'-dsvg')

%% panel C: combined set for when Test strategy is deterministic

% combine cases where True is deterministic and conditional: we don't know what True is!
combined_Test_deterministic = [P.true_deterministic_test_deterministic(:); P.true_cond_deterministic_test_deterministic(:)];

% removing cases where Test is logically equivalent to True strategy in one condition
combined_Test_deterministic_without_equivalence = [P.true_deterministic_test_deterministic(:); P.true_cond_deterministic_test_deterministic_no_equivalence(:)];

% plot
plot_cdfs(combined_Test_deterministic,combined_Test_deterministic_without_equivalence,n);

print([exportpath 'Test_deterministic_P_match_N_' num2str(n)],'-dpng')
print([exportpath 'Test_deterministic_P_match_N_' num2str(n)'],'-dsvg')

%% panel D: Test conditional, True deterministic 
plot_cdfs(P.true_deterministic_test_cond_deterministic(:),P.true_deterministic_test_cond_deterministic_no_equivalence,n);
print([exportpath 'Test_conditional_True_deterministic_P_match_N_' num2str(n)],'-dpng')
print([exportpath 'Test_conditional_True_deterministic_P_match_N_' num2str(n)'],'-dsvg')

%% panel E: Test conditional, True conditional

plot_cdfs(P.true_cond_deterministic_test_cond_deterministic(:),P.true_cond_deterministic_test_cond_deterministic_no_equivalence(:),n);
print([exportpath 'Test_conditional_True_conditional_P_match_N_' num2str(n)],'-dpng')
print([exportpath 'Test_conditional_True_conditional_P_match_N_' num2str(n)'],'-dsvg')

%% panel F: combined set for when Test strategy is conditional

% combine cases where True is deterministic and conditional: we don't know what True is!
combined_Test_conditional = [P.true_deterministic_test_cond_deterministic(:); P.true_cond_deterministic_test_cond_deterministic(:)];

% removing cases where Test is logically equivalent to True strategy in one condition
combined_Test_conditional_without_equivalence = [P.true_deterministic_test_cond_deterministic_no_equivalence(:); P.true_cond_deterministic_test_cond_deterministic_no_equivalence(:)];

% plot
plot_cdfs(combined_Test_conditional,combined_Test_conditional_without_equivalence,n);

print([exportpath 'Test_conditional_P_match_N_' num2str(n)'],'-dpng')
print([exportpath 'Test_conditional_P_match_N_' num2str(n)'],'-dsvg')

%% common function for plotting panels
function plot_cdfs(p_for_all_cases,p_for_cases_with_omissions,n)

run figure_properties.m

% compute ECDFs
[F_cdf_all,X_cdf_all] = ecdf(p_for_all_cases);
if ~isempty(p_for_cases_with_omissions)
    [F_cdf_omitted,X_cdf_omitted] = ecdf(p_for_cases_with_omissions);
end

figure('Units', 'centimeters', 'PaperPositionMode', 'auto','Position',[10 15 figsize.square]); 
line([1-1/n,1-1/n],[0,1],'Color',colourmaps.true_colour,'LineStyle','--'); hold on
stairs(X_cdf_all,F_cdf_all,'Color',colourmaps.cdf_all) ; hold on
if ~isempty(p_for_cases_with_omissions)
    stairs(X_cdf_omitted,F_cdf_omitted,'Color',colourmaps.cdf_no_equivalence)    
end
line([1,1],[0,1],'Color',colourmaps.cdf_true)
ylabel('Probability')
xlabel('P(match)')

FormatFig_For_Export(gcf,fontsize,fontname,widths.axis);


end


