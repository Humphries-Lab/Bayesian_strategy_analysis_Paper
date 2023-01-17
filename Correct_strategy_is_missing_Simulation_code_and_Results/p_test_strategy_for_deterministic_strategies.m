%% script to compute p(test=1|strategy=true) for deterministic True strategies (conditional or not)
% i.e. what is the expected probability that the tested strategy will be a
% success on each trial, given the true strategy being used by the agent
% where test != true. 
% 
% (And therefore what p(strategy_test | choices(1..t)) will be found by the
% algorithm at asymptote)
% 
% Answers the question: if true strategy is missing, how large do we expect
% the p(strategy_test) to be, and so how interpret whether it is true or
% not?.
%
% Two possible true strategies: deterministic, conditionally
% deterministic
%
% Two possible test strategies: deterministic, conditionally deterministic
%
% The set of possible probabilities for deterministic tests are always:
% {0, 1/n , 1-1/n, 1}
% mutually exclusive; randomly chose 1 option; randomly chose any but 1
% option; logically equivalent
%
% The set of possible probabilities for conditionally deterministic tests are always:
% {0, 1/n , 1-1/n, 1} x p(condition is met)
% mutually exclusive to true; randomly chose 1 option; randomly chose any but 1
% option; logically equivalent
%
% To compare the resulting probabilities, we omit all cases where the Test
% strategy is logically equivalent.
%
% Mark Humphries 2/11/2022

clearvars; close all;

n = 3;  % number of options in the choice task

p_test_1_give_true_x_deterministic = [0, 1/n, 1-1/n, 1];    % possible probability values for deterministic test
% p_test_1_give_true_x_cond_deterministic = [1/n^2, (n-1)/n^2, (n^2 - 2*n + 1)/n^2];  % possible probability values for conditional deterministic test

%% 1. true = deterministic; test = deterministic
% not allowed: p(true=1) = p(test=1|true=1) = 1; as are logically equivalent

% set of possible probabilities for conditional probability
probabilities.p_test_1_given_true_1 = p_test_1_give_true_x_deterministic;  
probabilities.p_true_1 = 1;
probabilities.p_test_1_given_true_0 = 0;
probabilities.p_true_0 = 0;
probabilities.p_test_1_given_true_null = 0;
probabilities.p_true_null = 0;

% calculate probabilities
P.true_deterministic_test_deterministic = compute_p_test_is_1_given_true_strategy(probabilities);

% remove dis-allowed cases
P.true_deterministic_test_deterministic(probabilities.p_true_1==1,probabilities.p_test_1_given_true_1==1)= nan;


[Table_proportions.true_deterministic_test_deterministic,True_probabilities.true_deterministic_test_deterministic] = ...
        proportion_test_less_than_true(probabilities,P.true_deterministic_test_deterministic(:));

% plot
plot_discrete_distribution(P.true_deterministic_test_deterministic(~isnan(P.true_deterministic_test_deterministic)));
title('true: deterministic; test; deterministic')

figure
ecdf(P.true_deterministic_test_deterministic(:));  hold on
stem(True_probabilities.true_deterministic_test_deterministic, ones(numel(True_probabilities.true_deterministic_test_deterministic),1), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'none', 'LineWidth', 1.5);
title('true: deterministic; test: conditional deterministic')
ylabel('Probability')
xlabel('P(test=1)')
exportPPTfig(gcf,'ECDF_for_true_deterministic_and_test_deterministic.png',[pwd '\'],[10 15 10 10])


%% 2. true = deterministic; test = conditional deterministic

% set of possible probabilities for conditional probability
probabilities.p_true_1 = 1;
probabilities.p_true_0 = 0;
probabilities.p_true_null = 0;

% p_condition_met_on_trial_t_minus_1 = [1/n, 1-1/n];      % model: choice is random variable
p_condition_met_on_trial_t_minus_1 = linspace(0,1,10);  % condition met is random variable, determined by True strategy and/or task
p_consistent_with_Test_on_trial_t = [0 1/n 1-1/n 1];      
probabilities.p_test_1_given_true_1 = p_condition_met_on_trial_t_minus_1 .* p_consistent_with_Test_on_trial_t'; probabilities.p_test_1_given_true_1 = probabilities.p_test_1_given_true_1(:)';
probabilities.p_test_1_given_true_null = 0;
probabilities.p_test_1_given_true_0 = 0;

% compute all probabilities
P.true_deterministic_test_cond_deterministic = compute_p_test_is_1_given_true_strategy(probabilities);

% remove dis-allowed cases: when Test is idential to True
P.true_deterministic_test_cond_deterministic(probabilities.p_true_1==1,probabilities.p_test_1_given_true_1==1)= nan;

[Table_proportions.true_deterministic_test_cond_deterministic,True_probabilities.true_deterministic_test_cond_deterministic] = ...
    proportion_test_less_than_true(probabilities,P.true_deterministic_test_cond_deterministic(:));

% plot
figure
ecdf(P.true_deterministic_test_cond_deterministic(:));  hold on
stem(True_probabilities.true_deterministic_test_cond_deterministic, ones(numel(True_probabilities.true_deterministic_test_cond_deterministic),1), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'none', 'LineWidth', 1.5);

title('true: deterministic; test: conditional deterministic')
ylabel('Probability')
xlabel('P(test=1)')
exportPPTfig(gcf,'ECDF_for_true_deterministic_and_test_cond.png',[pwd '\'],[10 15 10 10])

% without logical equivalence
p_consistent_with_Test_on_trial_t = [0 1/n 1-1/n];      
probabilities.p_test_1_given_true_1 = p_condition_met_on_trial_t_minus_1 .* p_consistent_with_Test_on_trial_t'; probabilities.p_test_1_given_true_1 = probabilities.p_test_1_given_true_1(:)';

% compute all probabilities
P.true_deterministic_test_cond_deterministic_no_equivalence = compute_p_test_is_1_given_true_strategy(probabilities);

% remove dis-allowed cases
P.true_deterministic_test_cond_deterministic_no_equivalence(probabilities.p_true_1==1,probabilities.p_test_1_given_true_1==1)= nan;

[Table_proportions.true_deterministic_test_cond_deterministic_no_equivalence,True_probabilities.true_deterministic_test_cond_deterministic_no_equivalence] = ...
    proportion_test_less_than_true(probabilities,P.true_deterministic_test_cond_deterministic_no_equivalence(:));

figure
ecdf(P.true_deterministic_test_cond_deterministic_no_equivalence(:));  hold on
stem(True_probabilities.true_deterministic_test_cond_deterministic_no_equivalence, ones(numel(True_probabilities.true_deterministic_test_cond_deterministic_no_equivalence),1), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'none', 'LineWidth', 1.5);

title('true: deterministic; test: conditional deterministic: no equivalence')
ylabel('Probability')
xlabel('P(test=1)')
exportPPTfig(gcf,'ECDF_for_true_deterministic_and_test_cond_no_equivalence.png',[pwd '\'],[10 15 10 10])

%% 3. true = conditionally deterministic; test = deterministic

% not allowed: p(test=1|true=1) = p(test=1|true=null) = 1
% not allowed: p(true=1) = p(test=1|true=1) = 1
probabilities.p_true_1 = linspace(0.1,1,10);       % range of possible probabilites of meeting condition to execute True explore strategy
probabilities.p_true_null = 1 - probabilities.p_true_1;
probabilities.p_true_0 = zeros(1,numel(probabilities.p_true_1));
probabilities.p_test_1_given_true_1 = p_test_1_give_true_x_deterministic;  
probabilities.p_test_1_given_true_0 = 0;
probabilities.p_test_1_given_true_null = probabilities.p_test_1_given_true_1;

% all combinations
P.true_cond_deterministic_test_deterministic = compute_p_test_is_1_given_true_strategy(probabilities);

% strip combinations that are not allowed
P.true_cond_deterministic_test_deterministic(:,probabilities.p_test_1_given_true_1==1,probabilities.p_test_1_given_true_null==1)= nan;
P.true_cond_deterministic_test_deterministic(probabilities.p_true_1==1,probabilities.p_test_1_given_true_1==1,:)= nan;

[Table_proportions.true_cond_deterministic_test_deterministic,True_probabilities.true_cond_deterministic_test_deterministic] = ...
        proportion_test_less_than_true(probabilities,P.true_cond_deterministic_test_deterministic(:));

figure
ecdf(P.true_cond_deterministic_test_deterministic(:));  hold on
stem(True_probabilities.true_cond_deterministic_test_deterministic, ones(numel(True_probabilities.true_cond_deterministic_test_deterministic),1), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'none', 'LineWidth', 1.5);
title('true: conditional deterministic; test; deterministic')
xlabel('Probability of test=1')
exportPPTfig(gcf,'ECDF_for_true_cond_and_test_deterministic.png',[pwd '\'],[10 15 10 10])


% strip combinations where either of the p(test=1|true=X) are exactly
% matched to strategy (e.g. True or A)
P.true_cond_deterministic_test_deterministic_no_equivalence = P.true_cond_deterministic_test_deterministic;
P.true_cond_deterministic_test_deterministic_no_equivalence(:,probabilities.p_test_1_given_true_1==1,:) = nan;
P.true_cond_deterministic_test_deterministic_no_equivalence(:,:,probabilities.p_test_1_given_true_null==1) = nan;

[Table_proportions.true_cond_deterministic_test_deterministic_no_equivalence,True_probabilities.true_cond_deterministic_test_deterministic_no_equivalence] = ...
        proportion_test_less_than_true(probabilities,P.true_cond_deterministic_test_deterministic_no_equivalence(:));

figure
ecdf(P.true_cond_deterministic_test_deterministic_no_equivalence(:));  hold on
stem(True_probabilities.true_cond_deterministic_test_deterministic_no_equivalence, ones(numel(True_probabilities.true_cond_deterministic_test_deterministic_no_equivalence),1), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'none', 'LineWidth', 1.5);
title('true: conditional deterministic; test; deterministic; no equivalence')
xlabel('Probability of test=1')
exportPPTfig(gcf,'ECDF_for_true_cond_and_test_deterministic_no_equivalence.png',[pwd '\'],[10 15 10 10])


%% 4. true = conditionally deterministic; test = conditionally deterministic
% not allowed: p(test=1|true=1) = p(test=1|true=null) = 1
probabilities.p_true_1 = linspace(0.1,1,10);       % range of possible probabilites of true
probabilities.p_true_0 = zeros(1,numel(probabilities.p_true_1));
probabilities.p_true_null = 1 - probabilities.p_true_1;

p_condition_met_on_trial_t_minus_1 = linspace(0,1,10);  % condition met is random variable, determined by True strategy and/or task
p_consistent_with_Test_on_trial_t = p_test_1_give_true_x_deterministic; 
probabilities.p_test_1_given_true_1 = p_condition_met_on_trial_t_minus_1 .* p_consistent_with_Test_on_trial_t'; probabilities.p_test_1_given_true_1 = probabilities.p_test_1_given_true_1(:)';
probabilities.p_test_1_given_true_0 = 0;
probabilities.p_test_1_given_true_null = probabilities.p_test_1_given_true_1;

P.true_cond_deterministic_test_cond_deterministic = compute_p_test_is_1_given_true_strategy(probabilities);
% strip conditions that are not allowed
P.true_cond_deterministic_test_cond_deterministic(:,probabilities.p_test_1_given_true_1==1,probabilities.p_test_1_given_true_null==1)= nan;
P.true_cond_deterministic_test_cond_deterministic(probabilities.p_true_1==1,probabilities.p_test_1_given_true_1==1,:)= nan;

[Table_proportions.true_cond_deterministic_test_cond_deterministic,True_probabilities.true_cond_deterministic_test_cond_deterministic] = ...
        proportion_test_less_than_true(probabilities,P.true_cond_deterministic_test_cond_deterministic(:));

figure
ecdf(P.true_cond_deterministic_test_cond_deterministic(:)); hold on
stem(True_probabilities.true_cond_deterministic_test_cond_deterministic, ones(numel(True_probabilities.true_cond_deterministic_test_cond_deterministic),1), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'none', 'LineWidth', 1.5);
title('true: conditional deterministic; test; conditionally deterministic')
xlabel('Probability of test=1')
ylabel('Cumulative probability')
exportPPTfig(gcf,'ECDF_for_true_cond_and_test_cond.png',[pwd '\'],[10 15 10 10])

% now with no logical equivalence between True and Test in either condition
p_consistent_with_Test_on_trial_t = [0 1/n 1-1/n]; 
probabilities.p_test_1_given_true_1 = p_condition_met_on_trial_t_minus_1 .* p_consistent_with_Test_on_trial_t'; probabilities.p_test_1_given_true_1 = probabilities.p_test_1_given_true_1(:)';
probabilities.p_test_1_given_true_null = probabilities.p_test_1_given_true_1;

P.true_cond_deterministic_test_cond_deterministic_no_equivalence = compute_p_test_is_1_given_true_strategy(probabilities);
% strip conditions that are not allowed
P.true_cond_deterministic_test_cond_deterministic_no_equivalence(:,probabilities.p_test_1_given_true_1==1,probabilities.p_test_1_given_true_null==1)= nan;
P.true_cond_deterministic_test_cond_deterministic_no_equivalence(probabilities.p_true_1==1,probabilities.p_test_1_given_true_1==1,:)= nan;

[Table_proportions.true_cond_deterministic_test_cond_deterministic_no_equivalence,True_probabilities.true_cond_deterministic_test_cond_deterministic_no_equivalence] = ...
        proportion_test_less_than_true(probabilities,P.true_cond_deterministic_test_cond_deterministic_no_equivalence(:));

figure
ecdf(P.true_cond_deterministic_test_cond_deterministic_no_equivalence(:)); hold on
stem(True_probabilities.true_cond_deterministic_test_cond_deterministic_no_equivalence, ones(numel(True_probabilities.true_cond_deterministic_test_cond_deterministic_no_equivalence),1), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'none', 'LineWidth', 1.5);
title('true: conditional deterministic; test; conditionally deterministic: no equivalence')
xlabel('Probability of test=1')
ylabel('Cumulative probability')
exportPPTfig(gcf,'ECDF_for_true_cond_and_test_cond_no_equivalence.png',[pwd '\'],[10 15 10 10])


%% save all combinations
save(['Results_for_Missing_Strategies_True_is_Deterministic_N_' num2str(n)],'Table_proportions','P','True_probabilities','n')

