%% script to compute p(test=1|strategy=true) 
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
% Four possible true strategies: deterministic, conditionally
% deterministic, stochastic, conditionally stochastic
%
% Two possible test strategies: deterministic, conditionally deterministic
% (because stochastic cases are found by MAP match to them)
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
%
% [2/11/2022]: work stopped on this script, as getting too unwieldy. Split
% off deterministic-only variants into separate script. So not all
% solutions below are correct!
%
% Mark Humphries 18/10/2022

clearvars; close all;

n = 2;  % number of options in the choice task

p_test_1_give_true_x_deterministic = [0, 1/n, 1-1/n, 1];    % possible probability values for deterministic test
p_test_1_give_true_x_cond_deterministic = [1/n^2, (n-1)/n^2, (n^2 - 2*n + 1)/n^2];  % possible probability values for conditional deterministic test

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


%% 2. true = stochastic; test = deterministic

% no restriction: p(test=1|true=1) = 1- p(test=1|true=0) [when meet condition,
% Test strategy can either match true strategy choice, or it cannot; not
% both!]

probabilities.p_true_1 = linspace(1/n,0.9,5);       % range of possible true stochastic strategies
probabilities.p_true_0 = 1 - probabilities.p_true_1;
probabilities.p_true_null = zeros(1,numel(probabilities.p_true_1));

probabilities.p_test_1_given_true_1 = p_test_1_give_true_x_deterministic; 
probabilities.p_test_1_given_true_0 = probabilities.p_test_1_given_true_1;  
probabilities.p_test_1_given_true_null = 0;

P.true_stochastic_test_deterministic = compute_p_test_is_1_given_true_strategy(probabilities);
[Table_proportions.true_stochastic_test_deterministic,True_probabilities.true_stochastic_test_deterministic] = ...
    proportion_test_less_than_true(probabilities,P.true_stochastic_test_deterministic(:));

figure
ecdf(P.true_stochastic_test_deterministic(:));  hold on
stem(True_probabilities.true_stochastic_test_deterministic, ones(numel(True_probabilities.true_stochastic_test_deterministic),1), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'none', 'LineWidth', 1.5);
title('true: stochastic; test; deterministic')
xlabel('Probability of test=1')

%% 3. true = conditionally deterministic; test = deterministic

% not allowed: p(test=1|true=1) = p(test=1|true=null) = 1
% not allowed: p(true=1) = p(test=1|true=1) = 1; as are logically equivalent

probabilities.p_true_1 = [1/n 1-1/n 1];       % range of possible probabilites of true; assuming random selection models, or sequence (wins)
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

%% 4. true = conditional stochastic; test = deterministic

% extra set of conditional probabilities here, for p(execute | condition
% met) [ whereas this  = 1 for determinstic conditional]

% restriction: p(test=1|true=1) = 1- p(test=1|true=0) [when meet condition,
% Test strategy can either match true stratgy choice, or it cannot; not
% both!]

% not allowed: p(test=1|true=0) = p(test=1|true=null) = 1
 
p_execute = 0.5:0.1:0.9; % probability of executing conditional strategy when condition is met
p_condition_met = [1/n, 1-1/n, 1];
% these sum to 1 (all possible outcomes for P(true=X)
probabilities.p_true_1 = p_execute .* p_condition_met'; probabilities.p_true_1 = probabilities.p_true_1(:)'; % make row vector
probabilities.p_true_0 = (1-p_execute) .* p_condition_met'; probabilities.p_true_0 = probabilities.p_true_0(:)';
probabilities.p_true_null = repmat(1 - p_condition_met',1,numel(p_execute)); probabilities.p_true_null = probabilities.p_true_null(:)';

probabilities.p_test_1_given_true_1 = p_test_1_give_true_x_deterministic;
probabilities.p_test_1_given_true_0 = 1 - p_test_1_give_true_x_deterministic;  % invert probability: not assessed in current code!
probabilities.p_test_1_given_true_null = p_test_1_give_true_x_deterministic;

P.true_cond_stochastic_test_deterministic = compute_p_test_is_1_given_true_strategy(probabilities);
% strip conditions that are not allowed
P.true_cond_stochastic_test_deterministic(:,:,probabilities.p_test_1_given_true_0==1,probabilities.p_test_1_given_true_null==1)= nan;

[Table_proportions.true_cond_stochastic_test_deterministic, True_probabilities.true_cond_stochastic_test_deterministic] = ...
        proportion_test_less_than_true(probabilities,P.true_cond_stochastic_test_deterministic(:));

figure
ecdf(P.true_cond_stochastic_test_deterministic(:));  hold on
stem(True_probabilities.true_cond_stochastic_test_deterministic, ones(numel(True_probabilities.true_cond_stochastic_test_deterministic),1), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'none', 'LineWidth', 1.5);
title('true: conditional stochastic; test; deterministic')
xlabel('Probability of test=1')

%% 5. true = deterministic; test = conditional deterministic

% set of possible probabilities for conditional probability
probabilities.p_true_1 = 1;
probabilities.p_true_0 = 0;
probabilities.p_true_null = 0;

probabilities.p_test_1_given_true_1 = p_test_1_give_true_x_cond_deterministic;  
probabilities.p_test_1_given_true_null = 0;
probabilities.p_test_1_given_true_0 = 0;

P.true_deterministic_test_cond_deterministic = compute_p_test_is_1_given_true_strategy(probabilities);
[Table_proportions.true_deterministic_test_cond_deterministic,True_probabilities.true_deterministic_test_cond_deterministic] = ...
    proportion_test_less_than_true(probabilities,P.true_deterministic_test_cond_deterministic(:));

% plot
plot_discrete_distribution(P.true_deterministic_test_cond_deterministic);
title('true: deterministic; test: conditional deterministic')
ylabel('Probability')
xlabel('P(test=1)')
exportPPTfig(gcf,'ECDF_for_true_cond_and_test_cond.png',[pwd '\'],[10 15 10 10])


%% 6. true = stochastic; test = conditional deterministic
probabilities.p_true_1 = linspace(1/n,0.9,5);       % range of possible true stochastic strategies
probabilities.p_true_0 = 1 - probabilities.p_true_1;
probabilities.p_true_null = zeros(1,numel(probabilities.p_true_1));
probabilities.p_test_1_given_true_1 = p_test_1_give_true_x_cond_deterministic;  
probabilities.p_test_1_given_true_0 = 1 - probabilities.p_test_1_given_true_1;
probabilities.p_test_1_given_true_null = 0;

P.true_stochastic_test_cond_deterministic = compute_p_test_is_1_given_true_strategy(probabilities);
[Table_proportions.true_stochastic_test_cond_deterministic,True_probabilities.true_stochastic_test_cond_deterministic] = ...
        proportion_test_less_than_true(probabilities,P.true_stochastic_test_cond_deterministic(:));

figure
ecdf(P.true_stochastic_test_cond_deterministic(:));  hold on
stem(True_probabilities.true_stochastic_test_cond_deterministic, ones(numel(True_probabilities.true_stochastic_test_cond_deterministic),1), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'none', 'LineWidth', 1.5);
title('true: stochastic; test; conditionally deterministic')
xlabel('Probability of test=1')

%% 7. true = conditionally deterministic; test = conditionally deterministic
% not allowed: p(test=1|true=1) = p(test=1|true=null) = 1
probabilities.p_true_1 = [1/n 1-1/n 1];       % range of possible probabilites of true
probabilities.p_true_0 = zeros(1,numel(probabilities.p_true_1));
probabilities.p_true_null = 1 - probabilities.p_true_1;

probabilities.p_test_1_given_true_1 = p_test_1_give_true_x_cond_deterministic; 
probabilities.p_test_1_given_true_0 = 0;
probabilities.p_test_1_given_true_null = probabilities.p_test_1_given_true_1;

P.true_cond_deterministic_test_cond_deterministic = compute_p_test_is_1_given_true_strategy(probabilities);
% strip conditions that are not allowed
P.true_cond_deterministic_test_deterministic(:,probabilities.p_test_1_given_true_1==1,probabilities.p_test_1_given_true_null==1)= nan;

[Table_proportions.true_cond_deterministic_test_cond_deterministic,True_probabilities.true_cond_deterministic_test_cond_deterministic] = ...
        proportion_test_less_than_true(probabilities,P.true_cond_deterministic_test_cond_deterministic(:));

figure
ecdf(P.true_cond_deterministic_test_cond_deterministic(:)); hold on
stem(True_probabilities.true_cond_deterministic_test_cond_deterministic, ones(numel(True_probabilities.true_cond_deterministic_test_cond_deterministic),1), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'none', 'LineWidth', 1.5);
title('true: conditional deterministic; test; conditionally deterministic')
xlabel('Probability of test=1')
ylabel('Cumulative probability')
exportPPTfig(gcf,'ECDF_for_true_cond_and_test_cond.png',[pwd '\'],[10 15 10 10])

    
%% 8. true = conditional stochastic; test = conditional deterministic

% extra set of conditional probabilities here, for p(execute | condition
% met) [ whereas this  = 1 for determinstic conditional]

% restriction: p(test=1|true=1) = 1- p(test=1|true=0) [when meet condition,
% Test strategy can either match true stratgy choice, or it cannot; not
% both!]

% not allowed: p(test=1|true=0) = p(test=1|true=null) = 1


p_execute = 0.5:0.1:0.9; % probability of executing conditional strategy when condition is met
p_condition_met = [1/n, 1-1/n, 1];
probabilities.p_true_1 = p_execute .* p_condition_met'; probabilities.p_true_1 = probabilities.p_true_1(:)'; % make row vector
probabilities.p_true_0 = (1-p_execute) .* p_condition_met'; probabilities.p_true_0 = probabilities.p_true_0(:)';
probabilities.p_true_null = repmat(1 - p_condition_met',1,numel(p_execute)); probabilities.p_true_null = probabilities.p_true_null(:)';
probabilities.p_test_1_given_true_1 = p_test_1_give_true_x_cond_deterministic; 
probabilities.p_test_1_given_true_0 = 1 - p_test_1_give_true_x_cond_deterministic;
probabilities.p_test_1_given_true_null = p_test_1_give_true_x_cond_deterministic;

P.true_cond_stochastic_test_cond_deterministic = compute_p_test_is_1_given_true_strategy(probabilities);
% strip conditions that are not allowed
P.true_cond_stochastic_test_cond_deterministic(:,:,probabilities.p_test_1_given_true_0==1,probabilities.p_test_1_given_true_null==1)= nan;

[Table_proportions.true_cond_stochastic_test_cond_deterministic,True_probabilities.true_cond_stochastic_test_cond_deterministic] = ...
        proportion_test_less_than_true(probabilities,P.true_cond_stochastic_test_cond_deterministic(:));

figure
ecdf(P.true_cond_stochastic_test_cond_deterministic(:));  hold on
stem(True_probabilities.true_cond_stochastic_test_cond_deterministic, ones(numel(True_probabilities.true_cond_stochastic_test_cond_deterministic),1), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'none', 'LineWidth', 1.5);
title('true: conditional stochastic; test; conditional deterministic')
xlabel('Probability of test=1')


