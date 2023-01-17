function [Table_proportions,true_probabilities] = proportion_test_less_than_true(probabilities,p_test_1_given_true)

% function to compute how many test strategy probabilities are lower than
% the corresponding true strategy detection probability. Returns Table with columns
% (p(true), proportions); and array of unique probabilities for True
% strategy
%
% Mark Humphries 18/10/2022

% calculate expected p(strategy|choices) for true strategy
% which is: p(true=1) / p(true=1)+p(true=0);
% for deterministic strategies, this is p(true=1); (as the denominator = 1)
% but for conditional strategies it isn't! Because p(true=null) does not
% contribute to the p(strategy|choices)

true_probabilities = unique(probabilities.p_true_1 ./ (probabilities.p_true_1 + probabilities.p_true_0));

n_true_probabilities = numel(true_probabilities);
n_tests = numel(p_test_1_given_true);

proportions = zeros(n_true_probabilities,1);
for iTrue = 1:n_true_probabilities
    proportions(iTrue) = sum(p_test_1_given_true < true_probabilities(iTrue)) / n_tests;
end

% make column vector
if size(true_probabilities,2) > 1 
    true_probabilities = true_probabilities'; end

Table_proportions = table(true_probabilities,proportions,'VariableNames',{'P_True_is_1','Proportion_test_strategies_less_than_true'});