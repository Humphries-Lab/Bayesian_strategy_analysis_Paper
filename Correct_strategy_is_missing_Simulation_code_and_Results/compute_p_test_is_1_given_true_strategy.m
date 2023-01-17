function p = compute_p_test_is_1_given_true_strategy(probabilities)

% function to compute the probability given each conditional probability:
% p(strategy_test = 1 | strategy_true) = 
%   p(strategy_test = 1 | true = 1)p(true=1) +
%   p(strategy_test = 1 | true = 0)p(true=0) +
%   p(strategy_test = 1 | true = null)p(true=null) +
%
% Pass "probabilities" as struct with fields: 
% Each field is a ROW vector
%
% This function computes all possible combinations of probabilities of p(test=1|true=X).
% Keeps each probability for p(true=X) as summing to 1
% Thus, returned p has full dimensions:
% #true conditions x #p(test=1|true=1) x #p(test=1|true=0) x
% #p(test=1|true=null)
% 
% Any case where p(test=1|true=X) is singleton is removed;
%
% Mark Humphries 18/10/22


% individual product terms - as many rows as p(true) are specified
probs_true_1 = probabilities.p_test_1_given_true_1 .* probabilities.p_true_1';
probs_true_0 = probabilities.p_test_1_given_true_0 .* probabilities.p_true_0';
probs_true_null = probabilities.p_test_1_given_true_null .* probabilities.p_true_null';

% keyboard

% all summations along rows: for each set of p(true=1), p(true=0), and
% p(true=null), find all combinations
p = zeros(numel(probabilities.p_true_1),numel(probabilities.p_test_1_given_true_1),numel(probabilities.p_test_1_given_true_0),numel(probabilities.p_test_1_given_true_null));
for iTrueSet = 1:numel(probabilities.p_true_1)
    % for each row
    for i1 = 1:numel(probabilities.p_test_1_given_true_1)
        for i0 = 1:numel(probabilities.p_test_1_given_true_0)
            for iNull = 1:numel(probabilities.p_test_1_given_true_null)
                p(iTrueSet,i1,i0,iNull) = probs_true_1(iTrueSet,i1) + probs_true_0(iTrueSet,i0) + probs_true_null(iTrueSet,iNull);
            end
        end
    end
end
p = squeeze(p); % get rid of singleton dimensions
    