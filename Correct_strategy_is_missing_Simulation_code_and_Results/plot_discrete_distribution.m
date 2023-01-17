function plot_discrete_distribution(discrete_values)

% pass an array of numbers that take discrete values
% and plot bar chart only at those values, into current figure

values = unique(discrete_values);
for iValue = 1:numel(values)
    frequency_of_values(iValue) = sum(values(iValue) == discrete_values);
end
figure
bar(values, frequency_of_values);
title('true: deterministic; test; deterministic')
set(gca,'XTick',values)
