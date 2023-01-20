function statistic = Summaries_of_Beta_distribution(alpha,beta,type,varargin)

% SUMMARIES_OF_BETA_DISTRIBUTION compute the requested summary statistic from the Beta PDF
% E = SUMMARIES_OF_BETA_DISTRIBUTION(a,b,T) given the Beta(a,b) parameters, computes the
% requested summary statistic T, returned as a scalar E
%   T:
%       'MAP': maximum a posteriori
%       'Mean': mean of beta distribution
%       'Var': variance
%       'Precision': 1/variance
%       'Percentile': Xth percentile of distribution, giving the Xth credible interval 
%       [define X as 4th argument: 0.95 = 95th percentile etc]; E will be 2-element array
%       [lower upper].
%
%  NOTE: if either of (a,b) are NaN, then the summary statistic T is
%  returned as a NaN
%
% 12/10/16: Initial version
% 31/03/22: testing MAP explicit calculation; and tidying up.
% 06/04/22: handle NaN values for alpha and beta
% Mark Humphries 

if isnan(alpha) || isnan(beta)
    statistic = nan;
else

    switch type
        case 'MAP'
            x = 0:0.001:1;
            y = betapdf(x,alpha,beta);
            [~,iy] = sort(y,'descend'); % sort values, then max will be first row of indices
            statistic = x(iy(1));  % maximum a posterioi estimate

        case 'Mean'
            statistic = alpha ./ (alpha + beta);
        case {'Var','Precision'}
            statistic = (alpha*beta) / ((alpha+beta).^2 .* (alpha+beta+1));
            if strcmp(type,'Precision')
                statistic = 1./statistic;
            end
        case 'Percentile'
            Prct = varargin{1};
            delta = 1-Prct;
            P = [delta/2 Prct + delta/2];
            statistic = betainv(P,alpha,beta);
        otherwise
            error([type ' is unknown estimator type'])
    end
end
