function [R1,R2,D,CI,DIST] = compare_corr_ind(varargin)
%COMPARE_CORR_IND Compare correlations between two independent groups using
%the percentile bootstrap.
%
% Uses the method described in Wilcox, R.R. (2012)
% Introduction to robust estimation and hypothesis testing.
% Academic Press, San Diego, CA.
% pages 467-468.
%
% For each group, generate 599 independent bootstrap samples. For each
% bootstrap sample, compute the correlation coefficient. Save the
% difference between the coefficients for group 1 and group 2. Use this
% distribution of bootstrap differences to define a 95% confidence interval for
% the difference. For Pearson's correlation, the bounds of the confidence interval
% are adjusted depending on sample size. For other techniques, there is no
% adjustement, i.e. the bounds are values 15th and 585th.
%
% FORMAT: [R1,R2,D,CI,DIST] = compare_correlations(data1, data2, method)
%
% INPUT: data1 = data from group 1 - matrix n1*2
%        data2 = data from group 2 - matrix n2*2
%        method = type of correlation to use 'Pearson' 'Spearman' 'bendcorr'
%        'Skipped_P' or 'Skipped_S' (for skipped corr Pearson or Spearman)
%
% OUTPUT: R1 = correlation for group 1
%         R2 = correlation for group 2
%         D = the observed difference
%         CI = 95% of the observed difference
%         DIST = distribution of bootstrap differences
%
% Example: two groups, one with a weak correlation, 
%          one with a strong correlation
%  rng(21)
%  n1 = 100;
%  mu1 = 1;
%  sigma1 = 2;
%  r1 = 0.1;
%  data1 = (mu1 + sigma1*randn(n1,2)) * chol([1 r1; r1 1]);
%  n2 = 70;
%  mu2 = 3;
%  sigma2 = 1;
%  r2 = 0.5;
%  data2 = (mu2 + sigma2*randn(n2,2)) * chol([1 r2; r2 1]);
%  [R1,R2,D,CI,DIST] = compare_corr_ind(data1, data2, 'Spearman');
%
% See also COMPARE_CORR_DEP

% Copyright 2013-2018 Corr_toolbox team
% -------------------------------------
% Cyril Pernet v1 08 April 2013
% Guillaume Rousselet v2 2018-06-21

%% input checks

if nargin<3
    error('not enough input arguments');
end

names{1} = 'Pearson';
names{2} = 'Spearman';
names{3} = 'bendcorr';
names{4} = 'Skipped_P';
names{5} = 'Skipped_S';

if nargin == 3
    type = 2;
else
    type = varargin{4};
end

if type == 1
    if sum((size(varargin{1})==size(varargin{2}))) ~=2
        error('for paired data, data1 and data2 must have the same size')
    end
end
data1 = varargin{1};
data2 = varargin{2};

method = varargin{3};
if isempty(cell2mat(strfind(names,method)))
    error('unknown method - choose among: ''Pearson'' ''Spearman'' ''bendcorr'' ''Skipped_P'' ''Skipped_S''');
end

%% bootstrap
nboot = 600;
low = round((5/100*nboot)/2);
high = nboot - low;

% boostrap data
n1 = size(data1,1); table1 = randi(n1,n1,599);
X1 = data1(:,1); X1 = X1(table1);
Y1 = data1(:,2); Y1 = Y1(table1);

n2 = size(data2,1); table2 = randi(n2,n2,599);
X2 = data2(:,1); X2 = X2(table2);
Y2 = data2(:,2); Y2 = Y2(table2);

if strcmp(method,names{1})
    R1 = Pearson(data1(:,1),data1(:,2),0);
    R2 = Pearson(data2(:,1),data2(:,2),0);
    D = R1 - R2;
    boot_r1 = Pearson(X1,Y1,0);
    boot_r2 = Pearson(X2,Y2,0);
    
    % adjust percentile following Wilcox 2012
    N  = n1+n2;
    if N<40
        low = 7; high = 593;
    elseif N>=40 && N<80
        low = 8; high = 592;
    elseif N>=80 && N<180
        low = 11; high = 588;
    elseif N>=180 && N<250
        low = 14; high = 585;
    elseif N>=250
        low = 15; high = 584;
    end
    
elseif strcmp(method,names{2})
    R1 = Spearman(data1(:,1),data1(:,2),0);
    R2 = Spearman(data2(:,1),data2(:,2),0);
    D = R1 - R2;
    boot_r1 = Spearman(X1,Y1,0);
    boot_r2 = Spearman(X2,Y2,0);
    
elseif strcmp(method,names{3})
    R1 = bendcorr(data1(:,1),data1(:,2),0);
    R2 = bendcorr(data2(:,1),data2(:,2),0);
    D = R1 - R2;
    boot_r1 = bendcorr(X1,Y1,0);
    boot_r2 = bendcorr(X2,Y2,0);
    
else
    if strcmp(method,names{4}), estimator = 'Pearson';
    else
        estimator = 'Spearman';
    end
    R1 = skipped_correlation(data1(:,1),data1(:,2),0,estimator);
    R2 = skipped_correlation(data2(:,1),data2(:,2),0,estimator);
    D =  R1 - R2;
    boot_r1 = skipped_correlation(X1,Y1,0,estimator);
    boot_r2 = skipped_correlation(X2,Y2,0,estimator);
end

d = sort(boot_r1 - boot_r2);
CI = [d(low) d(high)];
DIST = boot_r1 - boot_r2;




