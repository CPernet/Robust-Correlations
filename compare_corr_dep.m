function [R1,R2,D,CI,DIST] = compare_corr_dep(varargin)
%COMPARE_CORR_DEP Compare correlations between two dependent groups using
%the percentile bootstrap
%
% Two cases are considered: overlapping and non-overlapping correlations.
% In the overlapping case, 3 measurements are made in the same participants
% and we compute correlations between variable 1 and variables 2 and 3. For
% instance, variable 1 could be a behavioural measurement and variables 2
% and 3 activity from 2 brain areas.
% In the non-overlapping case, 2 measurements are made in the same
% participants on 2 separate occasions, and for each occasion a correlation
% is computed between variables. For instance, this could correspond to
% measurements made before and after a training procedure.
%
% The percentile bootstrap does not give satisfactory results when comparing Pearsons' correlations.
% For this special case, there is R code available from
% https://dornsife.usc.edu/labs/rwilcox/software/. Download the file
% Rallfun-vXX.txt. For the overlapping case, use the function TWOpov.
% For the non-overlapping case, use the function TWOpNOV.
%
% References:
% Wilcox, R.R. (2009) 
% Comparing Pearson Correlations: Dealing with Heteroscedasticity and Nonnormality. 
% Commun Stat-Simul C, 38, 2220-2234.
%
% Wilcox, R.R. (2016) 
% Comparing dependent robust correlations. 
% Brit J Math Stat Psy, 69, 215-224.
%
% FORMAT: [R1,R2,D,CI,DIST] = compare_correlations(data, method)
%
% INPUT: data = n*3 matrix in the overlap case (column 1 is correlated with columns 2 and 3);
%               n*4 matrix in the non-overlap case (column 1 is correlated
%               with column 2 and column 3 is correlated with column 4)
%        method = type of correlation to use 'Spearman' 'bendcorr'
%        'Skipped_P' or 'Skipped_S' (for skipped corr Pearson or Spearman)
%
% OUTPUT: R1 = correlation 1
%         R2 = correlation 2
%         D = the observed difference
%         CI = 95% of the observed difference
%         DIST = distribution of bootstrap differences
%
% See also COMPARE_CORR_IND

% Copyright 2013-2018 Corr_toolbox team
% -------------------------------------
% Cyril Pernet v1 08 April 2013
% Guillaume Rousselet v2 2018-06-21

%% input checks

if nargin<2
    error('not enough input arguments');
end

names{1} = 'Spearman';
names{2} = 'bendcorr';
names{3} = 'Skipped_P';
names{4} = 'Skipped_S';

data = varargin{1};
if size(data,2) == 3
    type = 1;
elseif size(data,2) == 4
    type = 2;
else
    error('data should have 3 or 4 columns')
end

method = varargin{2};
if isempty(cell2mat(strfind(names,method)))
    error('unknown method - choose among: ''Spearman'' ''bendcorr'' ''Skipped_P'' ''Skipped_S''');
end

%% bootstrap
nboot = 600;
low = round((5/100*nboot)/2);
high = nboot - low;

switch type
    case {1}  % overlap case
        % boostrap data
        table = randi(size(data,1), size(data,1),599);
        X1 = data(:,1); X1 = X1(table);
        Y1 = data(:,2); Y1 = Y1(table);
        Y2 = data(:,3); Y2 = Y2(table);
                   
        if strcmp(method,names{1})
            R1 = Spearman(data(:,1),data(:,2),0);
            R2 = Spearman(data(:,1),data(:,3),0);
            D = R1 - R2;
            boot_r1 = Spearman(X1,Y1,0);
            boot_r2 = Spearman(X1,Y2,0);
            
        elseif strcmp(method,names{2})
            R1 = bendcorr(data(:,1),data(:,2),0);
            R2 = bendcorr(data(:,1),data(:,3),0);
            D = R1 - R2;
            boot_r1 = bendcorr(X1,Y1,0);
            boot_r2 = bendcorr(X1,Y2,0);
            
        else
            if strcmp(method,names{3}), estimator = 'Pearson';
            else
                estimator = 'Spearman'; 
            end
            R1 = skipped_correlation(data(:,1),data(:,2),0,estimator);
            R2 = skipped_correlation(data(:,1),data(:,3),0,estimator);
            D = R1 - R2;
            boot_r1 = skipped_correlation(X1,Y1,0,estimator);
            boot_r2 = skipped_correlation(X1,Y2,0,estimator);
        end
        
    case {2}  % non-overlap case
        % boostrap data
        n1 = size(data,1); table1 = randi(n1,n1,599);
        X1 = data(:,1); X1 = X1(table1);
        Y1 = data(:,2); Y1 = Y1(table1);
        X2 = data(:,3); X2 = X2(table1);
        Y2 = data(:,4); Y2 = Y2(table1);
                 
        if strcmp(method,names{1})
            R1 = Spearman(data(:,1),data(:,2),0);
            R2 = Spearman(data(:,3),data(:,4),0);
            D = R1 - R2;
            boot_r1 = Spearman(X1,Y1,0);
            boot_r2 = Spearman(X2,Y2,0);
        elseif strcmp(method,names{2})
            R1 = bendcorr(data(:,1),data(:,2),0);
            R2 = bendcorr(data(:,3),data(:,4),0);
            D = R1 - R2;
            boot_r1 = bendcorr(X1,Y1,0);
            boot_r2 = bendcorr(X2,Y2,0);
        else
            if strcmp(method,names{3}), estimator = 'Pearson';
            else
                estimator = 'Spearman'; 
            end
            R1 = skipped_correlation(data(:,1),data(:,2),0,estimator);
            R2 = skipped_correlation(data(:,3),data(:,4),0,estimator);
            D = R1 - R2;
            boot_r1 = skipped_correlation(X1,Y1,0,estimator);
            boot_r2 = skipped_correlation(X2,Y2,0,estimator);
        end
end

d = sort(boot_r1-boot_r2);
CI = [d(low) d(high)];
DIST = boot_r1 - boot_r2;




