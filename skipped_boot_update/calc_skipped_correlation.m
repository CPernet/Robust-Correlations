function [r,t,h,outid,a,b,coef, column]=calc_skipped_correlation(x,y, hypothesis)

% performs a robust correlation using pearson/spearman correlation on
% data cleaned up for bivariate outliers - that is after finding the
% central point in the distribution using the mid covariance determinant,
% orthogonal distances are computed to this point, and any data outside the
% bound defined by the idealf estimator of the interquartile range is removed.
% 
% FORMAT:
%          [r,t,h] = skipped_correlation(X,Y);
%          [r,t,h] = skipped_correlation(X,Y,fig_flag);
%          [r,t,h,outid,hboot,CI] = skipped_correlation(X,Y,fig_flag);
%
% INPUTS:  X and Y are 2 vectors or matrices, in the latter case,
%          correlations are computed column-wise 
%          fig_flag (1/0) indicates to plot the data or not
%
% OUTPUTS:
%          r is the pearson/spearman correlation 
%          t is the T value associated to the skipped correlation
%          h is the hypothesis of no association at alpha = 5% 
%          outid is the index of bivariate outliers
% 
%          optional:
%
%          hboot 1/0 declares the test significant based on CI (h depends on t)
%          CI is the robust confidence interval computed by bootstrapping the 
%          cleaned-up data set and taking the .95 centile values
%
% This code rely on the mid covariance determinant as implemented in LIBRA
% - Verboven, S., Hubert, M. (2005), LIBRA: a MATLAB Library for Robust Analysis,
% Chemometrics and Intelligent Laboratory Systems, 75, 127-136.
% - Rousseeuw, P.J. (1984), "Least Median of Squares Regression,"
% Journal of the American Statistical Association, Vol. 79, pp. 871-881.
% 
% The quantile of observations whose covariance is minimized is 
% floor((n+size(X,2)*2+1)/2)),
% i.e. ((number of observations + number of variables*2)+1) / 2, 
% thus for a correlation this is floor(n/2 + 5/2).
%git 
% See also MCDCOV, IDEALF.

% Cyril Pernet & Guillaume Rousselet, v1 - April 2012
% ---------------------------------------------------
%  Copyright (C) Corr_toolbox 2012

%% compute
gval = sqrt(chi2inv(0.975,2)); % in fact depends on size(X,2) but here always = 2
[n,p] = size(x);
for column = 1:p
    if p>1

        fprintf('skipped correlation: processing pair %g \n',column); 

    end
    
    % get the centre of the bivariate distributionsi 
    X = [x(:,column) y(:,column)];

    result=mcdcov(X,'cor',1,'plots',0,'h',floor((n+size(X,2)*2+1)/2));

    center = result.center;
    
     
    % orthogonal projection to the lines joining the center
    % followed by outlier detection using mad median rule
    
    vec=1:n;
    for i=1:n % for each row
        dis=NaN(n,1);
        B = (X(i,:)-center)';
        BB = B.^2;
        bot = sum(BB);
        if bot~=0
            for j=1:n
                A = (X(j,:)-center)';
                dis(j)= norm(A'*B/bot.*B); 
            end
            % MAD median rule
            %[outliers,value] = madmedianrule(dis,2);
            %record{i} = dis > (median(dis)+gval.*value);
            % IQR rule
            [ql,qu]=idealf(dis);
            record{i} = (dis > median(dis)+gval.*(qu-ql)) ; % + (dis < median(dis)-gval.*(qu-ql));
        end
    end
    
    try
        flag = nan(n,1);
        flag = sum(cell2mat(record),2); % if any point is flagged
        
    catch ME  % this can happen to have an empty cell so loop
        flag = nan(n,size(record,2));
        index = 1;
        for s=1:size(record,2)
          if ~isempty(record{s})
              flag(:,index) = record{s};
              index = index+1;
          end
        end
        flag(:,index:end) = [];
        flag = sum(flag,2);
    end
            
    if sum(flag)==0
        outid{column}=[];
    else
        flag=(flag>=1);
        outid{column}=vec(flag);
    end
    keep=vec(~flag);
    
    %% Pearson/Spearman correlation
    
    if  p == 1 % in the special case of a single test Pearson is valid too
        a{column} = x(keep);
        b{column} = y(keep);
        
        rp = sum(detrend(a{column},'constant').*detrend(b{column},'constant')) ./ ...
            (sum(detrend(a{column},'constant').^2).*sum(detrend(b{column},'constant').^2)).^(1/2);
        tp = rp*sqrt((n-2)/(1-rp.^2));
        r.Pearson = rp; t.Pearson = tp;
        
        coef.Pearson = pinv([a{column} ones(length(a{column}),1)]) * b{column};
        
        xrank = tiedrank(a{column},0); yrank = tiedrank(b{column},0);
        rs = sum(detrend(xrank,'constant').*detrend(yrank,'constant')) ./ ...
            (sum(detrend(xrank,'constant').^2).*sum(detrend(yrank,'constant').^2)).^(1/2);
        ts = rs*sqrt((n-2)/(1-rs.^2));
        r.Spearman = rs; t.Spearman = ts;
        
        % get regression lines for Spearman
        coef.Spearman = pinv([xrank ones(length(a{column}),1)])*yrank;
            
    
    else % multiple tests, only use Spearman to control type 1 error
        a{column} = x(keep,column); xrank = tiedrank(a{column},0); 
        b{column} = y(keep,column); yrank = tiedrank(b{column},0);
        r(column) = sum(detrend(xrank,'constant').*detrend(yrank,'constant')) ./ ...
            (sum(detrend(xrank,'constant').^2).*sum(detrend(yrank,'constant').^2)).^(1/2);
        t(column) = r(column)*sqrt((n-2)/(1-r(column).^2));
        coef.Spearman{column} =  pinv([xrank ones(length(a{column}),1)])*yrank; % get regression lines for spearman
    end
end

%% get h

% the default test of 0 correlation is for alpha = 5%

c = 6.947 / n + 2.3197; % valid for n between 10 and 200
if p == 1
    h.Pearson = abs(tp) >= c;
    h.Spearman = abs(ts) >= c;
else
    h= abs(t) >= c;
end

%% adjustement for multiple testing using the .95 quantile of Tmax
if p>1 && p<=10
    switch hypothesis
        
        case 1 % Hypothesis of 0 correlation between all pairs
            
            if p == 2;  q = 5.333*n^-1 + 2.374;    end
            if p == 3;  q = 8.8*n^-1 + 2.78;       end
            if p == 4;  q = 25.67*n^-1.2 + 3.03;   end
            if p == 5;  q = 32.83*n^-1.2 + 3.208;  end
            if p == 6;  q = 51.53*n^-1.3 + 3.372;  end
            if p == 7;  q = 75.02*n^-1.4 + 3.502;  end
            if p == 8;  q = 111.34*n^-1.5 + 3.722; end
            if p == 9;  q = 123.16*n^-1.5 + 3.825; end
            if p == 10; q = 126.72*n^-1.5 + 3.943; end
            
        case 2 % Hypothesis of 0 correlation between x1 and all y
            
            if p == 2;  q = 5.333*n^-1 + 2.374;   end
            if p == 3;  q = 8.811*n^-1 + 2.54;    end
            if p == 4;  q = 14.89*n^-1.2 + 2.666; end
            if p == 5;  q = 20.59*n^-1.2 + 2.920; end
            if p == 6;  q = 51.01*n^-1.5 + 2.999; end
            if p == 7;  q = 52.15*n^-1.5 + 3.097; end
            if p == 8;  q = 59.13*n^-1.5 + 3.258; end
            if p == 9;  q = 64.93*n^-1.5 + 3.286; end
            if p == 10; q = 58.5*n^-1.5 + 3.414;  end
    end
    
   h = abs(t) >= q; 
end

    
end