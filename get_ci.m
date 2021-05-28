function [CI,z0] = get_ci(X,Y,varargin)

% routine to return bootstrapped CI, typically called by a correlation function
%
% FORMAT CI = get_ci(X,Y,'correlation','name','method','name','alpha',value)
%        for instance CI = get_ci(X,Y,'correlation', Pearson, 'method', 'bca')
%
% INPUTS X and Y are the data used for a given correlation
%        'correlation' can be 'Pearson','Spearman','Bend','Skipped Pearson','Skipped Spearman'
%        'method','Percentile','BC','BCa', 'all'
%        'alpha', is the alpha value for a 1-alpha CI
%
% OUTPUT CI is a [2 p] or [6 p] matrix with the lower/upper bounds in rows 
%                replicated for the Percentile bootstrap,BC and BCa if 'all'
%                is selected as a method. The p columns are for the p [XY] 
%                pairs used as input, and corrected for multiple comparisons.
%       z0 is the bootstrap bias
%
% Dr Cyril Pernet - University of Edinburgh
% -----------------------------------------
% Copyright (C) Corr_toolbox 2020

%% input checks
if nargin < 2
    error('two inputs requested');
else
    [X,Y] = check_corr_XY(X,Y);
end

% defaults
nboot  = 1000;
alphav = 5/100;

% check options
for v=1:2:size(varargin,2)
    if strcmpi(varargin{v},'alpha')
        alphav = varargin{v+1};
        if ~isnumeric(alphav)
            alphav = str2double(alphav);
        end
    elseif strcmpi(varargin{v},'correlation')
        correlation = varargin{v+1};
    elseif strcmpi(varargin{v},'method')
        method = varargin{v+1};
    end
end
[n,p]  = size(X);
alphav = alphav / p; % Bonferroni corrected
zalpha = icdf('Normal',alphav/2,0,1);

%% compute the correlation
if strcmpi(correlation,'Pearson')
    r = sum(detrend(X,'constant').*detrend(Y,'constant')) ./ ...
        (sum(detrend(X,'constant').^2).*sum(detrend(Y,'constant').^2)).^(1/2);
elseif strcmpi(correlation,'Spearman')
    Xrank = tiedrank(X,0);
    Yrank = tiedrank(Y,0);
    r = sum(detrend(Xrank,'constant').*detrend(Yrank,'constant')) ./ ...
        (sum(detrend(Xrank,'constant').^2).*sum(detrend(Yrank,'constant').^2)).^(1/2);
elseif strcmpi(correlation,'Kendall')
    [Xrank,Xadj] = tiedrank(X,1);
    [Yrank,Yadj] = tiedrank(Y,1);
    n0 = (n*(n-1))/2;
    t1 = Xadj(1); n1 = sum((t1*(t1-1))/2);
    t2 = Yadj(1); n2 = sum((t2*(t2-1))/2);
    K = 0;
    for k = 1:n-1
        K = K + sum(sign(Xrank(k)-Xrank(k+1:n)).*sign(Yrank(k)-Yrank(k+1:n)));
    end
    r = K / sqrt((n0-n1)*(n0-n2));
end

%% get the bootstrapped statistics
% make sure the bias is not infinit! 
run_boot = 1;

while run_boot == 1
    
    % make a resampling boot_table with at least 3 unique pairs
    for B=nboot:-1:1
        go = 0;
        while go == 0
            tmp = randi(n,n,1);
            if length(unique(tmp))>=3
                boot_table(:,B) = tmp;
                go = 1;
            end
        end
    end
    
    % compute nboot values
    for B=nboot:-1:1
        if strcmpi(correlation,'Pearson')
            rb(B,:) = sum(detrend(X(boot_table(:,B),:),'constant').*detrend(Y(boot_table(:,B),:),'constant')) ./ ...
                (sum(detrend(X(boot_table(:,B),:),'constant').^2).*sum(detrend(Y(boot_table(:,B),:),'constant').^2)).^(1/2);
        elseif strcmpi(correlation,'Spearman')
            xrank = tiedrank(X(boot_table(:,B),:),0); 
            yrank = tiedrank(Y(boot_table(:,B),:),0);
            rb(B,:) = sum(detrend(xrank,'constant').*detrend(yrank,'constant')) ./ ...
                (sum(detrend(xrank,'constant').^2).*sum(detrend(yrank,'constant').^2)).^(1/2);
        elseif strcmpi(correlation,'Kendall')
            xrank = tiedrank(X(boot_table(:,B),:),1);
            yrank = tiedrank(Y(boot_table(:,B),:),1);
            K = 0;
            for k = 1:n-1
                K = K + sum(sign(xrank(k)-xrank(k+1:n)).*sign(yrank(k)-yrank(k+1:n)));
            end
            rb(B,:) = K / sqrt((n0-n1)*(n0-n2)); % reuse the same n0,n1,n2 as with X,Y
        end
    end
    rb = sort(rb);
    
    % bias value zo (include half of the bootstrap values that are
    % tied with the original sample value into z0)
    z0 = icdf('Normal',mean(rb<r)+mean(rb==r,1)/2,0,1);
    if ~isinf(z0)
        run_boot = 0; 
    else
        disp('oh infnite bootrapped bias in get_ci, recomputing ... ')
    end
end


%% Compute the CI

if strcmpi(method,'percentile') || strcmpi(method,'all')
    low  = round((alphav*nboot)/2);
    high = nboot - low;
    CI = cibounds(rb,low,high);
end

if strcmpi(method,'bias corrected') || strcmpi(method,'all')
    low  = floor(nboot*cdf('Normal',(2*z0+(zalpha)),0,1)); 
    high = ceil(nboot*cdf('Normal',(2*z0-(zalpha)),0,1));
    if strcmpi(method,'all')
        CI = [CI;NaN(2,p)];
        CI([3 4],:) = cibounds(rb,low,high);
    else
        CI = cibounds(rb,low,high);
    end
end

if strcmpi(method,'bias corrected accelerated') || strcmpi(method,'all')
    % compute acceleration factor: DiCiccio and Efron (1996)
    for ss=n:-1:1
        index = 1:n; index(ss) = []; % jackknife
        if strcmpi(correlation,'Pearson')
            jr(ss,:) = Pearson(X(index,:),Y(index,:),'figure','off');
        elseif  strcmpi(correlation,'Spearman')
            jr(ss,:) = Spearman(X(index,:),Y(index,:),'figure','off');
        elseif  strcmpi(correlation,'Kendall')
            jr(ss,:) = Kendall(X(index,:),Y(index,:),'figure','off');
        end
    end
    skew = skewness(mean(jr)-jr) /sqrt(n); % adjusted skewness
    acc  = skew/6;                         % acceleration
    low  = floor(nboot*cdf('Normal',(z0 +(z0+zalpha)./(1-acc.*(z0+zalpha))),0,1));
    high = ceil(nboot*cdf('Normal',(z0 +(z0-zalpha)./(1-acc.*(z0-zalpha))),0,1));
    if strcmpi(method,'all')
        CI = [CI;NaN(2,p)];
        CI([5 6],:) = cibounds(rb,low,high);
    else
        CI = cibounds(rb,low,high);
    end
end

% cleanup
CI(CI<-1) = -1;
CI(CI>1) = 1;

function CI = cibounds(rb,low,high)

% check inputs
[n,p] = size(rb);

if isnan(low) || low == 0
    low = 1; high = p;
end

if size(low,2) == 1
    low  = repmat(low,1,p);
    high = repmat(high,1,p);
end

% make CI 
CI = repmat([NaN NaN]',1,p);
for c=p:-1:1
    CI(:,c)  = [rb(low(c),c) ; rb(high(c),c)]';
end
