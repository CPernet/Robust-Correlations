function [r,t,pval,CI,hboot] = Pearson(X,Y,varargin)

% Computes the Pearson correlation along with the bias corrected and accelerated
% percentile bootstrap Confidence Interval. Note there are no good method for CI but
% the BCa gives a relatively good probability coverage under the null.
%
% FORMAT:  [r,t,p]    = Pearson(X,Y,options)
%          [r,t,p,CI] = Pearson(X,Y,'figure','on','alpha',0.05)
%
% INPUTS:  X and Y are 2 vectors or matrices (correlations are computed column-wise)
%                  if X and Y are a vector and a matrix, the vector is
%                  replicated to match the matrix size
%          options are 'figure', if X and Y are vectors, this is 'on' by
%                                default, if X and Y are matrices this 'off' by default
%                      'alpha', the alpha level to use for the confidence
%                               interval (default is 5% Bonferroni adjusted)
%
% OUTPUTS: r is the Pearson correlation
%          t is the associated t value
%          pval is the corresponding p value
%          CI is the percentile bootstrap confidence interval
%
% If X and Y are matrices of size [n p], p correlations are computed
% consequently, the CI are adjusted at a level alpha/p (Bonferonni
% correction)
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
method = 'bias corrected accelerated';
if size(X,2) == 1 % and therefore size(Y,2) = 1 as well
    figflag = 'on';
else
    figflag = 'off';
end

% check options
for v=1:2:size(varargin,2)
    if strcmpi(varargin{v},'figure')
        figflag = varargin{v+1};
    elseif strcmpi(varargin{v},'alpha')
        alphav = varargin{v+1};
        if ~isnumeric(alphav)
            alphav = str2double(alphav);
        end
    elseif strcmpi(varargin{v},'method') % undocumented since only BCa ~works
        method = varargin{v+1};
    end
end
[n,p]  = size(X);

%% Basic Pearson

r    = sum(detrend(X,'constant').*detrend(Y,'constant')) ./ ...
    (sum(detrend(X,'constant').^2).*sum(detrend(Y,'constant').^2)).^(1/2);
t    = r.*sqrt((n-2)./(1-r.^2));
pval = 2*tcdf(-abs(t),n-2);

%% compute the CI
if nargout > 3
    
    % get bootstrapped r
    NB = round(1.2*nboot); % always do a few more to avoid NaNs
    table = randi(n,n,NB);
    for B=NB:-1:1
        rb(B,:) = sum(detrend(X(table(:,B),:),'constant').*detrend(Y(table(:,B),:),'constant')) ./ ...
            (sum(detrend(X(table(:,B),:),'constant').^2).*sum(detrend(Y(table(:,B),:),'constant').^2)).^(1/2);
        if strcmpi(figflag ,'on')
            for c=size(X,2):-1:1
                b              = pinv([X(table(:,B),c) ones(n,1)])*Y(table(:,B),c);
                slope(B,c)     = b(1);
            end
        end
    end
    
    rb(isnan(rb)) = [];                % remove possible NaN
    keep          = randi(NB,nboot,1); % choose nboot out of NB
    rb            = rb(keep);          % retain nboot sample as intended
    [rb,index]    = sort(rb,1);        % keep same index as r for slopes
    if strcmpi(figflag ,'on')
        slope = slope(keep);
        slope = slope(index,:);
    end
    
    % get the CI
    nanval = sum(isnan(rb)); % should be 0, but JIC
    nboot  = nboot-nanval;   % adjust valid number of bootstraps
    alphav = alphav / p;     % Bonferroni corrected
    
    if strcmpi(method,'percentile')
        low    = round((alphav*nboot)/2);
        if low == 0
            error('adjusted CI cannot be computed, too many tests for the number of observations')
        else
            high = nboot - low;
        end
    else 
        % bias value zo (include half of the bootstrap values that are
        % tied with the original sample value into z0)
        z0     = icdf('Normal',mean(rb<r)+mean(rb==r,1)/2,0,1);
        zalpha = icdf('Normal',alphav/2,0,1);
       
        if strcmpi(method,'bias corrected')
            low  = floor(nboot*cdf('Normal',(2*z0+(zalpha)),0,1));
            high = ceil(nboot*cdf('Normal',(2*z0-(zalpha)),0,1));
        elseif strcmpi(method,'bias corrected accelerated')
            % compute acceleration factor: DiCiccio and Efron (1996)
            for ss=n:-1:1
                index = 1:n; index(ss) = []; % jackknife
                jr(ss,:) = Pearson(X(index,:),Y(index,:),'figure','off');
            end
            skew = skewness(mean(jr)-jr) /sqrt(n); % adjusted skewness
            acc  = skew/6;                         % acceleration
            low  = floor(nboot*cdf('Normal',(z0 +(z0+zalpha)./(1-acc.*(z0+zalpha))),0,1));
            high = ceil(nboot*cdf('Normal',(z0 +(z0-zalpha)./(1-acc.*(z0-zalpha))),0,1));
        end
    end
     
    % CI and h
    for c=p:-1:1
        CI(:,c)  = [rb(low(c),c) ; rb(high(c),c)];
        hboot(c) = (CI(1,c) > 0) + (CI(2,c)< 0);
        if strcmpi(figflag ,'on')
            CIslope(1,c)   = slope(low(c),c);
            CIslope(2,c)   = slope(high(c),c);
            b              = pinv([X(:,c) ones(n,1)])*Y(:,c);
            CIintercept(c) = b(2);
        end
    end
end

%% plots
if strcmpi(figflag ,'on')
    for f=1:length(r)
        figure('Name',sprintf('Pearson correlation X%g Y%g',f,f));
        set(gcf,'Color','w');
        
        if nargout>3
            subplot(1,2,1);
            M = sprintf('r=%g \n %g%%CI [%.2f %.2f]',r(f),(1-alphav)*100,CI(1,f),CI(2,f));
        else
            M = sprintf('r=%g \n p=%g',r(f),pval(f));
        end
        
        scatter(X(:,f),Y(:,f),100,'filled'); grid on; box on;
        xlabel('X','FontSize',12); ylabel('Y','FontSize',12);
        title(M,'FontSize',14); h=lsline; set(h,'Color','r','LineWidth',4);
        
        if nargout>3 % if bootstrap done plot CI
            y1 = refline(CIslope(1,f),CIintercept(c)); set(y1,'Color','r');
            y2 = refline(CIslope(2,f),CIintercept(c)); set(y2,'Color','r');
            y1 = get(y1); y2 = get(y2);
            xpoints=[y1.XData(1):y1.XData(2),y2.XData(2):-1:y2.XData(1)];
            step1 = y1.YData(2)-y1.YData(1); step1 = step1 / (y1.XData(2)-y1.XData(1));
            step2 = y2.YData(2)-y2.YData(1); step2 = step2 / (y2.XData(2)-y2.XData(1));
            filled=[y1.YData(1):step1:y1.YData(2),y2.YData(2):-step2:y2.YData(1)];
            hold on; fillhandle=fill(xpoints,filled,[1 0 0]);
            set(fillhandle,'EdgeColor',[1 0 0],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
            
            subplot(1,2,2); k = round(1 + log2(nboot(c)));
            MV = histogram(rb(:,f),k); MV = max(MV.Values); grid on;
            title('Bootstrapped correlations','FontSize',14); hold on
            xlabel('boot correlations','FontSize',12);ylabel('frequency','FontSize',12)
            plot(repmat(CI(1,f),MV,1),1:MV,'r','LineWidth',4);
            plot(repmat(CI(2,f),MV,1),1:MV,'r','LineWidth',4);
            axis tight; colormap([.4 .4 1]); box on;
        end
    end
end





