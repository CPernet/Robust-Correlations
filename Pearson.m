function [r,t,pval,CI,hboot] = Pearson(X,Y,varargin)

% compute the Pearson correlation along with the bootstrap CI
%
% FORMAT:  [r,t,p]          = Pearson(X,Y,options)
%          [r,t,p]          = Pearson(X,Y,'figure','on','alpha',0.05)
%          [r,t,p,CI,hboot] = Pearson(X,Y,'figure','on','alpha',0.05)
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
%
%          optional: CI is the percentile bootstrap confidence interval
%                    hboot 1/0 declares the test significant based on CI 
%
% If X and Y are matrices of size [n p], p correlations are computed
% consequently, the CI are adjusted at a level alpha/p (Bonferonni
% correction) and hboot is based on these adjusted CI (pval remains
% uncorrected)
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
alphav = 5/100;
if size(X,2) == 1 % and therefore size(Y,2) = 1 as well
    figflag = 'on';
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
    end
end
[n,p] = size(X);

%% Basic Pearson

r    = sum(detrend(X,'constant').*detrend(Y,'constant')) ./ ...
       (sum(detrend(X,'constant').^2).*sum(detrend(Y,'constant').^2)).^(1/2);
t    = r.*sqrt((n-2)./(1-r.^2));
pval = 2*tcdf(-abs(t),n-2);

%% compute the CI
if nargout > 3
    % adjust boot parameters
    if p == 1
        nboot = 599;
        % adjust percentiles following Wilcox
        if n < 40
            low = 7 ; high = 593;
        elseif n >= 40 && n < 80
            low = 8 ; high = 592;
        elseif n >= 80 && n < 180
            low = 11 ; high = 588;
        elseif n >= 180 && n < 250
            low = 14 ; high = 585;
        elseif n >= 250
            low = 15 ; high = 584;
        end
        
    else
        nboot  = 1000;
        alphav = alphav / p;
        low = round((alphav*nboot)/2);
        if low == 0
            error('adjusted CI cannot be computed, too many tests for the number of observations')
        else
            high = nboot - low;
        end
    end
    
    % compute hboot and CI
    table = randi(n,n,nboot);
    for B=nboot:-1:1
        rb(B,:) = sum(detrend(X(table(:,B),:),'constant').*detrend(Y(table(:,B),:),'constant')) ./ ...
            (sum(detrend(X(table(:,B),:),'constant').^2).*sum(detrend(Y(table(:,B),:),'constant').^2)).^(1/2);
        if strcmpi(figflag ,'on')
            for c=size(X,2):-1:1
                b              = pinv([X(table(:,B),c) ones(n,1)])*Y(table(:,B),c);
                slope(B,c)     = b(1);
                intercept(B,c) = b(2,:);
            end
        end
    end
    
    [rb,index] = sort(rb,1); % keep same index as r for slopes
    intercept  = intercept(index);
    slope      = slope(index);
    
    % CI and h
    adj_nboot = nboot - sum(isnan(rb));
    adj_low   = round((alphav*adj_nboot)/2);
    adj_high  = adj_nboot - adj_low;
    
    for c=size(X,2):-1:1
        CI(:,c)              = [rb(adj_low(c),c) ; rb(adj_high(c),c)];
        hboot(c)             = (rb(adj_low(c),c) > 0) + (rb(adj_high(c),c) < 0);
        if strcmpi(figflag ,'on')
            CIslope(:,c)     = [slope(adj_low(c),c) ; slope(adj_high(c),c)];
            CIintercept(:,c) = [intercept(adj_low(c),c) ; intercept(adj_high(c),c)];
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
            M = sprintf('Pearson corr r=%g \n %g%%CI [%g %g]',r(f),(1-alphav)*100,CI(1,f),CI(2,f));
        else
            M = sprintf('Pearson corr r=%g \n p=%g',r(f),pval(f));
        end
        
        scatter(X(:,f),Y(:,f),100,'filled'); grid on
        xlabel('X','FontSize',14); ylabel('Y','FontSize',14);
        title(M,'FontSize',16);
        h=lsline; set(h,'Color','r','LineWidth',4);
        box on;set(gca,'Fontsize',14)
        
        if nargout>3 % if bootstrap done plot CI
            y1 = refline(CIslope(1),CIintercept(1)); set(y1,'Color','r');
            y2 = refline(CIslope(2),CIintercept(2)); set(y2,'Color','r');
            y1 = get(y1); y2 = get(y2);
            xpoints=[[y1.XData(1):y1.XData(2)],[y2.XData(2):-1:y2.XData(1)]];
            step1 = y1.YData(2)-y1.YData(1); step1 = step1 / (y1.XData(2)-y1.XData(1));
            step2 = y2.YData(2)-y2.YData(1); step2 = step2 / (y2.XData(2)-y2.XData(1));
            filled=[[y1.YData(1):step1:y1.YData(2)],[y2.YData(2):-step2:y2.YData(1)]];
            hold on; fillhandle=fill(xpoints,filled,[1 0 0]);
            set(fillhandle,'EdgeColor',[1 0 0],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
            
            subplot(1,2,2); k = round(1 + log2(length(rb))); hist(rb,k); grid on;
            title({'Bootstrapped correlations';['h=',num2str(hboot)]},'FontSize',16); hold on
            xlabel('boot correlations','FontSize',14);ylabel('frequency','FontSize',14)
            plot(repmat(CI(1),max(hist(rb,k)),1),[1:max(hist(rb,k))],'r','LineWidth',4);
            plot(repmat(CI(2),max(hist(rb,k)),1),[1:max(hist(rb,k))],'r','LineWidth',4);
            axis tight; colormap([.4 .4 1])
            box on;set(gca,'Fontsize',14)
        end
    end
end






