function [r,t,pval,CI] = Kendall(X,Y,varargin)

% Computes the Kendall correlation with its bootstrap CI.
%
% FORMAT:  [r,t,p] = Kendall(X,Y)
%          [r,t,p] = Kendall(X,Y,fig_flag,level)
%          [r,t,p,hboot,CI] = Kendall(X,Y,fig_flag,level)
%
% INPUTS:  X and Y are 2 vectors or matrices, in the latter case,
%          correlations are computed column-wise 
%          fig_flag indicates to plot (1 - default) the data or not (0)
%          level is the desired alpha level (5/100 is the default)
%
% OUTPUTS: r is the Spearman correlation
%          t is the associated t value
%          pval is the corresponding p value
%          hboot 1/0 declares the test significant based on CI
%          CI is the percentile bootstrap confidence interval
%
% If X and Y are matrices of size [n p], p correlations are computed
% and the CIs are adjusted at the alpha/p level (Bonferonni
% correction); hboot is based on these adjusted CIs but pval remains
% uncorrected.
%
% This function requires the tiedrank.m function from the matlab stat toolbox. 
%
% See also TIEDRANK.

% Cyril Pernet v1
% ---------------------------------
%  Copyright (C) Corr_toolbox 2012

%% input checks
if nargin < 2
    error('two inputs requested');
else
    [X,Y] = check_corr_XY(X,Y);
end

% defaults
nboot              = 1000;
alphav             = 5/100;
heteroscedasticity = 'unspecified'; % undocumented, on/off forces to return a CI (never a NaN)
if size(X,2) == 1 % and therefore size(Y,2) = 1 as well
    figflag = 'on';
else
    figflag = 'off';
end

% check options
for v=1:2:size(varargin,2)
    if strcmpi(varargin{v},'figure')
        figflag = varargin{v+1};
    elseif strcmpi(varargin{v}(1:6),'hetero') 
        heteroscedasticity = varargin{v+1};
        if ~strcmpi(varargin{v},'heteroscedasticity')
            fprintf('taking argument in ''%s'' as heteroscedasticity \n',varargin{v})
        end
    elseif strcmpi(varargin{v},'alpha')
        alphav = varargin{v+1};
        if ~isnumeric(alphav)
            alphav = str2double(alphav);
        end
    end
end

[n,p]  = size(X);
alphav = alphav./ p;


%% Compute 
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

% stats
if nargout > 1 || strcmpi(figflag ,'on')
    zr      = 0.5 * log((1+r)./(1-r));
    zalpha  = icdf('Normal',alphav/2,0,1);
    t       = r.*(sqrt(n-2)) ./ sqrt((1-r.^2));
    pval    = 2*tcdf(-abs(t),n-2);
    S       = 1/sqrt(n-3); % assumed standard error
    if p==1
        [t(2),pval(2),~,S_hc4]     = get_hc4stats(r,zr,zscore(X,0,1),zscore(Y,0,1));
    else
        [t(2,:),pval(2,:),~,S_hc4] = get_hc4stats(r,zr,zscore(X,0,1),zscore(Y,0,1));
    end
    
    for column = p:-1:1
        if all(pval(:,column)<alphav) || all(pval(:,column)>alphav) || strcmpi(heteroscedasticity,'off')
            CI(:,column) = tanh([zr(:,column)-abs(zalpha)*S ; zr(:,column)+abs(zalpha)*S]);
        else
            CI(:,column) = [NaN NaN]';
        end
    end
    
    %% undocumented - do a percentile t for the hc4 CI (kinda work for n>80)
    if strcmpi(heteroscedasticity,'on')
        
        for B=nboot:-1:1
            % make a resampling boot_table with enough unique pairs
            go = 0;
            while go == 0
                tmp = randi(n,n,1);
                if length(unique(tmp))>=6
                    boot_table(:,B) = tmp;
                    go = 1;
                end
            end
        end
        
        % resample and compute
        ibot = round(alphav*nboot/2)+1;
        itop = nboot-ibot+2;
        for column = p:-1:1
            zX           = X(:,p);
            zX           = zscore(tiedrank(zX(boot_table),1),0,1);
            zY           = Y(:,p);
            zY           = zscore(tiedrank(zY(boot_table),1),0,1);
            [~,~,B,S]    = get_hc4stats(r(column),zr(column),zX,zY); % all bootstraped betas and S
            v            = sort((B-r(column))./sqrt(S));
            CI(:,column) = [r(column)-v(itop)*sqrt(S_hc4(column)) r(column)-v(ibot)*sqrt(S_hc4(column))]';
        end
    end
end

%% figure
if strcmpi(figflag ,'on')
    if ~exist('boot_table','var')
        boot_table = randi(n,n,nboot);
    end
    for B=nboot:-1:1
        rX = tiedrank(X(boot_table(:,B),:),1);
        rY = tiedrank(Y(boot_table(:,B),:),1);
        t1 = Xadj(1); n1 = sum((t1*(t1-1))/2);
        t2 = Yadj(1); n2 = sum((t2*(t2-1))/2);
        K = 0;
        for k = 1:n-1
            K = K + sum(sign(rX(k)-rX(k+1:n)).*sign(rY(k)-rY(k+1:n)));
        end
        rb(B,:) = K / sqrt((n0-n1)*(n0-n2));
    end
end

if strcmpi(figflag ,'on')
    for f=1:length(r)
        figure('Name',sprintf('Kendall correlation X%g Y%g',f,f));
        set(gcf,'Color','w'); subplot(1,2,1);
        scatter(X(:,f),Y(:,f),100,'filled'); grid on; box on;
        xlabel('X','FontSize',12); ylabel('Y','FontSize',12);
        M = sprintf('r=%g \n %g%%CI [%.2f %.2f]',r(f),(1-alphav)*100,CI(1,f),CI(2,f));
        title(M,'FontSize',14); h=lsline; set(h,'Color','r','LineWidth',4);
        
        subplot(1,2,2); k = round(1 + log2(nboot));
        MV = histogram(rb(:,f),k); MV = max(MV.Values); grid on;
        title(sprintf('Bootstrapped correlations \n median=%g',median(rb)),'FontSize',14);
        xlabel('boot correlations','FontSize',12);ylabel('frequency','FontSize',12)
        hold on; plot(repmat(CI(1,f),MV,1),1:MV,'r','LineWidth',4);
        axis tight; colormap([.4 .4 1]); box on;
        plot(median(rb),MV/2,'ko','LIneWidth',3)
        
        if all(~isnan(CI(1:2,f))) % plot CI
            plot(repmat(CI(2,f),MV,1),1:MV,'r','LineWidth',4);
            subplot(1,2,1); hold on
            betas = pinv([X(:,f) ones(n,1)])*Y(:,f);
            transform = betas(1)/r(f);
            y1 = refline(CI(1,f)*transform,betas(2)); set(y1,'Color','r');
            y2 = refline(CI(2,f)*transform,betas(2)); set(y2,'Color','r');
            y1 = get(y1); y2 = get(y2);
            xpoints=[y1.XData(1):y1.XData(2),y2.XData(2):-1:y2.XData(1)];
            step1 = y1.YData(2)-y1.YData(1); step1 = step1 / (y1.XData(2)-y1.XData(1));
            step2 = y2.YData(2)-y2.YData(1); step2 = step2 / (y2.XData(2)-y2.XData(1));
            filled=[y1.YData(1):step1:y1.YData(2),y2.YData(2):-step2:y2.YData(1)];
            hold on; fillhandle=fill(xpoints,filled,[1 0 0]);
            set(fillhandle,'EdgeColor',[1 0 0],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
        end
    end
end

