function [r,t,p,CI,H,pH] = bendcorr(X,Y,varargin)

% Computes the percentage bend correlation along with the bootstrap CI.
% When H0 is rejected (say pval<0.05), it is reasonable to conclude that X
% and Y are dependent, but the nature of that dependence is unclear. If fact,
% in case of heteroscedasticity, even for rho = 0, H0 can be rejected because
% the t value relies on a wrong estimate of the variance, i.e. the test 
% statisitics t, tests if X and Y are independent rather than testing if 
% rho =0. Therefore, if pval is below the alpha level (H0 rejected) it should
% be compared to a robust t-test version that correct the covariance. Here, 
% we return in t(2) and pval(2) such values using the HC4 estimator of the 
% standard error (Godfrey, 2006). If that second p-value is also below alpha,
% CI are retuned using the standard z-transform. If pval is above alpha, but
% the robust estimate is below, the discrepency is likely due to heteroscedasticity
% and no good CI can be derived (returns NaN).
%
% FORMAT:  [r,t,p,CI]      = bendcorr(X,Y,'beta',0.2)
%          [r,t,p,CI]      = bendcorr(X,Y,'beta',0.2,'figure','on','alpha',0.05)
%          [r,t,p,CI,H,pH] = bendcorr(X,Y,'beta',0.2,'figure','on','alpha',0.05)
%
% INPUTS:  X and Y are 2 vectors or matrices (correlations are computed column-wise)
%                  if X and Y are a vector and a matrix, the vector is
%                  replicated to match the matrix size
%          options are 'beta', the amount of trimming: 0 <= beta <= 0.5
%                      (beta a.k.a. the bending constant for omega - default = 0.2)
%                       'figure', if X and Y are vectors, this is 'on' by
%                                default, if X and Y are matrices this 'off' by default
%                      'alpha', the alpha level to use for the confidence interval
%          If X and Y are matrices of size [n p], p correlations are computed and
%          the CI are adjusted at a level alpha/p (Bonferonni correction)
%
% OUTPUTS: r is the percentage bend correlation
%          t are the t-values using the standard equation and using the hc4 variance estimate
%          pval are the p values corresponding to t-tests
%          CI is 1-alpha confidence interval (only if all(pvals<alpha) or all(pvals>alpha))
%     for multiple comparisons
%          H is the measure of association between all pairs
%          pH is the p value for an omnibus test of independence between all pairs 
%
%
% If 'figure' is on, the data scatter plot is shown with the least
% square fit of the data and the confidence intervals from the r value. As a
% sanity check, the boostrapped correlations are also shown reporting the
% median of bootstrapped value.
%
% References:
% Godfrey (2006). Tests for regression models with heteroscedasticity of
%                 unknown form. Comp Stat & Data Analysis, 50, 2715-2733 
% Wilcox (2017). Introduction to Robust Estimation and Hypothesis Testing.
%                4th Ed. Acedemic Press
%
% This function requires the tiedrank.m function from the matlab stat toolbox. 
% See also TIEDRANK.
%
% Dr Cyril Pernet - University of Edinburgh
% Dr Guillaume Rousselet - University of Glasgow
% -----------------------------------------
% Copyright (C) Corr_toolbox 2020

%% input checks
if nargin < 2
    error('two inputs requested');
else
    [X,Y] = check_corr_XY(X,Y);
end

% defaults
nboot              = 1000;
beta               = 0.2;
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
    elseif strcmpi(varargin{v},'beta')
        beta = varargin{v+1};
        if beta>1
            beta = beta/100;
        end
    end
end

[n,p]  = size(X);
alphav = alphav./ p;
if beta<0 || beta>.5
    error('beta must be between 0 and 50%')
end

%% compute
% --------
if nargout > 4
    [r,t,p,XX,YY,H,pH] = bend_compute(X,beta,heteroscedasticity);
else
    [r,t,p,XX,YY] = bend_compute(X,beta,heteroscedasticity);
end

if nargout > 3 
    % bootstrap
    % -----------
    nboot = 1000;
    low = round((alphav*nboot)/2);
    if low == 0
        error('adjusted CI cannot be computed, too many tests for the number of observations')
    else
        high = nboot - low;
    end
     
    % make a resampling boot_table with enough unique pairs
    for B=nboot:-1:1
        go = 0;
        while go == 0
            tmp = randi(n,n,1);
            if length(unique(tmp))>=6
                boot_table(:,B) = tmp;
                go = 1;
            end
        end
    end
    
    % get bootrapped bend correlation
    for B=nboot:-1:1
        tmp = X(boot_table(:,B),:);
        rb(B,:) = bend_compute(tmp,beta);
    end
    rb = sort(rb);

    % CI and h
    for c=size(X,2)/2:-1:1
        CI(:,c) = [rb(low(c),c) ; rb(high(c),c)];
        hboot(c) = (rb(low(c),c) > 0) + (rb(high(c),c) < 0);
    end
end
 
%% plot
% -----
if strcmpi(figflag ,'on')
    for f=1:length(r)
        figure('Name',sprintf('Percentage bend correlation X%g Y%g',f,f));
        set(gcf,'Color','w'); subplot(1,2,1);
        scatter(X(XX{1},1),X(XX{1},2),110,'r','LineWidth',3);
        scatter(X(YY{1},1),X(YY{1},2),110,'g','LineWidth',3);
        scatter(X(intersect(XX{1},YY{1}),1),X(intersect(XX{1},YY{1}),2),110,'k','LineWidth',3);
        xlabel('X','FontSize',12); ylabel('Y','FontSize',12);grid on; box on; 
        M = sprintf('r=%g \n %g%%CI [%.2f %.2f]',r(f),(1-alphav)*100,CI(1,f),CI(2,f));
        title(M,'FontSize',14); h=lsline; set(h,'Color','r','LineWidth',4);
        
        subplot(1,2,2); k = round(1 + log2(nboot));
        MV = histogram(rb(:,f),k); MV = max(MV.Values); grid on;
        title(sprintf('Bootstrapped correlations \n median=%g h=%g',median(rb),hboot(f)),'FontSize',14);
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
end

function [r,t,p,XX,YY,H,pH] = bend_compute(X,beta,heteroscedasticity)

H= []; pH = [];

%% Medians and absolute deviation from the medians
% ---------------------------------------------
M = repmat(median(X),size(X,1),1);
W = sort(abs(X-M),1); 
 
% limits
% -------
m     = floor((1-beta)*size(X,1));
omega = W(m,:); 

%% Compute the correlation
% ------------------------
P           = (X-M)./ repmat(omega,size(X,1),1); 
P(isnan(P)) = 0; 
P(isinf(P)) = 0; % correct if omega = 0
comb        = [(1:size(X,2)/2)',((1:size(X,2)/2)+size(X,2)/2)']; % all pairs of columns

for j = size(comb,1):-1:1
    
    % column 1
    psi          = P(:,comb(j,1)); 
    i1           = length(psi(psi<-1)); 
    i2           = length(psi(psi>1)); 
    sx           = X(:,comb(j,1)); 
    sx(psi<(-1)) = 0; 
    sx(psi>1)    = 0; 
    pbos         = (sum(sx)+ omega(comb(j,1))*(i2-i1)) / (size(X,1)-i1-i2); 
    a            = (X(:,comb(j,1))-pbos)./repmat(omega(comb(j,1)),size(X,1),1); 
        
    % column 2
    psi          = P(:,comb(j,2));
    i1           = length(psi(psi<-1));
    i2           = length(psi(psi>1));
    sx           = X(:,comb(j,2));
    sx(psi<(-1)) = 0;
    sx(psi>1)    = 0;
    pbos         = (sum(sx)+ omega(comb(j,2))*(i2-i1)) / (size(X,1)-i1-i2);
    b            = (X(:,comb(j,2))-pbos)./repmat(omega(comb(j,2)),size(X,1),1);
    
     % return values of a,b to plot 
     XX{j} = union(find(a <= -1),find(a >= 1)); 
     YY{j} = union(find(b <= -1),find(b >= 1));
     
     % bend
     a(a<=-1) = -1; a(a>=1) = 1; 
     b(b<=-1) = -1; b(b>=1) = 1; 
     
     % get r, t and p
     r(j) = sum(a.*b)/sqrt(sum(a.^2)*sum(b.^2));
     if strcmpi(heteroscedasticity,'on')
         [~,~,~,S] = get_hc4stats(r,zr,zscore(X,0,1),zscore(Y,0,1));
     else
        S = sqrt((size(X,1) - 2)/(1 - r(j).^2)); 
     end
     t(j) = r(j)*S;
     p(j) = 2*(1 - tcdf(abs(t(j)),size(X,1)-2));
end

if size(X,2) > 2 && nargout > 5
    bv = 48*(size(X,1)-2.5).^2;
    for j=length(comb):-1:1
        c(j) = sqrt((size(X,1)-2.5)*log(1+t(j)^2/(size(X,1)-2))); 
        z(j) = c(j) + (c(j)^3+3*c(j))/bv - ( (4*c(j).^7+33*c(j).^5+240*c(j)^3+855*c(j)) / (10*bv.^2+8*bv*c(j).^4+1000*bv) ); 
    end
    H =  sum(z.^2);
    pH= 1- cdf('chi2',H,(size(X,2)*(size(X,2)-1))/2);
end

end