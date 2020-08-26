function [r,z,pval,CI] = Kendall(X,Y,varargin)

% Computes the Kendall correlation along with the alpha percent CI using
% the Z-transform but a critical value based on the Student t-distribution
% and modified estimate of the standard error for a corr<0.8, which gives
% a relativey good probability coverage under normality assumptions.
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
% FORMAT:  [r,t,p,CI] = Kendall(X,Y)
%          [r,t,p,CI] = Kendall(X,Y,'figure','on','alpha',0.05)
%
% INPUTS:  X and Y are 2 vectors or matrices (correlations are computed column-wise)
%                  if X and Y are a vector and a matrix, the vector is
%                  replicated to match the matrix size
%          options are 'figure', if X and Y are vectors, this is 'on' by
%                                default, if X and Y are matrices this 'off' by default
%                      'alpha', the alpha level to use for the confidence interval
%          If X and Y are matrices of size [n p], p correlations are computed and
%          the CI are adjusted at a level alpha/p (Bonferonni correction)
%
% OUTPUTS: r is the Kendall correlation
%          t are the t-values using the standard equation and using the hc4 variance estimate
%          pval are the p values corresponding to t-tests
%          CI is 1-alpha confidence interval (only if all(pvals<alpha) or all(pvals>alpha))
%
% If 'figure' is on, the data scatter plot is shown with the least
% square fit of the data and the confidence intervals from the r value. As a
% sanity check, the boostrapped correlations are also shown reporting the
% median of bootstrapped value.
%
% References:
% Kendall, M.G. (1955) Rank Correlation Methods. Griffin, London.
% Fieller, Hartley & Pearson, (1957) Tests for rank correlation
%                                    coefficients: Biometrika, 44, 470-481
% Godfrey (2006). Tests for regression models with heteroscedasticity of
%                 unknown form. Comp Stat & Data Analysis, 50, 2715-2733
% Wilcox (2017). Introduction to Robust Estimation and Hypothesis Testing.
%                4th Ed. Acedemic Press
%
% This function requires the tiedrank.m function from the matlab stat toolbox.
% See also TIEDRANK.
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


%% Compute in matrix form
T        = [X Y];
[i1, i2] = find(tril(ones(n, 'uint8'), -1)); % all indices lower n*n matrix
tau      = sign(T(i2, :) - T(i1, :));        % all pair diffences for all columns
tau      = tau' * tau;                       % between and within sums
temp     = diag(tau);                        % sum of within
Adj      = sqrt(temp * temp');
tau      = tau ./ Adj;                       % squared form
r        = diag(tau(p+1:end,1:p))';          % keep tau between X1..p and Y1..p
Adj      = diag(Adj(p+1:end,1:p))';          % reused in the bootstrap

% stats
if nargout > 1 || strcmpi(figflag ,'on')
    K       = r.*Adj;
    v0      = n*(n-1)*(2*n+5);
    zalpha  = icdf('Normal',alphav/2,0,1);

    for column = p:-1:1
        U            = unique(X(:,column));
        out          = arrayfun(@(x) find(X(:,column)==x),U,'UniformOutput',false);
        Len          = cellfun(@(c) length(c), out, 'UniformOutput', true);
        t1           = Len(Len ~= 1); % number of tied values within each X repeats
        vt1          = sum(t1.*(t1-1).*(2.*t1+5));
        U            = unique(Y(:,column));
        out          = arrayfun(@(x) find(Y(:,column)==x),U,'UniformOutput',false);
        Len          = cellfun(@(c) length(c), out, 'UniformOutput', true);
        t2           = Len(Len ~= 1); % number of tied values within each Y repeats
        vt2          = sum(t2.*(t2-1).*(2.*t2+5));
        v1           = sum(t1.*(t1-1))*sum(t2.*(t2-1)./(2*n*(n-1)));
        v2           = sum(t1.*(t1-1).*(t1-2))*sum(t2.*(t2-1).*(t2-2)./(9*n*(n-1)*(n-2)));
        v            = (v0-vt1-vt2)/18+v1+v2;
        z(column)    = K(column) / sqrt(v);
        tmp          = cdf('Normal',z(column),0,1);
        pval(column) = min(tmp,1-tmp);

        if nargout > 3
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
                ibot         = round(alphav*nboot/2)+1;
                itop         = nboot-ibot+2;
                zX           = X(:,column);
                zX           = zscore(tiedrank(zX(boot_table),0),0,1);
                zY           = Y(:,column);
                zY           = zscore(tiedrank(zY(boot_table),0),0,1);
                [~,~,B,S]    = get_hc4stats(r(column),zX,zY); % all bootstraped betas and S
                v            = sort((B-r(column))./sqrt(S));
                CI(:,column) = [r(column)-v(itop)*sqrt(S_hc4(column)) r(column)-v(ibot)*sqrt(S_hc4(column))]';
            else
                if abs(r(column)) < 0.8
                    S   = 0.437/sqrt(n-4); % Fieller, Hartley & Pearson, 1957
                else
                    S   = 1/sqrt(n-3); % Fisher 1925
                end
                CI(:,column) = tanh([z(column)-abs(zalpha)*S ; z(column)+abs(zalpha)*S]);
            end
        end
    end
end

%% figure
if strcmpi(figflag ,'on')
    if ~exist('boot_table','var')
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
    end

    for column=p:-1:1
        tmp1         = X(:,column);
        tmp2         = Y(:,column);
        T            = [tmp1(boot_table) tmp2(boot_table)]; % all boostraps at once
        tau          = sign(T(i2, :) - T(i1, :));           % all pair diffences for all columns
        tau          = tau' * tau;                          % between and within sums
        tau          = diag(tau(nboot+1:end,1:nboot))';
        rb(:,column) = tau ./ Adj(column);                  % use the original adjustement for ties
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
        title(sprintf('Bootstrapped correlations \n median=%g',median(rb(:,f))),'FontSize',14);
        xlabel('boot correlations','FontSize',12);ylabel('frequency','FontSize',12)
        axis tight; colormap([.4 .4 1]); box on;
        hold on; plot(median(rb(:,f)),MV/2,'ko','LIneWidth',3)

        if all(~isnan(CI(1:2,f))) % plot CI
            plot(repmat(CI(1,f),MV,1),1:MV,'r','LineWidth',4);
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
