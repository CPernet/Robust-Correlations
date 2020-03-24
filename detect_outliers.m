function outliers = detect_outliers(X,Y)

% find univariate and multivariate outliers based on the MAD
%
% FORMAT outliers = detect_outliers(X,Y)
%
% INPUTS X and Y are two vectors of the same length
%
% OUTPUT outliers is a structure
%        outliers.univariate.X retuns indices of univariate outliers in X
%        outliers.univariate.Y retuns indices of univariate outliers in Y
%        outliers.bivariate returns bivariate outliers 
%
% for univariate outliers simply looks if the data points are located above
% the 75th quantile of the MADN
% for multivariate outliers, computes all the distances to the center of
% the cloud and check if the distances are above the 75th quantile of the
% MADN of distances - see Wilcox 2012 and the skip correlation function
%
% Cyril Pernet v1 23/07/2012
% ---------------------------------
%  Copyright (C) Corr_toolbox 2012



%% check inputs

if nargin <2 
    error('not enough input arguments');
elseif nargin > 3
    error('too many input arguments');
end

% transpose if x or y are not in column
if size(X,1) == 1 && size(X,2) > 1; x = x'; end
if size(Y,1) == 1 && size(Y,2) > 1; y = y'; end

if numel(size(X))>2 || numel(size(Y))>2
    error('only taking vectors as input')
end

if numel(X) ~= numel(Y)
    error('vector must be of the same length')
end

%% univariate outliers

[n,p] = size(X);
if n > 9
    Xoutliers = madmedianrule(X,1);
    Youtliers = madmedianrule(Y,1);
else
    Xoutliers = madmedianrule(X,2);
    Youtliers = madmedianrule(Y,2);
end

% figure
figure('Name','Outlier detection'); set(gcf,'Color','w');
if sum([Xoutliers;Youtliers]) > 0
    subplot(1,2,1); scatter(X,Y,99,'o','filled'); grid on
    xlabel('X','Fontsize',12); ylabel('Y','Fontsize',12);
    title('Data scatter plot','Fontsize',16); hold on
    common = intersect(find(Xoutliers),find(Youtliers));
    if ~isempty(common)
        scatter(X(common),Y(common),100,'k','filled')
        if common ~= find(Xoutliers)
            scatter(X(find(Xoutliers)),Y(find(Youtliers)),100,'r','filled')
            scatter(X(find(Youtliers)),Y(find(Youtliers)),100,'g','filled')
            legend('data','outliers in X','outliers in Y','outliers in X and in Y')
        else
            legend('data','outliers in X and in Y')
        end
        outliers.univariate.X = Xoutliers;
        outliers.univariate.Y = Youtliers;
    elseif ~isempty(find(Xoutliers)) && ~isempty(find(Youtliers))
        scatter(X(find(Xoutliers)),Y(find(Xoutliers)),100,'r','filled')
        scatter(X(find(Youtliers)),Y(find(Youtliers)),100,'g','filled')
        legend('data','outliers in X','outliers in Y');
        outliers.univariate.X = Xoutliers;
        outliers.univariate.Y = Youtliers;
    elseif ~isempty(find(Xoutliers)) && isempty(find(Youtliers))
        scatter(X(find(Xoutliers)),Y(find(Xoutliers)),100,'r','filled')
        legend('data','outliers in X')
        outliers.univariate.X = Xoutliers;
        outliers.univariate.Y = 'none';
    elseif ~isempty(find(Youtliers)) && isempty(find(Xoutliers))
        scatter(X(find(Youtliers)),Y(find(Youtliers)),100,'g','filled')
        legend('data','outliers in Y')
        outliers.univariate.X = 'none';
        outliers.univariate.Y = Youtliers;
    end
    disp(' '); 
    fprintf('%g outliers found in X \n',sum(Xoutliers))
    fprintf('%g outliers found in Y \n',sum(Youtliers))
else
    subplot(1,2,1); scatter(X,Y,99,'o','fill'); grid on
    xlabel('x','Fontsize',12); ylabel('y','Fontsize',12);
    title('Data scatter plot - no outliers','Fontsize',16);
    fprintf('no outlier found in X or Y',sum(Xoutliers))
    outliers.univariate = 'none';
end
axis([min(X)-10/100*min(X) max(X)+10/100*max(X) min(Y) max(Y)])


%% Multivariate outlier

% get the centre of the bivariate distributions
X = [X Y]; gval = sqrt(chi2inv(0.975,rank(X))); 
result=mcdcov(X,'cor',1,'plots',0,'h',floor((n+size(X,2)*2+1)/2));
center = result.center;
flag = zeros(n,1);
     
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
        [out,value] = madmedianrule(dis,2);
        record{i} = dis > (median(dis)+gval.*value);
    end
end

flag = sum(cell2mat(record),2); % if any point is flagged

if sum(flag)==0
    outid=[];
else
    flag=(flag>=1);
    outid=vec(flag);
end

% figure
subplot(1,2,2); 
a = X(vec(~flag),1);
b = X(vec(~flag),2);
scatter(a,b,100,'b','fill');
grid on; hold on;
hh = lsline; set(hh,'Color','r','LineWidth',4);
xlabel('X','Fontsize',12); ylabel('Y','Fontsize',12);
title('Bivariate outliers','Fontsize',16);
scatter(X(outid,1),X(outid,2),100,'r','filled');
try
    A = min(min(a, X(outid,1))); B = max(max(a, X(outid,1)));
    C = min(min(b, X(outid,2))); D = max(max(b, X(outid,2)));
    axis([A-10/100*A B+10/100*B C D])
catch
    axis([min(X(:,1))-10/100*min(X(:,1)) max(X(:,1))+10/100*max(X(:,1)) min(X(:,2)) max(X(:,2))])
end
if isempty(outid)
    fprintf('no bivariate outliers detected');
    outliers.bivariate = 'none';
else
    fprintf('pair(s) %g found as bivariate outliers \n',outid);
    tmp = zeros(n,1); tmp(outid) = 1;
    outliers.bivariate = tmp;
end

