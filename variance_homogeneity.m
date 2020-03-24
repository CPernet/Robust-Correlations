function [h,CI] = variance_homogeneity(x,y,condition)

% use percentile bootstrap to compare variances using 95% CI
%
% FORMAT [h,CI] = variance_homogeneity(x,y)
%        [h,CI] = variance_homogeneity(x,y,condition)
%
% INPUTS  x and y are two vectors of the same length
%         condition if 1 (default) indiceates x and y need to be conditioned on each other
%
% OUTPUTS h indicates if the data have the same variance (0) or not (1)
%         CI is the 95% confidence interval of the difference between variances
%
% Cyril Pernet v1
% ---------------------------------
%  Copyright (C) Corr_toolbox 2012
  
if nargin == 2
    condition = 1;
end

if size(x)~=size(y)
    error('X and Y must have the same size')
end


% computes
nboot = 600;
nm = length(x);
if nm < 40; l=6; u=593;
elseif nm >= 40 && nm < 80; l=7; u=592;
elseif nm >= 80 && nm < 180; l=10; u=589;
elseif nm >= 180 && nm < 250; l=13; u=586;
elseif nm >= 250; l=15; u=584; end
 
% boostrap
table = randi(nm,nm,nboot);
for B=1:nboot
    % resample and standardize
    a = x(table(:,B)); a = (a-mean(a)) / sum((a-mean(a)).^2);
    b = y(table(:,B)); b = (b-mean(b)) / sum((b-mean(b)).^2);
    if condition == 1
        [values,variances]=conditional(a,b);
        Diff(B) = variances(1) - variances(2);
    else
        Diff(B) = var(a) - var(b);
    end
end

Diff = sort(Diff);
CI = [Diff(l+1) Diff(u)];
if sum(isnan(Diff)) ~=0
    adj_nboot = nboot - sum(isnan(Diff));
    adj_l = round((5/100*adj_nboot)/2);
    adj_u = adj_nboot - adj_l;
    CI = [Diff(adj_l+1) Diff(adj_u)];
end


% plot
figure('Name','Heteroscedasticity test');
set(gcf,'Color','w');
k = round(1 + log2(nboot));
[n,p] = hist(Diff,k); 
bar(p,n,1,'FaceColor',[0.5 0.5 1]);
grid on; axis tight; 
ylabel('frequency','Fontsize',12); hold on
plot(repmat(CI(1),max(hist(Diff,k)),1),[1:max(hist(Diff,k))],'r','LineWidth',4);
plot(repmat(CI(2),max(hist(Diff,k)),1),[1:max(hist(Diff,k))],'r','LineWidth',4);
xlabel('differences of conditional variances between X and Y','Fontsize',12);
if CI(1) < 0 && CI(2) > 0
    h = 0;
    if condition == 1
        mytitle = sprintf('Test on conditional variances: \n data are homoscedastic 95%% CI [%g %g]',CI(1),CI(2));
    else
        mytitle = sprintf('Test on variances: \n data are homoscedastic 95%% CI [%g %g]',CI(1),CI(2));
    end
else
    h = 1;
    if condition == 1
        mytitle = sprintf('Test on conditional variances: \n data are heteroscedastic 95%% CI [%g %g]',CI(1),CI(2));
    else
        mytitle = sprintf('Test on variances: \n data are heteroscedastic 95%% CI [%g %g]',CI(1),CI(2));
    end
end
title(mytitle,'Fontsize',14)
    


