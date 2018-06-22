% Copyright 2018 Corr_toolbox team
% Guillaume Rousselet v1 2018-06-21

%% =======================================================================
%% COMPARE_CORR_IND
%% =======================================================================
%% GENERATE DATA
% Two variables from 2 independent groups: the variables are correlated
% independently in each group and the correlations compared.
rng(21)
n1 = 100;
mu1 = 1;
sigma1 = 2;
r1 = 0.1;
data1 = (mu1 + sigma1*randn(n1,2)) * chol([1 r1; r1 1]);
n2 = 70;
mu2 = 3;
sigma2 = 1;
r2 = 0.5;
data2 = (mu2 + sigma2*randn(n2,2)) * chol([1 r2; r2 1]);
[R1, R2, D, CI, DIST] = compare_corr_ind(data1, data2, 'Spearman');
% [R1, R2, D, CI, DIST] = compare_corr_ind(data1, data2, 'Pearson');
% [R1, R2, D, CI, DIST] = compare_corr_ind(data1, data2, 'bendcorr');
% [R1, R2, D, CI, DIST] = compare_corr_ind(data1, data2, 'Skipped_P');
% [R1, R2, D, CI, DIST] = compare_corr_ind(data1, data2, 'Skipped_S');

%% ILLUSTRATE 2 GROUPS
figure('Color','white', 'NumberTitle', 'off')
subplot(1,2,1)
scatter(data1(:,1),data1(:,2),'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',.5,'MarkerFaceAlpha',.5)
axis square
set(gca, 'FontSize', 14)
title('Group 1','FontSize',24,'FontWeight','bold')
xlabel('Variable 1','FontSize',20,'FontWeight','bold')
ylabel('Variable 2','FontSize',20,'FontWeight','bold')

subplot(1,2,2)
scatter(data2(:,1),data2(:,2),'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',.5,'MarkerFaceAlpha',.5)
axis square
set(gca, 'FontSize', 14)
title('Group 2','FontSize',24,'FontWeight','bold')
xlabel('Variable 1','FontSize',20,'FontWeight','bold')
ylabel('Variable 2','FontSize',20,'FontWeight','bold')

%% RESULTS
fprintf('Correlation in group 1 = %.2f\n',R1)
fprintf('Correlation in group 2 = %.2f\n',R2)
fprintf('Correlation difference = %.2f [%.2f, %.2f]\n',D, CI(1), CI(2))

% Spearman(data1(:,1),data1(:,2))
% Spearman(data2(:,1),data2(:,2))

%% ILLUSTRATE BOOTSTRAP DIFFERENCES

figure('Color','white', 'NumberTitle', 'off')
hold on
histogram(DIST, 20)
ylim = get(gca,'ylim');
plot([D D],ylim, 'Color', [0.9100 0.4100 0.1700], 'LineWidth', 4)
plot([CI(1) CI(1)],ylim, 'Color', [0.9100 0.4100 0.1700], 'LineWidth', 2)
plot([CI(2) CI(2)],ylim, 'Color', [0.9100 0.4100 0.1700], 'LineWidth', 2)
set(gca, 'FontSize', 14)
title('Bootstrap differences','FontSize',24,'FontWeight','bold')
hold off

%% =======================================================================
%% COMPARE_CORR_DEP: OVERLAPPING CASE
%% =======================================================================
%% GENERATE DATA
% Three variables from 1 group: variable 1 is correlated with variable 2
% and variable 3 and the correlations are compared.
rng(2)
Np = 100; % sample size / number of participants
a = randn(Np,1);
b = a * .1 + rand(Np,1);
c = a * .4 + rand(Np,1);
data = [a b c];
[R1, R2, D, CI, DIST] = compare_corr_dep(data, 'Spearman');
% [R1, R2, D, CI, DIST] = compare_corr_dep(data, 'bendcorr');
% [R1, R2, D, CI, DIST] = compare_corr_dep(data, 'Skipped_P');
% [R1, R2, D, CI, DIST] = compare_corr_dep(data, 'Skipped_S');

%% ILLUSTRATE 2 CORRELATIONS
figure('Color','white', 'NumberTitle', 'off')
subplot(1,2,1)
scatter(data(:,1),data(:,2),'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',.5,'MarkerFaceAlpha',.5)
axis square
set(gca, 'FontSize', 14)
title('Correlation 1','FontSize',24,'FontWeight','bold')
xlabel('Variable 1','FontSize',20,'FontWeight','bold')
ylabel('Variable 2','FontSize',20,'FontWeight','bold')

subplot(1,2,2)
scatter(data(:,1),data(:,3),'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',.5,'MarkerFaceAlpha',.5)
axis square
set(gca, 'FontSize', 14)
title('Correlation 2','FontSize',24,'FontWeight','bold')
xlabel('Variable 1','FontSize',20,'FontWeight','bold')
ylabel('Variable 3','FontSize',20,'FontWeight','bold')

%% RESULTS
fprintf('Correlation 1 = %.2f\n',R1)
fprintf('Correlation 2 = %.2f\n',R2)
fprintf('Correlation difference = %.2f [%.2f, %.2f]\n',D, CI(1), CI(2))

% Spearman(data1(:,1),data1(:,2))
% Spearman(data2(:,1),data2(:,2))

%% ILLUSTRATE BOOTSTRAP DIFFERENCES

figure('Color','white', 'NumberTitle', 'off')
hold on
histogram(DIST, 20)
ylim = get(gca,'ylim');
plot([D D],ylim, 'Color', [0.9100 0.4100 0.1700], 'LineWidth', 4)
plot([CI(1) CI(1)],ylim, 'Color', [0.9100 0.4100 0.1700], 'LineWidth', 2)
plot([CI(2) CI(2)],ylim, 'Color', [0.9100 0.4100 0.1700], 'LineWidth', 2)
set(gca, 'FontSize', 14)
title('Bootstrap differences','FontSize',24,'FontWeight','bold')
hold off

%% =======================================================================
%% COMPARE_CORR_DEP: NON-OVERLAPPING CASE
%% =======================================================================
%% GENERATE DATA
% Two variables from 1 group but measured in 2 occasions: variable 1 is correlated with variable 2
% independently in each occasion and the correlations are compared.
rng(2)
Np = 100; % sample size / number of participants
a1 = randn(Np,1);
b1 = a1 * .1 + rand(Np,1);
a2 = a1 + rand(Np,1);
b2 = a2 * .2 + rand(Np,1);
data = [a1 b1 a2 b2];
[R1, R2, D, CI, DIST] = compare_corr_dep(data, 'Spearman');
% [R1, R2, D, CI, DIST] = compare_corr_dep(data, 'bendcorr');
% [R1, R2, D, CI, DIST] = compare_corr_dep(data, 'Skipped_P');
% [R1, R2, D, CI, DIST] = compare_corr_dep(data, 'Skipped_S');

%% ILLUSTRATE 2 CORRELATIONS
figure('Color','white', 'NumberTitle', 'off')
subplot(1,2,1)
scatter(data(:,1),data(:,2),'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',.5,'MarkerFaceAlpha',.5)
axis square
set(gca, 'FontSize', 14)
title('Session 1','FontSize',24,'FontWeight','bold')
xlabel('Variable 1','FontSize',20,'FontWeight','bold')
ylabel('Variable 2','FontSize',20,'FontWeight','bold')

subplot(1,2,2)
scatter(data(:,3),data(:,4),'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',.5,'MarkerFaceAlpha',.5)
axis square
set(gca, 'FontSize', 14)
title('Session 2','FontSize',24,'FontWeight','bold')
xlabel('Variable 1','FontSize',20,'FontWeight','bold')
ylabel('Variable 2','FontSize',20,'FontWeight','bold')

%% RESULTS
fprintf('Correlation in session 1 = %.2f\n',R1)
fprintf('Correlation in session 2 = %.2f\n',R2)
fprintf('Correlation difference = %.2f [%.2f, %.2f]\n',D, CI(1), CI(2))

% Spearman(data1(:,1),data1(:,2))
% Spearman(data2(:,1),data2(:,2))

%% ILLUSTRATE BOOTSTRAP DIFFERENCES

figure('Color','white', 'NumberTitle', 'off')
hold on
histogram(DIST, 20)
ylim = get(gca,'ylim');
plot([D D],ylim, 'Color', [0.9100 0.4100 0.1700], 'LineWidth', 4)
plot([CI(1) CI(1)],ylim, 'Color', [0.9100 0.4100 0.1700], 'LineWidth', 2)
plot([CI(2) CI(2)],ylim, 'Color', [0.9100 0.4100 0.1700], 'LineWidth', 2)
set(gca, 'FontSize', 14)
title('Bootstrap differences','FontSize',24,'FontWeight','bold')
hold off