function corr_normplot(x,y)

% plots for the bivariate data set defined by [x,y] 
%
% FORMAT corr_normplot(x,y)
%
% INPUTS x and y are two vectors of the same length
%
% Cyril Pernet v1 20/06/2012
% -----------------------------
% Copyright (C) Corr_toolbox 2012

if size(x)~=size(y)
    error('X and Y must have the same size')
end

[r c] = size(x);
if r == 1 && c > 1
    x = x'; 
    y = y';
elseif r > 1 && c > 1
    error('X and Y must be 2 vectors, more than 1 column/row detected')
end


figure('Name','Histrograms and scatter plot')
set(gcf,'Color','w');

% 1st univariate histogram
subplot(3,5,2:3);
[nu,x1,h1,xp,yp]=univar(x); 
bar(x1,nu/(length(x)*h1),1,'FaceColor',[0.5 0.5 1]);
v = max(yp) + 0.02*max(yp);
grid on; axis([min(x)-1/10*min(x) max(x)+1/10*max(x) 0 v]); hold on
plot(xp,yp,'r','LineWidth',3); title('Density histogram for X','Fontsize',14); 
ylabel('Freq.','FontSize',12); xlabel('X.','FontSize',12)

% 2nd univariate histogram
subplot(3,5,[6 11]);
[nu,x2,h2,xp,yp]=univar(y);
bar(x2,nu/(length(y)*h2),1,'FaceColor',[0.5 0.5 1]);
v = max(yp) + 0.02*max(yp);
grid on; axis([min(y)-1/10*min(y) max(y)+1/10*max(y) 0 v]); hold on
plot(xp,yp,'r','LineWidth',3); view(-90,90) 
title('Density histogram for Y','Fontsize',14); 
ylabel('Freq.','FontSize',12); xlabel('Y.','FontSize',12); 
drawnow

% scatter plot
subplot(3,5,[7 8 12 13]);
scatter(x,y,100,'filled'); grid on
xlabel('x','Fontsize',12); ylabel('y','Fontsize',12);
axis([min(x)-1/10*min(x) max(x)+1/10*max(x) min(y)-1/10*min(y) max(y)+1/10*max(y)])
title('Scatter plot','Fontsize',14);
drawnow

% joint histogram
subplot(3,5,[9 10 14 15]);
k = round(1 + log2(length(x)));
hist3([x y],[k k],'FaceAlpha',.65);
xlabel('X'); ylabel('Y'); title('Bivariate histogram','Fontsize',14)
set(gcf,'renderer','opengl');
drawnow
try
    set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
end



