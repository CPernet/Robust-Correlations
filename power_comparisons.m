%% let's compare power curves
clear
load power_1
pow   = sum(pvals<.05,4);% translated in error rate
pow   = binofit(pow(:),1000); % avg and CI
pow   = reshape(pow,7,7,20); % reshape
powho =(pow(:,:,10:-1:1) + pow(:,:,11:20))./2;
pvaho = pvals;
load power_2
pow   = sum(pvals<.05,4);% translated in error rate
pow   = binofit(pow(:),1000); % avg and CI
pow   = reshape(pow,7,7,20); % reshape
powhe =(pow(:,:,10:-1:1) + pow(:,:,11:20))./2;
pvahe = pvals;
load power_3
pow   = sum(pvals<.05,5);% translated in error rate
pow   = binofit(pow(:),1000); % avg and CI
pow   = reshape(pow,7,9,7,20); % reshape
powou = squeeze(mean(pow,2)); % average over distances
powou = (powou(:,:,10:-1:1) + powou(:,:,11:20))./2;
pvaou = pvals;
clear pow pvals


%% what is the difference due to distributons
% avgho = mean(mean(powho,3),2);
% avghe = mean(mean(powhe,3),2);
% avgou = mean(mean(powou,3),2);

M = [];
for test = 1:7
    A = squeeze(powho(test,:,:));
    B = squeeze(powhe(test,:,:));
    C = squeeze(powou(test,:,:));
    M = [M A(:)-B(:) A(:)-C(:)];
end
[h,CI] = rst_1ttest(M,'mean');

% also worth checking the effect of outlier distance
pow   = sum(pvaou<.05,5);% translated in error rate
pow   = binofit(pow(:),1000); % avg and CI
pow   = reshape(pow,7,9,7,20); % reshape
powou = (pow(:,:,:,10:-1:1) + pow(:,:,:,11:20))./2;

% figure
idx = 1;
C = flipud(rst_colour_maps(7));
distance = -4:4;
Name = {'Pearson','HC4','Spearman','HC4','Kendall','30% bend','HC4'};
figure('Name','Effect of distance')
for test = 1:7
    for d=1:9
        subplot(7,9,idx); hold on
        for ss=1:7
            plot(rho(11:20),squeeze(powou(test,d,ss,:)),'LineWidth',2,'Color',C(ss,:));
            axis([0.09 1.01 -0.01 1.01])
        end
        if test ==1; title(sprintf('distance=%g',distance(d))); end
        if any(idx == [1:9:91]); ylabel(sprintf('%s\npower',Name{test})); end
        if test ==7; xlabel('pop rho'); end
        grid on; box on; idx = idx+1;
    end
end


%% what is the difference due to variance estimator
figure('Name','Power curve differences')
ci = NaN(3,3,2,7,10);
idx = 1;
for dist = 1:3
    for test = [1 3 6]
        
        if dist == 1 && test == 1
            tmp = squeeze(pvaho([1 2],:,:,:));
            A = squeeze(powho(1,:,:));
            B = squeeze(powho(2,:,:));
        elseif dist == 1 && test == 3
            tmp = squeeze(pvaho([3 4],:,:,:));
            A = squeeze(powho(3,:,:));
            B = squeeze(powho(4,:,:));
        elseif dist == 1 && test == 6
            tmp = squeeze(pvaho([6 7],:,:,:));
            A = squeeze(powho(6,:,:));
            B = squeeze(powho(7,:,:));
        elseif dist == 2 && test == 1
            tmp = squeeze(pvahe([1 2],:,:,:));
            A = squeeze(powhe(1,:,:));
            B = squeeze(powhe(2,:,:));
        elseif dist == 2 && test == 3
            tmp = squeeze(pvahe([3 4],:,:,:));
            A = squeeze(powhe(3,:,:));
            B = squeeze(powhe(4,:,:));
        elseif dist == 2 && test == 6
            tmp = squeeze(pvahe([6 7],:,:,:));
            A = squeeze(powhe(6,:,:));
            B = squeeze(powhe(7,:,:));
        elseif dist == 3 && test == 1
            tmp = squeeze(pvaou([1 2],:,:,:,:));
            A = squeeze(mean(powou(1,:,:,:),2));
            B = squeeze(mean(powou(2,:,:,:),2));
        elseif dist == 3 && test == 3
            tmp = squeeze(pvaou([3 4],:,:,:,:));
            A = squeeze(mean(powou(3,:,:,:),2));
            B = squeeze(mean(powou(4,:,:,:),2));
        elseif dist == 3 && test == 6
            tmp = squeeze(pvaou([6 7],:,:,:,:));
            A = squeeze(mean(powou(6,:,:,:),2));
            B = squeeze(mean(powou(7,:,:,:),2));
        end
        
        subplot(3,3,idx); hold on
        for boot = 1000:-1:1
            if dist == 3
                pow   = sum(tmp(:,:,:,:,randi(1000,1000,1))<.05,5);
                pow   = binofit(pow(:),1000); % avg and CI
                pow   = squeeze(mean(reshape(pow,2,9,7,20),2));  % reshape
            else
                pow   = sum(tmp(:,:,:,randi(1000,1000,1))<.05,4);
                pow   = binofit(pow(:),1000); % avg and CI
                pow   = reshape(pow,2,7,20);  % reshape
            end
            powb(:,:,:,boot) =(pow(:,:,10:-1:1) + pow(:,:,11:20))./2;
        end
        powd = sort(squeeze(powb(1,:,:,:) - powb(2,:,:,:)),3);
        ci(dist,test,1,:,:) = powd(:,:,5);
        ci(dist,test,2,:,:) = powd(:,:,995);
        for ss=1:7
            plot(rho(11:20),(A(ss,:)-B(ss,:))','LineWidth',2,'Color',C(ss,:)); hold on
            xpoints = [powd(ss,:,5) fliplr(powd(ss,:,995))];
            fillhandle=fill([rho(11:20) rho(20:-1:11)],xpoints,C(ss,:));
            set(fillhandle,'LineWidth',2,'EdgeColor',C(ss,:),'FaceAlpha',0.2,'EdgeAlpha',0.2);%set edge color
        end
        
        grid on; box on; idx = idx+1; axis([0.09 1.01 -0.005 0.4])
        if dist == 1
            if test == 1
                title('Pearson - HC4')
            elseif test == 3
                title('Spearman - HC4')
            else
                title('30% Bend - HC4')
            end
        end
        
        if test == 1
            if dist == 1
                ylabel(sprintf('Normal data\npower difference'))
            elseif dist == 2
                ylabel(sprintf('Heteroscedastic data\npower difference'))
            else
                ylabel(sprintf('Normal+outlier data\npower difference'))
            end
        end
    end
end
ci(:,[2 4 5],:,:,:) = [];
