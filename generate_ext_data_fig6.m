addpath('functions')

%% Panel A

clearvars; close all;

load('data/sub-01/emg/sub-01_task-trapezoid25percentmvc_run-01_emg.mat')

subj_coding=1:24;
subj_nr=1;

[b,a]=butter(3,10/(2e3/2),'low'); % force filtering

[r,lags]=xcorr(-signal.data(69,:)',-signal.data(75,:));
find_lag=lags(r==max(r));

% rearrange spike trains according to recruitment threshold
% use average of first 10 firings
rte_avg=10;
i=1;

tmp=[];
for j=1:size(edition.Distimeclean{i},2)
    tmp=[tmp median(edition.Distimeclean{1,i}{1,j}(1:rte_avg))];
end

[~,indx]=sort(tmp,'ascend');

edition.Distimeclean{i}=edition.Distimeclean{i}(indx);
edition.Pulsetrainclean{i}=edition.Pulsetrain{i}(indx,:);

tmp=0;
fs=2e3;
i=1;
figure;set(gcf,'units','points','position',[563,489,508,196]);
cmap=turbo(size(edition.Distimeclean{i},2));

for j=1:size(edition.Distimeclean{i},2)
    locs_new=interpolate_spike_train(edition.Distimeclean{1}{j},edition.Pulsetrainclean{1}(j,:));
    locs_new=locs_new-round(find_lag*20);
    for k=1:size(locs_new,2)
        hold on;
        plot([locs_new(k) locs_new(k)]./(fs*20),[j-0.5 j+0.5],'Color',cmap(j,:),'LineWidth',0.5);
        hold off;
    end
end

xlim([6 57])
ylim([0 50]);
yticks([0:10:50]);
ylim_tmp=get(gca,'ylim');

force_sig=filtfilt(b,a,signal.path);
hold on;
plot(linspace(0,size(force_sig,2)/2e3,size(force_sig,2)),ylim_tmp(2)*force_sig/max(force_sig),'Color',[150,150,150]/255,'LineWidth',1);
hold off;

set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',12);

xlabel('Time (s)');
ylabel('Motoneuron #');

g=gcf;
g.Renderer='painters';

%% Panels B-C

ind=41;

figure;set(gcf,'units','points','position',[563,489,508,196]);
locs_new=interpolate_spike_train(edition.Distimeclean{1}{ind},edition.Pulsetrainclean{1}(ind,:));
locs_new=locs_new-round(find_lag*20);
force_sig=filtfilt(b,a,signal.path);
hold on;
plot(linspace(0,size(force_sig,2)/2e3,size(force_sig,2)),100./(50/(20+20))*force_sig/max(force_sig),'Color',[150,150,150]/255,'LineWidth',1);
plot_broken_axis(locs_new(1:end-1)/(40e3),1000./(diff(locs_new)/(40)),50,[20 50],[50 100],20,[6 57],cmap(ind,:),0,0,0);
hold off;
xlim([6 57]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',12);
xlabel('Time (s)');
ylabel('Firing rate (Hz)');

figure;set(gcf,'units','points','position',[563,489,508,196]);
hold on;
plot(linspace(0,size(edition.Pulsetrainclean{1},2)/2e3,size(edition.Pulsetrainclean{1},2)),edition.Pulsetrainclean{1}(ind,:))
hold off;
xlim([6 57]);
ylim([-0.5 1.5]);
set(gca,'YTickLabel',[]);
set(gca,'Visible','off');
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',12);

figure;set(gcf,'units','points','position',[563,489,508,196]);
hold on;
plot(linspace(0,size(edition.Pulsetrainclean{1},2)/2e3,size(edition.Pulsetrainclean{1},2)),edition.Pulsetrainclean{1}(ind,:))
hold off;
xlim([28.2 29.2]);
ylim([-0.3 0.8]);
set(gca,'YTickLabel',[]);
set(gca,'Visible','off');
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',12);

%% Panels D-E

ind=43;

figure;set(gcf,'units','points','position',[563,489,508,196]);
locs_new=interpolate_spike_train(edition.Distimeclean{1}{ind},edition.Pulsetrainclean{1}(ind,:));
locs_new=locs_new-round(find_lag*20);
force_sig=filtfilt(b,a,signal.path);
hold on;
plot(linspace(0,size(force_sig,2)/2e3,size(force_sig,2)),100./(50/(20+20))*force_sig/max(force_sig),'Color',[150,150,150]/255,'LineWidth',1);
plot_broken_axis(locs_new(1:end-1)/(40e3),1000./(diff(locs_new)/(40)),50,[20 50],[50 100],20,[6 57],cmap(ind,:),0,0,0);
hold off;
xlim([6 57]);
ylim_tmp=get(gca,'ylim');
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',12);
xlabel('Time (s)');
ylabel('Firing rate (Hz)');

figure;set(gcf,'units','points','position',[563,489,508,196]);
hold on;
plot(linspace(0,size(edition.Pulsetrainclean{1},2)/2e3,size(edition.Pulsetrainclean{1},2)),edition.Pulsetrainclean{1}(ind,:))
hold off;
xlim([6 57]);
ylim([-0.5 1.5]);
set(gca,'YTickLabel',[]);
set(gca,'Visible','off');
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',12);

figure;set(gcf,'units','points','position',[563,489,508,196]);
hold on;
plot(linspace(0,size(edition.Pulsetrainclean{1},2)/2e3,size(edition.Pulsetrainclean{1},2)),edition.Pulsetrainclean{1}(ind,:))
hold off;
xlim([42.75 43.75]);
ylim([-0.3 0.8]);
set(gca,'YTickLabel',[]);
set(gca,'Visible','off');
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',12);