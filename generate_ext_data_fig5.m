addpath('functions')

%% Panel A

clearvars; close all;

load('data/sub-02/emg/sub-02_task-trapezoid05percentmvc_run-01_emg.mat')

subj_coding=1:24;
subj_nr=13;

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
figure;set(gcf,'units','points','position',[345,586,627,168]);
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
ylim([0 30]);
yticks([0:10:30]);
ylim_tmp=get(gca,'ylim');

force_sig=filtfilt(b,a,signal.ref_signal);
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

%% Panel B

MUind=[2 8 20 3 11 16];

for i=1:3
    ind=MUind(i);

    figure;set(gcf,'units','points','position',[345,489,310,265]);
    locs_new=interpolate_spike_train(edition.Distimeclean{1}{ind},edition.Pulsetrainclean{1}(ind,:));
    locs_new=locs_new-round(find_lag*20);
    force_sig=filtfilt(b,a,signal.ref_signal);
    hold on;
    plot(linspace(0,size(force_sig,2)/2e3,size(force_sig,2)),30*force_sig/max(force_sig),'Color',[150,150,150]/255,'LineWidth',1);
    scatter(locs_new(1:end-1)/(40e3),1000./(diff(locs_new)/(40)),'MarkerFaceColor',cmap(ind,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75);
    hold off;
    xlim([6 57]);
    set(gca,'TickDir','out');
    set(gcf,'color','w');
    set(gca,'FontSize',20);
    xlabel('Time (s)');
    ylabel('Firing rate (Hz)');
end

%% Panel C

MUind=[2 8 12 3 11 16];

for i=4:6
    ind=MUind(i);

    figure;set(gcf,'units','points','position',[345,489,310,265]);
    locs_new=interpolate_spike_train(edition.Distimeclean{1}{ind},edition.Pulsetrainclean{1}(ind,:));
    locs_new=locs_new-round(find_lag*20);
    force_sig=filtfilt(b,a,signal.ref_signal);
    hold on;
    plot(linspace(0,size(force_sig,2)/2e3,size(force_sig,2)),400./(200/(20+20))*force_sig/max(force_sig),'Color',[150,150,150]/255,'LineWidth',1);
    plot_broken_axis(locs_new(1:end-1)/(40e3),1000./(diff(locs_new)/(40)),50,[20 200],[200 400],20,[6 57],cmap(ind,:),0,0,0);
    hold off;
    xlim([6 57]);
    set(gca,'TickDir','out');
    set(gcf,'color','w');
    set(gca,'FontSize',20);
    xlabel('Time (s)');
    ylabel('Firing rate (Hz)');
end