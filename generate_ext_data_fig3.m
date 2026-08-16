addpath('functions')

%% Panel A

load('data/sub-12/emg/sub-12_task-sinusoid1hz_run-01_emg.mat')

subj_coding=1:24;
subj_nr=12;

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

fs=2e3;
i=1;

cmap=turbo(size(edition.Distimeclean{i},2));

MUind=[1 24];

figure;set(gcf,'units','points','position',[591,383,411,325]);
ind=MUind(1);
locs_new=interpolate_spike_train(edition.Distimeclean{1}{ind},edition.Pulsetrainclean{1}(ind,:));
locs_new=locs_new-round(find_lag*20);
force_sig=filtfilt(b,a,signal.ref_signal);
hold on;
plot(linspace(0,size(force_sig,2)/2e3,size(force_sig,2)),30*force_sig/max(force_sig),'Color',[150,150,150]/255,'LineWidth',1);
scatter(locs_new(1:end-1)/(40e3),1000./(diff(locs_new)/(40)),'MarkerFaceColor',cmap(ind,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75);
hold off;
ylim([0 30]);
xlim([4.3 44.3]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',20);
xlabel('Time (s)');
ylabel('Firing rate (Hz)');
g=gcf;
g.Renderer='painters';

%% Panel B

figure;set(gcf,'units','points','position',[591,383,411,325]);
ind=MUind(1);
locs_new=interpolate_spike_train(edition.Distimeclean{1}{ind},edition.Pulsetrainclean{1}(ind,:));
locs_new=locs_new-round(find_lag*20);
force_sig=filtfilt(b,a,signal.ref_signal);
hold on;
plot(linspace(0,size(force_sig,2)/2e3,size(force_sig,2)),30*force_sig/max(force_sig),'Color',[150,150,150]/255,'LineWidth',1);
plot(locs_new(1:end-1)/(40e3),1000./(diff(locs_new)/(40)),':','Color',cmap(ind,:));
scatter(locs_new(1:end-1)/(40e3),1000./(diff(locs_new)/(40)),'MarkerFaceColor',cmap(ind,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75);
hold off;
ylim([0 30]);
set(gca,'YTickLabel',[]);
h = gca;
h.XAxis.Visible = 'off';
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',20);
set(gca,'Visible','off')
xlim([24.7277 28.5781])

%% Panel C

figure;set(gcf,'units','points','position',[591,383,411,325]);
ind=MUind(1);
sig_tmp=edition.Pulsetrainclean{1}(ind,find_lag:end);
hold on;
plot(linspace(0,length(sig_tmp)/2e3,length(sig_tmp)),sig_tmp)
hold off;
xlim([24.7277 28.5781]);
ylim([-0.3 0.8]);
set(gca,'YTickLabel',[]);
set(gca,'Visible','off');
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);

%% Panel D

figure;set(gcf,'units','points','position',[591,383,411,325]);
ind=MUind(2);
locs_new=interpolate_spike_train(edition.Distimeclean{1}{ind},edition.Pulsetrainclean{1}(ind,:));
locs_new=locs_new-round(find_lag*20);
force_sig=filtfilt(b,a,signal.ref_signal);
hold on;
plot(linspace(0,size(force_sig,2)/2e3,size(force_sig,2)),200./(100/(30+30))*force_sig/max(force_sig),'Color',[150,150,150]/255,'LineWidth',1);
plot_broken_axis(locs_new(1:end-1)/(40e3),1000./(diff(locs_new)/(40)),50,[30 100],[100 200],30,[4.3 44.3],cmap(ind,:),0,0);
hold off;
xlim([4.3 44.3]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',20);
xlabel('Time (s)');
ylabel('Firing rate (Hz)')
g=gcf;
g.Renderer='painters';

%% Panel E

figure;set(gcf,'units','points','position',[345,489,406,265]);
ind=MUind(2);
locs_new=interpolate_spike_train(edition.Distimeclean{1}{ind},edition.Pulsetrainclean{1}(ind,:));
locs_new=locs_new-round(find_lag*20);
force_sig=filtfilt(b,a,signal.ref_signal);
hold on;
plot(linspace(0,size(force_sig,2)/2e3,size(force_sig,2)),200./(100/(30+30))*force_sig/max(force_sig),'Color',[150,150,150]/255,'LineWidth',1);
[x,y]=plot_broken_axis(locs_new(1:end-1)/(40e3),1000./(diff(locs_new)/(40)),50,[30 100],[100 200],30,[4.3 44.3],cmap(ind,:));
hold off;
hold on;
plot(x,y,':o','Color',cmap(ind,:),'LineWidth',1);
hold off;
set(gca,'YTickLabel',[]);
xlim([24.7277 28.5781])
h = gca;
h.XAxis.Visible = 'off';
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',20);
set(gca,'Visible','off')
g=gcf;
g.Renderer='painters';

%% Panel F

figure;set(gcf,'units','points','position',[345,489,406,265]);
ind=MUind(2);
sig_tmp=edition.Pulsetrainclean{1}(ind,find_lag:end);
hold on;
plot(linspace(0,length(sig_tmp)/2e3,length(sig_tmp)),sig_tmp)
hold off;
xlim([24.7277 28.5781]);
ylim([-0.3 0.8]);
set(gca,'YTickLabel',[]);
set(gca,'Visible','off');
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);

g=gcf;
g.Renderer='painters';