addpath('functions');

%% Panel B

clearvars; close all;

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
edition.Pulsetrainclean{i}=edition.Pulsetrainclean{i}(indx,:);

tmp=0;
fs=2e3;
i=1;
figure;set(gcf,'units','points','position',[179,465,739,282]);
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

xlim([4.3 44.3]);
ylim([0 30]);
yticks([0:10:30]);
ylim_tmp=get(gca,'ylim');
ylabel('Motoneuron #')

force_sig=filtfilt(b,a,signal.ref_signal);
hold on;
plot(linspace(0,size(force_sig,2)/2e3,size(force_sig,2)),ylim_tmp(2)*force_sig/max(force_sig),'Color',[150,150,150]/255,'LineWidth',1);
hold off;

set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',16);

xlim([4.3 44.3]);
ylim([0 30]);
yticks([0:10:30]);
ylim_tmp=get(gca,'ylim');

g=gcf;
g.Renderer='painters';

%% Panels D-G
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
ylabel('Firing rate (Hz)');
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',16);
g=gcf;
g.Renderer='painters';

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
ylabel('Firing rate (Hz)')
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',16);
g=gcf;
g.Renderer='painters';

% Zoom in phase advance
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
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',20);
xlim([24.7277 28.5781])
g=gcf;
g.Renderer='painters';

figure;set(gcf,'units','points','position',[591,383,411,325]);
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
xlim([24.7277 28.5781])
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',20);
g=gcf;
g.Renderer='painters';

%% Panel C

load_data=1;

num_subjects=14;
doublet_thresh=100; % Hz

if load_data==0
    num_MUs_vec=zeros(1,num_subjects);
    prop_doublets=zeros(num_subjects,2);

    for i=1:num_subjects
        disp(i)
        if i<10
            tmp=['data/sub-0',num2str(i),'/emg/sub-0',num2str(i),'_task-sinusoid2hz_run-01_emg.mat'];
        else
            tmp=['data/sub-',num2str(i),'/emg/sub-',num2str(i),'_task-sinusoid2hz_run-01_emg.mat'];
        end
        load(tmp)
        num_MUs_vec(i)=size(edition.Distimeclean{1},2);
        doublet_iter=0;
        for j=1:size(edition.Distimeclean{1},2)
            fr=1000./(diff(edition.Distimeclean{1}{j})/2);
            doublet_iter=doublet_iter+any(fr>=50);
        end
        prop_doublets(i,:)=[doublet_iter/num_MUs_vec(i) (num_MUs_vec(i)-doublet_iter)/num_MUs_vec(i)];
    end
    save('data/figure-files/proportion_initial_doublets.mat','prop_doublets','num_MUs_vec')
else
    load('data/figure-files/proportion_initial_doublets.mat');
end

cmap=turbo(27);
MUind=[1 24];

figure;set(gcf,'units','points','position',[472,285,399,475]);

yyaxis left
hold on;
f=bar(100.*prop_doublets,'stacked','BarWidth', 1);
hold off;
f(1).FaceColor = cmap(MUind(1),:);
f(1).FaceAlpha = 0.7;
f(2).FaceColor = cmap(MUind(2),:);
f(2).FaceAlpha = 0.7;
set(gca,'XTick',[]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
xlim([0.5 14.5]);
ylim([0 100]);
yticks(0:25:100);
ax = gca;
ax.YColor = cmap(2,:);

subj_nr_vec=1:14;

col=[247,165,30]/255;

yyaxis right
hold on;
plot(subj_nr_vec,num_MUs_vec,'-','Color',col,'MarkerFaceColor',col,'MarkerSize',12,'MarkerEdgeColor',col,'LineWidth',1.5);
plot(subj_nr_vec,num_MUs_vec,'o','Color',col,'MarkerFaceColor',col,'MarkerSize',12,'MarkerEdgeColor',col,'LineWidth',1.5);

hold off;

ax = gca;
ax.YColor = col;
ylim([0 40]);
yticks(0:10:40);

g=gcf;
g.Renderer='painters';
