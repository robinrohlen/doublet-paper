addpath('functions')

%% Panel A

clearvars; close all;

load('data/sub-13/emg/sub-13_task-trapezoid15percentmvc_run-01_emg.mat')

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

ind=26;

figure;set(gcf,'units','points','position',[563,489,508,196]);
locs_new=interpolate_spike_train(edition.Distimeclean{1}{ind},edition.Pulsetrainclean{1}(ind,:));
locs_new=locs_new-round(find_lag*20);
force_sig=filtfilt(b,a,signal.ref_signal);
hold on;
plot(linspace(0,size(force_sig,2)/2e3,size(force_sig,2)),202./(100/(20+20))*force_sig/max(force_sig),'Color',[150,150,150]/255,'LineWidth',1);
plot_broken_axis(locs_new(1:end-1)/(40e3),1000./(diff(locs_new)/(40)),50,[20 100],[100 201.01],20,[6 57],cmap(ind,:),0,0,0);
hold off;
xlim([6 57]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',12);
xlabel('Time (s)');
ylabel('Firing rate (Hz)');

%% Panel C

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

figure;set(gcf,'units','points','position',[345,489,270,265]);
hold on;
plot(linspace(0,size(edition.Pulsetrainclean{1},2)/2e3,size(edition.Pulsetrainclean{1},2)),edition.Pulsetrainclean{1}(ind,:))
hold off;
xlim([23.7 24.7]);
ylim([-0.3 0.8]);
set(gca,'YTickLabel',[]);
set(gca,'Visible','off');
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',12);

figure;set(gcf,'units','points','position',[345,489,270,265]);
hold on;
plot(linspace(0,size(edition.Pulsetrainclean{1},2)/2e3,size(edition.Pulsetrainclean{1},2)),edition.Pulsetrainclean{1}(ind,:))
hold off;
xlim([43.1 44.1]);
ylim([-0.3 0.8]);
set(gca,'YTickLabel',[]);
set(gca,'Visible','off');
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',12);

%% Panel D

num_subjects=14;
mvc_levels1={'05','10','15','20','25'};

isi_cov_mean_singlet=[];
isi_cov_mean_doublet=[];

for subj_nr=1:num_subjects

    for file_num=1:size(mvc_levels1,2)
        if subj_nr<10
            tmp=['data/sub-0',num2str(subj_nr),'/emg/sub-0',num2str(subj_nr),'_task-trapezoid',mvc_levels1{file_num},'percentmvc_run-01_emg.mat'];
        else
            tmp=['data/sub-',num2str(subj_nr),'/emg/sub-',num2str(subj_nr),'_task-trapezoid',mvc_levels1{file_num},'percentmvc_run-01_emg.mat'];
        end
        load(tmp)

        for k=1:size(edition.Distimeclean{1},2)

            locs=edition.Distimeclean{1}{k};
            fr=1000./(diff(locs)/2);
            indx=find(fr>=100); % find doublets
            if length(indx)>=1
                isi_cov_mean_doublet=cat(2,isi_cov_mean_doublet,std(1000./(diff(locs)/2))./mean(1000./(diff(locs)/2)));
            else
                isi_cov_mean_singlet=cat(2,isi_cov_mean_singlet,std(1000./(diff(locs)/2))./mean(1000./(diff(locs)/2)));
            end

        end
    end
end

figure;set(gcf,'units','points','position',[293,430,413,312]);
histogram(100.*isi_cov_mean_singlet,40);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',12);
xlabel('ISI CoV (%)');
ylabel('Number of motoneurons');
title('Motoneurons with only singlet spikes (N=955)','FontWeight','normal')
xlim([0 300]);
g=gcf;
g.Renderer='painters';

disp([num2str(round(mean(100.*isi_cov_mean_singlet),1)),' +/- ',num2str(round(std(100.*isi_cov_mean_singlet),1))])

%% Panel E

figure;set(gcf,'units','points','position',[293,430,413,312]);
histogram(100.*isi_cov_mean_doublet,100);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',12);
xlabel('ISI CoV (%)');
ylabel('Number of motoneurons');
title('Motoneurons with doublets (N=189)','FontWeight','normal')
xlim([0 300]);
g=gcf;
g.Renderer='painters';

disp([num2str(round(mean(100.*isi_cov_mean_doublet),1)),' +/- ',num2str(round(std(100.*isi_cov_mean_doublet),1))])
