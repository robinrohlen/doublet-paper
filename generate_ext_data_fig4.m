addpath('functions')

%% Panel A

clearvars; close all;

load('data/sub-13/emg/sub-13_task-trapezoid20percentmvc_run-01_emg.mat')

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

i=1;

cmap=turbo(size(edition.Distimeclean{i},2));

% Fit lines to firing rates
locs=edition.Distimeclean{1};
pulsetrain=edition.Pulsetrain{1};
fs=signal.fsamp;

signal_len = size(signal.data,2);

t_sig_s=cell(1,size(locs,2));
fr_sig_s=cell(1,size(locs,2));
t_sig_d=cell(1,size(locs,2));
fr_sig_d=cell(1,size(locs,2));
doublet_ind=zeros(1,size(locs,2));

for i=1:size(locs,2)
    locs_i=interpolate_spike_train(locs{i},pulsetrain(i,:));
    t=locs_i(1:end-1)/(fs*20);
    fr=1./(diff(locs_i)/(fs*20));

    fr_s=fr(find(fr<50));
    t_s=t(find(fr<50));
    fr_d=fr(find(fr>=100));
    t_d=t(find(fr>=100));

    locs_s=cumsum(fs./fr_s)+t_s(1)*fs;
    if ~isempty(fr_d)
        locs_d=cumsum(fs./fr_d)+t_d(1)*fs;
        st_d=zeros(1,signal_len);
        st_d(round(locs_d))=1;
        if length(fr_d)>=3
            doublet_ind(i)=1;

            tmp_ind=find(diff(t_d)<=1.5);
            tmp_ind=unique([1 tmp_ind]);
            try
                tmp_ind=[tmp_ind tmp_ind(end)+1];
                fr_d_tmp=fr_d(tmp_ind);
                t_d_tmp=t_d(tmp_ind);
            catch
                fr_d_tmp=fr_d(tmp_ind);
                t_d_tmp=t_d(tmp_ind);
            end
            fr_d_tmp=fr_d_tmp(find(fr_d_tmp<=300));
            t_d_tmp=t_d_tmp(find(fr_d_tmp<=300));
            f = @(b,t_d_tmp) b(1).*exp(b(2).*t_d_tmp)+b(3);
            B = fminsearch(@(b) norm(fr_d_tmp - f(b,t_d_tmp)), [1000 -0.1 150]);

            % add interp with spline
            t_sig_d{i}=linspace(t_d_tmp(1),t_d_tmp(end),size(t_d_tmp,2)*10);
            fr_sig_d{i}=interp1(t_d_tmp,f(B,t_d_tmp),t_sig_d{i},'spline');
        end
    end

    st_s=zeros(1,signal_len);
    st_s(round(locs_s))=1;

    t_sig_s{i}=linspace(0,signal_len/fs,signal_len);
    tmp_ind=find(diff(t_s)<=1);
    t_sig_s{i}=t_sig_s{i}(round(fs*t_s(tmp_ind(1))):round(fs*t_s(tmp_ind(end))));
    fr_sig_s{i}=conv(st_s,hann(2*fs),'same'); % todo: modify to remove baseline
    fr_sig_s{i}=fr_sig_s{i}(round(fs*t_s(tmp_ind(1))):round(fs*t_s(tmp_ind(end))));
end

iter_s=0;
iter_d=0;

MUind=[1 3 16 26];

ind=MUind(3);

figure;set(gcf,'units','points','position',[345,489,406,265]);
locs_new=interpolate_spike_train(edition.Distimeclean{1}{ind},edition.Pulsetrainclean{1}(ind,:));
locs_new=locs_new-round(find_lag*20);
force_sig=filtfilt(b,a,signal.ref_signal);
hold on;
plot(linspace(0,size(force_sig,2)/2e3,size(force_sig,2)),200./(100/(20+20))*force_sig/max(force_sig),'Color',[150,150,150]/255,'LineWidth',1);
plot_broken_axis(locs_new(1:end-1)/(40e3),1000./(diff(locs_new)/(40)),50,[20 100],[100 200],20,[6 57],cmap(ind,:),0,0);
plot_broken_axis(t_sig_s{ind}-2,fr_sig_s{ind},50,[20 100],[100 200],20,[6 57],cmap(ind,:),0,0,1);
plot_broken_axis(t_sig_d{ind}-2,fr_sig_d{ind},50,[20 100],[100 200],20,[6 57],cmap(ind,:),0,0,1);
hold off;
xlim([6 57]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
xlabel('Time (s)');
ylabel('Firing rate (Hz)');

%% Panel B
figure;set(gcf,'units','points','position',[345,489,406,265]);
hold on;
plot(linspace(0,size(edition.Pulsetrainclean{1},2)/2e3,size(edition.Pulsetrainclean{1},2)),edition.Pulsetrainclean{1}(ind,:))
hold off;
xlim([6 57]);
ylim([-0.5 1.5]);
set(gca,'YTickLabel',[]);
set(gca,'Visible','off');
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);

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
set(gca,'FontSize',18);

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
set(gca,'FontSize',18);

%% Panel C
figure;set(gcf,'units','points','position',[345,489,406,265]);
ind=MUind(1);
locs_new=interpolate_spike_train(edition.Distimeclean{1}{ind},edition.Pulsetrainclean{1}(ind,:));
locs_new=locs_new-round(find_lag*20);
force_sig=filtfilt(b,a,signal.ref_signal);
hold on;
plot(linspace(0,size(force_sig,2)/2e3,size(force_sig,2)),30*force_sig/max(force_sig),'Color',[150,150,150]/255,'LineWidth',1);
scatter(locs_new(1:end-1)/(40e3),1000./(diff(locs_new)/(40)),'MarkerFaceColor',cmap(ind,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75);
plot(t_sig_s{ind}-2,fr_sig_s{ind},'Color',cmap(ind,:),'LineWidth',3);
hold off;
ylim([0 30]);
xlim([6 57]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
xlabel('Time (s)');
ylabel('Firing rate (Hz)');

%% Panel D
figure;set(gcf,'units','points','position',[345,489,406,265]);
hold on;
plot(linspace(0,size(edition.Pulsetrainclean{1},2)/2e3,size(edition.Pulsetrainclean{1},2)),edition.Pulsetrainclean{1}(ind,:))
hold off;
xlim([6 57]);
ylim([-0.5 1.5]);
set(gca,'YTickLabel',[]);
set(gca,'Visible','off');
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);

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
set(gca,'FontSize',18);

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
set(gca,'FontSize',18);