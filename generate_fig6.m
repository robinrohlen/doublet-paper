addpath('functions');

%% Panel B

clearvars; close all;

load('data/sub-13/emg/sub-13_task-trapezoid20percentmvc_run-01_emg.mat')

subj_coding=1:24;

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
figure;set(gcf,'units','points','position',[340,468,460,282]);
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
set(gca,'FontSize',18);

g=gcf;
g.Renderer='painters';

%% Panel H
locs=edition.Distimeclean{1};
pulsetrain=edition.Pulsetrainclean{1};
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
    fr_sig_s{i}=conv(st_s,hann(2*fs),'same');
    fr_sig_s{i}=fr_sig_s{i}(round(fs*t_s(tmp_ind(1))):round(fs*t_s(tmp_ind(end))));
end

iter_s=0;
iter_d=0;

figure;set(gcf,'units','points','position',[345,538,291,216]);
hold on;
plot(linspace(0,size(force_sig,2)/2e3,size(force_sig,2)),30*force_sig/max(force_sig),'Color',[150,150,150]/255,'LineWidth',1);
hold off;

for i=1:size(t_sig_s,2)
    if doublet_ind(i)==0
        iter_s=iter_s+1;
        hold on;
        plot(t_sig_s{i}-2,fr_sig_s{i},'Color',cmap(i,:),'LineWidth',3);
        hold off;
        xlim([6 57]);
        ylim([0 30]);
        set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',18);
    end
end
ylabel('Firing rate (Hz)');

g=gcf;
g.Renderer='painters';

%% Panel D

MUind=[1 3 16 26];

figure;set(gcf,'units','points','position',[345,538,291,216]);
ind=MUind(3);
locs_new=interpolate_spike_train(edition.Distimeclean{1}{ind},edition.Pulsetrainclean{1}(ind,:));
locs_new=locs_new-round(find_lag*20);
force_sig=filtfilt(b,a,signal.ref_signal);
hold on;
plot(linspace(0,size(force_sig,2)/2e3,size(force_sig,2)),200./(100/(20+20))*force_sig/max(force_sig),'Color',[150,150,150]/255,'LineWidth',1);
plot_broken_axis(locs_new(1:end-1)/(40e3),1000./(diff(locs_new)/(40)),50,[20 100],[100 200],20,[6 57],cmap(ind,:),0,0);
plot_broken_axis(t_sig_s{ind}-2,fr_sig_s{ind},3,[20 100],[100 200],20,[6 57],cmap(ind,:),0,0,1);
plot_broken_axis(t_sig_d{ind}-2,fr_sig_d{ind},3,[20 100],[100 200],20,[6 57],cmap(ind,:),0,0,1);
hold off;
xlim([6 57]);
ylabel('Firing rate (Hz)')
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);

g=gcf;
g.Renderer='painters';

%% Panel E

figure;set(gcf,'units','points','position',[345,538,291,216]);
hold on;
plot(linspace(0,signal_len/fs,signal_len),80*force_sig./max(force_sig),'LineWidth',1,'Color',[150,150,150]/255);
hold off;

for i=1:size(t_sig_s,2)
    if doublet_ind(i)==1
        iter_d=iter_d+1;
        hold on;
        plot_broken_axis(t_sig_s{i}-2,fr_sig_s{i},3,[20 100],[100 200],20,[6 57],cmap(i,:),0,0,1);
        plot_broken_axis(t_sig_d{i}-2,fr_sig_d{i},3,[20 100],[100 200],20,[6 57],cmap(i,:),0,0,1);
        hold off;
        xlim([6 57]);
        set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',18);
    end
end
ylabel('Firing rate (Hz)');

g=gcf;
g.Renderer='painters';

%% Panel F

figure;set(gcf,'units','points','position',[345,538,291,216]);
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
ylabel('Firing rate (Hz)')
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);

g=gcf;
g.Renderer='painters';

%% Panel G

figure;set(gcf,'units','points','position',[345,538,291,216]);
ind=MUind(2);
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
ylim_tmp=get(gca,'ylim');
ylabel('Firing rate (Hz)')
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);

g=gcf;
g.Renderer='painters';

%% Panel I

clearvars; close all;

load_data=1;

if load_data==0

    num_doublets=5; % requirement: at least 5 doublets in a spike train
    options = optimset('Display','off');

    % Set params
    fs=2e3; % Sample rate
    win=-50:50; % +/- 50 samples (~ +/- 25 ms) for triggering
    t = 1e3*[0:length(win)-1]/fs; % time vector for figures
    doublet_crit=100; % Hz

    mvc_levels1={'05','10','15','20','25'};
    mvc_levels2={'02','05','10','20'};

    singlet_fr=[];
    doublet_fr=[];

    iter=0;
    for subj_nr=1:24
        disp(['Processing subject ',num2str(subj_nr),'...']);

        singlet_force_rt=[];
        doublet_force_rt=[];

        if subj_nr<=14 % experiment 1
            for file_num=1:size(mvc_levels1,2)
                if subj_nr<10
                    tmp=['data/sub-0',num2str(subj_nr),'/emg/sub-0',num2str(subj_nr),'_task-trapezoid',mvc_levels1{file_num},'percentmvc_run-01_emg.mat'];
                else
                    tmp=['data/sub-',num2str(subj_nr),'/emg/sub-',num2str(subj_nr),'_task-trapezoid',mvc_levels1{file_num},'percentmvc_run-01_emg.mat'];
                end
                load(tmp)

                locs=edition.Distimeclean{1};
                pulsetrain=edition.Pulsetrain{1};
                fs=signal.fsamp;

                signal_len = size(signal.data,2);

                clearvars signal edition parameters

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

                        if length(fr_d)>=num_doublets
                            iter=iter+1;
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
                            [B,fval,exitflag,output] = fminsearch(@(b) norm(fr_d_tmp - f(b,t_d_tmp)), [1000 -0.1 150],options);

                            if exitflag<=0
                                iter=iter-1;
                                continue;
                            end

                            % add interp with spline
                            t_sig_d{iter}=linspace(t_d_tmp(1),t_d_tmp(end),size(t_d_tmp,2)*10);
                            fr_sig_d{iter}=interp1(t_d_tmp,f(B,t_d_tmp),t_sig_d{iter},'spline');
                        end
                    end

                    st_s=zeros(1,signal_len);
                    st_s(round(locs_s))=1;
                end
            end
        else % experiment 2
            for file_num=1:size(mvc_levels2,2)
                tmp=['data/sub-',num2str(subj_nr),'/emg/sub-',num2str(subj_nr),'_task-trapezoid',mvc_levels2{file_num},'percentmvc_run-01_emg.mat'];
                load(tmp)

                locs=edition.Distimeclean{1};
                pulsetrain=edition.Pulsetrain{1};
                fs=signal.fsamp;

                signal_len = size(signal.data,2);

                clearvars signal edition parameters

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

                        if length(fr_d)>=num_doublets
                            iter=iter+1;
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
                            [B,fval,exitflag,output] = fminsearch(@(b) norm(fr_d_tmp - f(b,t_d_tmp)), [1000 -0.1 150],options);

                            if exitflag<=0
                                iter=iter-1;
                                continue;
                            end

                            % add interp with spline
                            t_sig_d{iter}=linspace(t_d_tmp(1),t_d_tmp(end),size(t_d_tmp,2)*10);
                            fr_sig_d{iter}=interp1(t_d_tmp,f(B,t_d_tmp),t_sig_d{iter},'spline');
                        end
                    end

                    st_s=zeros(1,signal_len);
                    st_s(round(locs_s))=1;
                end
            end
        end
    end
    save('data/figure-files/adaptation_curves.mat','t_sig_d','fr_sig_d')
else
    load('data/figure-files/adaptation_curves.mat');
end

cmap=turbo(26);

figure;set(gcf,'units','points','position',[345,538,291,216]);
hold on;
for i=1:size(fr_sig_d,2)
    if ~isempty(fr_sig_d{i})
        p1=plot(t_sig_d{i}-t_sig_d{i}(1),fr_sig_d{i}./fr_sig_d{i}(1),'Color',[150,150,150]/255);
    end
end

t_vec=5:0.1:25;
curve_vals=cell(1,length(t_vec));
mean_curve=zeros(1,length(t_vec));

for t=1:length(t_vec)
    curve_vals{t}=[];
    for i=1:size(fr_sig_d,2)
        if ~isempty(fr_sig_d{i})
            tmp=t_sig_d{i}-t_sig_d{i}(1);
            curve_vals{t}=[curve_vals{t} fr_sig_d{i}(find(tmp>=t_vec(t)-5 & tmp<t_vec(t)+5))./fr_sig_d{i}(1)];
        end
    end
    mean_curve(t)=mean(curve_vals{t});
end

hold on;
p3=plot(t_sig_d{53}-t_sig_d{53}(1),fr_sig_d{53}./fr_sig_d{53}(1),'Color',cmap(16,:),'LineWidth',6);
p2=plot([0 t_vec],[1 mean_curve],'k','LineWidth',6);
hold off;

xlim([0 40]);
ylim([0.5 1.0]);
yticks(0.5:0.1:1);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);

ylabel('Firing rate (n.u.)');

g=gcf;
g.Renderer='painters';

%% Panel J
clearvars; close all;

cmap=turbo(27);
MUind=[1 24];

load_data=0;

if load_data==0
    doublet_all=[];
    singlet_all=[];
    num_MUs_vec=[];
    prop_doublets=[];

    for subj_nr=1:24
        disp(['Processing subject ',num2str(subj_nr),'...']);

        singlet_force_rt=[];
        doublet_force_rt=[];

        if subj_nr<=14 % experiment 1
            if subj_nr<10
                tmp=['data/sub-0',num2str(subj_nr),'/emg/sub-0',num2str(subj_nr),'_task-trapezoid_emg_trapz_unique.mat'];
            else
                tmp=['data/sub-',num2str(subj_nr),'/emg/sub-',num2str(subj_nr),'_task-trapezoid_emg_trapz_unique.mat'];
            end
            load(tmp)
        else % experiment 2
            tmp=['data/sub-',num2str(subj_nr),'/emg/sub-',num2str(subj_nr),'_task-trapezoid_emg_trapz_unique.mat'];
            load(tmp)
        end

        doublet_all = [doublet_all doublet_force_rt];
        singlet_all = [singlet_all singlet_force_rt];
        num_MUs_vec = [num_MUs_vec size(muap_unique,2)];
        prop_doublets = cat(1,prop_doublets,[sum(doublet_unique==1)/size(muap_unique,2) 1-sum(doublet_unique==1)/size(muap_unique,2)]);
    end

    singlet_all=singlet_all(~isnan(singlet_all));
    singlet_all=singlet_all(isfinite(singlet_all));
    singlet_all=singlet_all(singlet_all<=50);

    save('data/figure-files/unique_mus.mat','doublet_all','singlet_all','prop_doublets','num_MUs_vec');
else
    load('data/figure-files/unique_mus.mat');
end

figure;set(gcf,'units','points','position',[345,538,291,216]);
hold on;
vp=violinplot([doublet_all,singlet_all],[1*ones(size(doublet_all)) 2*ones(size(singlet_all))],'BoxWidth',0.04,'ShowData',false);
plot([0.5 2.5],[1 1],'k:');
hold off;
xlim([0.5 2.5]);
ylim([0 5]);
yticks(0:1:5);
xticklabels({'Doublets','Singlets'})
ylabel('Relative force (% MVC)');
set(gcf,'color','w');
set(gca,'FontSize',12);
set(gca,'TickDir','out');

vp(1).ScatterPlot.MarkerFaceColor=cmap(MUind(1),:);
vp(2).ScatterPlot.MarkerFaceColor=cmap(MUind(2),:);
vp(1).ScatterPlot.MarkerFaceAlpha=0.5;
vp(2).ScatterPlot.MarkerFaceAlpha=0.5;
vp(1).ViolinPlot.EdgeAlpha=0;
vp(2).ViolinPlot.EdgeAlpha=0;
vp(1).MedianPlot.SizeData=100;
vp(2).MedianPlot.SizeData=100;
vp(1).MedianPlot.MarkerEdgeColor=[50 50 50]/255;
vp(2).MedianPlot.MarkerEdgeColor=[50 50 50]/255;
vp(1).ViolinPlot.FaceAlpha=0;
vp(2).ViolinPlot.FaceAlpha=0;

vp(1).BoxColor=cmap(MUind(1),:);
vp(2).BoxColor=cmap(MUind(2),:);
vp(1).BoxWidth=0.15;
vp(2).BoxWidth=0.15;
ax = gca;
box(ax,'off');

g=gcf;
g.Renderer='painters';

%% Panel K
clearvars; close all;

load_data=0;

if load_data==0
    doublet_all=[];
    singlet_all=[];
    num_MUs_vec=[];
    prop_doublets=[];

    for subj_nr=1:24
        disp(['Processing subject ',num2str(subj_nr),'...']);

        singlet_force_rt=[];
        doublet_force_rt=[];

        if subj_nr<=14 % experiment 1
            if subj_nr<10
                tmp=['data/sub-0',num2str(subj_nr),'/emg/sub-0',num2str(subj_nr),'_task-trapezoid_emg_trapz_unique.mat'];
            else
                tmp=['data/sub-',num2str(subj_nr),'/emg/sub-',num2str(subj_nr),'_task-trapezoid_emg_trapz_unique.mat'];
            end
            load(tmp)
        else % experiment 2
            tmp=['data/sub-',num2str(subj_nr),'/emg/sub-',num2str(subj_nr),'_task-trapezoid_emg_trapz_unique.mat'];
            load(tmp)
        end

        doublet_all = [doublet_all doublet_force_rt];
        singlet_all = [singlet_all singlet_force_rt];
        num_MUs_vec = [num_MUs_vec size(muap_unique,2)];
        prop_doublets = cat(1,prop_doublets,[sum(doublet_unique==1)/size(muap_unique,2) 1-sum(doublet_unique==1)/size(muap_unique,2)]);
    end

    singlet_all=singlet_all(~isnan(singlet_all));
    singlet_all=singlet_all(isfinite(singlet_all));
    singlet_all=singlet_all(singlet_all<=50);

    save('data/figure-files/unique_mus.mat','doublet_all','singlet_all','prop_doublets','num_MUs_vec');
else
    load('data/figure-files/unique_mus.mat');
end

cmap=turbo(27);
MUind=[1 24];

results = tost_bootstrap(doublet_all, 1, 0.05, 0.05, @median);

figure;set(gcf,'units','points','position',[345,538,291,216]);
hold on;
vp=violinplot([doublet_all,singlet_all],[1*ones(size(doublet_all)) 2*ones(size(singlet_all))],'BoxWidth',0.04);
plot([0.5 2.5],[1 1],'k:');
hold off;
xlim([0.5 2.5]);
ylim([0 10]);
xticklabels({'Doublets','Singlets'})
set(gcf,'color','w');
set(gca,'FontSize',18);
set(gca,'TickDir','out');

vp(1).ScatterPlot.MarkerFaceColor=cmap(MUind(1),:);
vp(2).ScatterPlot.MarkerFaceColor=cmap(MUind(2),:);
vp(1).ScatterPlot.MarkerFaceAlpha=0.5;
vp(2).ScatterPlot.MarkerFaceAlpha=0.5;
vp(1).ViolinPlot.EdgeAlpha=0;
vp(2).ViolinPlot.EdgeAlpha=0;
vp(1).MedianPlot.SizeData=100;
vp(2).MedianPlot.SizeData=100;
vp(1).MedianPlot.MarkerEdgeColor=[50 50 50]/255;
vp(2).MedianPlot.MarkerEdgeColor=[50 50 50]/255;
vp(1).ViolinPlot.FaceAlpha=0;
vp(2).ViolinPlot.FaceAlpha=0;

vp(1).BoxColor=cmap(MUind(1),:);
vp(2).BoxColor=cmap(MUind(2),:);
vp(1).BoxWidth=0.15;
vp(2).BoxWidth=0.15;
vp(1).BoxPlot.EdgeColor=[1 1 1];

ax = gca;
box(ax,'off');

ylim([0.9 1.1]);
xlim([0.5 1.5]);

hold on;
patch([1.35 1.4 1.4 1.35],[0.95 0.95 1.05 1.05],cmap(MUind(1),:),'EdgeColor','none','FaceAlpha',0.5);
plot([1 1.4],[0.95 0.95],'--','Color','r');
plot([1 1.4],[1.05 1.05],'--','Color','r');
plot([1.375 1.375],[1.028333, 1.039538],'r-','LineWidth',4);
plot(1.375,1.032717,'o','MarkerFaceColor','w','MarkerEdgeColor','r','MarkerSize',6);
hold off;

g=gcf;
g.Renderer='painters';

%% Panel C

clearvars; close all;

load_data=0;

if load_data==0
    doublet_all=[];
    singlet_all=[];
    num_MUs_vec=[];
    prop_doublets=[];

    for subj_nr=1:24
        disp(['Processing subject ',num2str(subj_nr),'...']);

        singlet_force_rt=[];
        doublet_force_rt=[];

        if subj_nr<=14 % experiment 1
            if subj_nr<10
                tmp=['data/sub-0',num2str(subj_nr),'/emg/sub-0',num2str(subj_nr),'_task-trapezoid_emg_trapz_unique.mat'];
            else
                tmp=['data/sub-',num2str(subj_nr),'/emg/sub-',num2str(subj_nr),'_task-trapezoid_emg_trapz_unique.mat'];
            end
            load(tmp)
        else % experiment 2
            tmp=['data/sub-',num2str(subj_nr),'/emg/sub-',num2str(subj_nr),'_task-trapezoid_emg_trapz_unique.mat'];
            load(tmp)
        end

        doublet_all = [doublet_all doublet_force_rt];
        singlet_all = [singlet_all singlet_force_rt];
        num_MUs_vec = [num_MUs_vec size(muap_unique,2)];
        prop_doublets = cat(1,prop_doublets,[sum(doublet_unique==1)/size(muap_unique,2) 1-sum(doublet_unique==1)/size(muap_unique,2)]);
    end

    singlet_all=singlet_all(~isnan(singlet_all));
    singlet_all=singlet_all(isfinite(singlet_all));
    singlet_all=singlet_all(singlet_all<=50);

    save('data/figure-files/unique_mus.mat','doublet_all','singlet_all','prop_doublets','num_MUs_vec');
else
    load('data/figure-files/unique_mus.mat');
end

cmap=turbo(27);
MUind=[1 24];

figure;set(gcf,'units','points','position',[472,285,399,475]);%[345,489,406,265]);

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
set(gca,'XTickLabel',[]);
set(gcf,'color','w');
set(gca,'FontSize',18);
xlim([0.5 24.5]);
ylim([0 100]);
yticks(0:25:100);
ax = gca;
ax.YColor = cmap(2,:);

subj_nr_vec=1:24;

col=[247,165,30]/255;

yyaxis right
hold on;
plot(subj_nr_vec,num_MUs_vec,'-o','Color',col,'MarkerFaceColor',col,'MarkerSize',12,'MarkerEdgeColor',col,'LineWidth',1.5);
hold off;

ax = gca;
ax.YColor = col;
ylim([0 80]);
yticks(0:20:80);

g=gcf;
g.Renderer='painters';

%% Panel L
clearvars; close all;

load_data=1;

if load_data==0
    % Set params
    fs=2e3; % Sample rate
    win=-50:50; % +/- 50 samples (~ +/- 25 ms) for triggering
    t = 1e3*[0:length(win)-1]/fs; % time vector for figures
    doublet_crit=100; % Hz

    mvc_levels1={'05','10','15','20','25'};
    mvc_levels2={'02','05','10','20'};

    singlet_fr=[];
    doublet_fr=[];

    for subj_nr=1:24
        disp(['Processing subject ',num2str(subj_nr),'...']);

        singlet_force_rt=[];
        doublet_force_rt=[];

        if subj_nr<=14
            for file_num=1:size(mvc_levels1,2)
                if subj_nr<10
                    tmp=['data/sub-0',num2str(subj_nr),'/emg/sub-0',num2str(subj_nr),'_task-trapezoid',mvc_levels1{file_num},'percentmvc_run-01_emg.mat'];
                else
                    tmp=['data/sub-',num2str(subj_nr),'/emg/sub-',num2str(subj_nr),'_task-trapezoid',mvc_levels1{file_num},'percentmvc_run-01_emg.mat'];
                end
                load(tmp)

                iter=0;
                % Loop through each MU to extract its MUAP
                for MUnum=1:size(edition.Distimeclean{1},2)
                    iter=iter+1;

                    % Extract the discharges used for triggering
                    locs=edition.Distimeclean{1}{MUnum};

                    % Find doublets
                    fr=1000./(diff(locs)/2);
                    doublet_ind=find(fr >= doublet_crit);

                    if length(doublet_ind)>=2
                        doublet(MUnum)=1;

                        if doublet_ind(end)==size(fr,2)
                            doublet_ind(end)=[];
                        end
                        doublet_fr=[doublet_fr mean(fr(doublet_ind+1))];

                        % Remove doublets and post-doublet ISIs
                        locs(sort([doublet_ind doublet_ind+1]))=[];
                        fr=1000./(diff(locs)/2);
                        if ~isempty(locs)
                            singlet_fr=[singlet_fr mean(fr)];
                        else
                            singlet_fr=[singlet_fr 0];
                        end
                    end
                end
            end
        else
            for file_num=1:size(mvc_levels2,2)
                tmp=['data/sub-',num2str(subj_nr),'/emg/sub-',num2str(subj_nr),'_task-trapezoid',mvc_levels2{file_num},'percentmvc_run-01_emg.mat'];
                load(tmp)

                iter=0;
                % Loop through each MU to extract its MUAP
                for MUnum=1:size(edition.Distimeclean{1},2)
                    iter=iter+1;

                    % Extract the discharges used for triggering
                    locs=edition.Distimeclean{1}{MUnum};

                    % Find doublets
                    fr=1000./(diff(locs)/2);
                    doublet_ind=find(fr >= doublet_crit);

                    if length(doublet_ind)>=2
                        doublet(MUnum)=1;

                        if doublet_ind(end)==size(fr,2)
                            doublet_ind(end)=[];
                        end
                        doublet_fr=[doublet_fr mean(fr(doublet_ind+1))];

                        % Remove doublets and post-doublet ISIs
                        locs(sort([doublet_ind doublet_ind+1]))=[];
                        fr=1000./(diff(locs)/2);
                        if ~isempty(locs)
                            singlet_fr=[singlet_fr mean(fr)];
                        else
                            singlet_fr=[singlet_fr 0];
                        end
                    end
                end
            end
        end
    end
    save('data/figure-files/post_firingrates.mat','singlet_fr','doublet_fr')
else
    load('data/figure-files/post_firingrates.mat');
end

fr=[singlet_fr; doublet_fr];
fr=fr(:,find(fr(1,:)<20));
fr=fr(:,find(fr(2,:)<20));
mean_fr=median(fr,2)';

figure;set(gcf,'units','points','position',[345,538,291,216]);

hold on;
p1=plot(1:2,flip(fr),'-o','Color',[150,150,150]/255,'MarkerFaceColor',[150,150,150]/255,'MarkerSize',6);
p2=plot(1:2,flip(mean_fr),'-o','Color',[0,0,0]/255,'MarkerFaceColor',[0,0,0]/255,'LineWidth',6,'MarkerSize',6);
hold off;
xlim([0.5 2.5]);
ylim([0 20]);
xticks(1:2);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
xticklabels({'doublet','singlet'})
ylabel('Firing rate (Hz)');

g=gcf;
g.Renderer='painters';

fr=[singlet_fr; doublet_fr];
[p,h,stats]=signrank(diff(fr));

disp(['The paired median difference of singlet and doublet firing rates was significantly greater than 0 (Wilcoxon signed-rank test: W = ',num2str(stats.signedrank),', p < ',num2str(round(p,3)),').']);
