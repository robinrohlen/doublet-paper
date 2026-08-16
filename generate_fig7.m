addpath('functions')

%% Panel B

clearvars; close all;

% Choose between MU #1 or #11
% MUind=1;
MUind=11;

load('data/figure-files/twitch_population_locs.mat')

rec_ind=1:2:8;

subj=10; % 10 th subject in experiment 2 (#24 in total)
rec=4;

load('data/sub-24/emg/sub-24_task-trapezoid20percentmvc_run-01_emg.mat');
load('data/sub-24/us/sub-24_task-trapezoid20percentmvc_run-01_sync_locs.mat');
load('data/sub-24/us/sub-24_task-trapezoid20percentmvc_run-01_time.mat');
load('data/sub-24/us/sub-24_task-trapezoid20percentmvc_run-01_timevelocity.mat');
load('data/sub-24/us/sub-24_task-trapezoid20percentmvc_run-01_velocity.mat');

v_2d_all=v_2d_all(1:39,:,:);

locs=edition.Distimeclean{1}{MUind};

us_trig_tmp(1)=us_trig_locs(rec_ind(4))-emg_trig_locs(rec_ind(4))+1;
us_trig_tmp(2)=us_trig_locs(rec_ind(4)+1)-emg_trig_locs(rec_ind(4))+1;

locs=locs(find(locs>=us_trig_tmp(1) & locs<=us_trig_tmp(2))) - us_trig_tmp(1);
locs=round(1000*locs/2042.5);

fr=1000./diff(locs);
locs_d=save_doublet_locs{subj}{rec}{MUind};
locs_s=save_singlet_locs{subj}{rec}{MUind};

win_us=-100:100;
FR=1e03;
numCoords=10;

% STA
for i=1:2
    if i==1 % doublet
        locs_tmp=locs_d;
    elseif i==2 % singlet
        locs_tmp=locs_s;
    end

    TMP_MAT = zeros(size(v_2d_all,1),size(v_2d_all,2),length(win_us),length(locs_tmp));

    for tmpind = 1:length(locs_tmp)
        tmptmp = squeeze(v_2d_all(:,:,locs_tmp(tmpind)+win_us));
        TMP_MAT = cat(4,TMP_MAT,tmptmp);
    end

    MEAN_TMP_MAT=mean(TMP_MAT,4);
    STD_TMP_MAT=std(TMP_MAT,[],4);
    clearvars TMP_MAT

    direction=sum(MEAN_TMP_MAT(:,:,1:floor(length(win_us)/2))./STD_TMP_MAT(:,:,1:floor(length(win_us)/2)),3)-sum(MEAN_TMP_MAT(:,:,ceil(length(win_us)/2):length(win_us)-1)./STD_TMP_MAT(:,:,ceil(length(win_us)/2):length(win_us)-1),3);
    MU_intensity=-sign(direction).*sum((MEAN_TMP_MAT.^2)./STD_TMP_MAT,3);

    yc_sta_pos=zeros(numCoords,1);
    xc_sta_pos=zeros(numCoords,1);
    yc_sta_neg=zeros(numCoords,1);
    xc_sta_neg=zeros(numCoords,1);
    maxValsPos=maxk(MU_intensity(:),numCoords);
    maxValsNeg=mink(MU_intensity(:),numCoords);
    for indtmp=1:numCoords
        [yc_sta_pos(indtmp),xc_sta_pos(indtmp)]=find(MU_intensity==maxValsPos(indtmp));
        [yc_sta_neg(indtmp),xc_sta_neg(indtmp)]=find(MU_intensity==maxValsNeg(indtmp));
    end

    medianSTAtwitchPos=zeros(numCoords,length(win_us));
    medianSTAtwitchNeg=zeros(numCoords,length(win_us));
    for ind=1:numCoords,medianSTAtwitchPos(ind,:)=MEAN_TMP_MAT(yc_sta_pos(ind),xc_sta_pos(ind),:);end
    for ind=1:numCoords,medianSTAtwitchNeg(ind,:)=MEAN_TMP_MAT(yc_sta_neg(ind),xc_sta_neg(ind),:);end
    medianSTAtwitchPos=median(medianSTAtwitchPos);
    medianSTAtwitchNeg=median(medianSTAtwitchNeg);

    if i==1 % doublet
        twitch_d=medianSTAtwitchPos;
        sta_map_d=MU_intensity;
    elseif i==2 % singlet
        twitch_s=medianSTAtwitchPos;
        sta_map_s=MU_intensity;
    end
end

x_rect=38.4.*[min(xc_sta_pos) max(xc_sta_pos)]./128;y_rect=20.*[min(yc_sta_pos) max(yc_sta_pos)]./39;

unf_tet_sig=zeros(numCoords,size(v_2d_all,3));
for i=1:numCoords
    unf_tet_sig(i,:)=squeeze(v_2d_all(yc_sta_pos(i),xc_sta_pos(i),:));
end

% muap est
win_emg=-50:50;

sig=notchsignals(signal.data(1:64,:),signal.fsamp);
sig=bandpassingals(sig,signal.fsamp);
% Loop through each MU to extract its muap_d

% Extract the discharges used for triggering
locs_d1=edition.Distimeclean{1}{MUind};
locs_s1=edition.Distimeclean{1}{MUind};

% Remove doublets in averaging
doublet_ind=find(1000./(diff(locs_s1)/2) >= 100);
if ~isempty(doublet_ind)
    locs_s1(sort([doublet_ind doublet_ind+1]))=[];
end

no_doublet_ind=find(1000./(diff(locs_d1)/2) < 100);
if ~isempty(no_doublet_ind)
    locs_d1(sort(no_doublet_ind))=[];
end

% Pre-define muap_d cell
muap_d=zeros(64,length(win_emg),length(locs_d1));
muap_nd=zeros(64,length(win_emg),length(locs_s1));

% Extract the muap_d through STA for each channel
for ch=1:size(muap_d,1)
    iterTrig=0;
    for trig=1:length(locs_d1)
        muap_d(ch,:,trig)=sig(ch,locs_d1(trig)+win_emg);
    end
end

for ch=1:size(muap_nd,1)
    iterTrig=0;
    for trig=1:length(locs_s1)
        muap_nd(ch,:,trig)=sig(ch,locs_s1(trig)+win_emg);
    end
end

tmp_nd=muap_grid(mean(muap_nd,3));
tmp_nd_diff=diff(tmp_nd,1);
tmp_d=muap_grid(mean(muap_d,3));
tmp_d_diff=diff(tmp_d,1);

muap_nd=squeeze(tmp_nd_diff(3,3,:));
muap_d=squeeze(tmp_d_diff(3,3,:));

muap_nd_interp=-interp1(linspace(-50,50,length(win_emg)),muap_nd,linspace(-50,50,5*length(win_emg)),'spline');
muap_d_interp=-interp1(linspace(-50,50,length(win_emg)),muap_d,linspace(-50,50,5*length(win_emg)),'spline');

cmap_map = crameri('lajolla');
cmap_map = flip(cmap_map);

cmap_lines=[92 123 207;222 166 90;220 110 75;68 150 136;152 102 171;196 120 138]/255;

figure;set(gcf,'units','points','position',[572,524,178,137]);
imagesc([0 38.4],[0 20],sta_map_d)
cb1=colorbar;
cb1.FontSize=12;
cb1.Ticks=[-3 7];
colormap(cmap_map);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',12);

figure;set(gcf,'units','points','position',[572,524,178,137]);
imagesc([0 38.4],[0 20],sta_map_s)
cb1=colorbar;
cb1.FontSize=12;
cb1.Ticks=[-1 1];
colormap(cmap_map);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',12);

%% Panel C

unf_tet_sig=mean(unf_tet_sig,1);

figure;set(gcf,'units','points','position',[622,397,363,271]);

us_trig_tmp(1)=us_trig_locs(rec_ind(rec))-emg_trig_locs(rec_ind(rec))+1;
us_trig_tmp(2)=us_trig_locs(rec_ind(rec)+1)-emg_trig_locs(rec_ind(rec))+1;

locs1=interpolate_spike_train(edition.Distimeclean{1}{11},edition.Pulsetrainclean{1}(11,:));%edition.Distimeclean{1}{11};
locs1=locs1/20;
locs1=locs1(find(locs1>=us_trig_tmp(1) & locs1<=us_trig_tmp(2))) - us_trig_tmp(1);
locs1=1000*locs1/2042.5;

locs2=interpolate_spike_train(edition.Distimeclean{1}{1},edition.Pulsetrainclean{1}(1,:));%edition.Distimeclean{1}{1};
locs2=locs2/20;
locs2=locs2(find(locs2>=us_trig_tmp(1) & locs2<=us_trig_tmp(2))) - us_trig_tmp(1);
locs2=1000*locs2/2042.5;

plot_broken_axis(locs1(1:end-1)/1e3,1000./(diff(locs1)),100,[20 100],[100 200],20,[0 size(signal.data,2)/2042.5],cmap_lines(3,:),0,0);
plot_broken_axis(locs2(1:end-1)/1e3,1000./(diff(locs2)),100,[20 100],[100 200],20,[0 size(signal.data,2)/2042.5],cmap_lines(4,:),0,0);

save_ylim=ylim;
hold on;
h=plot(linspace(0,30,length(unf_tet_sig)),0.85*diff(save_ylim)*(unf_tet_sig-min(unf_tet_sig))/max(unf_tet_sig-min(unf_tet_sig)),'Color',[120,120,120]/255,'LineWidth',0.5);
p1=plot(-100,-100,'o','MarkerFaceColor',cmap_lines(3,:),'MarkerEdgeColor','k','MarkerSize',10);%,'MarkerEdgeColor','k','MarkerFaceAlpha',1.0,'MarkerEdgeAlpha',1.0);
p2=plot(-100,-100,'o','MarkerFaceColor',cmap_lines(4,:),'MarkerEdgeColor','k','MarkerSize',10);%,'MarkerEdgeColor','k','MarkerFaceAlpha',1.0,'MarkerEdgeAlpha',1.0);
p3=plot(-100,-100,'-','Color',[120,120,120]/255,'LineWidth',4);
hold off;
xlim([0 20]);
uistack(h,'bottom');
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',12);

g=gcf;
g.Renderer='painters';

%% Panel D

figure;set(gcf,'units','points','position',[622,434,417,234]);
hold on;
plot(linspace(0,30,length(unf_tet_sig)),1.5*diff(save_ylim)*(unf_tet_sig-min(unf_tet_sig))/max(unf_tet_sig-min(unf_tet_sig)),'Color',[120,120,120]/255,'LineWidth',3);
plot(locs1(1:end-1)/1e3,1000./(diff(locs1)),':o','Color',cmap_lines(3,:),'LineWidth',1,'MarkerSize',10);
scatter(locs1(1:end-1)/1e3,1000./(diff(locs1)),100,'MarkerFaceColor',cmap_lines(3,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75);
plot(locs2(1:end-1)/1e3,1000./(diff(locs2)),':o','Color',cmap_lines(4,:),'LineWidth',1,'MarkerSize',10);
scatter(locs2(1:end-1)/1e3,1000./(diff(locs2)),100,'MarkerFaceColor',cmap_lines(4,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75);
hold off;
xlim([7.6819 8.3819]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',22);
g=gcf;
g.Renderer='painters';

%% Panel E
figure;set(gcf,'units','points','position',[622,434,417,234]);
hold on;
plot(linspace(0,30,length(unf_tet_sig)),1.5*diff(save_ylim)*(unf_tet_sig-min(unf_tet_sig))/max(unf_tet_sig-min(unf_tet_sig)),'Color',[120,120,120]/255,'LineWidth',3);
plot(locs1(1:end-1)/1e3,1000./(diff(locs1)),':o','Color',cmap_lines(3,:),'LineWidth',1,'MarkerSize',10);
scatter(locs1(1:end-1)/1e3,1000./(diff(locs1)),100,'MarkerFaceColor',cmap_lines(3,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75);
plot(locs2(1:end-1)/1e3,1000./(diff(locs2)),':o','Color',cmap_lines(4,:),'LineWidth',1,'MarkerSize',10);
scatter(locs2(1:end-1)/1e3,1000./(diff(locs2)),100,'MarkerFaceColor',cmap_lines(4,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75);
hold off;
xlim([15.4964 15.4964+0.7]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',22);
g=gcf;
g.Renderer='painters';

%% Panel F

figure;set(gcf,'units','points','position',[622,460,166,208]);
hold on;
p2=plot(linspace(-50/2,50/2,5*length(win_emg)),-muap_d_interp,'LineWidth',3,'Color',cmap_lines(1,:));
p1=plot(linspace(-50/2,50/2,5*length(win_emg)),-(muap_nd_interp+0.025),'LineWidth',2,'Color',cmap_lines(2,:));
hold off;
axis tight;
ylim([-530 470]);
if MUind==1
    xlim([-12.5 12.5]);
elseif MUind==11
    xlim([-17.5 7.5]);
end
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',20);
g=gcf;
g.Renderer='painters';

%% Panel G

figure;set(gcf,'units','points','position',[622,460,166,208]);
hold on;
p1=plot(win_us(76:end),twitch_d(76:end)-min(twitch_d(76:end)),'LineWidth',3,'Color',cmap_lines(1,:));
p2=plot(win_us(76:end),twitch_s(76:end)-min(twitch_s(76:end)),'LineWidth',3,'Color',cmap_lines(2,:));
p3=plot(win_us(76:end),2*twitch_s(76:end)-min(2*twitch_s(76:end)),':','LineWidth',3,'Color',cmap_lines(2,:));
hold off;
axis tight;
ylim([-0.2 0.4])
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',20);
g=gcf;
g.Renderer='painters';

%% Panel H

rng(0);

fs=1e3;
T=10; % s
isi_cov=0.10;

num_iter=100;
gain_range=0.5:0.01:1.5;

mean_force_norm=zeros(length(gain_range),num_iter);
mean_force=zeros(length(gain_range),num_iter);

for iter=1:num_iter
    % Randomly generate for each iteration
    fr_singlet=11.2+(13.6-11.2).*rand(1,1);
    fr_doublet=(0.56+(0.75-0.56).*rand(1,1))*fr_singlet;

    Tc=75+(90-75).*rand(1,1);
    Thr=(5/3)*Tc;
    [t,F]=RaikovaForceTwitch5p(1e3,0,0,Tc,Thr,500,1);

    for i=1:length(gain_range)

        gain=gain_range(i);

        num_firings=round(T*fr_singlet);
        locs_singlet=round(cumsum(1000/fr_singlet+isi_cov*(1000/fr_singlet)*randn(1,num_firings)));
        st=zeros(1,T*fs);
        st(locs_singlet)=1;
        force_singlet=conv(st,F);

        rng(iter);
        num_firings=round(T*fr_doublet);
        intra_doublet_isi=5*ones(1,num_firings);
        locs_tmp=1000/fr_doublet+isi_cov*(1000/fr_doublet)*randn(1,num_firings);

        tmp=0;
        for ind=1:num_firings
            tmp=tmp+1;
            locs_doublet(tmp)=locs_tmp(ind);

            tmp=tmp+1;
            locs_doublet(tmp)=intra_doublet_isi(ind);
        end

        locs_doublet=cumsum(round(locs_doublet));
        st=zeros(1,T*fs);
        st(locs_doublet)=1;
        force_doublet=gain*conv(st,F);

        mean_force_norm(i,iter)=mean(force_doublet)./length(locs_doublet)-mean(force_singlet)./length(locs_singlet);
        mean_force(i,iter)=mean(force_doublet)-mean(force_singlet);

        clearvars locs_doublet
    end
end

% generate figure
cmap=[92 123 207;220 110 75;68 150 136;222 166 90;152 102 171;196 120 138]/255;

nIter = size(mean_force,2);
mf_mean = mean(mean_force,2);
mf_sem  = std(mean_force,0,2) ./ sqrt(nIter);
mfn_mean = mean(mean_force_norm,2);
mfn_sem  = std(mean_force_norm,0,2) ./ sqrt(nIter);

figure; set(gcf,'units','points','position',[573,461,156,199]);

hold on;

% SEM shaded area
fill([gain_range fliplr(gain_range)], ...
     [mf_mean-mf_sem; flipud(mf_mean+mf_sem)], ...
     [150,150,150]/255, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

% Mean line
plot(gain_range, mf_mean, 'Color', cmap(1,:), 'LineWidth', 2);

plot(xlim,[0 0],'k:');
plot([0.79 0.79],[0 0],'o','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor','k','MarkerSize',10);
ylim([-1 1]);
yticks([-1 0 1]);
xticks(0.5:0.5:1.5);
set(gca,'TickDir','out','FontSize',12);
hold off;

hold on;

mfn_mean_plot = 1e3 * mfn_mean / 5;
mfn_sem_plot  = 1e3 * mfn_sem / 5;

% SEM shaded area
fill([gain_range fliplr(gain_range)], ...
     [mfn_mean_plot-mfn_sem_plot; flipud(mfn_mean_plot+mfn_sem_plot)], ...
     [150,150,150]/255, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

% Mean line
plot(gain_range, mfn_mean_plot, 'Color', cmap(2,:), 'LineWidth', 2);

plot([1.035 1.035],[0 0],'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k','MarkerSize',10);
ylim([-1 1]);
yticks([-1 0 1]);
xticks(0.5:0.5:1.5);
set(gca,'TickDir','out','FontSize',12);
hold off;
xlabel('Gain');
ylabel('Delta mean force (n.u.)');

set(gcf,'color','w');

g=gcf;
g.Renderer='painters';

%% Panel I

clearvars; close all;

rng(0);

fs=1e3;
T=10; % s
isi_cov=0.10;

num_iter=100;
gain_range=1.035;

cov_force_singlet=zeros(length(gain_range),num_iter);
cov_force_doublet=zeros(length(gain_range),num_iter);

for iter=1:num_iter
    % Randomly generate for each iteration
    fr_singlet=11.2+(13.6-11.2).*rand(1,1);
    fr_doublet=(0.56+(0.75-0.56).*rand(1,1))*fr_singlet;

    Tc=75+(90-75).*rand(1,1);
    Thr=(5/3)*Tc;
    [t,F]=RaikovaForceTwitch5p(1e3,0,0,Tc,Thr,500,1);

    for i=1:length(gain_range)

        gain=gain_range(i);

        num_firings=round(T*fr_singlet);
        locs_singlet=round(cumsum(1000/fr_singlet+isi_cov*(1000/fr_singlet)*randn(1,num_firings)));
        st=zeros(1,T*fs);
        st(locs_singlet)=1;
        force_singlet=conv(st,F);

        rng(iter);
        num_firings=round(T*fr_doublet);
        intra_doublet_isi=5*ones(1,num_firings);
        locs_tmp=1000/fr_doublet+isi_cov*(1000/fr_doublet)*randn(1,num_firings);

        tmp=0;
        for ind=1:num_firings
            tmp=tmp+1;
            locs_doublet(tmp)=locs_tmp(ind);

            tmp=tmp+1;
            locs_doublet(tmp)=intra_doublet_isi(ind);
        end

        locs_doublet=cumsum(round(locs_doublet));
        st=zeros(1,T*fs);
        st(locs_doublet)=1;
        force_doublet=gain*conv(st,F);

        cov_force_singlet(i,iter)=std(force_singlet(1001:9000))./mean(force_singlet(1001:9000));
        cov_force_doublet(i,iter)=std(force_doublet(1001:9000))./mean(force_doublet(1001:9000));

        clearvars locs_doublet
    end
end

cmap=[92 123 207;220 110 75;68 150 136;222 166 90;152 102 171;196 120 138]/255;

figure;set(gcf,'units','points','position',[564,524,156,218]);
hold on;
vp=violinplot(100.*[cov_force_singlet;cov_force_doublet]',[1*ones(length(cov_force_singlet)); 2*ones(length(cov_force_doublet))]','BoxWidth',0.04);
hold off;
xlim([0.5 2.5])
ylim([0 50]);
yticks(0:25:50);
xticklabels({'Singlets','Doublets'})
ylabel('Force CoV (%)');
set(gcf,'color','w');
set(gca,'FontSize',12);
set(gca,'TickDir','out');

vp(1).ScatterPlot.MarkerFaceColor=cmap(1,:);
vp(2).ScatterPlot.MarkerFaceColor=cmap(2,:);
vp(1).ScatterPlot.MarkerFaceAlpha=0.5;
vp(2).ScatterPlot.MarkerFaceAlpha=0.5;
vp(1).MedianPlot.SizeData=100;
vp(2).MedianPlot.SizeData=100;
vp(1).MedianPlot.MarkerEdgeColor=[50 50 50]/255;
vp(2).MedianPlot.MarkerEdgeColor=[50 50 50]/255;

ax = gca;
box(ax,'off');

g=gcf;
g.Renderer='painters';
