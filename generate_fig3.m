addpath('simulation model');
addpath('functions');

%% Panel B
dt=0.001;

cmap=[92 123 207;220 110 75;68 150 136;222 166 90;152 102 171;196 120 138]/255;

param = [0.03 0.30 -90 0.20 0.04 0.05 0.02 0.02 0.07 0.08 0.15 150 25];
[V0,~,~,t,I_output0] = doublet_3comp_model([0.00 1.70], 'soma', 'euler', param);
[V1,~,~,t,I_output1] = doublet_3comp_model([0.00 2.70], 'soma', 'euler', param);

figure;set(gcf,'units','points','position',[488,585,392,175]);
hold on;
plot(t,V0,'Color',cmap(1,:),'LineWidth',2);
hold off;
axis tight;
xlim([975 1100]);
ylim([-80 30]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

figure;set(gcf,'units','points','position',[487,585,392,108]);
hold on;
plot(t,I_output0.I_stim,'Color','k','LineWidth',2);
hold off;
axis tight;
xlim([975 1100]);
ylim([0.0 6.0]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

figure;set(gcf,'units','points','position',[488,585,392,175]);
hold on;
plot(t,V1,'Color',cmap(1,:),'LineWidth',2);
hold off;
axis tight;
xlim([975 1100]);
ylim([-80 30]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

figure;set(gcf,'units','points','position',[487,585,392,108]);
hold on;
plot(t,I_output1.I_stim,'Color','k','LineWidth',2);
hold off;
axis tight;
xlim([975 1100]);
ylim([0.0 6.0]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

%% Panel C
[~,locs1]=findpeaks(V0,'MinPeakHeight',-10,'MinPeakDistance',2/dt);
fr1=1000./(dt*diff(locs1));
[~,locs2]=findpeaks(V1,'MinPeakHeight',-10,'MinPeakDistance',2/dt);
fr2=1000./(dt*diff(locs2));

figure;set(gcf,'units','points','position',[591,490,249,218]);
hold on;
plot_broken_axis(locs2(1:end-1)/(1/dt),fr2,75,[20 100],[100 150],20,[900 2000],cmap(1,:),1,0);
hold off;
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

%% Panel D

clearvars;close all;

cmap=[92 123 207;220 110 75;68 150 136;222 166 90;152 102 171;196 120 138]/255;

load_data=1;

if load_data==0
    cmap=[92 123 207;220 110 75;68 150 136;222 166 90;152 102 171;196 120 138]/255;

    param = [0.03 0.30 -90 0.20 0.04 0.05 0.02 0.02 0.07 0.08 0.15 150 25];

    I=[0.20 1.35]; % 0.2 1.4 works
    ramp_length=[flip(diff(I)./linspace(0.00015,0.3,30))]; % 10 works
    dV_th_peak=0.010; % 10 mV/ms

    dt = 0.001;

    V_thr_first=zeros(1,length(ramp_length));
    V_thr_second=zeros(1,length(ramp_length));
    for i=1:length(ramp_length)
        disp([num2str(i),'/',num2str(length(ramp_length))])

        T=2500+1000+ramp_length(i);
        stim_length=T-2500;
        Nt = round(T/dt);

        I_stim = I(1)*ones(Nt,1);
        I_stim(round(2500/dt):round((2500+stim_length)/dt)) = I(2);
        I_stim((round(2500/dt)-1):round((2500+ramp_length(i))/dt)) = linspace(I(1),I(2),round(ramp_length(i)/dt)+2);

        % First spike
        [Vs,Vd,Va,t,I_output] = doublet_3comp_model(I_stim, 'soma', 'euler', param);

        dVs=gradient(Vs);

        [~,locs]=findpeaks(Vs,'MinPeakHeight',0);
        fr=1000./(dt*diff(locs));
        fr_idx=find(fr>=100);
        locs(fr_idx+1)=[];
        dVs_roi=dVs((1000/dt):locs(1));
        %dVs_roi=(dVs((1000/dt):locs(1))-dV_th_peak).^2;
        %[maxVal,maxInd]=max(dVs_roi);
        
        abs_dVs_roi=abs((dVs_roi-dV_th_peak));
        [maxVal,maxInd]=max(abs_dVs_roi);
        abs_dVs_roi=abs_dVs_roi((maxInd-2/dt):maxInd);
        [minVal,minInd]=min(abs_dVs_roi);

        V_thr_first(i)=Vs((1000/dt)+(maxInd-2/dt)+minInd);
        % V_thr_first(i)=Vs((1003/dt)+minInd);
       
        % Second spike
        dVs_roi=dVs(locs(1):locs(2));
        abs_dVs_roi=abs((dVs_roi-dV_th_peak));
        [maxVal,maxInd]=max(abs_dVs_roi);
        abs_dVs_roi=abs_dVs_roi((maxInd-2/dt):maxInd);
        [minVal,minInd]=min(abs_dVs_roi);
        
        V_thr_second(i)=Vs(locs(1)+(maxInd-2/dt)+minInd);
        Vs_save{i}=Vs;
    end
    save('data/figure-files/ramp_data_thr.mat','I','ramp_length','V_thr_first','V_thr_second');
else
    load('data/figure-files/ramp_data_thr.mat');
end

ramp_speed=diff(I)./ramp_length;

figure;set(gcf,'units','points','position',[549,402,317,278]);
yyaxis left;
hold on;
plot(ramp_speed(1:1:end),V_thr_first(1:1:end),'o','MarkerFaceColor',cmap(5,:),'Color',cmap(5,:),'MarkerSize',8);
hold off;
% ylim([-56 -53]);
% yticks([-56 -55 -54 -53]);
ylim([-58 -56]);
yticks([-58 -57 -56]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);

yyaxis right;
hold on;
plot(ramp_speed(1:1:end),V_thr_second(1:1:end),'o','MarkerFaceColor',cmap(4,:),'Color',cmap(4,:),'MarkerSize',8);
hold off;
% ylim([-56 -53]);
% yticks([-56 -55 -54 -53]);
ylim([-58 -56]);
yticks([-58 -57 -56]);
xlim([0 0.3]);
xticks(0.0:0.1:0.3);
hold on;
plot(xlim,[1.0 1.0],':','Color',cmap(4,:));
hold off;

set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
ax = gca;
ax.YAxis(1).Color = cmap(5,:);
ax.YAxis(2).Color = cmap(4,:);
g=gcf;
g.Renderer='painters';

%% Panel E
clearvars;close all;

cmap=[92 123 207;220 110 75;68 150 136;222 166 90;152 102 171;196 120 138]/255;

load_data=1;

if load_data==0
    cmap=[92 123 207;220 110 75;68 150 136;222 166 90;152 102 171;196 120 138]/255;

    param = [0.03 0.30 -90 0.20 0.04 0.05 0.02 0.02 0.07 0.08 0.15 150 25];

    I=[0.0 2.7];
    ramp_length=[flip(diff(I)./linspace(0.00015,0.3,30))]; % step in ms

    dt = 0.001;
    T=3500;
    stim_length=T-2500;
    Nt = round(T/dt);

    fr_vs_stim=zeros(1,length(ramp_length));
    isi_ratio=zeros(1,length(ramp_length));
    for i=1:length(ramp_length)
        disp([num2str(i),'/',num2str(length(ramp_length))])
        I_stim = I(1)*ones(Nt,1);
        I_stim(round(2500/dt):round((2500+stim_length)/dt)) = I(2);
        I_stim((round(2500/dt)-1):round((2500+ramp_length(i))/dt)) = linspace(I(1),I(2),round(ramp_length(i)/dt)+2);

        [Vs,Vd,Va,t,I_output] = doublet_3comp_model(I_stim, 'soma', 'euler', param);

        [~,locs]=findpeaks(Vs,'MinPeakHeight',-10);
        fr=1000./(dt*diff(locs));

        fr_vs_stim(i)=fr(1);
        isi_ratio(i)=fr(2)/fr(3);
    end
    save('data/figure-files/ramp_data_doublet.mat','I','ramp_length','fr_vs_stim','isi_ratio');
else
    load('data/figure-files/ramp_data_doublet.mat');
end

ramp_speed=diff(I)./ramp_length;

figure;set(gcf,'units','points','position',[549,402,317,278]);
yyaxis left;
hold on;
plot(ramp_speed(1:1:end),fr_vs_stim(1:1:end),'o','MarkerFaceColor',cmap(5,:),'Color',cmap(5,:),'MarkerSize',8);
hold off;
ylim([0 200]);
yticks(0:50:200);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);

yyaxis right;
hold on;
plot(ramp_speed(1:1:end),isi_ratio(1:1:end),'o','MarkerFaceColor',cmap(4,:),'Color',cmap(4,:),'MarkerSize',8);
hold off;
ylim([0.8 1.2]);
yticks(0.8:0.1:1.2);
xlim([0 0.3]);
xticks(0.0:0.1:0.3);
hold on;
plot(xlim,[1.0 1.0],':','Color',cmap(4,:));
hold off;

set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
ax = gca;
ax.YAxis(1).Color = cmap(5,:);
ax.YAxis(2).Color = cmap(4,:);
g=gcf;
g.Renderer='painters';

%% Panel F

clearvars; close all;

load_data=1;

param = [0.03 0.30 -90 0.20 0.04 0.05 0.02 0.02 0.07 0.08 0.15 150 25];

if load_data==0

    dt=0.001;
    % Set
    I=linspace(0.2,0.3,11);
    I_step=1.00;
    comp_mV=30; % ms
    delta_nA=0.2;
    doublet_thresh=100; % ms

    % Find the holding current level that is stable
    for i=1:length(I)
        disp(i)
        [Vs,~,~,t,I_output] = doublet_3comp_model([I(i) I_step], 'soma', 'euler', param);
        [~,locs]=findpeaks(Vs(((length(Vs)/4):(length(Vs)/2))),'MinPeakHeight',0,'MinPeakDistance',2/dt);
        if ~isempty(locs)
            I_rb=I(i);
            break;
        end
    end

    [Vs,~,~,t,I_output] = doublet_3comp_model([I(1)-0.03 I(1)-0.03], 'soma', 'euler', param);
    V2=mean(Vs((length(Vs)-comp_mV/dt):length(Vs)));

    [Vs,~,~,t,I_output] = doublet_3comp_model([I(1)-0.03-delta_nA I(1)-0.03-delta_nA], 'soma', 'euler', param);
    V1=mean(Vs((length(Vs)-comp_mV/dt):length(Vs)));

    R_in=(V2-V1)/delta_nA;

    % Set Rin-normalised step currents
    mV_steps=5;
    mV_eq_steps=[20 80]; % down (holding) and up (step size)
    num_eq_steps=mV_eq_steps./mV_steps;

    nA_steps=mV_steps/R_in;

    hold_off_steps=flip((I_rb-nA_steps):-nA_steps:(I_rb-num_eq_steps(1)*nA_steps));
    curr_inj_steps=(I_rb+nA_steps):nA_steps:(I_rb+num_eq_steps(2)*nA_steps);

    doublet_mat=-1.*ones(length(hold_off_steps),length(curr_inj_steps));
    singlet_mat=-1.*ones(length(hold_off_steps),length(curr_inj_steps));
    fr_mat=zeros(length(hold_off_steps),length(curr_inj_steps));
    time2spike_mat=zeros(length(hold_off_steps),length(curr_inj_steps));
    AI_mat=zeros(length(hold_off_steps),length(curr_inj_steps));

    for i=1:length(hold_off_steps)
        for j=1:length(curr_inj_steps)
            disp([num2str(i),'/',num2str(length(hold_off_steps)),' ',num2str(j),'/',num2str(length(curr_inj_steps))]);
            [Vs,~,~,t,I_output] = doublet_3comp_model([hold_off_steps(i) curr_inj_steps(j)], 'soma', 'euler', param);
            [~,locs]=findpeaks(Vs,'MinPeakHeight',0,'MinPeakDistance',2/dt);
            if ~isempty(locs)
                fr=1000./(diff(locs)/(1/dt));
                if fr(1)>=doublet_thresh
                    % Find first location of singlet
                    first_singlet=find(fr<doublet_thresh);
                    first_singlet=first_singlet(1);
                    
                    % If triplet or more, remove

                    % Store number of doublets and singlets
                    doublet_mat(i,j)=sum(fr>=doublet_thresh)-sum(diff(find(fr>=doublet_thresh))==1)-sum(diff(find(fr>=doublet_thresh))>2);
                    singlet_mat(i,j)=sum(fr<doublet_thresh)-doublet_mat(i,j);

                    % Store inst. freq. of doublet
                    fr_mat(i,j)=fr(1);

                    % Store time to peak
                    time2spike_mat(i,j)=t(locs(1))-1e3;

                    % Store adaptation index (AI)
                    fr=fr(fr>=doublet_thresh);
                    if length(fr)>=4
                        AI_mat(i,j)=(1e3/fr(end-1))/(1e3/fr(2));
                    end
                else
                    doublet_mat(i,j)=0;
                end
            end
        end
    end

    save('data/figure-files/holding_step_data_initial.mat','hold_off_steps','I_rb','R_in','curr_inj_steps','doublet_mat','singlet_mat','fr_mat','time2spike_mat','AI_mat','num_eq_steps');
else
    load('data/figure-files/holding_step_data_initial.mat');
end

mV_steps=5;
mV_eq_steps=[20 80]; % down (holding) and up (step size)
num_eq_steps=mV_eq_steps./mV_steps;

cmap=hot;

s1=figure;set(gcf,'units','points','position',[587,512,290,231]);
imagesc((hold_off_steps-I_rb),(curr_inj_steps-I_rb),doublet_mat');
cb1=colorbar;
cb1.TickLength=0;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);

colormap(s1,cmap([1 128],:));

%% Panel G
clearvars; close all;

cmap=[92 123 207;220 110 75;68 150 136;222 166 90;152 102 171;196 120 138]/255;

param = [0.03 0.30 -90 0.20 0.04 0.05 0.02 0.02 0.07 0.08 0.15 150 25];

[Vs1,Vd1,Va1,t1,I_output1] = doublet_3comp_model([0.00 1.0], 'soma', 'euler', param);
[~,locs1]=findpeaks(Vs1,'MinPeakHeight',-10);locs1=locs1(1);

[Vs2,Vd2,Va2,t2,I_output2] = doublet_3comp_model([0.00 1.5], 'soma', 'euler', param);
[~,locs2]=findpeaks(Vs2,'MinPeakHeight',-10);locs2=locs2(1);

[Vs3,Vd3,Va3,t3,I_output3] = doublet_3comp_model([0.00 2.5], 'soma', 'euler', param);
[~,locs3]=findpeaks(Vs3,'MinPeakHeight',-10);locs3=locs3(1);

figure;set(gcf,'units','points','position',[659,549,150,233]);
hold on;
plot(t1(1:end-(locs1-locs3-1)),Vd1((locs1-locs3):end),'--','LineWidth',1.5,'Color',cmap(3,:));
plot(t2(1:end-(locs2-locs3-1)),Vd2((locs2-locs3):end),':','LineWidth',1.5,'Color',cmap(3,:));
plot(t3,Vd3,'LineWidth',1.5,'Color',cmap(3,:));
hold off;
xlim([995 1020]);ylim([-80 30]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

figure;set(gcf,'units','points','position',[659,549,150,233]);
hold on;
plot(t1(1:end-(locs1-locs3-1)),Vs1((locs1-locs3):end),'--','LineWidth',1.5,'Color',cmap(1,:));
plot(t2(1:end-(locs2-locs3-1)),Vs2((locs2-locs3):end),':','LineWidth',1.5,'Color',cmap(1,:));
plot(t3,Vs3,'LineWidth',1.5,'Color',cmap(1,:));
hold off;
xlim([995 1020]);ylim([-80 30]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

figure;set(gcf,'units','points','position',[659,549,150,233]);
hold on;
plot(t1(1:end-(locs1-locs3-1)),Va1((locs1-locs3):end),'--','LineWidth',1.5,'Color',cmap(2,:));
plot(t2(1:end-(locs2-locs3-1)),Va2((locs2-locs3):end),':','LineWidth',1.5,'Color',cmap(2,:));
plot(t3,Va3,'LineWidth',1.5,'Color',cmap(2,:));
hold off;
xlim([995 1020]);ylim([-80 30]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

figure;set(gcf,'units','points','position',[659,549,150,233]);
hold on;
plot(t1(1:end-(locs1-locs3-1)),I_output1.IN_d((locs1-locs3):end)+I_output1.IP_d((locs1-locs3):end)+I_output1.IT_d((locs1-locs3):end),'--','LineWidth',1.5,'Color',cmap(3,:));
plot(t2(1:end-(locs2-locs3-1)),I_output2.IN_d((locs2-locs3):end)+I_output2.IP_d((locs2-locs3):end)+I_output2.IT_d((locs2-locs3):end),':','LineWidth',1.5,'Color',cmap(3,:));
plot(t3,I_output3.IN_d+I_output3.IP_d+I_output3.IT_d,'LineWidth',1.5,'Color',cmap(3,:));
hold off;
xlim([995 1020]);
ylim([-4 0]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

figure;set(gcf,'units','points','position',[659,549,150,233]);
hold on;
plot(t1(1:end-(locs1-locs3-1)),I_output1.INaP_s((locs1-locs3):end),'--','LineWidth',1.5,'Color',cmap(1,:));
plot(t2(1:end-(locs2-locs3-1)),I_output2.INaP_s((locs2-locs3):end),':','LineWidth',1.5,'Color',cmap(1,:));
plot(t3,I_output3.INaP_s,'LineWidth',1.5,'Color',cmap(1,:));
hold off;
xlim([995 1020]);
ylim([-4 0]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

figure;set(gcf,'units','points','position',[659,549,150,233]);
hold on;
plot(t1(1:end-(locs1-locs3-1)),I_output1.INaP_ax((locs1-locs3):end),'--','LineWidth',1.5,'Color',cmap(2,:));
plot(t2(1:end-(locs2-locs3-1)),I_output2.INaP_ax((locs2-locs3):end),':','LineWidth',1.5,'Color',cmap(2,:));
plot(t3,I_output3.INaP_ax,'LineWidth',1.5,'Color',cmap(2,:));
hold off;
xlim([995 1020]);
ylim([-4 0]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

%% Panel H

figure;set(gcf,'units','points','position',[659,549,284,233]);
hold on;
plot(t3,Vs3,'LineWidth',1.5,'Color',cmap(1,:));
plot(t3,Va3,'LineWidth',1.5,'Color',cmap(2,:));
plot(t3,Vd3,'LineWidth',1.5,'Color',cmap(3,:));
hold off;
xlim([999 1014]);ylim([-80 30]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

figure;set(gcf,'units','points','position',[659,549,284,233]);
hold on;
plot(t3,I_output3.INaP_s,'LineWidth',1.5,'Color',cmap(1,:));
plot(t3,I_output3.INaP_ax,'LineWidth',1.5,'Color',cmap(2,:));
plot(t3,I_output3.IN_d+I_output3.IP_d+I_output3.IT_d,'LineWidth',1.5,'Color',cmap(3,:));
hold off;
xlim([999 1014]);
ylim([-4 0]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

%% Panel I
clearvars; close all;

dt=0.001;

cmap=[92 123 207;220 110 75;68 150 136;222 166 90;152 102 171;196 120 138]/255;

param = [0.03 0.30 -90 0.20 0.04 0.05 0.02 0.02 0.07 0.08 0.15 150 25];
[Vs1,Vd1,Va1,t1,I_output1] = doublet_3comp_model([0.0 2.1], 'soma', 'euler', param);
[~,locs1]=findpeaks(Vs1,'MinPeakHeight',-10);

param(1)=param(1)*2;
[Vs2,Vd2,Va2,t2,I_output2] = doublet_3comp_model([0.0 2.1], 'soma', 'euler', param);
[~,locs2]=findpeaks(Vs2,'MinPeakHeight',-10);

figure;set(gcf,'units','points','position',[659,549,150,233]);
hold on;
plot(t1(1:end-(locs1(1)-locs2(1)-1)),Vs1((locs1(1)-locs2(1)):end),'-','LineWidth',1.5,'Color',cmap(1,:));
plot(t2,Vs2,'--','LineWidth',1.5,'Color',cmap(1,:));
hold off;
xlim([997.5 1015]);ylim([-80 30]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

figure;set(gcf,'units','points','position',[659,549,150,233]);
hold on;
plot(t1(1:end-(locs1(1)-locs2(1)-1)),I_output1.INaP_s((locs1(1)-locs2(1)):end),'-','LineWidth',1.5,'Color',cmap(1,:));
plot(t2,I_output2.INaP_s,'--','LineWidth',1.5,'Color',cmap(1,:));
hold off;
xlim([997.5 1015]);
ylim([-8 0]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

figure;set(gcf,'units','points','position',[659,549,150,233]);
hold on;
plot(t1(1:end-(locs1(1)-locs2(1)-1)),I_output1.INaP_ax((locs1(1)-locs2(1)):end),'-','LineWidth',1.5,'Color',cmap(2,:));
plot(t2,I_output2.INaP_ax,'--','LineWidth',1.5,'Color',cmap(2,:));
hold off;
xlim([997.5 1015]);
ylim([-8 0]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

%% Panel J
clearvars; close all;

dt=0.001;

param = [0.03 0.30 -90 0.20 0.04 0.05 0.02 0.02 0.07 0.08 0.15 150 25];
[Vs1,Vd1,Va1,t1,I_output1] = doublet_3comp_model([0.0 2.4], 'soma', 'euler', param);
[~,locs1]=findpeaks(Vs1,'MinPeakHeight',0);
fr1=1000./(dt*diff(locs1));

param = [0.00 0.30 -90 0.20 0.04 0.05 0.02 0.02 0.07 0.08 0.15 150 25];
[Vs2,Vd2,Va2,t2,I_output2] = doublet_3comp_model([0.0 2.4], 'soma', 'euler', param);
[~,locs2]=findpeaks(Vs2,'MinPeakHeight',0);
fr2=1000./(dt*diff(locs2));

param = [0.03 0.30 -90 0.20 0.00 0.05 0.02 0.02 0.07 0.08 0.15 150 25];
[Vs3,Vd3,Va3,t3,I_output3] = doublet_3comp_model([0.0 2.4], 'soma', 'euler', param);
[~,locs3]=findpeaks(Vs3,'MinPeakHeight',0);
fr3=1000./(dt*diff(locs3));

% param = [0.03 0.30 0.10 -90 0.20 0.04 0.20 0.60 0.004 0.040 ...
%             0.05 0.02 0.02 0.00 0.30 0.025 0.20 0.07 0.0 0.00 0.00];
param = [0.03 0.30 -90 0.20 0.04 0.05 0.02 0.02 0.07 0.00 0.00 150 25];
[Vs4,Vd4,Va4,t4,I_output4] = doublet_3comp_model([0.0 2.4], 'soma', 'euler', param);
[~,locs4]=findpeaks(Vs4,'MinPeakHeight',0);
fr4=1000./(dt*diff(locs4));

cmap=[92 123 207;220 110 75;68 150 136;222 166 90;152 102 171;196 120 138]/255;

figure;set(gcf,'units','points','position',[648,589,270,193]);
hold on;
plot(t1,Vs1,'-','LineWidth',1.5,'Color',cmap(1,:));
plot(t2(1:end-(locs2(1)-locs1(1)-1)),Vs2((locs2(1)-locs1(1)):end),'LineWidth',1.5,'Color',cmap(4,:));
plot(t3(1:end-(locs3(1)-locs1(1)-1)),Vs3((locs3(1)-locs1(1)):end),'LineWidth',1.5,'Color',cmap(5,:));
plot(t4(1:end-(locs4(1)-locs1(1)-1)),Vs4((locs4(1)-locs1(1)):end),'-','LineWidth',1.5,'Color',cmap(6,:));
hold off;
xlim([997.5 1015]);ylim([-80 30]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

%% Panel K

clearvars; close all;

dt=0.001;

gSK_sweep=0.15:0.05:0.45;

Vs=zeros(length(gSK_sweep),2e6);
for i=1:length(gSK_sweep)
    disp(i)
    param = [0.03 gSK_sweep(i) -90 0.20 0.04 0.05 0.02 0.02 0.07 0.08 0.15 150 25];
    [Vs(i,:),Vd,Va,t,I_output] = doublet_3comp_model([0.00 2.70], 'soma', 'euler', param);
end

cmap=turbo(length(gSK_sweep));

figure;set(gcf,'units','points','position',[648,589,270,193]);
hold on;
for i=1:size(Vs,1)
    plot(t,Vs(i,:),'LineWidth',1.5,'Color',cmap(i,:));
end
hold off;
xlim([997.5 1015]);ylim([-80 30]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

%% Panel L
clearvars; close all;

dt=0.001;

param = zeros(1,13);

I_range=[1.0:0.5:4.0];

cmap=turbo;
cmap=cmap(round(linspace(1,256,length(I_range))),:);

Vs=zeros(length(I_range),2e6);
for i=1:length(I_range)
    disp(i)
    [Vs(i,:),~,~,t,I_output] = doublet_3comp_model([0.00 I_range(i)], 'soma', 'euler', param);
end

figure;set(gcf,'units','points','position',[549,402,317,278]);
hold on;
for i=1:size(Vs,1)
    plot(t,Vs(i,:),'Color',cmap(i,:),'LineWidth',2);
end
hold off;
axis tight;
xlim([995 1010]);
ylim([-80 30]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

figure;set(gcf,'units','points','position',[554,453,317,112]);
hold on;
plot(t,I_output.I_stim,'Color','k','LineWidth',2);
hold off;
axis tight;
xlim([995 1010]);
ylim([0.0 10.0]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';
