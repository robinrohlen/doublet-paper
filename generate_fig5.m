addpath('simulation model');
addpath('functions');

%% Panel B
clearvars; close all;

dt=0.001;

cmap=[92 123 207;220 110 75;68 150 136;222 166 90;152 102 171;196 120 138]/255;

param = [0.03 0.30 -85 0.20 0.04 0.05 0.02 0.02 0.07 0.08 0.15 150 25];

[Vs,Vd,Va,t,I_output] = doublet_3comp_model([-0.30 0.30], 'soma', 'euler', param);

figure;set(gcf,'units','points','position',[488,585,631,175]);
hold on;
plot(t,Vd,'Color',cmap(3,:),'LineWidth',2);
hold off;
axis tight;
xlim([975 1600]);
ylim([-85 30]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

figure;set(gcf,'units','points','position',[488,585,631,175]);
hold on;
plot(t,Vs,'Color',cmap(1,:),'LineWidth',2);
hold off;
axis tight;
xlim([975 1600]);
ylim([-85 30]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

figure;set(gcf,'units','points','position',[488,585,631,175]);
hold on;
plot(t,Va,'Color',cmap(2,:),'LineWidth',2);
hold off;
axis tight;
xlim([975 1600]);
ylim([-85 30]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

figure;set(gcf,'units','points','position',[488,585,631,175]);
hold on;
plot(t,I_output.I_stim,'Color','k','LineWidth',2);
hold off;
axis tight;
xlim([975 1600]);
ylim([-0.3 2.0]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
set(gca,'Visible','off')
g=gcf;
g.Renderer='painters';

%% Panel C

[~,locs]=findpeaks(Vs,'MinPeakHeight',-10,'MinPeakDistance',2/dt);
fr=1000./(dt*diff(locs));

figure;set(gcf,'units','points','position',[591,490,249,218]);%[648,484,169,298]);
hold on;
plot_broken_axis(locs(1:end-1)/(1/dt),fr,75,[10 100],[100 200],10,[900 2000],cmap(1,:),1,0);
hold off;
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

%% Panel D
clearvars;close all;

load_data=1;

dt=0.001;
gc_range=flip(0.10:0.01:0.25);

if load_data==0
    for i=1:length(gc_range)
        disp(i)
        param = [0.03 0.30 -85 gc_range(i) 0.04 0.05 0.02 0.02 0.07 0.08 0.15 150 25];
        [Vs,Vd,Va,t,I_output] = doublet_3comp_model([-0.30 0.30], 'soma', 'euler', param);

        [~,locs]=findpeaks(Vs,'MinPeakHeight',-10,'MinPeakDistance',2/dt);
        fr=1000./(dt*diff(locs));
        locs_low{i}=locs;locs_high{i}=locs;
        fr_low{i}=fr;fr_high{i}=fr;
        locs_high{i}(fr<100)=[];
        fr_high{i}(fr<100)=[];
        locs_low{i}(fr>=100)=[];
        fr_low{i}(fr>=100)=[];
    end

    save('data/figure-files/gc_sweep.mat','locs_low','locs_high','fr_high','fr_low','gc_range');
else
    load('data/figure-files/gc_sweep.mat');
end

cmap=flip(turbo(size(fr_high,2)));

figure;set(gcf,'units','points','position',[591,383,411,325]);
for i=1:size(fr_high,2)
    hold on;
    plot_broken_axis(locs_low{i}(1:end-1)*dt,fr_low{i},100,[10 100],[100 200],10,[900 2000],cmap(i,:),1,0,1);
    plot_broken_axis(locs_high{i}(1:end-1)*dt,fr_high{i},100,[10 100],[100 200],10,[900 2000],cmap(i,:),1,0,1);
    hold off;
    disp(num2str(gc_range(i)))
    % pause;
end
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',16);
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
g=gcf;
g.Renderer='painters';

%% Panel E
clearvars; close all;

dt=0.001;

cmap=[92 123 207;220 110 75;68 150 136;222 166 90;152 102 171;196 120 138]/255;

param = [0.03 0.30 -85 0.20 0.04 0.05 0.02 0.02 0.07 0.08 0.15 150 25];
[Vs_control,~,~,t1,I_output1] = doublet_3comp_model([0 0.30], 'soma', 'euler', param);
[~,locs1]=findpeaks(Vs_control,'MinPeakHeight',-10);
fr1=1000./(dt*diff(locs1));

param = [0.00 0.30 -85 0.20 0.04 0.05 0.02 0.02 0.07 0.08 0.15 150 25];
[Vs_sNaP_block,~,~,t2,I_output2] = doublet_3comp_model([0 0.30], 'soma', 'euler', param);
[~,locs2]=findpeaks(Vs_sNaP_block,'MinPeakHeight',-10);
fr2=1000./(dt*diff(locs2));

param = [0.03 0.30 -85 0.20 0.00 0.05 0.02 0.02 0.07 0.08 0.15 150 25];
[Vs_aNaP_block,~,~,t3,I_output3] = doublet_3comp_model([0 0.70], 'soma', 'euler', param);
[~,locs3]=findpeaks(Vs_aNaP_block,'MinPeakHeight',-10);
fr3=1000./(dt*diff(locs3));

param = [0.03 0.30 -85 0.20 0.04 0.05 0.00 0.00 0.07 0.00 0.00 150 25];
[Vs_dendrCa_block,~,~,t4,I_output4] = doublet_3comp_model([0 0.185], 'soma', 'euler', param);
[~,locs4]=findpeaks(Vs_dendrCa_block,'MinPeakHeight',-10);
fr4=1000./(dt*diff(locs4));

figure;set(gcf,'units','points','position',[659,549,150,233]);
hold on;
plot(t1,Vs_control,'-','LineWidth',1.5,'Color',cmap(1,:));
plot(t2(1:end-(locs2(1)-locs1(1)-1)),Vs_sNaP_block((locs2(1)-locs1(1)):end),'LineWidth',1.5,'Color',cmap(4,:));
plot(t3(1:end-(locs3(1)-locs1(1)-1)),Vs_aNaP_block((locs3(1)-locs1(1)):end),'LineWidth',1.5,'Color',cmap(5,:));
plot(t4(1:end-(locs4(1)-locs1(1)-1)),Vs_dendrCa_block((locs4(1)-locs1(1)):end),'-','LineWidth',1.5,'Color',cmap(6,:));
hold off;
xlim([1015 1035]);ylim([-80 30]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

figure;set(gcf,'units','points','position',[668,470,218,144]);
hold on;
plot(t1,I_output1.I_stim,'Color',cmap(1,:),'LineWidth',2);
plot(t2,I_output2.I_stim,'Color',cmap(4,:),'LineWidth',2);
plot(t3,I_output3.I_stim,'Color',cmap(5,:),'LineWidth',2);
plot(t4,I_output4.I_stim,'Color',cmap(6,:),'LineWidth',2);
hold off;
axis tight;
xlim([750 2000]);
ylim([0.0 1.0]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
set(gca,'Visible','off')
g=gcf;
g.Renderer='painters';

%% Panel F

figure;set(gcf,'units','points','position',[659,549,150,233]);
hold on;
plot(t1,I_output1.INaP_s,'-','LineWidth',1.5,'Color',cmap(1,:));
plot(t2(1:end-(locs2(1)-locs1(1)-1)),I_output2.INaP_s((locs2(1)-locs1(1)):end),'LineWidth',1.5,'Color',cmap(4,:));
plot(t3(1:end-(locs3(1)-locs1(1)-1)),I_output3.INaP_s((locs3(1)-locs1(1)):end),'LineWidth',1.5,'Color',cmap(5,:));
plot(t4(1:end-(locs4(1)-locs1(1)-1)),I_output4.INaP_s((locs4(1)-locs1(1)):end),'-','LineWidth',1.5,'Color',cmap(6,:));
hold off;
xlim([1015 1035]);ylim([-4 0]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

figure;set(gcf,'units','points','position',[659,549,150,233]);
hold on;
plot(t1,I_output1.INaP_ax,'-','LineWidth',1.5,'Color',cmap(1,:));
plot(t2(1:end-(locs2(1)-locs1(1)-1)),I_output2.INaP_ax((locs2(1)-locs1(1)):end),'LineWidth',1.5,'Color',cmap(4,:));
plot(t3(1:end-(locs3(1)-locs1(1)-1)),I_output3.INaP_ax((locs3(1)-locs1(1)):end),'LineWidth',1.5,'Color',cmap(5,:));
plot(t4(1:end-(locs4(1)-locs1(1)-1)),I_output4.INaP_ax((locs4(1)-locs1(1)):end),'-','LineWidth',1.5,'Color',cmap(6,:));
hold off;
xlim([1015 1035]);ylim([-4 0]);
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',18);
g=gcf;
g.Renderer='painters';

%% Panel G

clearvars; close all;

load_data=1;

if load_data==0

    param = [0.03 0.30 -85 0.20 0.04 0.05 0.02 0.02 0.07 0.08 0.15 150 25];

    I=0.30; % nA
    dt=0.001;
    T=10/dt;
    from_spike=2; % from the second doublet spike

    beta = 1:6;

    f_ss_soma=zeros(1,length(beta));
    tau_adapt_soma=zeros(1,length(beta));
    f_ss_axon=zeros(1,length(beta));
    tau_adapt_axon=zeros(1,length(beta));

    for i=1:length(beta)
        disp([num2str(i)]);

        beta_tmp=beta(i);
        param(12)=150*beta_tmp;
        param(13)=25;

        % Soma
        [Vs,~,~,t,I_output] = doublet_3comp_model([-0.30 I], 'soma', 'euler', param, T);
        dt=diff(t(1:2));

        [~,locs]=findpeaks(Vs,'MinPeakHeight',-20,'MinPeakDistance',2/dt);
        fr=1000./(dt*diff(locs));
        locs(length(locs))=[];

        locs=locs(fr>=100);
        fr=fr(fr>=100);

        % f_intra = your vector of intra-doublet (high-frequency) values, in order
        n = locs*dt-locs(1)*dt;%(1:numel(fr))';               % spike/doublet time

        % Model: f(n) = f_ss + (f0 - f_ss)*exp(-(n-1)/tau)
        % Parameters: p(1)=f_ss (plateau), p(2)=f0-f_ss (amplitude), p(3)=tau
        model = @(p,n) p(1) + p(2).*exp(-(n-1)./p(3));

        % Initial guesses from the data
        p0 = [fr(end), fr(from_spike)-fr(end), beta_tmp*150];   % [plateau, amplitude, tau]

        cost = @(p) sum((model(p,n(from_spike:end)) - fr(from_spike:end)).^2);
        p_fit = fminsearch(cost, p0);
        
        f_ss_soma(i) = p_fit(1);
        tau_adapt_soma(i) = p_fit(3);
        model_fit_soma{i}=p_fit;
        fr_save_soma{i}=fr;
        locs_save_soma{i}=n;
        INaP_save_soma{i}=I_output.INaP_s;

        % Axon
        param(12)=150;
        param(13)=25*beta_tmp;
        [Vs,~,~,t,I_output] = doublet_3comp_model([-0.30 I], 'soma', 'euler', param, T);
        dt=diff(t(1:2));

        [~,locs]=findpeaks(Vs,'MinPeakHeight',-20,'MinPeakDistance',2/dt);
        fr=1000./(dt*diff(locs));
        locs(length(locs))=[];

        locs=locs(fr>=100);
        fr=fr(fr>=100);

        % f_intra = your vector of intra-doublet (high-frequency) values, in order
        n = locs*dt-locs(1)*dt;%(1:numel(fr))';               % spike/doublet time

        % Model: f(n) = f_ss + (f0 - f_ss)*exp(-(n-1)/tau)
        % Parameters: p(1)=f_ss (plateau), p(2)=f0-f_ss (amplitude), p(3)=tau
        model = @(p,n) p(1) + p(2).*exp(-(n-1)./p(3));

        % Initial guesses from the data
        p0 = [fr(end), fr(from_spike)-fr(end), beta_tmp*25];   % [plateau, amplitude, tau]

        cost = @(p) sum((model(p,n(from_spike:end)) - fr(from_spike:end)).^2);
        p_fit = fminsearch(cost, p0);
        
        f_ss_axon(i) = p_fit(1);
        tau_adapt_axon(i) = p_fit(3);
        model_fit_axon{i}=p_fit;
        fr_save_axon{i}=fr;
        locs_save_axon{i}=n;
        INaP_save_axon{i}=I_output.INaP_s;

        % Plot
        % figure; plot(n, fr, 'o'); hold on;
        % plot(n(from_spike:end), model(p_fit,n(from_spike:end)), '-');
        % xlabel('doublet index'); ylabel('intra-doublet frequency (Hz)');
        % legend('data','exp fit');
    end

    save('data/figure-files/NaP_adaptation.mat','beta','from_spike','model','I','f_ss_soma','tau_adapt_soma','model_fit_soma','fr_save_soma','locs_save_soma','f_ss_axon','tau_adapt_axon');
else
    load('data/figure-files/NaP_adaptation.mat');
end

cmap=[92 123 207;220 110 75;68 150 136;222 166 90;152 102 171;196 120 138]/255;

figure;set(gcf,'units','points','position',[591,383,411,325]);
hold on;
plot(locs_save_soma{1},fr_save_soma{1},'o','Color','k','MarkerFaceColor',cmap(5,:),'MarkerSize',8);%,'LineWidth',2);
plot(locs_save_soma{1}(from_spike:end), model(model_fit_soma{1},locs_save_soma{1}(from_spike:end)), '-','Color','k','LineWidth',2);
plot(locs_save_soma{6},fr_save_soma{6},'o','Color','k','MarkerFaceColor',cmap(6,:),'MarkerSize',8);%,'LineWidth',2);
plot(locs_save_soma{6}(from_spike:end), model(model_fit_soma{6},locs_save_soma{6}(from_spike:end)), '-','Color','k','LineWidth',2);
hold off;
axis tight;
ylim([130 190]);
yticks([130:20:190]);
xlim([-200 7500]);
ylabel('Intra-doublet firing rate (Hz)');
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',16);
g=gcf;
g.Renderer='painters';

%% Panel H

clearvars; close all;

try
    load('data/figure-files/NaP_adaptation.mat');
catch
    error('Cannot load NaP_adaptation.mat');
end

cmap=[92 123 207;220 110 75;68 150 136;222 166 90;152 102 171;196 120 138]/255;

figure;set(gcf,'units','points','position',[591,383,411,325]);
yyaxis left
hold on;
plot(beta*150,f_ss_soma,'-o','Color',cmap(5,:),'MarkerFaceColor',cmap(5,:),'LineWidth',6,'MarkerSize',10);
plot(beta*150,f_ss_axon,':d','Color',cmap(5,:),'MarkerFaceColor',cmap(5,:),'LineWidth',6,'MarkerSize',10);
hold off;
xlim([0 1000]);
xticks([0:500:1000]);
ylim([140 190]);
yticks([140:25:190]);
xlabel('NaP inactivation time (ms)');
ylabel('Steady-state frequency (Hz)');
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',16);

yyaxis right
hold on;
plot(beta*150,tau_adapt_soma,'-o','Color',cmap(6,:),'MarkerFaceColor',cmap(6,:),'LineWidth',6,'MarkerSize',10);
plot(beta*150,tau_adapt_axon,':d','Color',cmap(6,:),'MarkerFaceColor',cmap(6,:),'LineWidth',6,'MarkerSize',10);
hold off;
xlim([0 1000]);
xticks([0:500:1000]);
ylim([0 1000]);
yticks([0:500:1000]);
ylabel('Adaptation time constant (ms)');
set(gca,'TickDir','out');
set(gcf,'color','w');
set(gca,'FontSize',16);
ax = gca;
ax.YAxis(1).Color = cmap(5,:);
ax.YAxis(2).Color = cmap(6,:);
g=gcf;
g.Renderer='painters';

%% Panel I
clearvars; close all;

load_data=1;

param = [0.03 0.30 -85 0.20 0.04 0.05 0.02 0.02 0.07 0.08 0.15 150 25];

if load_data==0

    dt=0.001;
    % Set
    I=linspace(0.1,0.2,11);
    I_step=0.80;
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
    mV_steps=1;
    mV_eq_steps=[20 20]; % down (holding) and up (step size)
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

    save('data/figure-files/holding_step_data_repetitive.mat','hold_off_steps','I_rb','R_in','curr_inj_steps','doublet_mat','singlet_mat','fr_mat','time2spike_mat','AI_mat','num_eq_steps');
else
    load('data/figure-files/holding_step_data_repetitive.mat');
end

mV_steps=1;
mV_eq_steps=[20 20]; % down (holding) and up (step size)
num_eq_steps=mV_eq_steps./mV_steps;

cmap=hot;

s1=figure;set(gcf,'units','points','position',[587,512,290,231]);
imagesc((hold_off_steps-I_rb),(curr_inj_steps-I_rb),doublet_mat');
xticks([-0.7 -0.5 -0.3 -0.1]);
yticks([0.1:0.2:0.7]);
cb1=colorbar;
cb1.TickLength=0;
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
colormap(s1,cmap(linspace(1,256,6),:));
