addpath('functions')

%% Panels A-C

rng(100);

cmap=[92 123 207;220 110 75;68 150 136;222 166 90;152 102 171;196 120 138]/255;

fs=1e3;
T=5; % s
isi_cov=0.10;

twitch_gain=1.035;
num_MUs=[1 100 200];

for n=1:length(num_MUs)
    tot_force_singlet=zeros(1,(T+1)*1e3);
    tot_force_doublet=zeros(1,(T+1)*1e3);
    tot_force_halfpop=zeros(1,(T+1)*1e3);

    for i=1:num_MUs(n)
        disp(i);
        % Randomly generate for each iteration
        fr_singlet=11.2+(13.6-11.2).*rand(1,1);
        fr_doublet=(0.56+(0.75-0.56).*rand(1,1))*fr_singlet;

        Tc=75+(90-75).*rand(1,1);
        Thr=(5/3)*Tc;
        [t,F]=RaikovaForceTwitch5p(1e3,0,0,Tc,Thr,500,1);

        num_firings=round(T*fr_singlet);
        start_spike_time=1+(5*(1000/fr_singlet)-1).*rand(1,1);
        locs_singlet=round(cumsum([start_spike_time 1000/fr_singlet+isi_cov*(1000/fr_singlet)*randn(1,num_firings-1)]));
        st=zeros(1,T*fs);
        st(locs_singlet)=1;
        force_singlet=conv(st,F);

        num_firings=round(T*fr_doublet);
        intra_doublet_isi=5*ones(1,num_firings);
        locs_tmp=[start_spike_time 1000/fr_doublet+isi_cov*(1000/fr_doublet)*randn(1,num_firings-1)];

        tmp=0;
        for ind=1:num_firings
            tmp=tmp+1;
            locs_doublet(tmp)=locs_tmp(ind);

            tmp=tmp+1;
            locs_doublet(tmp)=intra_doublet_isi(ind);
        end

        locs_doublet=cumsum(round(locs_doublet));
        locs_doublet=locs_doublet(locs_doublet<=locs_singlet(end));
        st=zeros(1,T*fs);
        st(locs_doublet)=1;
        force_doublet=twitch_gain*conv(st,F);

        clearvars locs_doublet locs_singlet

        % singlet summation
        try
            tot_force_singlet=tot_force_singlet+force_singlet(1:length(tot_force_singlet));
        catch
            tot_force_singlet=tot_force_singlet+[force_singlet zeros(1,length(tot_force_singlet)-length(force_singlet))];
        end

        % doublet summation
        try
            tot_force_doublet=tot_force_doublet+force_doublet(1:length(tot_force_doublet));
        catch
            tot_force_doublet=tot_force_doublet+[force_doublet zeros(1,length(tot_force_doublet)-length(force_doublet))];
        end

        % half singlets + half doublets
        if num_MUs(n)>1
            try
                if i<=num_MUs(n)/2
                    tot_force_halfpop=tot_force_halfpop+force_singlet(1:length(tot_force_singlet));
                else
                    tot_force_halfpop=tot_force_halfpop+force_doublet(1:length(tot_force_doublet));
                end
            catch
                if i<=num_MUs(n)/2
                    tot_force_halfpop=tot_force_halfpop+[force_singlet zeros(1,length(tot_force_singlet)-length(force_singlet))];
                else
                    tot_force_halfpop=tot_force_halfpop+[force_doublet zeros(1,length(tot_force_doublet)-length(force_doublet))];
                end
            end
        end
    end

    figure;set(gcf,'units','points','position',[345,489,310,265]);

    hold on;
    plot(linspace(0,length(tot_force_singlet)/1e3,length(tot_force_singlet)),tot_force_doublet,'Color',cmap(2,:),'LineWidth',2);
    plot(linspace(0,length(tot_force_singlet)/1e3,length(tot_force_singlet)),tot_force_singlet,'Color',cmap(1,:),'LineWidth',2);
    hold off;

    if num_MUs(n)>1
        hold on;
        plot(linspace(0,length(tot_force_singlet)/1e3,length(tot_force_singlet)),tot_force_halfpop,'Color',cmap(3,:),'LineWidth',2);
        hold off;
    end

    if n==2
        l=legend('Singlet spiking','Doublet spiking','50% singlet + 50% doublet spiking','Location','south');
        l.Box='off';
    end

    xlabel('Time (s)');
    ylabel('Force (n.u.)');
    set(gca,'TickDir','out');
    set(gcf,'color','w');
    set(gca,'FontSize',12);

    if n==1
        title('N=1 motor unit forces','FontWeight','normal');
    elseif n==2
        title('Sum of N=100 motor unit forces','FontWeight','normal');
    else
        title('Sum of N=200 motor unit forces','FontWeight','normal');
    end
    g=gcf;
    g.Renderer='painters';
end
