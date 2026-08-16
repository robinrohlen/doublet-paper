addpath('functions');

clearvars; close all;

% Set params
num_spikes_rt=3; % x spikes within one second
isi_crit=250; % ms
doublet_crit=100; % Hz

mvc_levels1={'05','10','15','20','25'};
mvc_levels2={'02','05','10','20'};

singlet_fr=[];
doublet_fr=[];

iter=0;
for subj_nr=15:24
    disp(['Processing subject ',num2str(subj_nr),'...']);

    clearvars muap_unique doublet_unique rt_unique

    singlet_force_rt=[];
    doublet_force_rt=[];

    % Experiment 1
    % Loop through every subject
    if subj_nr<=14
        muap_all=cell(1,size(mvc_levels1,2));
        RT_all=cell(1,size(mvc_levels1,2));
        doublet_all=cell(1,size(mvc_levels1,2));
        % Loop through every trapezoid recording
        for file_num=1:size(mvc_levels1,2)
            if subj_nr<10
                tmp=['/Volumes/Expansion/Reproducing doublet paper figures/data/sub-0',num2str(subj_nr),'/emg/sub-0',num2str(subj_nr),'_task-trapezoid',mvc_levels1{file_num},'percentmvc_run-01_emg.mat'];
            else
                tmp=['/Volumes/Expansion/Reproducing doublet paper figures/data/sub-',num2str(subj_nr),'/emg/sub-',num2str(subj_nr),'_task-trapezoid',mvc_levels1{file_num},'percentmvc_run-01_emg.mat'];
            end
            load(tmp)

            fs=signal.fsamp;
            [b,a]=butter(3,10/(2e3/2),'low');
            win=-50:50; % +/- 50 samples (~ +/- 25 ms) for triggering
            t = 1e3*[0:length(win)-1]/fs; % time vector for figures

            % Pre-define muap cell
            muap=cell(1,size(edition.Distimeclean{1},2));
            RT=zeros(1,size(edition.Distimeclean{1},2));
            doublet=zeros(1,size(edition.Distimeclean{1},2));

            % Filter raw EMG signals (50 Hz notch and 20-500 Hz bandpass)
            sig=notchsignals(signal.data(1:64,:),signal.fsamp);
            sig=bandpassingals(sig,signal.fsamp);

            % Find zeros due to data loss
            find_zeros=find(signal.data(65,:)==0);

            % Remove zeros in data
            signal.data(:,find_zeros)=[];

            try
                force=filtfilt(b,a,signal.data(72,:));
                force=force-force(1);
                force=force./max(force);
                force=(max(signal.path(10001:110000))-signal.path(1)).*force./max(force)+signal.path(1);
                force=abs(force);
            catch
                force=filtfilt(b,a,signal.data(72,:));
                force=force-force(1);
                force=force./max(force);
                force=file_num.*5.*force./mean(force(round(length(force)/2)+(-20000:20000)));
                force=abs(force);
            end

            % find lag between Ssq+ and SyncStation
            [r,lags]=xcorr(-signal.data(75,:)',-signal.data(69,:));
            find_lag=lags(r==max(r));

            % Adjust force based on sync
            if find_lag<=0
                force=[zeros(1,abs(find_lag)) force(1:(length(force)-abs(find_lag)))];
            else
                force=[force((find_lag+1):end) zeros(1,find_lag)];
            end

            iter=0;
            % Loop through each MU to extract its MUAP
            for MUnum=1:size(muap,2)
                iter=iter+1;

                % Extract the discharges used for triggering
                locs=edition.Distimeclean{1}{MUnum};

                % Estimate RT
                try
                    for i=0:(length(locs)-num_spikes_rt)
                        locs_tmp=locs((1:num_spikes_rt+1)+i);
                        if mean(round(1e3*diff(locs((1:num_spikes_rt+1)+i))/fs))<=isi_crit
                            break;
                        end
                    end
                catch
                    i=0;
                end

                RT(MUnum)=mean(force(locs((1:num_spikes_rt+1)+i)));

                % Set up the muap matrix (4 grids x 64 channels = 256 channels)
                muap{MUnum}=zeros(64,length(win));

                % Find doublets
                doublet_ind=find(1000./(diff(locs)/2) >= doublet_crit);

                if length(doublet_ind)>=2
                    doublet(MUnum)=1;
                    % AT WHAT FORCE LEVELS DID THE DOUBLETS APPEAR
                    doublet_force_rt=[doublet_force_rt force(locs(doublet_ind))./force(locs(doublet_ind(1)))];
                    if any(force(locs(doublet_ind))./RT(MUnum)<0.6)
                        %pause;
                    end
                end

                % Remove doublets in averaging
                if ~isempty(doublet_ind)
                    locs(sort([doublet_ind doublet_ind+1]))=[];
                    if ~isempty(locs)
                        singlet_force_rt=[singlet_force_rt force(locs)./force(locs(1))];
                    end
                end

                % If only doublets
                if isempty(locs) || length(locs)<10
                    locs=edition.Distimeclean{1}{MUnum};
                end

                % Extract the MUAP through STA for each channel
                for ch=1:size(muap{MUnum},1)
                    iterTrig=0;
                    for trig=3:length(locs)-3
                        muap{MUnum}(ch,:)=muap{MUnum}(ch,:)+sig(ch,locs(trig)+win);
                        iterTrig=iterTrig+1;
                    end
                    muap{MUnum}(ch,:)=muap{MUnum}(ch,:)./iterTrig;
                end
            end
            RT_all{file_num}=RT;
            muap_all{file_num}=muap;
            doublet_all{file_num}=doublet;
        end
        disp('Let''s find unique MUs...');
        % FIND UNIQUE MUAPS
        plot_match=0;
        xcorr_cutoff=0.65;

        matched_MUs=cell(1,size(mvc_levels1,2));

        for tmp=1:(size(mvc_levels1,2)-1)

            % Matched MUs
            matched_MUs{tmp}=zeros(size(muap_all{tmp},2),size(mvc_levels1,2));
            matched_MUs{tmp}(:,tmp)=1:size(muap_all{tmp},2);

            for file_num=(tmp+1):size(mvc_levels1,2)
                if subj_nr<10
                    tmpfile=['/Volumes/Expansion/Reproducing doublet paper figures/data/sub-0',num2str(subj_nr),'/emg/sub-0',num2str(subj_nr),'_task-trapezoid',mvc_levels1{file_num},'percentmvc_run-01_emg.mat'];
                else
                    tmpfile=['/Volumes/Expansion/Reproducing doublet paper figures/data/sub-',num2str(subj_nr),'/emg/sub-',num2str(subj_nr),'_task-trapezoid',mvc_levels1{file_num},'percentmvc_run-01_emg.mat'];
                end
                load(tmpfile)

                signal.data(1:64,:)=notchsignals(signal.data(1:64,:),signal.fsamp);
                signal.data(1:64,:)=bandpassingals(signal.data(1:64,:),signal.fsamp);

                R=16;
                eSIG = extension(signal.data(1:64,:),R);
                [wSIG, whitening_matrix] = whitening(eSIG,'ZCA');

                % Loop through each spatio-temporal MUAP from previous recording
                for MU=1:size(muap_all{tmp},2)

                    % Compute MU filter
                    w = muap_all{tmp}{MU};
                    w = extension(w,R);
                    w = whitening_matrix * w;

                    % Reconstruction on current recording
                    sig=w'*wSIG;

                    % Select the source with highest skewness
                    save_skew=zeros(1,size(sig,1));
                    for ind=1:size(sig,1)
                        save_skew(ind)=skewness(sig(ind,:));
                    end
                    [~,maxInd]=max(save_skew);
                    w = w(:,maxInd);
                    w = w./norm(w);

                    % Reconstruction
                    sig=w'*wSIG;

                    [Pulsetrain, Dischargetimes, sil] = calcSIL(wSIG, w, signal.fsamp);

                    r_vec=zeros(1,size(edition.Pulsetrainclean{1},1));
                    for i=1:size(edition.Pulsetrainclean{1},1)
                        [r,lags]=xcorr(Pulsetrain,edition.Pulsetrainclean{1}(i,:),2*R,'normalized');
                        r_vec(i)=max(r);
                    end
                    [~,indx]=find(r_vec==max(r_vec));

                    % Plot
                    if plot_match==1
                        figure(1);
                        title(['Subject ',num2str(subj_nr),' file num ',num2str(file_num), ' (MU ',num2str(MU),')']);
                        subplot(2,1,1);
                        plot(r_vec,'-o');ylim([0 1]);hold on;plot(xlim,[xcorr_cutoff xcorr_cutoff],'k:');hold off;
                        subplot(2,1,2);
                        plot(Pulsetrain);hold on;plot(edition.Pulsetrainclean{1}(indx,:));hold off;
                        pause;
                    end

                    % If above a certain threshold, keep it, otherwise []
                    if max(r_vec) >= xcorr_cutoff
                        matched_MUs{tmp}(MU,file_num)=indx;
                    end
                end
            end
        end

        % Find unique MUs
        tmptmp{1}=1:size(muap_all{1},2);
        tmptmp{2}=1:size(muap_all{2},2);
        tmptmp{3}=1:size(muap_all{3},2);
        tmptmp{4}=1:size(muap_all{4},2);
        tmptmp{5}=1:size(muap_all{5},2);

        cntr=0;
        for i=1:size(matched_MUs,2)
            for j=tmptmp{i}

                % count number of unique muaps
                cntr=cntr+1;

                % find matched mus
                if ~isempty(matched_MUs{i})
                    mu_nums=matched_MUs{i}(j,:);

                    % ignore zeros
                    tmp_ind=find(mu_nums~=0);

                    % take muap from latest recording having a match
                    selected_muap_rec=tmp_ind(end);

                    % save unique muap
                    muap_unique{cntr}=muap_all{selected_muap_rec}{mu_nums(selected_muap_rec)};

                    doublet_tmp=[];
                    rt_tmp=[];
                    for k=tmp_ind
                        % remove indexes for that unique muap
                        tmptmp{k}(find(tmptmp{k}==mu_nums(k)))=[];
                        doublet_tmp=[doublet_tmp doublet_all{k}(mu_nums(k))];
                        rt_tmp=[rt_tmp RT_all{k}(mu_nums(k))];
                    end

                    % save if doublet
                    doublet_unique(cntr)=any(doublet_tmp);

                    % save mean rt
                    rt_unique(cntr)=mean(rt_tmp);
                end
            end
        end

        doublet_unique((size(muap_unique,2)+1):(size(muap_unique,2)+size(tmptmp{size(mvc_levels1,2)},2)))=doublet_all{size(mvc_levels1,2)}(tmptmp{size(mvc_levels1,2)});
        rt_unique((size(muap_unique,2)+1):(size(muap_unique,2)+size(tmptmp{size(mvc_levels1,2)},2)))=RT_all{size(mvc_levels1,2)}(tmptmp{size(mvc_levels1,2)});

        if ~isempty(tmptmp{size(mvc_levels1,2)})
            muap_unique((size(muap_unique,2)+1):(size(muap_unique,2)+size(tmptmp{size(mvc_levels1,2)},2)))=muap_all{size(mvc_levels1,2)}(tmptmp{size(mvc_levels1,2)});
        end

        disp('Saving...');
        if subj_nr<10
            tmp=['/Volumes/Expansion/Reproducing doublet paper figures/data/sub-0',num2str(subj_nr),'/emg/sub-0',num2str(subj_nr),'_task-trapezoid_emg_trapz_unique.mat'];
        else
            tmp=['/Volumes/Expansion/Reproducing doublet paper figures/data/sub-',num2str(subj_nr),'/emg/sub-',num2str(subj_nr),'_task-trapezoid_emg_trapz_unique.mat'];
        end
        save(tmp,'muap_unique','doublet_unique','rt_unique','doublet_force_rt','singlet_force_rt');

        % Experiment 2
        % Loop through every subject
    else
        muap_all=cell(1,size(mvc_levels2,2));
        RT_all=cell(1,size(mvc_levels2,2));
        doublet_all=cell(1,size(mvc_levels2,2));
        max_ind=[2 5 10 20];
        % Loop through every trapezoid recording
        for file_num=1:size(mvc_levels2,2)
            tmp=['/Volumes/Expansion/Reproducing doublet paper figures/data/sub-',num2str(subj_nr),'/emg/sub-',num2str(subj_nr),'_task-trapezoid',mvc_levels2{file_num},'percentmvc_run-01_emg.mat'];
            load(tmp)

            fs=2042.5; % actual sampling rate is 2042.5, although officially 2048 Hz
            [b,a]=butter(3,10/(fs/2),'low');
            win=-50:50; % +/- 50 samples (~ +/- 25 ms) for triggering
            t = 1e3*[0:length(win)-1]/fs; % time vector for figures

            % Pre-define muap cell
            muap=cell(1,size(edition.Distimeclean{1},2));
            RT=zeros(1,size(edition.Distimeclean{1},2));
            doublet=zeros(1,size(edition.Distimeclean{1},2));

            % Filter raw EMG signals (50 Hz notch and 20-500 Hz bandpass)
            sig=notchsignals(signal.data(1:64,:),signal.fsamp);
            sig=bandpassingals(sig,signal.fsamp);

            %signal.path=signal.ref_signal;
            force=filtfilt(b,a,signal.target);
            force=force-force(1);
            force=max_ind(file_num)*force./mean(force(round(length(force)/2)+(-20000:20000)));
            force(force<0)=0;

            iter=0;
            % Loop through each MU to extract its MUAP
            for MUnum=1:size(muap,2)
                iter=iter+1;

                % Extract the discharges used for triggering
                locs=edition.Distimeclean{1}{MUnum};

                % Estimate RT
                try
                    for i=0:(length(locs)-num_spikes_rt)
                        locs_tmp=locs((1:num_spikes_rt+1)+i);
                        if mean(round(1e3*diff(locs((1:num_spikes_rt+1)+i))/fs))<=isi_crit
                            break;
                        end
                    end
                catch
                    i=0;
                end

                RT(MUnum)=mean(force(locs((1:num_spikes_rt+1)+i)));

                % Set up the muap matrix (4 grids x 64 channels = 256 channels)
                muap{MUnum}=zeros(64,length(win));

                % Find doublets
                doublet_ind=find(1000./(diff(locs)/2) >= doublet_crit);

                if length(doublet_ind)>=2
                    doublet(MUnum)=1;
                    % AT WHAT FORCE LEVELS DID THE DOUBLETS APPEAR
                    doublet_force_rt=[doublet_force_rt force(locs(doublet_ind))./force(locs(doublet_ind(1)))];
                end

                % Remove doublets in averaging
                if ~isempty(doublet_ind)
                    locs(sort([doublet_ind doublet_ind+1]))=[];
                    if ~isempty(locs)
                        singlet_force_rt=[singlet_force_rt force(locs)./force(locs(1))];
                    end
                end

                % If only doublets
                if isempty(locs) || length(locs)<10
                    locs=edition.Distimeclean{1}{MUnum};
                end

                % Extract the MUAP through STA for each channel
                for ch=1:size(muap{MUnum},1)
                    iterTrig=0;
                    for trig=3:length(locs)-3
                        muap{MUnum}(ch,:)=muap{MUnum}(ch,:)+sig(ch,locs(trig)+win);
                        iterTrig=iterTrig+1;
                    end
                    muap{MUnum}(ch,:)=muap{MUnum}(ch,:)./iterTrig;
                end
            end
            RT_all{file_num}=RT;
            muap_all{file_num}=muap;
            doublet_all{file_num}=doublet;
        end
        disp('Let''s find unique MUs...');
        % FIND UNIQUE MUAPS
        plot_match=0;
        xcorr_cutoff=0.65;

        matched_MUs=cell(1,size(muap_all,2));

        for tmp=1:size(matched_MUs,2)

            % Matched MUs
            matched_MUs{tmp}=zeros(size(muap_all{tmp},2),size(matched_MUs,2));
            matched_MUs{tmp}(:,tmp)=1:size(muap_all{tmp},2);

            for file_num=(tmp+1):size(matched_MUs,2)
                tmpfile=['/Volumes/Expansion/Reproducing doublet paper figures/data/sub-',num2str(subj_nr),'/emg/sub-',num2str(subj_nr),'_task-trapezoid',mvc_levels2{file_num},'percentmvc_run-01_emg.mat'];
                load(tmpfile)

                signal.data(1:64,:)=notchsignals(signal.data(1:64,:),signal.fsamp);
                signal.data(1:64,:)=bandpassingals(signal.data(1:64,:),signal.fsamp);

                R=16;
                eSIG = extension(signal.data(1:64,:),R);
                [wSIG, whitening_matrix] = whitening(eSIG,'ZCA');

                % Loop through each spatio-temporal MUAP from previous recording
                for MU=1:size(muap_all{tmp},2)

                    % Compute MU filter
                    w = muap_all{tmp}{MU};
                    w = extension(w,R);
                    w = whitening_matrix * w;

                    % Reconstruction on current recording
                    sig=w'*wSIG;

                    % Select the source with highest skewness
                    save_skew=zeros(1,size(sig,1));
                    for ind=1:size(sig,1)
                        save_skew(ind)=skewness(sig(ind,:));
                    end
                    [~,maxInd]=max(save_skew);
                    w = w(:,maxInd);
                    w = w./norm(w);

                    % Reconstruction
                    sig=w'*wSIG;

                    [Pulsetrain, Dischargetimes, sil] = calcSIL(wSIG, w, signal.fsamp);

                    r_vec=zeros(1,size(edition.Pulsetrainclean{1},1));
                    for i=1:size(edition.Pulsetrainclean{1},1)
                        [r,lags]=xcorr(Pulsetrain,edition.Pulsetrainclean{1}(i,:),2*R,'normalized');
                        r_vec(i)=max(r);
                    end
                    [~,indx]=find(r_vec==max(r_vec));

                    % Plot
                    if plot_match==1
                        figure(1);
                        title(['Subject ',num2str(subj_nr),' file num ',num2str(file_num), ' (MU ',num2str(MU),')']);
                        subplot(2,1,1);
                        plot(r_vec,'-o');ylim([0 1]);hold on;plot(xlim,[xcorr_cutoff xcorr_cutoff],'k:');hold off;
                        subplot(2,1,2);
                        plot(Pulsetrain);hold on;plot(edition.Pulsetrainclean{1}(indx,:));hold off;
                        pause;
                    end

                    % If above a certain threshold, keep it, otherwise []
                    if max(r_vec) >= xcorr_cutoff
                        matched_MUs{tmp}(MU,file_num)=indx;
                    end
                end
            end
        end

        % Find unique MUs
        tmptmp{1}=1:size(muap_all{1},2);
        tmptmp{2}=1:size(muap_all{2},2);
        tmptmp{3}=1:size(muap_all{3},2);
        tmptmp{4}=1:size(muap_all{4},2);

        cntr=0;
        for i=1:size(matched_MUs,2)
            for j=tmptmp{i}

                % count number of unique muaps
                cntr=cntr+1;

                % find matched mus
                if ~isempty(matched_MUs{i})
                    mu_nums=matched_MUs{i}(j,:);

                    % ignore zeros
                    tmp_ind=find(mu_nums~=0);

                    % take muap from latest recording having a match
                    selected_muap_rec=tmp_ind(end);

                    % save unique muap
                    muap_unique{cntr}=muap_all{selected_muap_rec}{mu_nums(selected_muap_rec)};

                    doublet_tmp=[];
                    rt_tmp=[];
                    for k=tmp_ind
                        % remove indexes for that unique muap
                        tmptmp{k}(find(tmptmp{k}==mu_nums(k)))=[];
                        doublet_tmp=[doublet_tmp doublet_all{k}(mu_nums(k))];
                        rt_tmp=[rt_tmp RT_all{k}(mu_nums(k))];
                    end

                    % save if doublet
                    doublet_unique(cntr)=any(doublet_tmp);

                    % save mean rt
                    rt_unique(cntr)=mean(rt_tmp);
                end
            end
        end

        doublet_unique((size(muap_unique,2)+1):(size(muap_unique,2)+size(tmptmp{size(mvc_levels2,2)},2)))=doublet_all{size(mvc_levels2,2)}(tmptmp{size(mvc_levels2,2)});
        rt_unique((size(muap_unique,2)+1):(size(muap_unique,2)+size(tmptmp{size(mvc_levels2,2)},2)))=RT_all{4}(tmptmp{size(mvc_levels2,2)});

        if ~isempty(tmptmp{size(mvc_levels2,2)})
            muap_unique((size(muap_unique,2)+1):(size(muap_unique,2)+size(tmptmp{size(mvc_levels2,2)},2)))=muap_all{size(mvc_levels2,2)}(tmptmp{size(mvc_levels2,2)});
        end

        disp('Saving...');
        tmp=['/Volumes/Expansion/Reproducing doublet paper figures/data/sub-',num2str(subj_nr),'/emg/sub-',num2str(subj_nr),'_task-trapezoid_emg_trapz_unique.mat'];
        save(tmp,'muap_unique','doublet_unique','rt_unique','doublet_force_rt','singlet_force_rt');
    end
end