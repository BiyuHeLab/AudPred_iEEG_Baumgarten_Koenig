function TimelockData = NASTD_ECoG_TLA_CompTLdata_perTDPredp34...
    (preprocData_perTD)

TimelockData = [];
%1. Select trials per Predp34

label_Predp34 = [-1 0 1];% low, medium, high;
for i_Predp34 = 1:length(label_Predp34)
    
    %Determine current Predp34 trials
    TrialFilt_Predp34 = preprocData_perTD.behav.stim.predID == label_Predp34(i_Predp34); %filter to select trials for certain tone dur
    Index_Predp34 = find(TrialFilt_Predp34 == 1);
    
    %Determine min and amx number of sampels across trials
    Num_samples = [];
    for i_trial = 1:length(preprocData_perTD.trial)   
        Num_samples = [Num_samples; length(preprocData_perTD.trial{i_trial})];
    end
    max_samples = max(Num_samples);
    min_samples = min(Num_samples);
    
%2. Compute timelocked data across TD+Predp34 trials for each selected channel

    %2.1 For low freq
    proxy_TLA = preprocData_perTD;
    proxy_TLA.trial = preprocData_perTD.trial;
    
    cfg             = [];
    cfg.channel     = 'all';
    cfg.trials      = Index_Predp34;
    cfg.latency     = [proxy_TLA.time{1}(1) proxy_TLA.time{1}(min_samples)]; %restrict to TW present in all trials
    cfg.keeptrials  = 'yes';
    TimelockData{i_Predp34} = ft_timelockanalysis(cfg, proxy_TLA);  
    TimelockData{i_Predp34}.Trial_LF = TimelockData{i_Predp34}.trial;  
    rmfield(TimelockData{i_Predp34},'trial');
%     cfg.keeptrials  = 'no';
%     temp = ft_timelockanalysis(cfg, proxy_TLA);    
%     TimelockData{i_Predp34}.avg = temp.trial;
%     clear temp

%     %2.2 For normalized GammaAmp
%     proxy_TLA = preprocData_perTD;
%     proxy_TLA.trial = preprocData_perTD.GammaAmp;
%     
%     cfg             = [];
%     cfg.channel     = 'all';
%     cfg.trials      = Index_Predp34;
%     cfg.latency     = [proxy_TLA.time{1}(1) proxy_TLA.time{1}(min_samples)]; %restrict to TW present in all trials
%     cfg.keeptrials  = 'yes';
%     temp = ft_timelockanalysis(cfg, proxy_TLA);  
%     TimelockData{i_Predp34}.Trial_GammaAmp = temp.trial;
    
%     cfg.keeptrials  = 'no';
%     temp = ft_timelockanalysis(cfg, proxy_TLA);    
%     TimelockData{i_Predp34}.GammaAmp_avg = temp.trial;
     clear temp

    %2.3 For normalized LogGammaAmp
    proxy_TLA = preprocData_perTD;
    proxy_TLA.trial = preprocData_perTD.LogGammaAmp;
    
    cfg             = [];
    cfg.channel     = 'all';
    cfg.trials      = Index_Predp34;
    cfg.latency     = [proxy_TLA.time{1}(1) proxy_TLA.time{1}(min_samples)]; %restrict to TW present in all trials
    cfg.keeptrials  = 'yes';
    temp = ft_timelockanalysis(cfg, proxy_TLA);  
    TimelockData{i_Predp34}.Trial_LogGammaAmp = temp.trial;
    
%     cfg.keeptrials  = 'no';
%     temp = ft_timelockanalysis(cfg, proxy_TLA);    
%     TimelockData{i_Predp34}.LogGammaAmp_avg = temp.trial;
     clear temp


    %3. Manually compute Avg and STD across TD+Predp34 trials for each selected channel
%(for plotting)
    %3.1 For low freq
    proxy_SelTrialsperChan = NaN(length(Index_Predp34),min_samples);
    for i_chan = 1:length(preprocData_perTD.label)
        for i_trial = 1:length(Index_Predp34) 
            proxy_SelTrialsperChan(i_trial,1:min_samples) = ...
                preprocData_perTD.trial{Index_Predp34(i_trial)}(i_chan,1:min_samples);
        end  
        TimelockData{i_Predp34}.Avg_LF(i_chan,:) = mean(proxy_SelTrialsperChan);
        TimelockData{i_Predp34}.STD_LF(i_chan,:) = std(proxy_SelTrialsperChan);        
    end

%     %3.2 For normalized GammaAmp
%     proxy_SelTrialsperChan = NaN(length(Index_Predp34),min_samples);
%     for i_chan = 1:length(preprocData_perTD.label)
%         for i_trial = 1:length(Index_Predp34) 
%             proxy_SelTrialsperChan(i_trial,1:min_samples) = ...
%                 preprocData_perTD.GammaAmp{Index_Predp34(i_trial)}(i_chan,1:min_samples);
%         end  
%         TimelockData{i_Predp34}.Avg_GammaAmp(i_chan,:) = mean(proxy_SelTrialsperChan);
%         TimelockData{i_Predp34}.STD_GammaAmp(i_chan,:) = std(proxy_SelTrialsperChan);        
%     end   

    %3.3 For normalized LogGammaAmp
    proxy_SelTrialsperChan = NaN(length(Index_Predp34),min_samples);
    for i_chan = 1:length(preprocData_perTD.label)
        for i_trial = 1:length(Index_Predp34) 
            proxy_SelTrialsperChan(i_trial,1:min_samples) = ...
                preprocData_perTD.LogGammaAmp{Index_Predp34(i_trial)}(i_chan,1:min_samples);
        end  
        TimelockData{i_Predp34}.Avg_LogGammaAmp(i_chan,:) = mean(proxy_SelTrialsperChan);
        TimelockData{i_Predp34}.STD_LogGammaAmp(i_chan,:) = std(proxy_SelTrialsperChan);        
    end     
    
end

end