function NASTD_ECoG_TLA_PlotTrialAvg_BothTD_perChan...
    (sub, inputData, ChanIndex, ...
    preprocData_perTD, ...
    plot_STD, save_fig, paths_NASTD_ECoG)

% sub = 'NY688';
% inputData = 'LF';
% ChanIndex = 10; 

%% 1. Timelock across trials
TimelockData = [];
%1. Select trials per TD
ToneDur_text = {'0.2' '0.4'};

for i_TD = 1:length(ToneDur_text)   
    
    %Determine min and amx number of sampels across trials
    Num_samples = [];
    for i_trial = 1:length(preprocData_perTD{i_TD}.trial)   
        Num_samples = [Num_samples; length(preprocData_perTD{i_TD}.trial{i_trial})];
    end
    max_samples = max(Num_samples);
    min_samples = min(Num_samples);
    
%1.2. Compute timelocked data across TD+Predp34 trials for each selected channel

    %1.2.1 For low freq
    proxy_TLA = preprocData_perTD{i_TD};
    proxy_TLA.trial = preprocData_perTD{i_TD}.trial;
    
    cfg             = [];
    cfg.channel     = preprocData_perTD{i_TD}.label(ChanIndex);
    cfg.trials      = 'all';
    cfg.latency     = [proxy_TLA.time{1}(1) proxy_TLA.time{1}(min_samples)]; %restrict to TW present in all trials
    cfg.keeptrials  = 'yes';
    TimelockData{i_TD} = ft_timelockanalysis(cfg, proxy_TLA);  
    TimelockData{i_TD}.Trial_LF = TimelockData{i_TD}.trial;  
    rmfield(TimelockData{i_TD},'trial');

    %1.2.2 For normalized GammaAmp
    proxy_TLA = preprocData_perTD{i_TD};
    proxy_TLA.trial = preprocData_perTD{i_TD}.GammaAmp_Norm;
    
    cfg             = [];
    cfg.channel     = preprocData_perTD{i_TD}.label(ChanIndex);
    cfg.trials      = 'all';
    cfg.latency     = [proxy_TLA.time{1}(1) proxy_TLA.time{1}(min_samples)]; %restrict to TW present in all trials
    cfg.keeptrials  = 'yes';
    temp = ft_timelockanalysis(cfg, proxy_TLA);  
    TimelockData{i_TD}.Trial_GammaAmp = temp.trial;    
    clear temp

    %1.2.3 For normalized LogGammaAmp
    proxy_TLA = preprocData_perTD{i_TD};
    proxy_TLA.trial = preprocData_perTD{i_TD}.LogGammaAmp_Norm;
    
    cfg             = [];
    cfg.channel     = preprocData_perTD{i_TD}.label(ChanIndex);
    cfg.trials      = 'all';
    cfg.latency     = [proxy_TLA.time{1}(1) proxy_TLA.time{1}(min_samples)]; %restrict to TW present in all trials
    cfg.keeptrials  = 'yes';
    temp = ft_timelockanalysis(cfg, proxy_TLA);  
    TimelockData{i_TD}.Trial_LogGammaAmp = temp.trial;
    clear temp


    %1. 3. Manually compute Avg and STD across TD+Predp34 trials for each selected channel (for plotting)
    %1.3.1 For low freq
    proxy_SelTrialsperChan = NaN(1,min_samples);
        for i_trial = 1:length(preprocData_perTD{i_TD}.trial) 
            proxy_SelTrialsperChan(i_trial,1:min_samples) = ...
                preprocData_perTD{i_TD}.trial{i_trial}(ChanIndex,1:min_samples);
        end  
        TimelockData{i_TD}.Avg_LF(1,:) = mean(proxy_SelTrialsperChan);
        TimelockData{i_TD}.STD_LF(1,:) = std(proxy_SelTrialsperChan);        

    %3.2 For normalized GammaAmp
    proxy_SelTrialsperChan = NaN(1,min_samples);
        for i_trial = 1:length(preprocData_perTD{i_TD}.trial) 
            proxy_SelTrialsperChan(i_trial,1:min_samples) = ...
                preprocData_perTD{i_TD}.GammaAmp_Norm{i_trial}(ChanIndex,1:min_samples);
        end  
        TimelockData{i_TD}.Avg_GammaAmp(1,:) = mean(proxy_SelTrialsperChan);
        TimelockData{i_TD}.STD_GammaAmp(1,:) = std(proxy_SelTrialsperChan);        
       

    %3.3 For normalized LogGammaAmp
     proxy_SelTrialsperChan = NaN(1,min_samples);
       for i_trial = 1:length(preprocData_perTD{i_TD}.trial) 
            proxy_SelTrialsperChan(i_trial,1:min_samples) = ...
                preprocData_perTD{i_TD}.LogGammaAmp_Norm{i_trial}(ChanIndex,1:min_samples);
        end  
        TimelockData{i_TD}.Avg_LogGammaAmp(1,:) = mean(proxy_SelTrialsperChan);
        TimelockData{i_TD}.STD_LogGammaAmp(1,:) = std(proxy_SelTrialsperChan);        
        
    
end

%2 Plot
%2.1 Set up figure and plotting param
h = figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1])
hold on
subplot_position = 0;

plot_param.color = {[0, 0.28, 0.73],[1, 0.75, 0],[0.69, 0, 0.16]};
plot_param.linestyle = {'-','-','-'};
plot_param.lineThicknessToneMark = 0.5;
plot_param.Labels_ToneStart = {'T1','T2','T3','T4','T5','T6','T7','T8','T9','T10',...
    'T11','T12','T13','T14','T15','T16','T17','T18','T19','T20',...
    'T21','T22','T23','T24','T25','T26','T27','T28','T29','T30',...
    'T31','T32','T33(STABLE)','T34(RATED)','','RespDisp'};

ToneDur_text = {'0.2' '0.4'};
label_Predp34 = [-1 0 1];% low, medium, high;
  
%1.2 relabel input data
if strcmp(inputData,'LF')
    inputData_Avg = 'Avg_LF';
    inputData_STD = 'STD_LF';    
elseif strcmp(inputData,'GammaAmp')
    inputData_Avg = 'Avg_GammaAmp';
    inputData_STD = 'STD_GammaAmp';
elseif strcmp(inputData,'LogGammaAmp')
    inputData_Avg = 'Avg_LogGammaAmp';
    inputData_STD = 'STD_LogGammaAmp';
end

%2. Read out TP and samples corresponding to tone presentation and response window 
for i_tonedur = 1:length(preprocData_perTD)
    TP_StartTone{i_tonedur} = NaN(1,35);
    TP_StartTone{i_tonedur}(1) = 0;
    for i_tone = 1:35
        Dist = abs(preprocData_perTD{i_tonedur}.time{1} - ((str2num(ToneDur_text{i_tonedur})*i_tone)-str2num(ToneDur_text{i_tonedur})));
        minDist = min(Dist);
        i_minDist = find(Dist == minDist);
        TP_StartTone{i_tonedur}(i_tone) = preprocData_perTD{i_tonedur}.time{1}(i_minDist);
    end
    
    Dist = abs(preprocData_perTD{i_tonedur}.time{1} - (TP_StartTone{i_tonedur}(end)+0.4));%0.4 s until resp window appeared
    minDist = min(Dist);
    i_minDist = find(Dist == minDist);
    TP_StartTone{i_tonedur}(36) = preprocData_perTD{i_tonedur}.time{1}(i_minDist); 
    
    for i_tone = 1:length(TP_StartTone{i_tonedur})
        Sample_StartTone{i_tonedur}(i_tone) = find(preprocData_perTD{i_tonedur}.time{1} == TP_StartTone{i_tonedur}(i_tone));
    end       
end   

%2. Compute min/max vals across conditions and TD for common y-axis    
minAmp = 1000; %set high, so overwritten
maxAmp = -1000;
for i_tonedur = 1:length(TimelockData)
        minAmp_temp = min(TimelockData{i_tonedur}.(inputData_Avg)(1,Sample_StartTone{i_tonedur}(1):Sample_StartTone{i_tonedur}(end)));
        maxAmp_temp = max(TimelockData{i_tonedur}.(inputData_Avg)(1,Sample_StartTone{i_tonedur}(1):Sample_StartTone{i_tonedur}(end)));
        
        if minAmp_temp < minAmp
            minAmp = minAmp_temp;
        end
        if maxAmp_temp > maxAmp
            maxAmp = maxAmp_temp;
        end    
      
end

%3. Plot all Predp34 per TD in subpot and 1 subplot per TD
for i_tonedur = 1:length(TimelockData)

    subplot_position = subplot_position + 1;
    subplot(2,1,subplot_position);
    hold on;
    
       
        %GAvg per Predp34
        plot(1:length(TimelockData{i_tonedur}.(inputData_Avg)(1,:)),...
            TimelockData{i_tonedur}.(inputData_Avg)(1,:),...
            'color','b','LineWidth',2)
        ylim([minAmp maxAmp])

        %STD per Predp34
        if plot_STD == 1 
            
            shadedErrorBar(1:length(TimelockData{i_tonedur}.(inputData_Avg)(1,:)),...
                TimelockData{i_tonedur}.(inputData_Avg)(1,:), ...
                TimelockData{i_tonedur}.(inputData_STD)(1,:),...
                'lineprops',{'color','b'})
            
        ylim([minAmp*5 maxAmp*5])
            
        end
        
        hold on;
        
    
    %Vertical line for each tone/event
    %Create aray that has max/min deflection at sample of each tone/event
    Array_Tone_Pos = zeros(1,length(TimelockData{i_tonedur}.time));
    Array_Tone_Neg = zeros(1,length(TimelockData{i_tonedur}.time));
    
    for i_SampleToneStart = Sample_StartTone{i_tonedur}
        Array_Tone_Pos(i_SampleToneStart) = maxAmp;
        Array_Tone_Neg(i_SampleToneStart) = minAmp;        
    end
    
    plot(1:length(Array_Tone_Pos),...
        Array_Tone_Pos,...
        'k','LineWidth',plot_param.lineThicknessToneMark) %plot mark for each tone-start-sample
    plot(1:length(Array_Tone_Neg),...
        Array_Tone_Neg,...
        'k','LineWidth',plot_param.lineThicknessToneMark) %plot mark for each tone-start-sample
        
    xlim([1 length(TimelockData{i_tonedur}.(inputData_Avg)(1,:))])  

    %x-axis shows tones [tone number]
    set(gca,'FontSize',10,'XTickLabelRotation',45)
    set(gca,'xtick',Sample_StartTone{i_tonedur})
    x_label = plot_param.Labels_ToneStart;
    set(gca,'FontSize',10)
    xticklabels(x_label);
    
    
    title({['ToneDur: ' ToneDur_text{i_tonedur} ' s']});
end

    suptitle({[sub '; Chan: ' preprocData_perTD{i_tonedur}.label{ChanIndex} '; ' inputData],...
        ['ERP per TD(GAvg + STD across trials (n = ' num2str(size(preprocData_perTD{i_tonedur}.trial,1)) ')) timelocked to tone presentation']});

if save_fig == 1
    path_fig = ([paths_NASTD_ECoG.Fig_Timelocked_SingleSubs sub '/TD/' inputData '/']);
    mkdir(path_fig);
    filename     = [preprocData_perTD{i_tonedur}.label{ChanIndex} '_' inputData '_' sub '_TLdata.png'];
    figfile      = [path_fig filename];
    saveas(gcf, [figfile], 'png'); %save png version
    close
end


end
    
