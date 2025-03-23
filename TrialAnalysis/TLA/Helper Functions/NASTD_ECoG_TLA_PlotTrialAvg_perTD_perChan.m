function NASTD_ECoG_TLA_PlotTrialAvg_perTD_perChan...
    (sub, inputData, TD_text, ChanIndex, ...
    input_preprocData_perTD, ...
    plot_STD, save_fig, paths_NASTD_ECoG)


%Determine min and amx number of sampels across trials
Num_samples = [];
for i_trial = 1:length(input_preprocData_perTD.trial)
    Num_samples = [Num_samples; length(input_preprocData_perTD.trial{i_trial})];
end
max_samples = max(Num_samples);
min_samples = min(Num_samples);

%1.2. Compute timelocked data across TD+Predp34 trials for each selected channel

%1.2.1 For low freq
proxy_TLA = input_preprocData_perTD;
proxy_TLA.trial = input_preprocData_perTD.trial;

cfg             = [];
cfg.channel     = input_preprocData_perTD.label(ChanIndex);
cfg.trials      = 'all';
cfg.latency     = [proxy_TLA.time{1}(1) proxy_TLA.time{1}(min_samples)]; %restrict to TW present in all trials
cfg.keeptrials  = 'yes';
TimelockData = ft_timelockanalysis(cfg, proxy_TLA);
TimelockData.Trial_LF = TimelockData.trial;
rmfield(TimelockData,'trial');

%1.2.2 For normalized GammaAmp
proxy_TLA = input_preprocData_perTD;
proxy_TLA.trial = input_preprocData_perTD.GammaAmp_Norm;

cfg             = [];
cfg.channel     = input_preprocData_perTD.label(ChanIndex);
cfg.trials      = 'all';
cfg.latency     = [proxy_TLA.time{1}(1) proxy_TLA.time{1}(min_samples)]; %restrict to TW present in all trials
cfg.keeptrials  = 'yes';
temp = ft_timelockanalysis(cfg, proxy_TLA);
TimelockData.Trial_GammaAmp = temp.trial;
clear temp

%1.2.3 For normalized LogGammaAmp
proxy_TLA = input_preprocData_perTD;
proxy_TLA.trial = input_preprocData_perTD.LogGammaAmp_Norm;

cfg             = [];
cfg.channel     = input_preprocData_perTD.label(ChanIndex);
cfg.trials      = 'all';
cfg.latency     = [proxy_TLA.time{1}(1) proxy_TLA.time{1}(min_samples)]; %restrict to TW present in all trials
cfg.keeptrials  = 'yes';
temp = ft_timelockanalysis(cfg, proxy_TLA);
TimelockData.Trial_LogGammaAmp = temp.trial;
clear temp


%1. 3. Manually compute Avg and STD across TD+Predp34 trials for each selected channel (for plotting)
%1.3.1 For low freq
proxy_SelTrialsperChan = NaN(1,min_samples);
for i_trial = 1:length(input_preprocData_perTD.trial)
    proxy_SelTrialsperChan(i_trial,1:min_samples) = ...
        input_preprocData_perTD.trial{i_trial}(ChanIndex,1:min_samples);
end
TimelockData.Avg_LF(1,:) = mean(proxy_SelTrialsperChan);
TimelockData.STD_LF(1,:) = std(proxy_SelTrialsperChan);

%3.2 For normalized GammaAmp
proxy_SelTrialsperChan = NaN(1,min_samples);
for i_trial = 1:length(input_preprocData_perTD.trial)
    proxy_SelTrialsperChan(i_trial,1:min_samples) = ...
        input_preprocData_perTD.GammaAmp_Norm{i_trial}(ChanIndex,1:min_samples);
end
TimelockData.Avg_GammaAmp(1,:) = mean(proxy_SelTrialsperChan);
TimelockData.STD_GammaAmp(1,:) = std(proxy_SelTrialsperChan);


%3.3 For normalized LogGammaAmp
proxy_SelTrialsperChan = NaN(1,min_samples);
for i_trial = 1:length(input_preprocData_perTD.trial)
    proxy_SelTrialsperChan(i_trial,1:min_samples) = ...
        input_preprocData_perTD.LogGammaAmp_Norm{i_trial}(ChanIndex,1:min_samples);
end
TimelockData.Avg_LogGammaAmp(1,:) = mean(proxy_SelTrialsperChan);
TimelockData.STD_LogGammaAmp(1,:) = std(proxy_SelTrialsperChan);


%2 Plot
%2.1 Set up figure and plotting param
h = figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1])
hold on
subplot_position = 0;

plot_param.color = {[0, 0.28, 0.73],[1, 0.75, 0],[0.69, 0, 0.16]};
plot_param.linestyle = {'-','-','-'};
plot_param.lineThicknessToneMark = 0.5;
plot_param.Labels_ToneStart = {'T1','','','','','','','','','T10',...
    '','','','','','','','','','T20',...
    '','','','','','','','','','T30',...
    '','','T33','T34','',''};

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
TP_StartTone = NaN(1,35);
TP_StartTone(1) = 0;
for i_tone = 1:35
    Dist = abs(input_preprocData_perTD.time{1} - ((str2num(TD_text)*i_tone)-str2num(TD_text)));
    minDist = min(Dist);
    i_minDist = find(Dist == minDist);
    TP_StartTone(i_tone) = input_preprocData_perTD.time{1}(i_minDist);
end

Dist = abs(input_preprocData_perTD.time{1} - (TP_StartTone(end)+0.4));%0.4 s until resp window appeared
minDist = min(Dist);
i_minDist = find(Dist == minDist);
TP_StartTone(36) = input_preprocData_perTD.time{1}(i_minDist);

for i_tone = 1:length(TP_StartTone)
    Sample_StartTone(i_tone) = find(input_preprocData_perTD.time{1} == TP_StartTone(i_tone));
end

%2. Compute min/max vals across conditions and TD for common y-axis
minAmp = 1000; %set high, so overwritten
maxAmp = -1000;
minAmp_temp = min(TimelockData.(inputData_Avg)(1,Sample_StartTone(1):Sample_StartTone(end)));
maxAmp_temp = max(TimelockData.(inputData_Avg)(1,Sample_StartTone(1):Sample_StartTone(end)));

if minAmp_temp < minAmp
    minAmp = minAmp_temp;
end
if maxAmp_temp > maxAmp
    maxAmp = maxAmp_temp;
end


%3. Plot all Predp34 per TD in subpot and 1 subplot per TD

%GAvg per Predp34
if plot_STD == 0    

    % plot(1:length(TimelockData.(inputData_Avg)(1,:)),...
    %     TimelockData.(inputData_Avg)(1,:),...
    %     'color','b','LineWidth',2)
    plot(1:length(TimelockData.(inputData_Avg)(1,:)),...
        zscore(TimelockData.(inputData_Avg)(1,:)),...
        'color','b','LineWidth',2)
    ylim([-5 5])

%STD per Predp34
elseif plot_STD == 1    
    
    shadedErrorBar(1:length(TimelockData.(inputData_Avg)(1,:)),...
        TimelockData.(inputData_Avg)(1,:), ...
        TimelockData.(inputData_STD)(1,:),...
        'lineprops',{'color','b'})
    
    ylim([minAmp*5 maxAmp*5])
    
end

hold on;

%Vertical line for each tone/event
%Create aray that has max/min deflection at sample of each tone/event
Array_Tone_Pos = zeros(1,length(TimelockData.time));
Array_Tone_Neg = zeros(1,length(TimelockData.time));

for i_SampleToneStart = Sample_StartTone
    Array_Tone_Pos(i_SampleToneStart) = 5;
    Array_Tone_Neg(i_SampleToneStart) = -5;
end

plot(1:length(Array_Tone_Pos),...
    Array_Tone_Pos,...
    'k','LineWidth',plot_param.lineThicknessToneMark) %plot mark for each tone-start-sample
plot(1:length(Array_Tone_Neg),...
    Array_Tone_Neg,...
    'k','LineWidth',plot_param.lineThicknessToneMark) %plot mark for each tone-start-sample

xlim([1 length(TimelockData.(inputData_Avg)(1,:))])


%x-axis shows tones [tone number]
set(gca,'FontSize',10,'XTickLabelRotation',45)
set(gca,'xtick',Sample_StartTone)
x_label = plot_param.Labels_ToneStart;
set(gca,'FontSize',10)
xticklabels(x_label);

% suptitle({[sub '; Chan: ' input_preprocData_perTD.label{ChanIndex} '; ' inputData],['ToneDur: ' TD_text ' s']...
%     ['ERP per TD(GAvg + STD across trials (n = ' num2str(size(input_preprocData_perTD.trial,1)) ')) timelocked to tone presentation']});
 suptitle({[sub ' - ' input_preprocData_perTD.label{ChanIndex} ],['ToneDur: ' TD_text ' s']});


if save_fig == 1
    path_fig = ([paths_NASTD_ECoG.Fig_Timelocked_SingleSubs sub '/TD/' inputData '/']);
    mkdir(path_fig);
    filename     = [sub '_' inputData  input_preprocData_perTD.label{ChanIndex} '_' TD_text 'sTD.png'];
    figfile      = [path_fig filename];
    saveas(gcf, [figfile], 'png'); %save png version
    close
end

end

