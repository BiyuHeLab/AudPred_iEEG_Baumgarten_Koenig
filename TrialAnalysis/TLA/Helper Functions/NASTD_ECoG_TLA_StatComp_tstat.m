function TimelockStat_Predp34 = NASTD_ECoG_TLA_StatComp_tstat...
    (sub, i_chan, inputData, ComparedConditions, pval_Predp34comp, ...
    test_TimelockData_perTDPredp34, ...
    plot_fig, save_fig, paths_NASTD_ECoG)
%Aim: Statistical comparison of electrode-wise timelocked neural data.
%Contrast: different p*34 (i.e., prediction effect)
%TOI: each sample from T1:T34
%Statistical: t-test (independent)


%1.1 Label to-be-compared conditions
for i_conditions = 1:length(ComparedConditions)
    if ComparedConditions(i_conditions) == 1
        tempLabel_ComparedConditions{i_conditions} = 'lowPredp34 ';
    elseif ComparedConditions(i_conditions) == 2
        tempLabel_ComparedConditions{i_conditions} = 'medPredp34 ';
    elseif ComparedConditions(i_conditions) == 3
        tempLabel_ComparedConditions{i_conditions} = 'highPredp34 ';
    end
end
Label_ComparedConditions = strcat(tempLabel_ComparedConditions{1}, '-vs-', tempLabel_ComparedConditions{2});

%1.2 relabel input data
if strcmp(inputData,'LF')
    inputData_Trial = 'Trial_LF';
    inputData_Avg = 'Avg_LF';
    inputData_STD = 'STD_LF';
elseif strcmp(inputData,'GammaAmp')
    inputData_Trial = 'Trial_GammaAmp';
    inputData_Avg = 'Avg_GammaAmp';
    inputData_STD = 'STD_GammaAmp';
elseif strcmp(inputData,'LogGammaAmp')
    inputData_Trial = 'Trial_LogGammaAmp';
    inputData_Avg = 'Avg_LogGammaAmp';
    inputData_STD = 'STD_LogGammaAmp';
end

%2. Determine TP/samples for each tone start+end to define TOI
if length(test_TimelockData_perTDPredp34{1}.time) < 7000
    ToneDur_text = '0.2';
elseif length(test_TimelockData_perTDPredp34{1}.time) > 7000
    ToneDur_text = '0.4';
end

TP_Tone_StartStop = NaN(35,2);
TP_Tone_StartStop(1,1) = 0;
for i_tone = 2:35
    Dist = abs(test_TimelockData_perTDPredp34{1}.time - ((str2num(ToneDur_text)*i_tone)-str2num(ToneDur_text)));
    minDist = min(Dist);
    i_minDist = find(Dist == minDist);
    TP_Tone_StartStop(i_tone,1) = test_TimelockData_perTDPredp34{1}.time(i_minDist);
end
for i_tone = 1:34
    i_LastSampleTone = find(test_TimelockData_perTDPredp34{1}.time == TP_Tone_StartStop(i_tone+1,1));
    TP_Tone_StartStop(i_tone,2) = test_TimelockData_perTDPredp34{1}.time(i_LastSampleTone);
end
TP_Tone_StartStop(35,1) = TP_Tone_StartStop(34,2)+0.4;
TP_Tone_StartStop(35,2) = NaN;

%check if all tones are of equal length, if not then choose min length
%by deleting the last sample
minSeqLength_sec = min(TP_Tone_StartStop(1:34,2)-TP_Tone_StartStop(1:34,1));
for i_tone = 1:34
    if TP_Tone_StartStop(i_tone,2)-TP_Tone_StartStop(i_tone,1) > minSeqLength_sec
        TP_Tone_StartStop(i_tone,2) = ...
            test_TimelockData_perTDPredp34{1}.time(find(test_TimelockData_perTDPredp34{1}.time == ...
            TP_Tone_StartStop(i_tone,2))-1);
    end
end

Sample_Tone_StartStop = NaN(35,2);
for i_tone = 1:34
    Sample_Tone_StartStop(i_tone,1) = ...
        find(test_TimelockData_perTDPredp34{1}.time == TP_Tone_StartStop(i_tone,1));
    Sample_Tone_StartStop(i_tone,2) = ...
        find(TP_Tone_StartStop(i_tone,2) == test_TimelockData_perTDPredp34{1}.time);
end
    Sample_Tone_StartStop(35,1) = ...
        nearest(test_TimelockData_perTDPredp34{1}.time,TP_Tone_StartStop(35,1));


%3. Performt-test between two selected conditions with FT-timelockstat
cfg                     = [];

cfg.channel             = i_chan; %only for current channel
cfg.latency             = [TP_Tone_StartStop(1,1) TP_Tone_StartStop(35,1)]; %during sequence presentation
cfg.avgovertime         = 'no';
cfg.parameter           = inputData_Trial;

cfg.method              = 'montecarlo';
cfg.alpha               = pval_Predp34comp;
cfg.tail                = 0;        % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.correcttail         = 'prob';   %correct prob.subfield for 2sided test

cfg.correctm            = 'cluster';
cfg.clusteralpha        = 0.05;
cfg.clusterstatistic    = 'maxsum';
cfg.clustertail         = 0;

cfg.numrandomization    = 1000;

Num_trials              = size(test_TimelockData_perTDPredp34{1}.trial,1) + size(test_TimelockData_perTDPredp34{3}.trial,1);
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable (i.e., condition)

cfg.design(1,1:Num_trials)  = [ones(1,size(test_TimelockData_perTDPredp34{1}.trial,1)) 2*ones(1,size(test_TimelockData_perTDPredp34{3}.trial,1))];

cfg.statistic           = 'ft_statfun_indepsamplesT';

TimelockStat_Predp34 = ft_timelockstatistics(cfg, ...
    test_TimelockData_perTDPredp34{ComparedConditions(1)}, ...
    test_TimelockData_perTDPredp34{ComparedConditions(2)});

%3.3 Add subfield marking sign samples for plotting
TimelockStat_Predp34.plot_statp = ...
    zeros(1,size(test_TimelockData_perTDPredp34{1}.trial,3));

for i_sample = 1:length(TimelockStat_Predp34.time)
    %Find where sample from stat comparison is in TLA file
    stat_sample = TimelockStat_Predp34.time(i_sample) == ...
        test_TimelockData_perTDPredp34{1}.time;
    %
    if TimelockStat_Predp34.prob(i_sample) < cfg.alpha
        TimelockStat_Predp34.plot_statp(stat_sample) = 1;
    end
end

TimelockStat_Predp34.comparison = Label_ComparedConditions;
TimelockStat_Predp34.input = inputData;

%2.4 Display output
if sum(TimelockStat_Predp34.plot_statp) > 0
    Num_signSamples = length(find(TimelockStat_Predp34.plot_statp ~= 0));
    disp([num2str(Num_signSamples) ' sign. different samples found for electrode ' ...
        test_TimelockData_perTDPredp34{1}.label{i_chan} ' (' num2str(i_chan) ')']);
    TimelockStat_Predp34.SignDiff = 1;
    
    %Optional: Plot timelocked activity and sign. samples
    if plot_fig == 1
        max_amp1 = max(test_TimelockData_perTDPredp34{ComparedConditions(1)}.(inputData_Avg)(i_chan,:));
        max_amp2 = max(test_TimelockData_perTDPredp34{ComparedConditions(2)}.(inputData_Avg)(i_chan,:));
        max_amp = max([max_amp1 max_amp2]);
        min_amp1 = min(test_TimelockData_perTDPredp34{ComparedConditions(1)}.(inputData_Avg)(i_chan,:));
        min_amp2 = min(test_TimelockData_perTDPredp34{ComparedConditions(2)}.(inputData_Avg)(i_chan,:));
        min_amp = min([min_amp1 min_amp2]);
        
        figure;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        
        %         shadedErrorBar(1:length(test_TimelockData_perTDPredp34{ComparedConditions(1)}.(inputData_Avg)),...
        %             test_TimelockData_perTDPredp34{ComparedConditions(1)}.(inputData_Avg)(i_chan,:), ...
        %             test_TimelockData_perTDPredp34{ComparedConditions(1)}.(inputData_STD)(i_chan,:),...
        %             'lineprops',{'color','b'});
        %         hold on;
        plot(1:length(test_TimelockData_perTDPredp34{ComparedConditions(1)}.(inputData_Avg)), ...
            test_TimelockData_perTDPredp34{ComparedConditions(1)}.(inputData_Avg)(i_chan,:),'b','LineWidth',2)
        hold on;
        
        %         shadedErrorBar(1:length(test_TimelockData_perTDPredp34{ComparedConditions(2)}.(inputData_Avg)),...
        %             test_TimelockData_perTDPredp34{ComparedConditions(2)}.(inputData_Avg)(i_chan,:), ...
        %             test_TimelockData_perTDPredp34{ComparedConditions(2)}.(inputData_STD)(i_chan,:),...
        %             'lineprops',{'color','r'});
        %         hold on;
        plot(1:length(test_TimelockData_perTDPredp34{ComparedConditions(2)}.(inputData_Avg)), ...
            test_TimelockData_perTDPredp34{ComparedConditions(2)}.(inputData_Avg)(i_chan,:),'r','LineWidth',2)
        hold on;
        
        area(1:length(test_TimelockData_perTDPredp34{ComparedConditions(1)}.(inputData_Avg)),TimelockStat_Predp34.plot_statp*max_amp,...
            'basevalue',0,'FaceColor',[0.3, 0.3, 0.3],'FaceAlpha', 0.5,'LineStyle','none');
        hold on;
        area(1:length(test_TimelockData_perTDPredp34{ComparedConditions(1)}.(inputData_Avg)),TimelockStat_Predp34.plot_statp*min_amp,...
            'basevalue',0,'FaceColor',[0.3, 0.3, 0.3],'FaceAlpha', 0.5,'LineStyle','none');
        hold on;
        
        %Vertical line for each tone/event
        Array_Tone_Pos = zeros(1,length(test_TimelockData_perTDPredp34{1}.(inputData_Avg)));
        Array_Tone_Neg = zeros(1,length(test_TimelockData_perTDPredp34{1}.(inputData_Avg)));
        
        for i_SampleToneStart = [Sample_Tone_StartStop(1:34,1)' Sample_Tone_StartStop(34,2) Sample_Tone_StartStop(35,1)]
            Array_Tone_Pos(i_SampleToneStart) = max_amp;
            Array_Tone_Neg(i_SampleToneStart) = min_amp;
        end
        
        plot(1:length(Array_Tone_Pos),...
            Array_Tone_Pos,...
            'k','LineWidth',1) %plot mark for each tone-start-sample
        hold on;
        plot(1:length(Array_Tone_Neg),...
            Array_Tone_Neg,...
            'k','LineWidth',1) %plot mark for each tone-start-sample
        
        
        %x-axis shows tones [tone number]
        xlim([1  length(test_TimelockData_perTDPredp34{ComparedConditions(1)}.(inputData_Avg))])
        set(gca,'FontSize',10,'XTickLabelRotation',45)
        set(gca,'xtick',Sample_Tone_StartStop(:,1)')
        Labels_ToneStart = {'T1','T2','T3','T4','T5','T6','T7','T8','T9','T10',...
        'T11','T12','T13','T14','T15','T16','T17','T18','T19','T20',...
        'T21','T22','T23','T24','T25','T26','T27','T28','T29','T30',...
        'T31','T32','T33(STABLE)','T34(RATED)','', 'RespDisp'};
        x_label = Labels_ToneStart;
        set(gca,'FontSize',10)
        xticklabels(x_label);
        
        legend(tempLabel_ComparedConditions{1}, tempLabel_ComparedConditions{2},'Location','NorthWest')        
        
        title({[sub '; Chan: ' test_TimelockData_perTDPredp34{ComparedConditions(1)}.label{i_chan} ';ToneDur: ' ToneDur_text 's'], ...
            [inputData '; t-test (indep): ' Label_ComparedConditions]});
        
        if save_fig == 1
            path_fig = ([paths_NASTD_ECoG.Fig_Timelocked_SingleSubs sub '/Stat/tstat/' inputData '/']);
            mkdir(path_fig);
            filename     = [test_TimelockData_perTDPredp34{1}.label{i_chan} '_' inputData '_'  ToneDur_text 'sTD_' sub '_TLtstat.png'];
            figfile      = [path_fig filename];
            saveas(gcf, [figfile], 'png'); %save png version
            close
        end
        
        
    end
else
    disp(['No sign. different samples found for electrode ' ...
        test_TimelockData_perTDPredp34{1}.label{i_chan} ' (' num2str(i_chan) ')']);
    TimelockStat_Predp34.SignDiff = 0;
end

end