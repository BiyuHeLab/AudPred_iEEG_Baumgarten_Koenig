function NASTD_ECoG_TLA_PlotPredp34_perChan(sub, inputData, ChanIndex, TimelockData_perTDPredp34, plot_STD, save_fig, paths_NASTD_ECoG)

%1.1 Set up figure and plotting param
h = figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1])
hold on
subplot_position = 0;

plot_param.legend_label = {'predp34 = low','predp34 =  med','predp34 =  high'};
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
for i_tonedur = 1:length(TimelockData_perTDPredp34)

    TP_StartTone{i_tonedur} = NaN(1,35);
    TP_StartTone{i_tonedur}(1) = 0;
    for i_tone = 1:35
        Dist = abs(TimelockData_perTDPredp34{i_tonedur}{1}.time - ((str2num(ToneDur_text{i_tonedur})*i_tone)-str2num(ToneDur_text{i_tonedur})));
        minDist = min(Dist);
        i_minDist = find(Dist == minDist);
        TP_StartTone{i_tonedur}(i_tone) = TimelockData_perTDPredp34{i_tonedur}{1}.time(i_minDist);
    end
    
    Dist = abs(TimelockData_perTDPredp34{i_tonedur}{1}.time - (TP_StartTone{i_tonedur}(end)+0.4));%0.4 s until resp window appeared
    minDist = min(Dist);
    i_minDist = find(Dist == minDist);
    TP_StartTone{i_tonedur}(36) = TimelockData_perTDPredp34{i_tonedur}{1}.time(i_minDist); 
    
    for i_tone = 1:length(TP_StartTone{i_tonedur})
        Sample_StartTone{i_tonedur}(i_tone) = find(TimelockData_perTDPredp34{i_tonedur}{1}.time == TP_StartTone{i_tonedur}(i_tone));
    end   
    
end   
    
%2. Compute min/max vals across conditions and TD for common y-axis    
minAmp = 1000; %set high, so overwritten
maxAmp = -1000;
for i_tonedur = 1:length(TimelockData_perTDPredp34)
    for i_Predp34 = 1:length(label_Predp34)
        minAmp_temp = min(TimelockData_perTDPredp34{i_tonedur}{i_Predp34}.(inputData_Avg)(ChanIndex,Sample_StartTone{i_tonedur}(1):Sample_StartTone{i_tonedur}(end)));
        maxAmp_temp = max(TimelockData_perTDPredp34{i_tonedur}{i_Predp34}.(inputData_Avg)(ChanIndex,Sample_StartTone{i_tonedur}(1):Sample_StartTone{i_tonedur}(end)));
        
        if minAmp_temp < minAmp
            minAmp = minAmp_temp;
        end
        if maxAmp_temp > maxAmp
            maxAmp = maxAmp_temp;
        end    
    end  
end

%3. Plot all Predp34 per TD in subpot and 1 subplot per TD
for i_tonedur = 1:length(TimelockData_perTDPredp34)

    subplot_position = subplot_position + 1;
    subplot(2,1,subplot_position);
    hold on;
    
    for i_Predp34 = 1:length(label_Predp34)
        
        %GAvg per Predp34
        plot(1:length(TimelockData_perTDPredp34{i_tonedur}{i_Predp34}.(inputData_Avg)(ChanIndex,:)),...
            TimelockData_perTDPredp34{i_tonedur}{i_Predp34}.(inputData_Avg)(ChanIndex,:),...
            'color',plot_param.color{i_Predp34},'LineStyle',plot_param.linestyle{i_Predp34},'LineWidth',2)
        
        %STD per Predp34
        if plot_STD == 1
            
            plot_param.legend_label = {'predp34 = low','', 'predp34 =  med','', 'predp34 =  high',''};

            
            shadedErrorBar(1:length(TimelockData_perTDPredp34{i_tonedur}{i_Predp34}.(inputData_Avg)(ChanIndex,:)),...
                TimelockData_perTDPredp34{i_tonedur}{i_Predp34}.(inputData_Avg)(ChanIndex,:), ...
                TimelockData_perTDPredp34{i_tonedur}{i_Predp34}.(inputData_STD)(ChanIndex,:),...
                'lineprops',{'color',plot_param.color{i_Predp34}})
        end
        
        hold on;
        
    end
    
    %Vertical line for each tone/event
    %Create aray that has max/min deflection at sample of each tone/event
    Array_Tone_Pos = zeros(1,length(TimelockData_perTDPredp34{i_tonedur}{i_Predp34}.time));
    Array_Tone_Neg = zeros(1,length(TimelockData_perTDPredp34{i_tonedur}{i_Predp34}.time));
    
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
        
    xlim([1 length(TimelockData_perTDPredp34{i_tonedur}{i_Predp34}.(inputData_Avg)(ChanIndex,:))])  
    
    ylim([minAmp maxAmp])
    if plot_STD == 1
        ylim([minAmp*2 maxAmp*2])
    end
        

    
    %x-axis shows tones [tone number]
    set(gca,'FontSize',10,'XTickLabelRotation',45)
    set(gca,'xtick',Sample_StartTone{i_tonedur})
    x_label = plot_param.Labels_ToneStart;
    set(gca,'FontSize',10)
    xticklabels(x_label);   
     
    if i_tonedur == 1
        legend(plot_param.legend_label,'Location','NorthWest')
    end
    
    
    title({['ToneDur: ' ToneDur_text{i_tonedur} ' s']});
    suptitle({[sub '; Chan: ' TimelockData_perTDPredp34{i_tonedur}{i_Predp34}.label{ChanIndex} '; ' inputData],...
        ['ERP per p*34 (GAvg + STD across trials (n = ' num2str(size(TimelockData_perTDPredp34{i_tonedur}{i_Predp34}.trial,1)) ')) timelocked to tone presentation']});
    
end  


if save_fig == 1
    path_fig = ([paths_NASTD_ECoG.Fig_Timelocked_SingleSubs sub '/Predp34/' inputData '/']);
    mkdir(path_fig);
    filename     = [TimelockData_perTDPredp34{i_tonedur}{i_Predp34}.label{ChanIndex} '_' inputData '_' sub '_TLdata.png'];
    figfile      = [path_fig filename];
    saveas(gcf, [figfile], 'png'); %save png version
    close
end


end
    
