function NASTD_ECoG_StimCorr_ToneProc ...
    (sub, FuncInput_DataType, FuncInput_ToneDur_text,  ...
    param,...
    FuncInput_InputData, ...
    save_poststepFigs, paths_NASTD_ECoG)

%Aim: Compute single tone tracking in neural data to determine electrodes 
%that react to tones (inependent of sequence processing). 
%Compare N100 ERP against pre-tone baseline and select electrodes that are 
%sign. different from 0/baseline.

%% 0.1) Specify vars, paths, and setup fieldtrip
%Add base dir and own script dir
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));

path_dataoutput = [paths_NASTD_ECoG.ECoGdata_StimCorr '/' sub '/Data/'];
if (~exist(path_dataoutput, 'dir')); mkdir(path_dataoutput); end

path_fig = ([paths_NASTD_ECoG.ECoGdata_StimCorr ...
    '/' sub '/Figs/ToneProc/']);
if (~exist(path_fig, 'dir')); mkdir(path_fig); end
path_fig_TLdata = ([paths_NASTD_ECoG.ECoGdata_StimCorr ...
    '/' sub '/Figs/ToneProc/TLdata_WholeSeq/']);
if (~exist(path_fig_TLdata, 'dir')); mkdir(path_fig_TLdata); end


%% 1. Select current input data in FT-struct
FuncInput_InputData.trial = [];
FuncInput_InputData.trial = FuncInput_InputData.(FuncInput_DataType);

%% 2. Determine time points and samples for each tone, the P100, and all baselines (BL)
nTrials             = length(FuncInput_InputData.trial);
nSensors            = size(FuncInput_InputData.trial{1},1);
SampleFreq          = FuncInput_InputData.fsample;
ToneDur_Sec         = str2num(FuncInput_ToneDur_text);
nSamples_perTone    = ToneDur_Sec * SampleFreq;
SensorLabels        = FuncInput_InputData.label;

%2.2 Determine TP/samples for each tone start+end
TP_Tone_StartStop       = NaN(36,2);
TP_Tone_StartStop(1,1)  = 0; %p1 set as t = 0 in trial definition
for i_tone = 2:36
    Dist = abs(FuncInput_InputData.time{1} - ...
        ((str2num(FuncInput_ToneDur_text)*i_tone) - str2num(FuncInput_ToneDur_text)));
    minDist = min(Dist);
    i_minDist = find(Dist == minDist);
    TP_Tone_StartStop(i_tone,1) = FuncInput_InputData.time{1}(i_minDist);
end
for i_tone = 1:35
    i_LastSampleTone = find(FuncInput_InputData.time{1} == TP_Tone_StartStop(i_tone + 1,1));
    TP_Tone_StartStop(i_tone,2) = FuncInput_InputData.time{1}(i_LastSampleTone);
end
TP_Tone_StartStop = TP_Tone_StartStop(1:34,:);

%Check if all tones are of equal length
%(if not then choose min length by deleting the last sample of longer trials)
minSeqLength_sec = min(TP_Tone_StartStop(:,2)-TP_Tone_StartStop(:,1));
for i_tone = 1:34
    if (TP_Tone_StartStop(i_tone,2) - TP_Tone_StartStop(i_tone,1)) > minSeqLength_sec
        TP_Tone_StartStop(i_tone,2) = ...
            FuncInput_InputData.time{1}(...
            find(FuncInput_InputData.time{1} == TP_Tone_StartStop(i_tone,2))-1);
    end
end

%Determine samples corresponding to TP
Sample_Tone_StartStop = NaN(34,2);
for i_tone = 1:34
    Sample_Tone_StartStop(i_tone,1) = ...
        find(FuncInput_InputData.time{1} == TP_Tone_StartStop(i_tone,1));
    Sample_Tone_StartStop(i_tone,2) = ...
        find(TP_Tone_StartStop(i_tone,2) == FuncInput_InputData.time{1});
end

%Determine N100 start/stop samples relative to each tone TW
N100_samples = ...
    [round(param.N100_TW(1)*FuncInput_InputData.fsample) ...
    round(param.N100_TW(2)*FuncInput_InputData.fsample)];

%Determine pre-sequence and pre-tone baseline samples relative to each tone TW
BL_tone_samples = [1 round(param.BL_tone(2)*FuncInput_InputData.fsample)];


%% 3. Compute N100 ERPs based on pre-tone baseline (i.e., baseline normalization for every tone)

clear Data_perTone Data_perTone_BLperTone Data_avgTone_BLperTone Data_N100Amp_perTrialElec

%3.1 For each electrode and trial, cut neural data into time windows
%matching single tones
for i_elec = 1:nSensors
    for i_trial = 1:nTrials
        for i_tone = 1:length(Sample_Tone_StartStop) 
            
            Data_perTone{i_elec}(i_trial, i_tone, :) = ...
                FuncInput_InputData.trial{i_trial}(i_elec, ...
                Sample_Tone_StartStop(i_tone,1):Sample_Tone_StartStop(i_tone,2));
            
        end    
    end
end

%3.2 If input data is LP35Hz, square raw amplitude values to deal with
%positive/negative deflections. Take square root afterwards to return to
%original scale.
if strcmp(FuncInput_DataType,'LP35Hz')                            
    for i_elec = 1:nSensors
        for i_trial = 1:nTrials
            for i_tone = 1:length(Sample_Tone_StartStop)

                Data_perTone{i_elec}(i_trial, i_tone, :) = ...
                    sqrt(Data_perTone{i_elec}(i_trial, i_tone, :) .^2);
            end
        end
    end
end

%3.3 For each single tone, baseline-correct neural data with pre-tone baseline
for i_elec = 1:nSensors
    for i_trial = 1:nTrials
        for i_tone = 1:length(Sample_Tone_StartStop) 
            
            BL_perTone = [];
            BL_perTone = mean(Data_perTone{i_elec}...
                (i_trial, i_tone, BL_tone_samples(1):BL_tone_samples(end)));
            
            Data_perTone_BLperTone{i_elec}(i_trial, i_tone, :) = ...
                Data_perTone{i_elec}(i_trial, i_tone, :) - BL_perTone;
        end    
    end
end
 
%3.4 For each trial, average neural data across tones 
for i_elec = 1:nSensors
    for i_trial = 1:nTrials
        
        Data_avgTone_BLperTone{i_elec}(i_trial, :) = ...
            squeeze(mean(Data_perTone_BLperTone{i_elec}(i_trial, :, :)))';
        
    end
end

% %3.4.2 Plot ERFs averaged across tones for all trials for an exemplary electrode
% i_elec = randi(nSensors,1);
% h = figure;
% for i_trial = 1:nTrials
%     hold on;
%     plot((Data_avgTone_BLperTone{i_elec}(i_trial,:)))
% end
% hold on; plot([N100_samples(1), N100_samples(1)],...
%     [h.CurrentAxes.YLim(1), h.CurrentAxes.YLim(2)], 'k--')
% hold on; plot([N100_samples(2), N100_samples(2)],...
%     [h.CurrentAxes.YLim(1), h.CurrentAxes.YLim(2)], 'k--')
% hold on; plot(mean(Data_avgTone_BLperTone{i_elec}),'LineWidth',2,'Color',[0 0 1])
% 
% xlabel('Time [samples]')
% ylabel('Amplitude')
% title({'Neural data timelocked to individual tone onset', ...
%     ['(per trial, averaged across tones, pre-tone BL-corrected) for elec: ' ...
%     FuncInput_InputData.label{i_elec}]})

%3.4 Average amplitude values across N100 samples
for i_elec = 1:nSensors
    for i_trial = 1:nTrials
        Data_N100Amp_BLperTone_perTrialElec(i_trial,i_elec) = ...
            mean(Data_avgTone_BLperTone{i_elec}(i_trial, ...
            N100_samples(1) : N100_samples(2)));
    end
end

%3.5 Statistically compare across-trial N100 amplitudes against 0 for each electrode
for i_elec = 1:nSensors
    [~,p,~,stats] = ttest(Data_N100Amp_BLperTone_perTrialElec(:,i_elec),0,'tail','right');
    stat_N100Ampvs0_BLperTone.t(i_elec,1) = stats.tstat;
    stat_N100Ampvs0_BLperTone.p(i_elec,1) = p;    
end

%% 4. Save data
savefile = [path_dataoutput sub '_ToneProc_' FuncInput_DataType '_' ...
    FuncInput_ToneDur_text 'sTD.mat'];

save(savefile, ...
    'stat_N100Ampvs0_BLperTone', ...
    'SensorLabels', 'param', ...
    '-v7.3');

%% 5. Plot stats of N100 amplitudes on surface brain and highlight sign. electrodes
%Plot summary figure containing:
%1) t-values for N100 amplitudes vs. zero for all electrodes
%2) t-values for N100 amplitudes vs. zero for sign. electrodes only
%3) Table with anatomical labels for sign. electrodes

clear plot_struct

%5.0 Set up subplot structure
h = figure('visible','on'); %ensures figure doesn't pop during plotting
set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen

DimSubplot          = [2 2];
SubplotPosition     = [0 -0.05 0 0];
ColorbarPosition    = [0 0 0 0];
SizeFactor          = 3;
CounterSubplot      = 1;

sgtitle({...
    ['Neural data timelocked to each single tone (N100: ' ...
    num2str(param.N100_TW(1)) '-' num2str(param.N100_TW(2)) 's; ' ...
    'pre-tone baseline-correction, ' ...
    num2str(param.BL_tone(1)) '-' num2str(param.BL_tone(2)) 's)'], ...
    [sub ' - ' FuncInput_DataType ' - ' FuncInput_ToneDur_text ...
    ' - (RH elecs projected on LH)']}, ...
    'Interpreter','none')

%Find index of each electrode label in cleaned data .elec struct
%and read out elec-specific MNI coordinates for plotting
Index_Elecs2ChanPos = NaN(1,nSensors);
for i_elec = 1:length(FuncInput_InputData.label)
    Index_Elecs2ChanPos(1,i_elec) = ...
        find(strcmp(FuncInput_InputData.label{i_elec}, ...
        FuncInput_InputData.elec.label));
end

%Determine sign. electrodes and their full label
SignElecs.index = find(stat_N100Ampvs0_BLperTone.p < param.pval_plotting);
nSignElecs      = length(SignElecs.index);
for i_elec = 1:nSignElecs
    temp_elecindex = ...
        find(strcmp(FuncInput_InputData.label{SignElecs.index(i_elec)}, ...
        FuncInput_InputData.elec.label));
    SignElecs.anatlabel{i_elec,1} = ...
        [FuncInput_InputData.elec.label{temp_elecindex} ' ' ...
        FuncInput_InputData.elec.T1AnatLabel{temp_elecindex}];
end

plot_struct.coords = FuncInput_InputData.elec.chanpos(Index_Elecs2ChanPos,:);
%Project all electrodes on left hemisphere,
plot_struct.coords(:,1) = abs(plot_struct.coords(:,1)) * -1;

%5.1 Plot surface plot with correlation coefficient for identical sequences
plot_struct.dimord          = 'chan_time';
plot_struct.time            = 0;
plot_struct.sign_elecs      = logical(ones(1,nSensors));
plot_struct.clims           = [0 5]; 
plot_struct.chanSize        = ...
    ones(1,length(plot_struct.sign_elecs))*SizeFactor; %electrode size (arbitrary)
plot_struct.cmap            = 'jet';
plot_struct.textcolor_rgb   = [0 0 0];

for i_elec = 1:length(FuncInput_InputData.label)
    plot_struct.label{i_elec} = '';
end

plot_struct.avg         = stat_N100Ampvs0_BLperTone.t;

%Plot surface with highlighted sign electrodes
sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH...
    (plot_struct.coords, plot_struct.label, plot_struct.avg, plot_struct.sign_elecs,...
    plot_struct.chanSize, plot_struct.clims, plot_struct.cmap, plot_struct.textcolor_rgb, ...
    DimSubplot, CounterSubplot, SubplotPosition, ColorbarPosition, []);

sp_handle_surf{CounterSubplot} = sp_handle_surf_temp.L;
CounterSubplot = CounterSubplot + 1;
title(['all electrodes'],'FontSize',12)

%5.2 Plot surface plot with correlation coefficient for identical sequences
for i_signelecs = 1:nSignElecs    
    plot_struct.label{SignElecs.index(i_signelecs)} = ...
        FuncInput_InputData.label{SignElecs.index(i_signelecs)};
end

plot_struct.avg         = stat_N100Ampvs0_BLperTone.t;
plot_struct.sign_elecs  = logical(stat_N100Ampvs0_BLperTone.p < param.pval_plotting);
%Plot surface with highlighted sign electrodes
sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH...
    (plot_struct.coords, plot_struct.label, plot_struct.avg, plot_struct.sign_elecs,...
    plot_struct.chanSize, plot_struct.clims, plot_struct.cmap, plot_struct.textcolor_rgb, ...
    DimSubplot, CounterSubplot, SubplotPosition, ColorbarPosition, []);

sp_handle_surf{CounterSubplot} = sp_handle_surf_temp.L;
CounterSubplot = CounterSubplot + 1;
title(['p < ' num2str(param.pval_plotting) ' (uncorrected)'], 'FontSize',12)

%Add colorbar
ColorbarPosition = ...
    [sp_handle_surf{CounterSubplot-1}.Position(1)-0.05 ...
    sp_handle_surf{CounterSubplot-1}.Position(2)+0.05 0 0];
h = colorbar;
h.Position(1) = ColorbarPosition(1); %sets colorbar to the right
h.Position(2) = ColorbarPosition(2); %sets colorbar higher
h.Position(4) = h.Position(4)*0.8; %makes colorbar shorter
h.Label.String = {['t-value'],...
    ['(trial-wise N100 amplitude'], ...
    ['averaged across tones vs. 0)']};
h.FontSize = 14;
caxis(plot_struct.clims)

%5.3 Table with analysis & electrode information
subplot(DimSubplot(1), DimSubplot(2), 3)
textbox_info = {...
    [sub] ...
    ['# sign. elecs / all elecs: ' ...
    num2str(nSignElecs) ' / ' num2str(nSensors)] ...
    ''};
textbox_elec = [];

if nSignElecs <= 15 %one row
    for i_elec = 1:nSignElecs
        sel_elec = SignElecs.index(i_elec);
        textbox_elec{i_elec} = ...
            [FuncInput_InputData.label{sel_elec} ' = ' SignElecs.anatlabel{i_elec}];
    end
    textbox = [textbox_info textbox_elec];
    t = text(0, 0.5, 0, textbox, 'FontSize',12,'Interpreter','none');
else
    nTextColumns = ceil(nSignElecs/15);
    i_elec = 1;
    for i_columns = 1:nTextColumns
        while i_elec <= 15 * i_columns && i_elec <= nSignElecs
            sel_elec = SignElecs.index(i_elec);
            place_elec = i_elec - (15 * (i_columns -1));
            textbox_elec{i_columns}{place_elec} = ...
                [FuncInput_InputData.label{sel_elec} ' = ' SignElecs.anatlabel{i_elec}];
            i_elec = i_elec +1;
        end
        if i_columns == 1
            textbox = [textbox_info textbox_elec{i_columns}];
            t = text(-0.3, 0.5, 0, textbox, ...
                'FontSize',12,'Interpreter','none');
        else
            textbox = [];
            textbox = textbox_elec{i_columns};
            subplot(DimSubplot(1), DimSubplot(2), 3)
            t = text(-0.3 + (0.6 * (i_columns -1)), 0.5, 0, textbox, ...
                'FontSize',12,'Interpreter','none');
        end
    end
end
set(gca,'visible','off')

%5.4 Save Figure
if save_poststepFigs == 1
    filename     = ['N100ToneProc_Surf1H_' sub '_' ...
        FuncInput_DataType '_' FuncInput_ToneDur_text 'sTD.png'];
    figfile      = [path_fig filename];
    saveas(gcf, figfile, 'png'); %save png version
    close all;
end

%% 6 Plot across-trial-averaged whole-sequence time series neural data for each sign. electrode
%Aim: To get an impression of the neural response across the whole sequence

clear plot_struct

%6.0 Set up subplot structure
NumPlots        = ceil(nSignElecs/16);
elec_counter    = 0;
i_subplot       = 1;

for i_plot = 1:NumPlots
    
    %6.1 Set up current subplot
    h = figure('visible','on'); %ensures figure doesn't pop during plotting
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    DimSubplot          = [4 4];
    SubplotPosition     = [0 -0.05 0 0];
    
    currelecs = (1+elec_counter):(16+elec_counter);
    if currelecs(end) > nSignElecs
        currelecs(end) = nSignElecs;
    end
    
    sgtitle({...
        ['Timelocked neural data averaged across trials for all electrodes with sign. tone proc effect (N100: ' ...
        num2str(param.N100_TW(1)) '-' num2str(param.N100_TW(2)) 's; ' ...
        'pre-tone baseline-correction, ' ...
        num2str(param.BL_tone(1)) '-' num2str(param.BL_tone(2)) 's)'], ...
        [sub ' - ' FuncInput_DataType ' - ' FuncInput_ToneDur_text ...
        ' - Elec: ' num2str(currelecs(end)) '/' num2str(nSignElecs) ' - ' ...
        'Plot: ' num2str(i_plot) '/' num2str(NumPlots)]}, ...
        'Interpreter','none')
    
    while i_subplot <= 16 * i_plot && i_subplot <= nSignElecs
        
        %6.2 Read out and restructure each trial for current electrode
        i_currelec = SignElecs.index(i_subplot);
        
        samples_pertrial = [];
        for i_trial = 1:length(FuncInput_InputData.trial)
            samples_pertrial = ...
                [samples_pertrial, length(FuncInput_InputData.trial{i_trial})];
        end
        temvar_alltrial_currelec = ...
            nan(length(FuncInput_InputData.trial),max(samples_pertrial));
        for i_trial = 1:length(FuncInput_InputData.trial)
            temvar_alltrial_currelec...
                (i_trial,1:length(FuncInput_InputData.trial{i_trial})) = ...
                FuncInput_InputData.trial{i_trial}(i_currelec,:);
        end
        
        %6.3 Plot across-trial average and std
        h = subplot(DimSubplot(1), DimSubplot(2), i_subplot - (16 * (i_plot -1)));
        plot(nanmean(temvar_alltrial_currelec), 'b-', 'LineWidth',1)
        %     hold on;
        %     shadedErrorBar([],nanmean(temvar_alltrial_currelec), ...
        %         nanstd(temvar_alltrial_currelec), ...
        %         'lineprops',{'color',[0 0 1]})
        
        %6.4 add demarkation lines for each tone start
        for i_tone = 1:length(Sample_Tone_StartStop)
            hold on;
            if i_tone == 1 || i_tone == 10 || i_tone == 20 || i_tone == 30
                LineColor_tones = [0 0 0];            
            else
                LineColor_tones = [0.8 0.8 0.8];
            end
            
            plot([Sample_Tone_StartStop(i_tone,1), Sample_Tone_StartStop(i_tone,1)], ...
                [h.YLim(1), h.YLim(2)], ...
                'Color', LineColor_tones, 'LineWidth', 0.25);
        end
        h = subplot(DimSubplot(1), DimSubplot(2), i_subplot - (16 * (i_plot -1)));
        plot(nanmean(temvar_alltrial_currelec), 'b-', 'LineWidth',1)
        
        xlim([Sample_Tone_StartStop(1,1)-100 Sample_Tone_StartStop(end,end)+100])
        xlabel('Time')
        ylabel('Amplitude')        
        title(['Elec: ' FuncInput_InputData.label{i_currelec} ...
            ' (p = ' num2str(round(stat_N100Ampvs0_BLperTone.p(i_currelec),4)) ')'])
        
        i_subplot = i_subplot + 1;
        elec_counter = elec_counter + 1;        
    end
    
    %6.5 Save Figure
    if save_poststepFigs == 1
        filename     = ['TLdata_WholeSeqAvgTrials_' sub '_' ...
            FuncInput_DataType '_' FuncInput_ToneDur_text 'sTD_' ...
            num2str(i_plot) '.png'];
        figfile      = [path_fig_TLdata filename];
        saveas(gcf, figfile, 'png'); %save png version
        close all;
    end
    
end

end