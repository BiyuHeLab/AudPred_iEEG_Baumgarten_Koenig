function NASTD_ECoG_Predict_CompPred_TW_Subs...
    (sub, FuncInput_DataType, FuncInput_ToneDur_text,  ...
    param,...
    FuncInput_InputData, ...
    plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG)

%Aim: Compute prediction effect (How is neural activity at tone 33
%modulated by the predicted pitch of tone 34 (predp34)?) and
%prediction error effect (How is neural activity at tone 34
%modulated by the difference between predicted (predp34) and presented
%(p34) pitch of tone 34?)

%Method: Linear regression of pitch values p1-32 on windowed neural activity during p33.

%% 0.1) Specify vars, paths, and setup fieldtrip
%Add base dir and own script dir
addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));

path_dataoutput = [paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
    sub '/Data/' param.Label_TW 'sTW/'];
if (~exist(path_dataoutput, 'dir')); mkdir(path_dataoutput); end
    
path_fig = ([paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
    sub '/Figs/ERFs/PredEffect/' FuncInput_DataType '/']);
if (~exist(path_fig, 'dir')); mkdir(path_fig); end

%% 1. Select current input data in FT-struct
FuncInput_InputData.trial = [];
FuncInput_InputData.trial = FuncInput_InputData.(FuncInput_DataType);

%% 2) Extract p33 and p34 per trial
%2.1 Define time window parameters depending on TD
nTrials     = length(FuncInput_InputData.trial);
nSensors    = size(FuncInput_InputData.trial{1},1);
SampleFreq  = FuncInput_InputData.fsample;
ToneDur_Sec  = str2num(FuncInput_ToneDur_text);
nSamples_perTone = ToneDur_Sec * SampleFreq;

%2.2 Determine TP/samples for each tone start+end
TP_Tone_StartStop = NaN(36,2);
TP_Tone_StartStop(1,1) = 0; %p1 set as t = 0 in trial definition
for i_tone = 2:36
    Dist = abs(FuncInput_InputData.time{1} - ((str2num(FuncInput_ToneDur_text)*i_tone) - str2num(FuncInput_ToneDur_text)));
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

%2.3) Initialize and fill 3D data arrays (nSens, nSamples per tone, nTrials)
Data_p33 = zeros(nSensors, ...
    length(Sample_Tone_StartStop(1,1):Sample_Tone_StartStop(1,2)), ...
    nTrials);
Data_p34 = Data_p33;
%copy trial wise info (i.e., all samples for selected tone for all channels and all trials into data arrays
for i_trial = 1:nTrials
    Data_p33(:, :, i_trial) = ...
        FuncInput_InputData.trial{i_trial}(:, ...
        Sample_Tone_StartStop(param.ToneIndex,1) : Sample_Tone_StartStop(param.ToneIndex,2));
    Data_p34(:, :, i_trial) = ...
        FuncInput_InputData.trial{i_trial}(:, ...
        Sample_Tone_StartStop(param.ToneIndex+1,1) : Sample_Tone_StartStop(param.ToneIndex+1,2));    
    Data_p1to33(:, :, i_trial) = ...
        FuncInput_InputData.trial{i_trial}(:, ...
        Sample_Tone_StartStop(1,1) : Sample_Tone_StartStop(param.ToneIndex,2));
end

%% 3) Extract Predp34 and Prediction Error
%For each trial, compute/read out
%1) Predp34 (discretized to next full tone)
%2) p34 (i.e., the really presented last tone)
%3) Simple Prediction error (i.e., the distance between the final (p34) and penultimate (p33( tone
%4) Complex Prediction error (i.e., the distance between presented (p34) and predicted (Predp34) tone

for i_k = 1:length(param.PredSeqrange) %for each tone for which predition is analyzed
    
    series_start_ind = param.ToneIndex - param.PredSeqrange(i_k) + 1; %beginning of tone sequence part used for prediction
    series_end_ind   = param.ToneIndex; %end of tone sequence part used for prediction
    
    for i_trial = 1:nTrials
        series = FuncInput_InputData.behav.stim.series_f{i_trial}...
            (series_start_ind : series_end_ind); %(non-log) selected sequence tones in Hz
        beta = FuncInput_InputData.behav.stim.beta(i_trial); %always 1.5
         
        LogFreq_p33(i_trial) = log( FuncInput_InputData.behav.stim.series_f{i_trial}(param.ToneIndex) ); %log(p33)
        LogFreq_p34(i_trial) = log( FuncInput_InputData.behav.stim.series_f{i_trial}(param.ToneIndex+1) ); %log(p34)
        LogFreq_predp34_discretized(i_trial) = FuncInput_InputData.behav.stim.logf_pred(i_trial); %log(Predp34) discretized
        
        LogFreq_ComplexPredError(i_k, i_trial) = ...
            LogFreq_p34(i_trial) - LogFreq_predp34_discretized(i_trial);        
        LogFreq_SimplePredError(i_k, i_trial) = ...
            LogFreq_p34(i_trial) - LogFreq_p33(i_trial);        
    end
    
end

%% 4) Compute avg ERF activity per time window for current tone and next tone
%4.1) Define parameters for window used to commpute neural response
win_size    = param.SamplesTW; 
s_per_win = (1/SampleFreq)*win_size;
win_overlap = 0;

%4.2 Define number, start and end sample of window per tone
windows = [1 win_size];
while windows(end,end) < min(nSamples_perTone)
    windows = [windows; windows(end,:) + (win_size - win_overlap)];
end
if windows(end,end) > min(nSamples_perTone) 
    windows(end,:) = [];
end

windows_inms = (windows / SampleFreq) * 1000;

%4.3) Compute avg activity per window for each sensor and trial
for i_win = 1:size(windows,1)
    
    ind_start = windows(i_win, 1);
    ind_end   = windows(i_win, 2);
    
    ERF_p33_win{i_win} = squeeze( mean( Data_p33(:, ind_start:ind_end, :), 2) );
    ERF_p34_win{i_win} = squeeze( mean( Data_p34(:, ind_start:ind_end, :), 2) );
    %Output: channel*trial matrix per window with 1 ERF values averaged
    %across all samples per window
end

%% Optional: 5) Plot neural activity during p33 as a funciton of grouped Predp34
if plot_poststepFigs == 1
    %Separate trials by Predp34
    label_Predp34 = [-1 0 1];% low, medium, high;
    for i_Predp34 = 1:length(label_Predp34)
        TrialFilt_Predp34 = FuncInput_InputData.behav.stim.predID == label_Predp34(i_Predp34);
        Index_Predp34{i_Predp34} = find(TrialFilt_Predp34 == 1);
    end
    
    %Compute AVG + STD per Predp34
    for i_Predp34 = 1:length(label_Predp34)
        Avgp33{i_Predp34} = mean(Data_p33(:,:,Index_Predp34{i_Predp34}),3);
        STDp33{i_Predp34} = std(Data_p33(:,:,Index_Predp34{i_Predp34}),1,3);
        
        Avgp1to33{i_Predp34} = mean(Data_p1to33(:,:,Index_Predp34{i_Predp34}),3);
        STDp1to33{i_Predp34} = std(Data_p1to33(:,:,Index_Predp34{i_Predp34}),1,3);
    end
    
    %Plot ERF for all selected elecs
    NumFig = ceil(nSensors/30); %20 elecs per plot
    plot_param.color = {[0, 0.4470, 0.7410, 0.7],[0, 0.75, 0.75, 0.7],[0.8500, 0.3250, 0.0980, 0.7]};
    
    %ERF for p33
    i_elec = 0;
    for i_fig = 1:NumFig
        h = figure;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        for i_subplot = 1:30
            if i_elec < nSensors
                subplot(5,6,i_subplot)
                i_elec = i_elec+1;
                
                plot(1:length(Avgp33{1}(i_elec,:)),Avgp33{1}(i_elec,:),...
                    'color',plot_param.color{1},'LineWidth', 2)
                hold on;
                plot(1:length(Avgp33{2}(i_elec,:)),Avgp33{2}(i_elec,:),...
                    'color',plot_param.color{2},'LineWidth', 2)
                hold on;
                plot(1:length(Avgp33{3}(i_elec,:)),Avgp33{3}(i_elec,:),...
                    'color',plot_param.color{3},'LineWidth', 2)
                title(['Elec: ' FuncInput_InputData.label{i_elec}])
                axis tight
                a = axis;
                xlabel('TW p33')
                %Add demarcation lines for TW (50ms)
                for i_win = 1:size(windows,1)
                    ind_start = windows(i_win, 1);
                    hold on;
                    plot([ind_start ind_start],[a(3) a(4)], 'k--')
                end
                
            end
        end
        sgtitle([sub ' - ' FuncInput_DataType ' - ' FuncInput_ToneDur_text ...
            'msTD - ERF p33 as a function of p*34 (' ...
            num2str(i_fig) '/' num2str(NumFig) ')'],'Interpreter','none')
        
        if save_poststepFigs == 1
            filename     = ['p33ERF_AllElecs_' sub '_' FuncInput_DataType '_' ...
                FuncInput_ToneDur_text 'sTD_' num2str(i_fig) '.png'];
            figfile      = [path_fig filename];
            saveas(gcf, [figfile], 'png'); %save png version
        end
        close
    end
    
    %ERF for p1-p33
    i_elec = 0;
    for i_fig = 1:NumFig
        h = figure;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        for i_subplot = 1:30
            if i_elec < nSensors
                subplot(5,6,i_subplot)
                i_elec = i_elec+1;
                
                plot(1:length(Avgp1to33{1}(i_elec,:)),Avgp1to33{1}(i_elec,:),...
                    'color',plot_param.color{1},'LineWidth', 2)
                hold on;
                plot(1:length(Avgp1to33{2}(i_elec,:)),Avgp1to33{2}(i_elec,:),...
                    'color',plot_param.color{2},'LineWidth', 2)
                hold on;
                plot(1:length(Avgp1to33{3}(i_elec,:)),Avgp1to33{3}(i_elec,:),...
                    'color',plot_param.color{3},'LineWidth', 2)
                title(['Elec: ' FuncInput_InputData.label{i_elec}])
                axis tight
                a = axis;
                xlabel('p1-p33')
                %Add demarcation lines for tones
                for i_tone = 1:size(Sample_Tone_StartStop,1)
                    ind_start = Sample_Tone_StartStop(i_tone, 1) - 512;
                    hold on;
                    plot([ind_start ind_start],[a(3) a(4)], 'Color',[0.5 0.5 0.5])
                end
            end
        end
        sgtitle([sub ' - ' FuncInput_DataType ' - ' FuncInput_ToneDur_text ...
            'msTD - ERF p1-p33 as a function of p*34 (' ...
            num2str(i_fig) '/' num2str(NumFig) ')'],'Interpreter','none')
        
        if save_poststepFigs == 1
           filename     = ['p1to33ERF_AllElecs_' sub '_' FuncInput_DataType '_' ...
                FuncInput_ToneDur_text 'sTD_' num2str(i_fig) '.png'];
            figfile      = [path_fig filename];
            saveas(gcf, [figfile], 'png'); %save png version
            close
        end
    end
end

%% 6) Compute measures of association between ERF windows and prediction / prediction error at each sensor

for i_k = 1:length(param.PredSeqrange)
    for i_win = 1:size(windows,1)
        for i_sensor = 1:nSensors
            
            %6.1 Prediction Effect
            dv = ERF_p33_win{i_win}(i_sensor,:)';  %p33 ERF (averaged across window) at this sensor for all trials
            iv_pred = LogFreq_predp34_discretized(i_k,:)'; %predicted p34 tone pitch for all trials            
            % linear regression between ERF (per channel, time window) and predicted tone pitch
            stats = regstats(dv, iv_pred, 'linear', 'tstat');            
            PredEffect.stats{i_k}{i_win}{i_sensor}  = stats;
            PredEffect.tval{i_k}{i_win}(i_sensor)   = stats.tstat.t(2);
            PredEffect.pval{i_k}{i_win}(i_sensor)   = stats.tstat.pval(2);
            PredEffect.beta{i_k}{i_win}(i_sensor)   = stats.tstat.beta(2);
            stats = [];
            
            %6.2 Simple Prediction Error Effect
            dv = ERF_p34_win{i_win}(i_sensor,:)';  %p34 ERF (averaged across window) at this sensor for all trials 
            iv_simplepred_error = abs( LogFreq_SimplePredError(i_k,:)' );           
            iv_p34 = LogFreq_p34';            
            % linear regression with 2 predictors (simple prediction error + p34)
            stats = regstats(dv, [iv_simplepred_error, iv_p34], 'linear', 'tstat');            
            SimplePredErrEffect.stats{i_k}{i_win}{i_sensor} = stats;
            SimplePredErrEffect.tval{i_k}{i_win}(i_sensor)  = stats.tstat.t(2);
            SimplePredErrEffect.pval{i_k}{i_win}(i_sensor)  = stats.tstat.pval(2);
            SimplePredErrEffect.beta{i_k}{i_win}(i_sensor)  = stats.tstat.beta(2);
            stats = [];
                        
            %6.3 Complex Prediction Error Effect
            dv = ERF_p34_win{i_win}(i_sensor,:)'; %p34 ERF (averaged across window) at this sensor for all trials  
            iv_pred_error = abs( LogFreq_ComplexPredError(i_k,:)' );           
            iv_p34 = LogFreq_p34'; %presented final tone (p34)            
            % linear regression with 2 predictors (complex prediction error + p34)
            stats = regstats(dv, [iv_pred_error, iv_p34], 'linear', 'tstat');            
            ComplexPredErrEffect.stats{i_k}{i_win}{i_sensor} = stats;
            ComplexPredErrEffect.tval{i_k}{i_win}(i_sensor)  = stats.tstat.t(2);
            ComplexPredErrEffect.pval{i_k}{i_win}(i_sensor)  = stats.tstat.pval(2);
            ComplexPredErrEffect.beta{i_k}{i_win}(i_sensor)  = stats.tstat.beta(2);
            stats = [];
            
        end
        
    end
end

%% 7) Save variables
param_LoadedData = param;
labels_loadedData = FuncInput_InputData.label;

savefile = [path_dataoutput sub '_PredEffects_' FuncInput_DataType '_' FuncInput_ToneDur_text 'sTD.mat'];
save(savefile, 'PredEffect', 'SimplePredErrEffect', 'ComplexPredErrEffect', ...
    'param_LoadedData', 'labels_loadedData', ...
    'ERF_p33_win', 'ERF_p34_win', ...
    'Index_Predp34', ...
    'LogFreq*');

end