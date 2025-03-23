function NASTD_ECoG_Predict_ComparePredErrEffects_TWvsSample...
    (sub, FuncInput_EffectType, FuncInput_DataType, FuncInput_ToneDur_text,  ...
    param,...
    FuncInput_InputData, ...
    save_poststepFigs, paths_NASTD_ECoG)

%Aims:
%Plot neural signals per predp34-p34 distance during p34 and highlight both
%sign. time windows determined via TW-approach and sign. samples determined
%via sample-based approach.

%1. Load in TW and sample base data
%% 0.1) Specify vars, paths, and setup fieldtrip

%Add base dir and own script dir
addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));

path_inputdata_TW = [paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
    sub '/Data/' param.Label_TW 'sTW/'];
path_inputdata_Sample = [paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
    sub '/Data/Samplewise/'];

path_fig = ([paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
    sub '/Figs/ERFs/' FuncInput_EffectType '/' FuncInput_DataType '/']);
if (~exist(path_fig, 'dir')); mkdir(path_fig); end

SampleFreq  = FuncInput_InputData.fsample;
nSensors    = size(FuncInput_InputData.trial{1},1);
nTrials     = length(FuncInput_InputData.trial);

%% 1. Load in data
PredErrEffectType_label = FuncInput_EffectType;
%1.1 TW data
load([path_inputdata_TW sub '_PredEffects_' FuncInput_DataType '_' FuncInput_ToneDur_text  'sTD.mat'], ...
    FuncInput_EffectType);
PredErrEffect_TW = eval(FuncInput_EffectType);

%1.2 Sample-wise data
load([path_inputdata_Sample sub '_PredEffectsCluster_' FuncInput_DataType '_' FuncInput_ToneDur_text 'sTD.mat'], ...
   FuncInput_EffectType , 'Data_p34');
PredErrEffect_Sample = eval(FuncInput_EffectType);
clear PredEffect SimplePredErrEffect ComplexPredErrEffect

%% 2. Separate trials per predp34-p34/p33-p34 pitch distances (predictor for lin. regr.)

%Determine pitch distances
for i_trial = 1:nTrials
    LogFreq_p33(i_trial) = log( FuncInput_InputData.behav.stim.series_f{i_trial}(param.ToneIndex) ); %log(p33)
    LogFreq_p34(i_trial) = log( FuncInput_InputData.behav.stim.series_f{i_trial}(param.ToneIndex+1) ); %log(p34)
    LogFreq_predp34_discretized(i_trial) = FuncInput_InputData.behav.stim.logf_pred(i_trial); %log(Predp34) discretized
    
    LogFreq_ComplexPredError(i_trial) = ...
        LogFreq_p34(i_trial) - LogFreq_predp34_discretized(i_trial);
    LogFreq_SimplePredError(i_trial) = ...
        LogFreq_p34(i_trial) - LogFreq_p33(i_trial);
end

%Determine unique pitch distances used
if strcmp(PredErrEffectType_label,'ComplexPredErrEffect')
    unique_predErrPitchDiff = unique(round(LogFreq_ComplexPredError,4));
    sortedunique_predErrPitchDiff = unique(sort(abs(unique_predErrPitchDiff)));   
elseif strcmp(PredErrEffectType_label,'SimplePredErrEffect')
    unique_predErrPitchDiff = unique(round(LogFreq_SimplePredError,4));
    sortedunique_predErrPitchDiff = unique(sort(abs(unique_predErrPitchDiff)));   
end

%Categorize into bins
% %Unique pitch distances (direction of pitch distance matters)
% NumBins = length(unique_predErrPitchDiff); 
% PitchperBin = length(unique_predErrPitchDiff)/NumBins;
% PredErr_PitchDiff = nan(NumBins, PitchperBin);
% counter = 0;
% for i_numbins = 1:NumBins
%     counter = counter + 1;
%     temp = unique_predErrPitchDiff...
%         (1 + (PitchperBin * (counter - 1)) : counter * PitchperBin);
%     PredErr_PitchDiff(i_numbins,1:length(temp)) = temp;
% end
% %Read out corresponding trial indices
% for i_numbins = 1:NumBins
%     Index_PredError{i_numbins} = [];
%     if strcmp(PredErrEffectType_label,'ComplexPredErrEffect')
%         for i_pitchperbin = 1:PitchperBin
%             Index_PredError{i_numbins} = [Index_PredError{i_numbins} ...
%                 find(round(LogFreq_ComplexPredError,4) == ...
%                 PredErr_PitchDiff(i_numbins,i_pitchperbin))];
%         end
%    elseif strcmp(PredErrEffectType_label,'SimplePredErrEffect')
%         for i_pitchperbin = 1:PitchperBin
%             Index_PredError{i_numbins} = [Index_PredError{i_numbins} ...
%                 find(round(LogFreq_SimplePredError,4) == ...
%                 PredErr_PitchDiff(i_numbins,i_pitchperbin))];
%         end       
%     end     
% end

%Unique absolute pitch distances (direction of pitch distance does not matter)
%Lin regr. also computed with absolute PredErr Input
NumBins = length(sortedunique_predErrPitchDiff); %Unique absolute pitch distances
PitchperBin = length(sortedunique_predErrPitchDiff)/NumBins;
PredErr_PitchDiff = nan(NumBins, PitchperBin);
counter = 0;
for i_numbins = 1:NumBins
    counter = counter + 1;
    temp = sortedunique_predErrPitchDiff...
        (1 + (PitchperBin * (counter - 1)) : counter * PitchperBin);
    PredErr_PitchDiff(i_numbins,1:length(temp)) = temp;
end
%Read out corresponding trial indices
for i_numbins = 1:NumBins
    Index_PredError{i_numbins} = [];
    if strcmp(PredErrEffectType_label,'ComplexPredErrEffect')
        for i_pitchperbin = 1:PitchperBin
            Index_PredError{i_numbins} = [Index_PredError{i_numbins} ...
                find(abs(round(LogFreq_ComplexPredError,4)) == ...
                PredErr_PitchDiff(i_numbins,i_pitchperbin))];
        end
   elseif strcmp(PredErrEffectType_label,'SimplePredErrEffect')
        for i_pitchperbin = 1:PitchperBin
            Index_PredError{i_numbins} = [Index_PredError{i_numbins} ...
                find(abs(round(LogFreq_SimplePredError,4)) == ...
                PredErr_PitchDiff(i_numbins,i_pitchperbin))];
        end       
    end     
end

%% 3. Define TWs
win_size    = param.SamplesTW; 
win_overlap = 0;
ToneDur_Sec  = str2num(FuncInput_ToneDur_text);
nSamples_perTone = ToneDur_Sec * SampleFreq;

windows = [1 win_size];
while windows(end,end) < min(nSamples_perTone)
    windows = [windows; windows(end,:) + (win_size - win_overlap)];
end
if windows(end,end) > min(nSamples_perTone) 
    windows(end,:) = [];
end

%% 4. Plot neural signal during p34
%Compute AVG + STD across trials per bin
for i_numbins = 1:NumBins
    Avgp34{i_numbins} = mean(Data_p34(:,:,Index_PredError{i_numbins}),3);
    STDp34{i_numbins} = std(Data_p34(:,:,Index_PredError{i_numbins}),1,3);
end

%Plot ERF for all selected elecs
NumFig = ceil(nSensors/30); %20 elecs per plot
plot_param.color = parula;

%ERF for p34
i_elec = 0;
for i_fig = 1:NumFig
    h = figure('visible','on'); %ensures figure doesn't pup during plotting
%     h = figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    for i_subplot = 1:30
        if i_elec < nSensors
            subplot(5,6,i_subplot)
            i_elec = i_elec+1;
            
            %Plot signal traces
            for i_numbins = 1:NumBins
                plot(1:length(Avgp34{i_numbins}(i_elec,:)),Avgp34{i_numbins}(i_elec,:),...
                    'color',plot_param.color(i_numbins*(floor(length(plot_param.color)/NumBins)),:), ...
                    'LineWidth', 2)
                hold on;
            end
            axis tight
            a = axis;
            xlabel('samples (p34)')
            
            %Highlight sign. samples/clusters
            maxY = [];
            minY = [];
            for i_numbins = 1:NumBins
                maxY = [maxY max(Avgp34{i_numbins}(i_elec,:))];
                minY = [minY min(Avgp34{i_numbins}(i_elec,:))];
            end
            maxY = max(maxY);
            minY = min(minY);
            
            if ~isnan(PredErrEffect_Sample.clusterstat{i_elec}.cluster_pval)
                for i_cluster = 1:length(PredErrEffect_Sample.clusterstat{i_elec}.cluster_pval)
                    if PredErrEffect_Sample.clusterstat{i_elec}.cluster_pval(i_cluster) < param.alpha
                        sign_samples = ...
                            find(PredErrEffect_Sample.clusterstat{i_elec}.cluster_timecourse ...
                            == i_cluster);
                        plot(sign_samples, ...
                            ones(length(sign_samples)) * max([minY maxY])+abs(max([minY maxY])*0.05),...
                            'linewidth',7,'Color',[175 175 175]./255)
                        text(sign_samples(1),max([minY maxY])+abs(max([minY maxY])*0.05), ...
                            [num2str(round(...
                            PredErrEffect_Sample.clusterstat{i_elec}.cluster_pval(i_cluster),3))], ...
                            'FontSize',6);
                    end
                end
            end
            
            %Highlight sign. TWs
            signTW = [];
            for i_win = 1:size(windows,1)
                %Add demarcation lines for TW (50ms)
                ind_TWstart = windows(i_win, 1);
                ind_TWtop = windows(i_win, 2);
                hold on;
                plot([ind_TWstart ind_TWstart],...
                     [min([minY maxY])-+abs(max([minY maxY])*0.35) ...
                    max([minY maxY])+abs(max([minY maxY])*0.35)], ...
                    'k--','LineWidth', 0.2)
                
                %Add shading
                if PredErrEffect_TW.pval{1}{i_win}(i_elec) < param.alpha
                    signTW = [signTW 1];
                    hold on;
                    plot(ind_TWstart:ind_TWtop, ...
                        ones(length(ind_TWstart:ind_TWtop)) * max([minY maxY])+abs(max([minY maxY])*0.35),...
                        'linewidth',7,'Color',[150 150 150]./255)
                    text(ind_TWstart, max([minY maxY])+abs(max([minY maxY])*0.35), ...
                        [num2str(round(PredErrEffect_TW.pval{1}{i_win}(i_elec),3))], ...
                        'FontSize',6);
                else
                    signTW = [signTW 0];                  
                end
            end
            
            %Scale y-axis to show sign. text if present
            if any(PredErrEffect_Sample.clusterstat{i_elec}.cluster_pval < param.alpha) || sum(signTW) > 0
                ylim([min([minY maxY])-+abs(max([minY maxY])*0.4) ...
                    max([minY maxY])+abs(max([minY maxY])*0.4)]); %With space for sign. info
            else
                ylim([min([minY maxY]) max([minY maxY])]); %Without space for sign. info
            end
            
            %Title
            title(['Elec: ' FuncInput_InputData.label{i_elec}])            
        end
    end
    
    %Add subtitle
    if strcmp(PredErrEffectType_label,'ComplexPredErrEffect')       
        sgtitle({[sub ' - ' FuncInput_DataType ' - ' FuncInput_ToneDur_text ...
            'msTD'], ['Complex Pred Effect (p34 ERF as a function of absolute p34-p*34 pitch difference, p < ' num2str(param.alpha) '), [' ...
            num2str(i_fig) '/' num2str(NumFig) ']']},'Interpreter','none')
        filename     = ['p34ERF_CPE_AllElecs_TWvsClusterStat_' sub '_' FuncInput_DataType '_' ...
            FuncInput_ToneDur_text 'sTD_' num2str(i_fig) '.png'];
    elseif strcmp(PredErrEffectType_label,'SimplePredErrEffect')      
        sgtitle({[sub ' - ' FuncInput_DataType ' - ' FuncInput_ToneDur_text ...
            'msTD'], ['Simple Pred Effect (p34 ERF as a function of absolute p33-p34 pitch difference, p < ' num2str(param.alpha) '), [' ...
            num2str(i_fig) '/' num2str(NumFig) ']']},'Interpreter','none')
        filename     = ['p34ERF_SPE_AllElecs_TWvsClusterStat_' sub '_' FuncInput_DataType '_' ...
            FuncInput_ToneDur_text 'sTD_' num2str(i_fig) '.png'];        
    end
    
    if save_poststepFigs == 1
        figfile      = [path_fig filename];
        saveas(gcf, [figfile], 'png'); %save png version
        close
    end
    
end

end