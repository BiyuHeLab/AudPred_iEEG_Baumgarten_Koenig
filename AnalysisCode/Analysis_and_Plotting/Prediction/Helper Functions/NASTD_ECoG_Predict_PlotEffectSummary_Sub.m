function NASTD_ECoG_Predict_PlotEffectSummary_Sub...
    (sub, FuncInput_EffectType, FuncInput_DataType, FuncInput_ToneDur_text,  ...
    param,...
    FuncInput_InputData, ...
    save_poststepFigs, paths_NASTD_ECoG)

%Aim: Plot summary of sample-wise prediction/prediction error results showing:
%1) Location of sign. electrodes on template brain
%2) Corresponding ERFs for electrodes showing a sign. effect
%3) Shading indicating time points with sign. effect

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')

%Add base dir and own script dir
addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

path_inputdata = [paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
    sub '/Data/Samplewise/'];

path_fig = ([paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
    sub '/Figs/SamplewiseSummary/' FuncInput_EffectType '/']);
if (~exist(path_fig, 'dir')); mkdir(path_fig); end

SampleFreq  = FuncInput_InputData.fsample;
nSensors    = size(FuncInput_InputData.trial{1},1);
nTrials     = length(FuncInput_InputData.trial);

%% 0.2) Determine subject-specific parameters
NASTD_ECoG_subjectinfo %load subject info file (var: si)
subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;

%% 1) Load effect data and clean MEG data (for electrode info)
load([path_inputdata sub '_PredEffectsCluster_' FuncInput_DataType '_' FuncInput_ToneDur_text 'sTD.mat'], ...
   FuncInput_EffectType , 'Data_p33', 'Data_p34', 'labels_loadedData');
CurrentEffect = eval(FuncInput_EffectType);

clear PredEffect SimplePredErrEffect ComplexPredErrEffect

% loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];%path to indiv preprocessed ECoG data
% tic
% disp(['Loading preprocessed data set for sub: ' sub])
% load(loadfile_ECoGpreprocdata);
% disp(['done loading in ' num2str(toc) ' sec'])

%% 3) find color limits for plotting based on cluster_statsum
lim_clusterstatsum = [0 0];
Index_signelecs = [];
Array_signelecs = zeros(1,nSensors);
Array_clusterstatsum = zeros(1,nSensors);
NumSginClusters = 0;

for i_elec = 1:nSensors
    dv_elec = CurrentEffect.clusterstat{i_elec}.cluster_statSum;   
    if min(dv_elec) < lim_clusterstatsum(1)
        lim_clusterstatsum(1) = min(dv_elec);
    end
    if max(dv_elec) > lim_clusterstatsum(2)
        lim_clusterstatsum(2) = max(dv_elec);
    end 
    
    if min(CurrentEffect.clusterstat{i_elec}.cluster_pval) < param.alpha
        Index_signelecs = [Index_signelecs i_elec];
        Array_signelecs(i_elec) = 1;
        index_maxclusterstat = ...
            find(abs(CurrentEffect.clusterstat{i_elec}.cluster_statSum) == ...
            max(abs(CurrentEffect.clusterstat{i_elec}.cluster_statSum)));
        Array_clusterstatsum(i_elec) = ...
            CurrentEffect.clusterstat{i_elec}.cluster_statSum(index_maxclusterstat);
        
        for i_cluster = 1:length(CurrentEffect.clusterstat{i_elec}.cluster_pval)
            if CurrentEffect.clusterstat{i_elec}.cluster_pval(i_cluster) < param.alpha
                NumSginClusters = NumSginClusters +1;
            end            
        end
    end
end

%% 4) Plot figure
%4.1 Set up subplot structure based on number of sign. electrodes
h = figure('visible','on'); %ensures figure doesn't pup during plotting
set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
%Add subtitle
sgtitle({['Summary sample-wise effects (cluster-corrected across samples): ' sub ]},'Interpreter','none')

DimSubplot = [ceil(length(Index_signelecs)/3)+1,3];

CounterSubplot = 2;
CounterSurfplot = 2;
SubplotPosition = [0 0 0 0];
ColorbarPosition = [0 0 0 0];
SizeFactor = 4;
    
%4.2 Plot information about current subject
switch FuncInput_EffectType
    case 'PredEffect'
        effect_text = {['Effect: Prediction'] ['(explain p33 activity by p*34)'] ['p < ' num2str(param.alpha)]};
    case 'SimplePredErrEffect'
        effect_text = {['Effect: Simple Prediction error'] ['(explain p34 activity by absolute p33-p34 difference)'] ['p < ' num2str(param.alpha)]};
    case 'ComplexPredErrEffect'
        effect_text = {['Effect: Complex Prediction error '] ['(explain p34 activity by absolute p*34-p34 difference)'] ['p < ' num2str(param.alpha)]};
end

subplot(ceil(length(Index_signelecs)/2)+1,3,1);

t = text(0,1,['Subject: ' sub]);
text(0,t.Position(2)*0.6,[effect_text],'Interpreter','none');
text(0,t.Position(2)*0.4,['Input data: ' FuncInput_DataType],'Interpreter','none');
text(0,t.Position(2)*0.2,['Tone Duration: ' FuncInput_ToneDur_text 'ms']);
text(0,t.Position(2)*0.1,['# sign. elecs / all elecs: ' ...
    num2str(length(Index_signelecs)) ' / ' num2str(nSensors)]);
text(0,t.Position(2)*0.0,['# sign. cluster: ' num2str(NumSginClusters)]);
axis off

%4.3 Plot surface plot (both hemispheres) in first row
%Set up data struct for plotting
clear plot_struct
plot_struct.dimord      = 'chan_time';
plot_struct.time        = 0;
plot_struct.label       = labels_loadedData;
plot_struct.avg         = Array_clusterstatsum';
plot_struct.sign_elecs  = logical(Array_signelecs)';

%Find index of each electrode label in cleaned data .elec struct and read
%out elec-specific coordinates
Index_Elecs2ChanPos = NaN(1,nSensors);
for i_elec = 1:length(labels_loadedData)
    Index_Elecs2ChanPos(1,i_elec) = ...
        find(strcmp(labels_loadedData{i_elec}, ...
        FuncInput_InputData.elec.label));
end
plot_struct.coords = FuncInput_InputData.elec.chanpos(Index_Elecs2ChanPos,:); %MNI coordinates for selected electrodes

plot_struct.clims = [max(abs(lim_clusterstatsum))*(-1) max(abs(lim_clusterstatsum))]; %free symmetric scaling 
plot_struct.chanSize = ...
    ones(1,length(plot_struct.avg))*SizeFactor; %electrode size (arbitrary)
plot_struct.cmap = 'jet';

%Plot surface with highlighted sign electrodes
sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_SubplotHemis...
    (plot_struct.coords, plot_struct.label, plot_struct.avg, plot_struct.sign_elecs,...
    plot_struct.chanSize, plot_struct.clims, plot_struct.cmap, ...
    DimSubplot, CounterSubplot, SubplotPosition, ColorbarPosition, []);
sp_handle_surf{1} = sp_handle_surf_temp.L;
sp_handle_surf{2} = sp_handle_surf_temp.R;

CounterSubplot = CounterSubplot + 2;

%Add colorbar
ColorbarPosition = [sp_handle_surf{2}.Position(1) - 0.02 sp_handle_surf{1}.Position(2)*1.02 0 0.75];
h = colorbar;
h.Position(1) = ColorbarPosition(1); %sets colorbar to the right
h.Position(2) = ColorbarPosition(2); %sets colorbar higher
h.Position(4) = h.Position(4)*ColorbarPosition(4); %makes colorbar shorter
h.Label.String = ['summed max cluster statistic'];
h.FontSize = 12;
caxis(plot_struct.clims)

%4.4 Add time-resolved signal traces for each sign. electrode

%Prediction effect: p33 signal as function of p*34
if strcmp(FuncInput_EffectType, 'PredEffect')
    
    %Separate trials into predp34-conditions
    label_Predp34 = [-1 0 1];% low, medium, high;
    for i_Predp34 = 1:length(label_Predp34)
        TrialFilt_Predp34 = FuncInput_InputData.behav.stim.predID == label_Predp34(i_Predp34);
        Index_Predp34{i_Predp34} = find(TrialFilt_Predp34 == 1);
    end
    clear TrialFilt_Predp34
    
    %Compute AVG + STD across trials per Predp34
    for i_Predp34 = 1:length(label_Predp34)
        Avgp33{i_Predp34} = mean(Data_p33(:,:,Index_Predp34{i_Predp34}),3);
        STDp33{i_Predp34} = std(Data_p33(:,:,Index_Predp34{i_Predp34}),1,3);
    end
    
    %Determine trace colors
    plot_param.color = ...
        {[0, 0.4470, 0.7410, 0.7], ...
        [0, 0.75, 0.75, 0.7],[0.8500, ...
        0.3250, 0.0980, 0.7]};
    
    %Plot ERF for sign  elecs
    for i_subplot = 1:length(Index_signelecs)
        subplot(DimSubplot(1),DimSubplot(2),CounterSubplot)
        CounterSubplot = CounterSubplot +1;
        i_elec = Index_signelecs(i_subplot);
        
        %Plot signal traces
        plot(1:length(Avgp33{1}(i_elec,:)),Avgp33{1}(i_elec,:),...
            'color',plot_param.color{1},'LineWidth', 2)
        hold on;
        plot(1:length(Avgp33{2}(i_elec,:)),Avgp33{2}(i_elec,:),...
            'color',plot_param.color{2},'LineWidth', 2)
        hold on;
        plot(1:length(Avgp33{3}(i_elec,:)),Avgp33{3}(i_elec,:),...
            'color',plot_param.color{3},'LineWidth', 2)
        axis tight
        
        %Change x-axis labeling to time [s]
        temp_xticks = xticks;
        temp_text = {};
        for i_xtick = 1:length(temp_xticks)
            temp_text = [temp_text num2str(round(temp_xticks(i_xtick)/SampleFreq,2))];
        end
        xticklabels(temp_text)
        xlabel('Time after p33-start [s]')
                
        grey  = [100 100 100]./255;

        %Highlight sign. samples/clusters
        maxY = max([max(Avgp33{1}(i_elec,:)), max(Avgp33{2}(i_elec,:)), max(Avgp33{3}(i_elec,:))]);
        minY = min([min(Avgp33{1}(i_elec,:)), min(Avgp33{2}(i_elec,:)), min(Avgp33{3}(i_elec,:))]);
        if ~isnan(CurrentEffect.clusterstat{i_elec}.cluster_pval)
            for i_cluster = 1:length(CurrentEffect.clusterstat{i_elec}.cluster_pval)
                if CurrentEffect.clusterstat{i_elec}.cluster_pval(i_cluster) < param.alpha
                    sign_samples = ...
                        CurrentEffect.clusterstat{i_elec}.cluster_timecourse ...
                        == i_cluster;
                    area(1:length(sign_samples), ...
                        sign_samples * (max([minY maxY])+abs(max([minY maxY])*0.05)),...
                        'basevalue',0,'FaceColor',grey,'FaceAlpha', 0.3, 'LineStyle','none');
                    if min([minY maxY]) < 0
                        area(1:length(sign_samples), ...
                            sign_samples * (min([minY maxY])-abs(min([minY maxY])*0.05)),...
                            'basevalue',0,'FaceColor',grey,'FaceAlpha', 0.3, 'LineStyle','none'); 
                    end
%                     plot(sign_samples, ...
%                         ones(length(sign_samples)) * max([minY maxY])+abs(max([minY maxY])*0.05),...
%                         'linewidth',7,'Color',[175 175 175]./255)
                    shading_startsample = find(sign_samples == 1);
                    text(shading_startsample(1),max([minY maxY])+abs(max([minY maxY])*0.025), ...
                        [' p = ' num2str(round(...
                        CurrentEffect.clusterstat{i_elec}.cluster_pval(i_cluster),3))], ...
                        'FontSize',10,'FontWeight', 'bold');
                end
            end
        end
        
        %Scale y-axis to show sign. text if present
        ylim([min([minY maxY])-+abs(max([minY maxY])*0.05) ...
            max([minY maxY])+abs(max([minY maxY])*0.05)]); %With space for sign. info
        
        %Title
        title(['Elec: ' FuncInput_InputData.label{i_elec}])
    end
        
%Prediction error effect
elseif strcmp(FuncInput_EffectType, 'SimplePredErrEffect') || strcmp(FuncInput_EffectType, 'ComplexPredErrEffect')

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
    if strcmp(FuncInput_EffectType,'ComplexPredErrEffect')
        unique_predErrPitchDiff = unique(round(LogFreq_ComplexPredError,4));
        sortedunique_predErrPitchDiff = unique(sort(abs(unique_predErrPitchDiff)));   
    elseif strcmp(FuncInput_EffectType,'SimplePredErrEffect')
        unique_predErrPitchDiff = unique(round(LogFreq_SimplePredError,4));
        sortedunique_predErrPitchDiff = unique(sort(abs(unique_predErrPitchDiff)));   
    end
    
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
        if strcmp(FuncInput_EffectType,'ComplexPredErrEffect')
            for i_pitchperbin = 1:PitchperBin
                Index_PredError{i_numbins} = [Index_PredError{i_numbins} ...
                    find(abs(round(LogFreq_ComplexPredError,4)) == ...
                    PredErr_PitchDiff(i_numbins,i_pitchperbin))];
            end
       elseif strcmp(FuncInput_EffectType,'SimplePredErrEffect')
            for i_pitchperbin = 1:PitchperBin
                Index_PredError{i_numbins} = [Index_PredError{i_numbins} ...
                    find(abs(round(LogFreq_SimplePredError,4)) == ...
                    PredErr_PitchDiff(i_numbins,i_pitchperbin))];
            end       
        end     
    end
    
    %Compute AVG + STD across trials per bin
    for i_numbins = 1:NumBins
        Avgp34{i_numbins} = mean(Data_p34(:,:,Index_PredError{i_numbins}),3);
        STDp34{i_numbins} = std(Data_p34(:,:,Index_PredError{i_numbins}),1,3);
    end
    
    %Plot ERF for all selected elecs
    plot_param.color = parula;
    
    for i_subplot = 1:length(Index_signelecs)
        subplot(DimSubplot(1),DimSubplot(2),CounterSubplot)
        CounterSubplot = CounterSubplot +1;
        i_elec = Index_signelecs(i_subplot);
        
        %Plot signal traces
        for i_numbins = 1:NumBins
            plot(1:length(Avgp34{i_numbins}(i_elec,:)),Avgp34{i_numbins}(i_elec,:),...
                'color',plot_param.color(i_numbins*(floor(length(plot_param.color)/NumBins)),:), ...
                'LineWidth', 2)
            hold on;
        end
        axis tight
        
        %Change x-axis labeling to time [s]
        temp_xticks = xticks;
        temp_text = {};
        for i_xtick = 1:length(temp_xticks)
            temp_text = [temp_text num2str(round(temp_xticks(i_xtick)/SampleFreq,2))];
        end
        xticklabels(temp_text)
        xlabel('Time after p34-start [s]')
                    
        %Highlight sign. samples/clusters
        maxY = [];
        minY = [];
        for i_numbins = 1:NumBins
            maxY = [maxY max(Avgp34{i_numbins}(i_elec,:))];
            minY = [minY min(Avgp34{i_numbins}(i_elec,:))];
        end
        maxY = max(maxY);
        minY = min(minY);
        grey  = [100 100 100]./255;

        if ~isnan(CurrentEffect.clusterstat{i_elec}.cluster_pval)
            for i_cluster = 1:length(CurrentEffect.clusterstat{i_elec}.cluster_pval)
                if CurrentEffect.clusterstat{i_elec}.cluster_pval(i_cluster) < param.alpha
                    hold on;
                    sign_samples = ...
                        CurrentEffect.clusterstat{i_elec}.cluster_timecourse ...
                        == i_cluster;
                    area(1:length(sign_samples), ...
                        sign_samples * (max([minY maxY])+abs(max([minY maxY])*0.05)),...
                        'basevalue',0,'FaceColor',grey,'FaceAlpha', 0.3, 'LineStyle','none');
                    if min([minY maxY]) < 0
                        area(1:length(sign_samples), ...
                            sign_samples * (min([minY maxY])-abs(min([minY maxY])*0.05)),...
                            'basevalue',0,'FaceColor',grey,'FaceAlpha', 0.3, 'LineStyle','none'); 
                    end
%                     plot(sign_samples, ...
%                         ones(length(sign_samples)) * max([minY maxY])+abs(max([minY maxY])*0.05),...
%                         'linewidth',7,'Color',[175 175 175]./255)
                    shading_startsample = find(sign_samples == 1);
                    text(shading_startsample(1),max([minY maxY])+abs(max([minY maxY])*0.025), ...
                        [' p = ' num2str(round(...
                        CurrentEffect.clusterstat{i_elec}.cluster_pval(i_cluster),3))], ...
                        'FontSize',10,'FontWeight', 'bold');
                end
            end
        end
        
        %Scale y-axis to show sign. text if present
        ylim([min([minY maxY])-+abs(max([minY maxY])*0.05) ...
            max([minY maxY])+abs(max([minY maxY])*0.05)]); %With space for sign. info
        
        %Title
        title(['Elec: ' FuncInput_InputData.label{i_elec}])
    end
end
            

 
%% 5. Save figure
if save_poststepFigs == 1
    filename     = [sub '_Summary_' FuncInput_EffectType '_' FuncInput_DataType ...
        '_' FuncInput_ToneDur_text 'msTD_p' num2str(param.alpha) '.png'];
    figfile      = [path_fig filename];
    saveas(gcf, figfile, 'png'); %save png version
    close all;
end

end