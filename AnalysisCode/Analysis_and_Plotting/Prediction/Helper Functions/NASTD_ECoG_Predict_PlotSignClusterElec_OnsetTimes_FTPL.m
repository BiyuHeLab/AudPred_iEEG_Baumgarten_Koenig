function Onset_times_ms = NASTD_ECoG_Predict_PlotSignClusterElec_OnsetTimes_FTPL...
    (subs, ...
    FuncInput_EffectType, FuncInput_DataType, FuncInput_ToneDur_text, ...
    param, ...
    save_poststepFigs, paths_NASTD_ECoG)

%Aim: Plot electrodes showing sign. cluster-corrected prediction effects
%per input data & TD condition aggregated over subjects and tone durations
%(TD) and projected on 1 Hemisphere.
%Based on sample-wise analysis, thus no time windows
%but instead entire analysis window.
%Plot Figure with 5 subfigures:
% - clusterstat
% - p-value derivative
% - onset
% - offset
% - duration

% set(0, 'DefaultFigureVisible', 'on');

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
%Add base dir and own script dir
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

if param.FDRcorrect == 1
    FDR_label = 'FDRcorr';
else
    FDR_label = 'uncorr';
end

path_fig = ([paths_NASTD_ECoG.ECoGdata_Prediction ...
    '/PredEffects/Allsub_n' num2str(length(subs)) ...
    '/Figs/PredEffects_Surf/' FuncInput_EffectType ...
    '/ClusterCorr/' param.ElecSelect '/' FDR_label '/']);
if (~exist(path_fig, 'dir')); mkdir(path_fig); end


%% 1) Load prediction effect data and aggregate relevant info across subjects
clear Param2plot
coords_allsub                       = [];
labels_allsub                       = [];
anatlabels_allsub                   = [];
Param2plot.all_subs.sub_index       = [];
Param2plot.all_subs.TD_index        = [];
usedElecs_chanposIndex              = [];

for i_sub = 1:length(subs)
    
    tic
    sub = subs{i_sub};
    disp(['-- Loading data for sub: ' sub ' --'])
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;
    
    for i_TD = 1:length(FuncInput_ToneDur_text)
        %Load prediction data and select current effect
        path_inputdata = [paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
            sub '/Data/Samplewise/'];
        load([path_inputdata sub '_PredEffectsFTPL_' ...
            FuncInput_DataType '_' ...
            FuncInput_ToneDur_text{i_TD} 'sTD_NoClusterCorr.mat'], ...
            FuncInput_EffectType , 'labels_loadedData');
        CurrentEffect = eval(FuncInput_EffectType);
        clear PredErrFTPL
        
        %Also load ECoG preproc data for channel labels and position
        loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];
        load(loadfile_ECoGpreprocdata);
        
        SampleFreq              = DataClean_AllTrials.fsample;
        nSensors_all            = size(CurrentEffect.stats,2);
        nSamples(i_TD)          = size(CurrentEffect.pval_FDR{1},2);
        
        %Load stimulus correlation data and select current effect
        if strcmp(param.ElecSelect, 'StimCorr')
            path_inputdata = [paths_NASTD_ECoG.ECoGdata_StimCorr '/' sub '/Data/'];
            load([path_inputdata sub '_StimCorr_' FuncInput_DataType '_' ...
                FuncInput_ToneDur_text{i_TD} 'sTD.mat'], ...
                'corr_ttest', 'SensorLabels');
            filt_signelecs_StimCorr = corr_ttest.p < 0.05;%only sign. stim corr elecs
        else
            filt_signelecs_StimCorr  = true(length(labels_loadedData),1);  %all elecs
        end
        ind_signelecs_StimCorr  = find(filt_signelecs_StimCorr);
        nSensors_sel            = sum(filt_signelecs_StimCorr);
        
        %Read electrode labels, coordinates, and anatomical labels for all
        %sign. StimCorr electrodes and aggregate them across subjects
        for i_elec = 1:length(ind_signelecs_StimCorr)
            usedElecs_chanposIndex(i_elec,1) = ...
                find(strcmp(labels_loadedData{ind_signelecs_StimCorr(i_elec)}, ...
                DataClean_AllTrials.elec.label));
        end
        coords_sub{i_sub, i_TD}   = ...
            DataClean_AllTrials.elec.chanpos(usedElecs_chanposIndex,:);
        coords_allsub       = [coords_allsub; coords_sub{i_sub, i_TD}];
        anatlabels_allsub   = ...
            [anatlabels_allsub; DataClean_AllTrials.elec.T1AnatLabel(usedElecs_chanposIndex)];
        for i_elec = 1:length(ind_signelecs_StimCorr)
            labels_allsub{end+1,1} = [labels_loadedData{ind_signelecs_StimCorr(i_elec)} ' ' sub];
        end
        
        %Store cluster information from all
        Param2plot.per_sub.pval{i_sub, i_TD}                  = nan(nSensors_sel,1);
        Param2plot.per_sub.pval_FDR{i_sub, i_TD}              = nan(nSensors_sel,1);
        Param2plot.per_sub.effect_onset{i_sub, i_TD}          = nan(nSensors_sel,1);
        
        for i_elec = 1:nSensors_sel
            %For each electrode, read out corresponding information
            sig_timepoints = find(CurrentEffect.pval_FDR{ind_signelecs_StimCorr(i_elec)} < 0.05);
            
            if ~isempty(sig_timepoints) %only take the information for the first significant timepoint
                i_t = sig_timepoints(1);
                Param2plot.per_sub.pval{i_sub, i_TD}(i_elec,1) = CurrentEffect.pval{ind_signelecs_StimCorr(i_elec)}(1,i_t);
                Param2plot.per_sub.pval_FDR{i_sub, i_TD}(i_elec,1) = CurrentEffect.pval_FDR{ind_signelecs_StimCorr(i_elec)}(1,i_t);
                Param2plot.per_sub.effect_onset{i_sub, i_TD}(i_elec,1) = i_t / nSamples(i_TD);
            end
        end
        
        %Restrict cluster data to current electrode selection (StimCorr or All)
        for i_elec = 1:nSensors_sel
            if filt_signelecs_StimCorr(ind_signelecs_StimCorr(i_elec)) == 0 %Restrict p-values to selected elecs (StimCorr or All)
                Param2plot.per_sub.pval{i_sub, i_TD}(i_elec, :) = NaN;
            end
        end
        
        %Compute derivative of p-value to estimate strength of effect
%         Param2plot.per_sub.pval_derivative{i_sub, i_TD} = ...
%             -(log10(Param2plot.per_sub.pval{i_sub, i_TD}));
        
        %Data aggregation over subjects
        %Array to differentiate subject entries
        Param2plot.all_subs.sub_index = ...
            [Param2plot.all_subs.sub_index; ...
            ones(length(coords_sub{i_sub, i_TD}),1)*i_sub];        
        %Array to differentiate TD entries for each electrode
        for i_elec = 1:nSensors_all
            Param2plot.all_subs.TD_index    = ...
                [Param2plot.all_subs.TD_index; ...
                FuncInput_ToneDur_text{i_TD} 's'];
        end
        
        %Clusterinfo aggregation over subjects (Dim: elec * cluster)
        if i_sub == 1 && i_TD == 1
            Param2plot.all_subs.pval            = ...
                Param2plot.per_sub.pval{i_sub, i_TD};
            Param2plot.all_subs.pval_FDR        = ...
                Param2plot.per_sub.pval_FDR{i_sub, i_TD};      
            Param2plot.all_subs.effect_onset    = ...
                Param2plot.per_sub.effect_onset{i_sub, i_TD};
        else
            Param2plot.all_subs.pval            = ...
                [Param2plot.all_subs.pval; ...
                Param2plot.per_sub.pval{i_sub, i_TD}];
            Param2plot.all_subs.pval_FDR        = ...
                [Param2plot.all_subs.pval_FDR; ...
                Param2plot.per_sub.pval_FDR{i_sub, i_TD}];
            Param2plot.all_subs.effect_onset    = ...
                [Param2plot.all_subs.effect_onset; ...
                Param2plot.per_sub.effect_onset{i_sub, i_TD}];
        end
        
        %Cleanup
        usedElecs_chanposIndex = [];
        clear CurrentEffect labels_loadedData corr_ttest DataClean_AllTrials temp*
    end
    disp([' -- done loading in ' num2str(toc) ' sec --'])
end

%% 2) Determine sign. electrodes
SignElecs.array         = any(Param2plot.all_subs.pval < param.pval_plotting,2); %1D filter denoting sign. elecs
SignElecs.index         = find(SignElecs.array);
SignElecs.num_elecs     = length(SignElecs.index);
SignElecs.num_cluster   = sum(sum(Param2plot.all_subs.pval < param.pval_plotting));

%Determine labels of sign. elecs and add number to subtitle
SignElecs.labels        = labels_allsub(SignElecs.index);
SignElecs.anatlabels    = anatlabels_allsub(SignElecs.index);
SignElecs.fulllabels    = [];
for i_elec = 1:length(SignElecs.labels)
    SignElecs.fulllabels{i_elec,1} = ...
        [SignElecs.labels{i_elec}, ' ', SignElecs.anatlabels{i_elec}];
end

if ~isempty(SignElecs.labels)
    sign_title = ...
        [num2str(SignElecs.num_elecs)...
        ' / ' num2str(length(labels_allsub)) ' sign. elecs']; %1 line
else
    sign_title = ...
        ['0 / ' num2str(length(labels_allsub)) ' sign. elecs']; %1 line
end

%Create electrode labels for plotting
for i_elec = 1:length(labels_allsub)%No electrode labels
    labels_plotting_empty{i_elec} = '';
end
counter_elecs           = 0;
labels_plotting_number  = labels_plotting_empty;
for i_elec = SignElecs.index'
    counter_elecs = counter_elecs +1;
    labels_plotting_number{i_elec} = num2str(counter_elecs);
end

for i_elec = 1:length(SignElecs.index)
%     SignElecs.onset_minpcluster(i_elec,1) = ...
%         Param2plot.all_subs.effect_onset(SignElecs.index(i_elec),1); %
%         Commented out by LK on Nov 20, 2024 and changed to the code below
    SignElecs.onset_mincluster(i_elec,1) = ...
        Param2plot.all_subs.effect_onset(SignElecs.index(i_elec), :);
end

% Add information about which Tone Duration and convert the effect_onset
% into ms
SignElecs.TD = Param2plot.all_subs.TD_index(SignElecs.index,:);
Onset_times_ms = nan(length(SignElecs.index),2); % electrodes w/ sig. clusters x 2 tone durations
for i_elec = 1:length(SignElecs.index)
    if strcmp(SignElecs.TD(i_elec, :), '0.2s')
        Onset_times_ms(i_elec,1) = SignElecs.onset_mincluster(i_elec,1) * 200;
    else
        Onset_times_ms(i_elec,2) = SignElecs.onset_mincluster(i_elec,1) * 400;
    end
end
Index_TD1 = SignElecs.index(~isnan(Onset_times_ms(:,1)));
OnsetTimes_TD1 = Onset_times_ms(~isnan(Onset_times_ms(:,1)),1);
Index_TD2 = SignElecs.index(isnan(Onset_times_ms(:,1)));
OnsetTimes_TD2 = Onset_times_ms(~isnan(Onset_times_ms(:,2)),2);
SignElecs_TD1 = false(1898, 1);
SignElecs_TD1(Index_TD1) = 1;
SignElecs_TD2 = false(1898, 1);
SignElecs_TD2(Index_TD2) = 1;

%% 3) Prepare & Plot Figure 1: Clusterstat & p-value derivative
if ~isempty(SignElecs.index)
    
    %3.1 Prepare figure
    h = figure('visible','on'); %ensures figure doesn't pop up during plotting
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    set(gcf,'renderer','opengl');
    
    DimSubplot      = [1,2];
    CounterSubplot  = 1;
    CounterSurfplot = 1;
    SizeFactor      = 3;
    SubplotPosition = [0 -0.1 0 0];
    
    %3.2 Adjust header/title
    %3.2 Adjust header/title
    switch FuncInput_EffectType
        case 'PredErrFTPL'
            effect_text = ['Effect: Complex Prediction error (explain p34 activity by FTPL rating); p < ' num2str(param.pval_plotting) ' ' FDR_label];
    end
    
    Fig_title = {['Group level (n = ' num2str(length(subs)) ')'] ...
        [effect_text] ...
        ['Input data: ' FuncInput_DataType ', Pooled TD']} ;
    sgtitle(Fig_title,'FontSize',20,'Interpreter','none')
    
    %3.2 Subplot 1: p-val FDR
    for i_elec = 1:length(Param2plot.all_subs.pval_FDR)
        PlotInput(i_elec,1) = max(Param2plot.all_subs.pval_FDR(i_elec,:));
        %select custer based on maximal clusterstat for each elctrode
    end
    
    % Compute derivative of pval FDR for visual comprehension 
    PlotInput = -(log10(PlotInput));
    
    %Colorlim
    dv_max                  = max(PlotInput);
    dv_min                  = min(PlotInput);
    clims                   = [1.3 3.3];
    Label_Colorbar          = 'p-val FDR';
    
    pvals_ticks = [0.05, 0.01, 0.001];
    logticks = -log10(pvals_ticks);
    caxis(clims); % Assuming p-values only positive, adjust if needed
    
    %Data struct
    dat.dimord          = 'chan_time';
    dat.time            = 0;
    dat.label           = labels_allsub;
    dat.avg             = PlotInput;
    dat.sign_elecs      = SignElecs.array;
    dat.textcolor_rgb   = [1 0 0];
    
    chanSize    = ones(1,length(dat.avg))*SizeFactor; %electrode size (arbitrary)
    cmap        = 'parula';
    
    %Project all electrodes on one hemisphere,
    coords_allsub(:,1) = abs(coords_allsub(:,1)) * -1;
    
    sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH...
        (coords_allsub, labels_plotting_number, dat.avg, SignElecs.array,...
        chanSize, clims, cmap, dat.textcolor_rgb, ...
        DimSubplot, CounterSubplot, SubplotPosition, [], []);
    
    sp_handle_surf{CounterSurfplot} = sp_handle_surf_temp.L;
    CounterSubplot                  = CounterSubplot + 1;
    
    %Colorbar
    ColorbarPosition    = [0.1 0.58 0 0.5]; %Left of surface
    h                   = colorbar;
    h.Position(1)       = ColorbarPosition(1); %sets colorbar to the right
    h.Position(2)       = ColorbarPosition(2); %sets colorbar higher
    h.Position(4)       = h.Position(4)*ColorbarPosition(4); %makes colorbar shorter
    h.Label.String      = Label_Colorbar;
    h.Ticks = logticks;
    h.TickLabels = cellstr(num2str(pvals_ticks', '%.3f'));
    h.FontSize          = 16;
    caxis(clims)
    
    %% 4 Save figure
    if save_poststepFigs == 1
        filename     = ['Surf1HSignElec' param.ElecSelect ...
            '_Allsubn' num2str(length(subs)) '_' ...
            FuncInput_EffectType '_' FuncInput_DataType ...
            '_AllTD_Stat' FDR_label '_ONSET_ms.png'];
        figfile      = [path_fig filename];
        saveas(gcf, figfile, 'png'); %save png version
        close all;
    end
    
end

end