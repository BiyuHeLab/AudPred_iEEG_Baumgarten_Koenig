function NASTD_ECoG_Predict_PlotFDRSignClusterElec_AllSubTD...
    (subs, ...
    FuncInput_EffectType, FuncInput_DataType, FuncInput_ToneDur_text, ...
    param, ...
    save_poststepFigs, paths_NASTD_ECoG)

%Aim: Plot electrodes showing sign. cluster-corrected prediction effects
%per input data & TD condition aggregated over subjects and tone durations 
%(TD) and projected on 1 Hemisphere. 
%Use only electrodes that show sign. Stim-Corr effect.

%Based on sample-wise analysis, thus no time windows in analysis.
%Time windows added for movie plotting

% set(0, 'DefaultFigureVisible', 'on');

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
%Add base dir and own script dir
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

path_fig = ([paths_NASTD_ECoG.ECoGdata_Prediction ...
    '/PredEffects/Allsub_n' num2str(length(subs)) ...
    '/Figs/PredEffects_Surf/' FuncInput_EffectType '/ClusterCorr_StimCorrElecOnly/']);
if (~exist(path_fig, 'dir')); mkdir(path_fig); end

%% 1) Load prediction effect data and aggregate relevant info across subjects
clear Param2plot
coords_allsub                       = [];
labels_allsub                       = [];
anatlabels_allsub                   = [];
Param2plot.all_subs.sub_index       = [];
Param2plot.all_subs.TD_index        = [];
Param2plot.all_subs.filt_StimCorr   = [];
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
        load([path_inputdata sub '_PredEffectsCluster_' ...
            FuncInput_DataType '_' ...
            FuncInput_ToneDur_text{i_TD} 'sTD.mat'], ...
            FuncInput_EffectType , 'labels_loadedData');
        CurrentEffect = eval(FuncInput_EffectType);
        clear PredEffect SimplePredErrEffect ComplexPredErrEffect        
        
        %Also load ECoG preproc data for channel labels and position
        loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];
        load(loadfile_ECoGpreprocdata);
        
        SampleFreq              = DataClean_AllTrials.fsample;
        nSensors_all            = size(CurrentEffect.stats,2);
        nSamples(i_TD)          = size(CurrentEffect.clusterstat{1}.cluster_timecourse,2);
        
        %Load stimulus correlation data and select current effect
        path_inputdata = [paths_NASTD_ECoG.ECoGdata_StimCorr '/' sub '/Data/'];
        load([path_inputdata sub '_StimCorr_' FuncInput_DataType '_' ...
            FuncInput_ToneDur_text{i_TD} 'sTD.mat'], ...
            'corr_ttest', 'SensorLabels'); 
        
        if strcmp(param.ElecSelect, 'StimCorr')
            filt_signelecs_StimCorr = corr_ttest.p < 0.05;%only sign. stim corr elecs
        else
            filt_signelecs_StimCorr  = true(length(corr_ttest.p),1);  %all elecs
        end
        ind_signelecs_StimCorr  = find(filt_signelecs_StimCorr);
        nSensors_sel            = sum(filt_signelecs_StimCorr);
     
        %Read electrode labels, coordinates, and anatomical labels for all 
        %sign. StimCorr electrodes and aggregate them across subjects
        for i_elec = 1:length(ind_signelecs_StimCorr)
            usedElecs_chanposIndex(i_elec,1) = ...
                find(strcmp(SensorLabels{ind_signelecs_StimCorr(i_elec)}, ...
                DataClean_AllTrials.elec.label));
        end
        coords_sub{i_sub, i_TD}   = ...
            DataClean_AllTrials.elec.chanpos(usedElecs_chanposIndex,:);    
        coords_allsub       = [coords_allsub; coords_sub{i_sub, i_TD}];
        anatlabels_allsub   = ...
            [anatlabels_allsub; DataClean_AllTrials.elec.T1AnatLabel(usedElecs_chanposIndex)]; 
        for i_elec = 1:length(ind_signelecs_StimCorr)
            labels_allsub{end+1,1} = [SensorLabels{ind_signelecs_StimCorr(i_elec)} ' ' sub];        
        end   
        
        %Store clusterstat and p-values for selected sign. StimCorr elecs
        Param2plot.per_sub.clusterstat{i_sub, i_TD}       = nan(10, nSensors_sel);
        Param2plot.per_sub.pval{i_sub, i_TD}              = nan(10, nSensors_sel);
        Param2plot.per_sub.effect_onset{i_sub, i_TD}      = nan(10, nSensors_sel);
        Param2plot.per_sub.effect_offset{i_sub, i_TD}     = nan(10, nSensors_sel);
        Param2plot.per_sub.effect_duration{i_sub, i_TD}   = nan(10, nSensors_sel);
        
        
        %Note: Add minimal p-value cluster selection - currently across all
        %clusters
        for i_elec = 1:nSensors_sel
            for i_cluster = 1:length(CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_pval)
                if CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_pval(i_cluster) < param.pval_plotting
                    Param2plot.per_sub.clusterstat{i_sub, i_TD}(i_cluster, i_elec) = ...
                        CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_statSum(i_cluster);
                    Param2plot.per_sub.pval{i_sub, i_TD}(i_cluster, i_elec) = ...
                        CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_pval(i_cluster);
                    Param2plot.per_sub.effect_onset{i_sub, i_TD}(i_cluster, i_elec) = ... %Relative onset to combine both TD
                        CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_samples{i_cluster}(1) / nSamples(i_TD);
                    Param2plot.per_sub.effect_offset{i_sub, i_TD}(i_cluster, i_elec) = ...%Relative offset to combine both TD
                        CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_samples{i_cluster}(end) / nSamples(i_TD);
                    Param2plot.per_sub.effect_duration{i_sub, i_TD}(i_cluster, i_elec) = ...%Relative duration to combine both TD
                        length(CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_samples{i_cluster}) / nSamples(i_TD);
                end
            end                
        end
        
        %Derivative of p-value to estimate strength of effect
        Param2plot.per_sub.pval_derivative{i_sub, i_TD} = ...
            -(log10(Param2plot.per_sub.pval{i_sub, i_TD}));        
        
        %Select cluster with min p-value and restrict analysis to elecs
        %showing a sign. StimCorr effect
        cluster_pval_perelec = nan(nSensors_sel,1);
        for i_elec = 1:nSensors_sel
            if any(~isnan(Param2plot.per_sub.pval{i_sub, i_TD}(:,i_elec))) %Read out min p-value
                cluster_pval_perelec(i_elec) = ...
                    min(Param2plot.per_sub.pval{i_sub, i_TD}(:,i_elec));
            end
            if filt_signelecs_StimCorr(ind_signelecs_StimCorr(i_elec)) == 0 %Restrict p-values to StimCorr sign elecs
                cluster_pval_perelec(i_elec) = NaN;
            end
        end
        
        %Derivative of p-value to estimate strength of effect
        Param2plot.per_sub.pval_derivative{i_sub, i_TD} = ...
            -(log10(cluster_pval_perelec));
        
        %Perform FDR_correction on cluster-p-values across all sign. StimCorr electrodes        
        Param2plot.per_sub.pval{i_sub, i_TD} = nan(nSensors_sel,1);
        if param.FDRcorrect == 1
            FDR_label = 'FDRcorr';
            [~,~,~,cluster_pvalFDR_perelec] = fdr_bh(cluster_pval_perelec, param.pval_FDR, 'pdep','no');
            index_signelecsFDR = cluster_pvalFDR_perelec < param.pval_plotting;
            if any(index_signelecsFDR)
                for i_signelec = index_signelecsFDR
                    Param2plot.per_sub.pval{i_sub, i_TD}(i_signelec) = ...
                        cluster_pvalFDR_perelec(i_signelec);
                end
            end
        else
            FDR_label = 'uncorr';
            Param2plot.per_sub.pval{i_sub, i_TD} = cluster_pval_perelec;
        end        
               
        %Array to differentiate subject and TD entries for each electrode
        Param2plot.all_subs.sub_index   = ...
            [Param2plot.all_subs.sub_index; ...
            ones(length(coords_sub{i_sub}),1)*i_sub];
        for i_elec = 1:nSensors_all
            Param2plot.all_subs.TD_index    = ...
                [Param2plot.all_subs.TD_index; ...
                FuncInput_ToneDur_text{i_TD} 's'];
        end
               
        %clusterstat and p-value aggregation over subjects
        if i_sub == 1 && i_TD == 1
            Param2plot.all_subs.clusterstat     = ...
                Param2plot.per_sub.clusterstat{i_sub, i_TD}';
            Param2plot.all_subs.pval            = ...
                Param2plot.per_sub.pval{i_sub, i_TD};
            Param2plot.all_subs.pval_derivative = ...
                Param2plot.per_sub.pval_derivative{i_sub, i_TD};
            Param2plot.all_subs.effect_onset    = ...
                Param2plot.per_sub.effect_onset{i_sub, i_TD}';
            Param2plot.all_subs.effect_offset   = ...
                Param2plot.per_sub.effect_offset{i_sub, i_TD}';
            Param2plot.all_subs.effect_duration = ...
                Param2plot.per_sub.effect_duration{i_sub, i_TD}';
        else
            Param2plot.all_subs.clusterstat     = ...
                [Param2plot.all_subs.clusterstat; Param2plot.per_sub.clusterstat{i_sub, i_TD}'];
            Param2plot.all_subs.pval            = ...
                [Param2plot.all_subs.pval; Param2plot.per_sub.pval{i_sub, i_TD}];
            Param2plot.all_subs.pval_derivative = ...
                [Param2plot.all_subs.pval_derivative; Param2plot.per_sub.pval_derivative{i_sub, i_TD}];
            Param2plot.all_subs.effect_onset    = ...
                [Param2plot.all_subs.effect_onset; Param2plot.per_sub.effect_onset{i_sub, i_TD}'];
            Param2plot.all_subs.effect_offset   = ...
                [Param2plot.all_subs.effect_offset; Param2plot.per_sub.effect_offset{i_sub, i_TD}'];
            Param2plot.all_subs.effect_duration = ...
                [Param2plot.all_subs.effect_duration; Param2plot.per_sub.effect_duration{i_sub, i_TD}'];
        end
        
        %Cleanup
        usedElecs_chanposIndex = [];
        clear CurrentEffect labels_loadedData corr_ttest DataClean_AllTrials
    end
    disp([' -- done loading in ' num2str(toc) ' sec --'])      
end

%% 2) Determine sign. electrodes
temp_indexsignelec = nan(length(Param2plot.all_subs.pval),1);
for i_elec = 1:length(Param2plot.all_subs.pval)
    temp_indexsignelec(i_elec) = min(Param2plot.all_subs.pval(i_elec,:))';
end
SignElecs.array         = temp_indexsignelec < param.pval_plotting;
SignElecs.index         = find(SignElecs.array);
SignElecs.num_elecs     = length(SignElecs.index);
SignElecs.num_cluster   = sum(sum(Param2plot.all_subs.pval < param.pval_plotting));

%Determine labels of sign. elecs and add number to subtitle
SignElecs.labels        = labels_allsub(SignElecs.index);
SignElecs.anatlabels    = anatlabels_allsub(SignElecs.index);
SignElecs.fulllabels    = [];
for i_elec = 1:length(SignElecs.labels)
    SignElecs.fulllabels{i_elec,1} = ...
        [SignElecs.labels{i_elec}, ' ', ...
        Param2plot.all_subs.TD_index(SignElecs.index(i_elec),:), ' ', ...
        SignElecs.anatlabels{i_elec}];
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

%% 3) Prepare & Plot Figure 1: Clusterstat & p-value derivative
%3.1 Prepare figure
h = figure('visible','on'); %ensures figure doesn't pop up during plotting
set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
set(gcf,'renderer','opengl');

DimSubplot      = [2,3];
CounterSubplot  = 1;
CounterSurfplot = 1;
SizeFactor      = 4;
SubplotPosition = [0 -0.1 0 0];

%3.2 Adjust header/title
switch FuncInput_EffectType
    case 'PredEffect'
        effect_text = ['Effect: Prediction (explain p33 activity by p*34); p < ' num2str(param.pval_plotting) '(' FDR_label ')'];
    case 'SimplePredErrEffect'
        effect_text = ['Effect: Simple Prediction error (explain p34 activity by absolute p33-p34 difference); p < ' num2str(param.pval_plotting) '(' FDR_label ')'];
    case 'ComplexPredErrEffect'
        effect_text = ['Effect: Complex Prediction error (explain p34 activity by absolute p*34-p34 difference); p < ' num2str(param.pval_plotting) '(' FDR_label ')'];
end

Fig_title = {['Group level (n = ' num2str(length(subs)) ')'] ...
    [effect_text] ...
    ['Input data: ' FuncInput_DataType ', Pooled TD']} ;
sgtitle(Fig_title,'FontSize',20,'Interpreter','none')

%3.2 Subplot 1: Clusterstat
for i_elec = 1:length(Param2plot.all_subs.clusterstat)
    PlotInput(i_elec,1) = max(Param2plot.all_subs.clusterstat(i_elec,:)); 
    %select custer based on maximal clusterstat for each elctrode
end

%Colorlim
dv_max                  = max(max(Param2plot.all_subs.clusterstat));
dv_min                  = min(min(Param2plot.all_subs.clusterstat));
dv_absmax               = max([abs(dv_max) abs(dv_min)]);
clims_maxclusterstat    = [-dv_absmax dv_absmax];
clims                   = clims_maxclusterstat;
Label_Colorbar          = 'max sum cluster statistic';

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
h.FontSize          = 16;
caxis(clims)

%3.3 Subplot 2: p-value derivative
for i_elec = 1:length(Param2plot.all_subs.pval_derivative)
    PlotInput(i_elec,1) = max(Param2plot.all_subs.pval_derivative(i_elec,:));
end

%Colorlim
clims = [0 3];
Label_Colorbar = '- log10 (cluster p-value) ';

%prepare struct
dat.avg         = PlotInput;

%Plot surface
sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH...
    (coords_allsub, labels_plotting_empty, dat.avg, SignElecs.array,...
    chanSize, clims, cmap, dat.textcolor_rgb, ...
    DimSubplot, CounterSubplot, SubplotPosition, [], []);
sp_handle_surf{CounterSurfplot} = sp_handle_surf_temp.L;
CounterSubplot                  = CounterSubplot + 1;

%Colorbar
ColorbarPosition    = [0.39 0.58 0 0.5]; %Left of surface
h = colorbar;
h.Position(1)       = ColorbarPosition(1); %sets colorbar to the right
h.Position(2)       = ColorbarPosition(2); %sets colorbar higher
h.Position(4)       = h.Position(4)*ColorbarPosition(4); %makes colorbar shorter
h.Label.String      = Label_Colorbar;
h.FontSize          = 16;
caxis(clims)

%3.4 Subplot 3: Table with analysis & electrode information
subplot(DimSubplot(1), DimSubplot(2), 3)
textbox_info1 = {...
    [num2str(length(subs)) ' subjects'] ...
    [num2str(SignElecs.num_elecs) ' / ' num2str(length(SignElecs.array)) ' sign. elec / all elec'] ...
    [num2str(SignElecs.num_cluster) ' sign. clusters total'] []};
textbox_info2 = {};
if length(SignElecs.index) < 30 %one row
    for i_elec = 1:length(SignElecs.index)
        sel_elec = SignElecs.index(i_elec);
        textbox_info2{i_elec} = ...
            [labels_plotting_number{sel_elec} ' = ' SignElecs.fulllabels{i_elec}] ;
    end    
    textbox = [textbox_info1 textbox_info2];
    t = text(0, 0.5, 0, textbox, 'FontSize',8,'Interpreter','none');    
else
    for i_elec = 1:25 %two rows
        sel_elec = SignElecs.index(i_elec);
        textbox_info2{i_elec} = ...
            [labels_plotting_number{sel_elec} ' = ' SignElecs.fulllabels{i_elec}] ;
    end   
    for i_elec = 26:length(SignElecs.index)
        sel_elec = SignElecs.index(i_elec);
        textbox_info3{i_elec} = ...
            [labels_plotting_number{sel_elec} ' = ' SignElecs.fulllabels{i_elec}] ;
    end 
    textbox = [textbox_info1 textbox_info2];    
    t = text(-0.3, 0.5, 0, textbox, 'FontSize',7,'Interpreter','none');    
    t2 = text(0.6, 1, 0, textbox_info3, 'FontSize',7,'Interpreter','none');    
end
set(gca,'visible','off')
CounterSubplot = CounterSubplot + 1;

%3.5 Subplot 4: Effect onset
for i_elec = 1:length(Param2plot.all_subs.effect_onset)
    PlotInput(i_elec,1) = min(Param2plot.all_subs.effect_onset(i_elec,:));
end
clims           = [0 1];
Label_Colorbar  = {'effect onset [rel. to tone]'}; 
dat.avg         = PlotInput;

sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH...
    (coords_allsub, labels_plotting_empty, dat.avg, SignElecs.array,...
    chanSize, clims, cmap, dat.textcolor_rgb, ...
    DimSubplot, CounterSubplot, SubplotPosition, [], []);
sp_handle_surf{CounterSurfplot} = sp_handle_surf_temp.L;
CounterSubplot                  = CounterSubplot + 1;

ColorbarPosition    = [0.1 0.2 0 0.5]; %Left of surface
h = colorbar;
h.Position(1)       = ColorbarPosition(1); %sets colorbar to the right
h.Position(2)       = ColorbarPosition(2); %sets colorbar higher
h.Position(4)       = h.Position(4)*ColorbarPosition(4); %makes colorbar shorter
h.Label.String      = Label_Colorbar;
h.FontSize          = 16;
h.Ticks             = [0 0.25 0.5 0.75 1];
h.TickLength        = 0.1;
caxis(clims)

%3.6 Subplot 5: Effect offset
for i_elec = 1:length(Param2plot.all_subs.effect_offset)
    PlotInput(i_elec,1) = max(Param2plot.all_subs.effect_offset(i_elec,:));
end
clims           = [0 1];
Label_Colorbar  = {'effect offset [rel. to tone]'};
dat.avg         = PlotInput;

sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH...
    (coords_allsub, labels_plotting_empty, dat.avg, SignElecs.array,...
    chanSize, clims, cmap, dat.textcolor_rgb, ...
    DimSubplot, CounterSubplot, SubplotPosition, [], []);
sp_handle_surf{CounterSurfplot} = sp_handle_surf_temp.L;
CounterSubplot                  = CounterSubplot + 1;

ColorbarPosition    = [0.39 0.2 0 0.5]; %Left of surface
h                   = colorbar;
h.Position(1)       = ColorbarPosition(1); %sets colorbar to the right
h.Position(2)       = ColorbarPosition(2); %sets colorbar higher
h.Position(4)       = h.Position(4)*ColorbarPosition(4); %makes colorbar shorter
h.Label.String      = Label_Colorbar;
h.FontSize          = 16;
h.Ticks             = [0 0.25 0.5 0.75 1];
h.TickLength        = 0.1;
caxis(clims)

%3.7 Subplot 6: Effect duration
for i_elec = 1:length(Param2plot.all_subs.effect_duration)
    PlotInput(i_elec,1) = nansum(Param2plot.all_subs.effect_duration(i_elec,:));
end
clims           = [0 1];
Label_Colorbar  = {'effect duration [rel to tone]'};
dat.avg         = PlotInput;

sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH...
    (coords_allsub, labels_plotting_empty, dat.avg, SignElecs.array,...
    chanSize, clims, cmap, dat.textcolor_rgb, ...
    DimSubplot, CounterSubplot, SubplotPosition, [], []);
sp_handle_surf{CounterSurfplot} = sp_handle_surf_temp.L;

ColorbarPosition    = [0.67 0.2 0 0.5]; %In between both hemispheres
h                   = colorbar;
h.Position(1)       = ColorbarPosition(1); %sets colorbar to the right
h.Position(2)       = ColorbarPosition(2); %sets colorbar higher
h.Position(4)       = h.Position(4)*ColorbarPosition(4); %makes colorbar shorter
h.Label.String      = Label_Colorbar;
h.FontSize = 16;
h.Ticks             = [0 0.25 0.5 0.75 1];
h.TickLength        = 0.1;
caxis(clims)

%% 4 Save figure
if save_poststepFigs == 1
    filename     = ['Surf1HSignElec' FDR_label '_Allsubn' num2str(length(subs)) '_' ...
        FuncInput_EffectType '_' FuncInput_DataType ...
        '_AllTD_Stat.png'];
    figfile      = [path_fig filename];
    saveas(gcf, figfile, 'png'); %save png version
    close all;
end

end