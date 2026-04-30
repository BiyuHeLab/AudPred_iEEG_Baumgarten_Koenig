function NASTD_ECoG_Predict_PlotSignClusterElec_SubplotParam_AllSub...
    (subs, ...
    FuncInput_EffectType, FuncInput_DataType, FuncInput_ToneDur_text, ...
    param, ...
    save_poststepFigs, paths_NASTD_ECoG)

%Aim: Plot electrodes showing sign. cluster-corrected prediction effects
%per input data & TD condition aggregated over subjects. Based on
%sample-wise analysis, thus no time windows but instead entire analysis
%window. 
%Plots 2 figures, 1 with clusterstat & p-value derviative, 
%1 with effect on-, offset, and duration

% set(0, 'DefaultFigureVisible', 'on');

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')

%Add base dir and own script dir
addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

path_fig = ([paths_NASTD_ECoG.ECoGdata_Prediction ...
    '/PredEffects/Allsub_n' num2str(length(subs)) ...
    '/Figs/PredEffects_Surf/' FuncInput_EffectType '/ClusterCorr/']);
if (~exist(path_fig, 'dir')); mkdir(path_fig); end

%% 1) Load prediction effect data and aggregate relevant info across subjects
coords_allsub           = [];
labels_allsub           = [];
anatlabels_allsub       = [];
sub_index               = [];
usedElecs_chanposIndex  = [];

for i_sub = 1:length(subs)
    
    tic
    disp(['Loading data for sub: ' sub])

    sub = subs{i_sub};
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;
    
    %Load effect
    path_inputdata = [paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
        sub '/Data/Samplewise/'];
    load([path_inputdata sub '_PredEffectsCluster_' FuncInput_DataType '_' FuncInput_ToneDur_text 'sTD.mat'], ...
        FuncInput_EffectType , 'labels_loadedData');
    CurrentEffect = eval(FuncInput_EffectType);
    clear PredEffect SimplePredErrEffect ComplexPredErrEffect
    
    %Also load ECoG preproc data for channel labels and position
    loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];%path to indiv preprocessed ECoG data
    load(loadfile_ECoGpreprocdata);
    
    %Determine basic parameters
    SampleFreq  = DataClean_AllTrials.fsample;
    nSensors    = size(CurrentEffect.stats,2);
    nSamples    = size(CurrentEffect.clusterstat{1}.cluster_timecourse,2);
    
    %Read electrode labels, coordinates, and anatomical labels for all 
    %analyzed electrodes and aggregate them across subjects
    for i_elec = 1:length(labels_loadedData)
        usedElecs_chanposIndex(1,i_elec) = ...
            find(strcmp(labels_loadedData{i_elec}, ...
            DataClean_AllTrials.elec.label));
    end
    coords_sub{i_sub}   = ...
        DataClean_AllTrials.elec.chanpos(usedElecs_chanposIndex,:);    
    coords_allsub       = [coords_allsub; coords_sub{i_sub}];
    anatlabels_allsub   = ...
        [anatlabels_allsub; DataClean_AllTrials.elec.T1AnatLabel(usedElecs_chanposIndex)]; 
    for i_elec = 1:length(usedElecs_chanposIndex)
        labels_allsub{end+1,1} = [labels_loadedData{i_elec} '_' sub];        
    end   
    
    %Store clusterstat and p-values from sign. clusters only
    Param2plot.per_sub.clusterstat{i_sub}       = nan(10, nSensors);
    Param2plot.per_sub.pval{i_sub}              = nan(10, nSensors);
    Param2plot.per_sub.effect_onset{i_sub}      = nan(10, nSensors);
    Param2plot.per_sub.effect_offset{i_sub}     = nan(10, nSensors);
    Param2plot.per_sub.effect_duration{i_sub}   = nan(10, nSensors);
    
    for i_elec = 1:nSensors
        for i_cluster = 1:length(CurrentEffect.clusterstat{i_elec}.cluster_pval)
            if CurrentEffect.clusterstat{i_elec}.cluster_pval(i_cluster) < param.pval_plotting
                Param2plot.per_sub.clusterstat{i_sub}(i_cluster, i_elec) = ...
                    CurrentEffect.clusterstat{i_elec}.cluster_statSum(i_cluster);
                Param2plot.per_sub.pval{i_sub}(i_cluster, i_elec) = ...
                    CurrentEffect.clusterstat{i_elec}.cluster_pval(i_cluster);            
                Param2plot.per_sub.effect_onset{i_sub}(i_cluster, i_elec) = ...
                    CurrentEffect.clusterstat{i_elec}.cluster_samples{i_cluster}(1);
                Param2plot.per_sub.effect_offset{i_sub}(i_cluster, i_elec) = ...
                    CurrentEffect.clusterstat{i_elec}.cluster_samples{i_cluster}(end);
                Param2plot.per_sub.effect_duration{i_sub}(i_cluster, i_elec) = ...
                    length(CurrentEffect.clusterstat{i_elec}.cluster_samples{i_cluster});
            end
        end
    end
    %Derivative of p-value to estimate effect
    Param2plot.per_sub.pval_derivative{i_sub} = -(log10(Param2plot.per_sub.pval{i_sub}));
    
    %Array to differentiate subject entries
    Param2plot.all_subs.sub_index = [sub_index; ones(length(coords_sub{i_sub}),1)*i_sub];
    
    %clusterstat and p-value aggregation over subjects
    if i_sub == 1
        Param2plot.all_subs.clusterstat     = Param2plot.per_sub.clusterstat{i_sub}';
        Param2plot.all_subs.pval            = Param2plot.per_sub.pval{i_sub}';
        Param2plot.all_subs.pval_derivative = Param2plot.per_sub.pval_derivative{i_sub}';
        Param2plot.all_subs.effect_onset    = Param2plot.per_sub.effect_onset{i_sub}';
        Param2plot.all_subs.effect_offset   = Param2plot.per_sub.effect_offset{i_sub}';
        Param2plot.all_subs.effect_duration = Param2plot.per_sub.effect_duration{i_sub}';
    else
        Param2plot.all_subs.clusterstat     = [Param2plot.all_subs.clusterstat; Param2plot.per_sub.clusterstat{i_sub}'];
        Param2plot.all_subs.pval            = [Param2plot.all_subs.pval; Param2plot.per_sub.pval{i_sub}'];
        Param2plot.all_subs.pval_derivative = [Param2plot.all_subs.pval_derivative; Param2plot.per_sub.pval_derivative{i_sub}'];
        Param2plot.all_subs.effect_onset    = [Param2plot.all_subs.effect_onset; Param2plot.per_sub.effect_onset{i_sub}'];
        Param2plot.all_subs.effect_offset   = [Param2plot.all_subs.effect_offset; Param2plot.per_sub.effect_offset{i_sub}'];
        Param2plot.all_subs.effect_duration = [Param2plot.all_subs.effect_duration; Param2plot.per_sub.effect_duration{i_sub}'];
    end
    
    %Cleanup
    usedElecs_chanposIndex = [];
    clear CurrentEffect labels_loadedData used_elecs_chanposIndex
    disp(['done loading in ' num2str(toc) ' sec'])    
end

%% 2) Determine sign. electrodes
temp_indexsignelec      = nan(length(Param2plot.all_subs.pval),1);
for i_elec = 1:length(Param2plot.all_subs.pval)
    temp_indexsignelec(i_elec) = min(Param2plot.all_subs.pval(i_elec,:))';
end
SignElecs.array         = temp_indexsignelec < param.pval_plotting;
SignElecs.index         = find(SignElecs.array);
SignElecs.num_elecs     = length(SignElecs.index);
SignElecs.num_cluster   = sum(sum(Param2plot.all_subs.pval < param.pval_plotting));

%Determine labels of sign. elecs
SignElec.labels        = labels_allsub(SignElecs.index);
SignElec.anatlabels     = anatlabels_allsub(SignElecs.index);
[SignElec.labels, SignElec.anatlabels]

if ~isempty(SignElec.labels)
    sign_title = ...
        [num2str(SignElecs.num_elecs)...
        ' / ' num2str(length(labels_allsub)) ' sign. elecs']; %1 line
else
    sign_title = ...
        ['0 / ' num2str(length(labels_allsub)) ' supra-thresh. elecs']; %1 line
end

%Electrode labels?
if param.plotlabels == 1 %Visible electrode labels
    for i_elec = 1:length(labels_allsub)
        labels_plotting{i_elec} = ...
            [labels_allsub{i_elec}, '_', anatlabels_allsub{i_elec}];
    end
else
    for i_elec = 1:length(labels_allsub)%No electrode labels
        labels_plotting{i_elec} = '';
    end
end

%% 3) Prepare & Plot Figure 1: Clusterstat & p-value derivative
%3.1 Prepare figure
h = figure('visible','on'); %ensures figure doesn't pop up during plotting
set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
set(gcf,'renderer','opengl');

DimSubplot = [2,2]; %For both hemispheres
CounterSubplot = 1;
CounterSurfplot = 1;
SizeFactor = 4;

%3.2 Adjust header/title
switch FuncInput_EffectType
    case 'PredEffect'
        effect_text = ['Effect: Prediction (explain p33 activity by p*34); p < ' num2str(param.pval_plotting)];
    case 'SimplePredErrEffect'
        effect_text = ['Effect: Simple Prediction error (explain p34 activity by absolute p33-p34 difference); p < ' num2str(param.pval_plotting)];
    case 'ComplexPredErrEffect'
        effect_text = ['Effect: Complex Prediction error (explain p34 activity by absolute p*34-p34 difference); p < ' num2str(param.pval_plotting)];
end

Fig_title = {['Group level (n = ' num2str(length(subs)) ')'] ...
    [effect_text] ...
    ['Input data: ' FuncInput_DataType ', ' FuncInput_ToneDur_text 'ms TD']} ;
sgtitle(Fig_title,'FontSize',20,'Interpreter','none')

%3.2 Subplot 1: Clusterstat
for i_elec = 1:length(Param2plot.all_subs.clusterstat)
    PlotInput(i_elec,1) = max(Param2plot.all_subs.clusterstat(i_elec,:));
end

%Colorlim
dv_max                  = max(max(Param2plot.all_subs.clusterstat));
dv_min                  = min(min(Param2plot.all_subs.clusterstat));
dv_absmax               = max([abs(dv_max) abs(dv_min)]);
clims_maxclusterstat    = [-dv_absmax dv_absmax];
clims                   = clims_maxclusterstat;
Label_Colorbar          = 'max sum cluster stat';

%Data struct
dat.dimord      = 'chan_time';
dat.time        = 0;
dat.label       = labels_allsub;
dat.avg         = PlotInput;
dat.sign_elecs  = SignElecs.array;

chanSize    = ones(1,length(dat.avg))*SizeFactor; %electrode size (arbitrary)
cmap        = 'parula';

SubplotPosition = [0 -0.1 0 0];

%Surface plot
% coords_allsub(:,1) = abs(coords_allsub(:,1)) * -1; 
% %To project all electrodes on one hemisphere, 
% %take the absolute value (R) and mutliply it with -1 (L)

sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_SubplotHemis...
    (coords_allsub, labels_plotting, dat.avg, SignElecs.array,...
    chanSize, clims, cmap, ...
    DimSubplot, CounterSubplot, SubplotPosition, [], []);

sp_handle_surf{CounterSurfplot}     = sp_handle_surf_temp.L;
sp_handle_surf{CounterSurfplot+1}   = sp_handle_surf_temp.R;
CounterSubplot = CounterSubplot + 2;

%Colorbar
ColorbarPosition = [0.53 0.65 0 0.5]; %In between both hemispheres
h = colorbar;
h.Position(1) = ColorbarPosition(1); %sets colorbar to the right
h.Position(2) = ColorbarPosition(2); %sets colorbar higher
h.Position(4) = h.Position(4)*ColorbarPosition(4); %makes colorbar shorter
h.Label.String = Label_Colorbar;
h.FontSize = 16;
caxis(clims)

%3.3 Subplot 2: p-value derivative
for i_elec = 1:length(Param2plot.all_subs.pval_derivative)
    PlotInput(i_elec,1) = max(Param2plot.all_subs.pval_derivative(i_elec,:));
end

%Colorlim
clims = [0 max(PlotInput)];
Label_Colorbar = '- log10 (cluster p-value) ';

%prepare struct
dat.dimord      = 'chan_time';
dat.time        = 0;
dat.label       = labels_allsub;
dat.avg         = PlotInput;
dat.sign_elecs  = SignElecs.array;

chanSize = ones(1,length(dat.avg))*SizeFactor; %electrode size (arbitrary)

SubplotPosition = [0 0 0 0];

%Plot surface
sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_SubplotHemis...
    (coords_allsub, labels_plotting, dat.avg, SignElecs.array,...
    chanSize, clims, cmap, ...
    DimSubplot, CounterSubplot, SubplotPosition, [], []);
sp_handle_surf{CounterSurfplot} = sp_handle_surf_temp.L;
sp_handle_surf{CounterSurfplot+1} = sp_handle_surf_temp.R;

%Colorbar
ColorbarPosition = [0.53 0.2 0 0.5]; %In between both hemispheres
h = colorbar;
h.Position(1) = ColorbarPosition(1); %sets colorbar to the right
h.Position(2) = ColorbarPosition(2); %sets colorbar higher
h.Position(4) = h.Position(4)*ColorbarPosition(4); %makes colorbar shorter
h.Label.String = Label_Colorbar;
h.FontSize = 16;
caxis(clims)

%3.4 Add details via text
subplot(DimSubplot(1), DimSubplot(2), 3)
textbox_text = {...
    [num2str(length(subs)) ' subjects'] ...
    [num2str(SignElecs.num_elecs) ' / ' num2str(length(SignElecs.array)) ' sign. elec / all elec'] ...
    [num2str(SignElecs.num_cluster) ' sign. clusters total']};
t = text(0, -120, 100, textbox_text, 'FontSize',16,'Interpreter','none');

%3.5 Save figure
if save_poststepFigs == 1
    filename     = ['SurfSignElec_Allsubn' num2str(length(subs)) '_' ...
        FuncInput_EffectType '_' FuncInput_DataType ...
        '_' FuncInput_ToneDur_text 'msTD_Stat.png'];
    figfile      = [path_fig filename];
    saveas(gcf, figfile, 'png'); %save png version
    close all;
end

%% 4) Prepare & Plot Figure 2: Onset, Offset & Duration
%4.1 Prepare figure
h = figure('visible','off'); %ensures figure doesn't pop up during plotting
set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
set(gcf,'renderer','opengl');

DimSubplot = [3,2]; %For both hemispheres
CounterSubplot = 1;
CounterSurfplot = 1;
SizeFactor = 4;

%4.2 Adjust header/title
switch FuncInput_EffectType
    case 'PredEffect'
        effect_text = ['Effect: Prediction (explain p33 activity by p*34); p < ' num2str(param.pval_plotting)];
    case 'SimplePredErrEffect'
        effect_text = ['Effect: Simple Prediction error (explain p34 activity by absolute p33-p34 difference); p < ' num2str(param.pval_plotting)];
    case 'ComplexPredErrEffect'
        effect_text = ['Effect: Complex Prediction error (explain p34 activity by absolute p*34-p34 difference); p < ' num2str(param.pval_plotting)];
end

Fig_title = {['Group level (n = ' num2str(length(subs)) ')'] ...
    [effect_text] ...
    ['Input data: ' FuncInput_DataType ', ' FuncInput_ToneDur_text 'ms TD']} ;
sgtitle(Fig_title,'FontSize',20,'Interpreter','none')

%4.3 Subplot 1: Onset
for i_elec = 1:length(Param2plot.all_subs.effect_onset)
    PlotInput(i_elec,1) = min(Param2plot.all_subs.effect_onset(i_elec,:));
end
clims = [1 nSamples];
Label_Colorbar = {['effect onset'] ['[ms]']};

dat.dimord      = 'chan_time';
dat.time        = 0;
dat.label       = labels_allsub;
dat.avg         = PlotInput;
dat.sign_elecs  = SignElecs.array;

chanSize = ones(1,length(dat.avg))*SizeFactor; %electrode size (arbitrary)
cmap = 'parula';

SubplotPosition = [0 -0.05 0 0];

sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_SubplotHemis...
    (coords_allsub, labels_plotting, dat.avg, SignElecs.array,...
    chanSize, clims, cmap, ...
    DimSubplot, CounterSubplot, SubplotPosition, [], []);
sp_handle_surf{CounterSurfplot} = sp_handle_surf_temp.L;
sp_handle_surf{CounterSurfplot+1} = sp_handle_surf_temp.R;

hold on;
CounterSubplot = CounterSubplot + 2;

ColorbarPosition = [0.53 0.75 0 0.5]; %In between both hemispheres
h = colorbar;
h.Position(1) = ColorbarPosition(1); %sets colorbar to the right
h.Position(2) = ColorbarPosition(2); %sets colorbar higher
h.Position(4) = h.Position(4)*ColorbarPosition(4); %makes colorbar shorter
h.Label.String = Label_Colorbar;
h.FontSize = 16;
for i_ticklabel = 1:length(h.TickLabels) %Change from samples to ms
    h.TickLabels{i_ticklabel} = ...
        num2str(round(str2num(h.TickLabels{i_ticklabel}) / SampleFreq,2));
end
h.TickLength = 0.1;
caxis(clims)

%4.4 Subplot 2: Offset
for i_elec = 1:length(Param2plot.all_subs.effect_offset)
    PlotInput(i_elec,1) = max(Param2plot.all_subs.effect_offset(i_elec,:));
end
clims = [1 nSamples];
Label_Colorbar = {['effect offset'] ['[ms]']};

dat.dimord      = 'chan_time';
dat.time        = 0;
dat.label       = labels_allsub;
dat.avg         = PlotInput;
dat.sign_elecs  = SignElecs.array;

chanSize = ones(1,length(dat.avg))*SizeFactor; %electrode size (arbitrary)

SubplotPosition = [0 -0.025 0 0];

sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_SubplotHemis...
    (coords_allsub, labels_plotting, dat.avg, SignElecs.array,...
    chanSize, clims, cmap, ...
    DimSubplot, CounterSubplot, SubplotPosition, [], []);
sp_handle_surf{CounterSurfplot} = sp_handle_surf_temp.L;
sp_handle_surf{CounterSurfplot+1} = sp_handle_surf_temp.R;

CounterSubplot = CounterSubplot + 2;

ColorbarPosition = [0.53 0.45 0 0.5]; %In between both hemispheres
h = colorbar;
h.Position(1) = ColorbarPosition(1); %sets colorbar to the right
h.Position(2) = ColorbarPosition(2); %sets colorbar higher
h.Position(4) = h.Position(4)*ColorbarPosition(4); %makes colorbar shorter
h.Label.String = Label_Colorbar;
h.FontSize = 16;
for i_ticklabel = 1:length(h.TickLabels) %Change from samples to ms
    h.TickLabels{i_ticklabel} = ...
        num2str(round(str2num(h.TickLabels{i_ticklabel}) / SampleFreq,2));
end
h.TickLength = 0.1;
caxis(clims)

%Subplot 3: Duration
for i_elec = 1:length(Param2plot.all_subs.effect_duration)
    PlotInput(i_elec,1) = nansum(Param2plot.all_subs.effect_duration(i_elec,:));
end
clims = [1 nSamples];
Label_Colorbar = {['effect duration'] ['[ms]']};

dat.dimord      = 'chan_time';
dat.time        = 0;
dat.label       = labels_allsub;
dat.avg         = PlotInput;
dat.sign_elecs  = SignElecs.array;

chanSize = ones(1,length(dat.avg))*SizeFactor; %electrode size (arbitrary)

SubplotPosition = [0 0 0 0];

sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_SubplotHemis...
    (coords_allsub, labels_plotting, dat.avg, SignElecs.array,...
    chanSize, clims, cmap, ...
    DimSubplot, CounterSubplot, SubplotPosition, [], []);
sp_handle_surf{CounterSurfplot} = sp_handle_surf_temp.L;
sp_handle_surf{CounterSurfplot+1} = sp_handle_surf_temp.R;

ColorbarPosition = [0.53 0.18 0 0.5]; %In between both hemispheres
h = colorbar;
h.Position(1) = ColorbarPosition(1); %sets colorbar to the right
h.Position(2) = ColorbarPosition(2); %sets colorbar higher
h.Position(4) = h.Position(4)*ColorbarPosition(4); %makes colorbar shorter
h.Label.String = Label_Colorbar;
h.FontSize = 16;
for i_ticklabel = 1:length(h.TickLabels) %Change from samples to ms
    h.TickLabels{i_ticklabel} = ...
        num2str(round(str2num(h.TickLabels{i_ticklabel}) / SampleFreq,2));
end
h.TickLength = 0.1;
caxis(clims)

%4.4 Add details via text
subplot(DimSubplot(1), DimSubplot(2), 5)
textbox_text = {...
    [num2str(length(subs)) ' subjects'] ...
    [num2str(SignElecs.num_elecs) ' / ' num2str(length(SignElecs.array)) ' sign. elec / all elec'] ...
    [num2str(SignElecs.num_cluster) ' sign. clusters total']};
t = text(0, 250, 200, textbox_text, 'FontSize',16,'Interpreter','none');

%4.5 Save Figure
if save_poststepFigs == 1
    filename     = ['SurfSignElec_Allsubn' num2str(length(subs)) '_' ...
        FuncInput_EffectType '_' FuncInput_DataType ...
        '_' FuncInput_ToneDur_text 'msTD_Dur.png'];
    figfile      = [path_fig filename];
    saveas(gcf, figfile, 'png'); %save png version
    close all;
end

end