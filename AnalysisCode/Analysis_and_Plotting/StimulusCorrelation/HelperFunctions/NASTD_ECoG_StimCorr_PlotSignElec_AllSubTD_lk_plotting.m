function NASTD_ECoG_StimCorr_PlotSignElec_AllSubTD ...
    (subs, ...
    FuncInput_DataType, FuncInput_ToneDur_text, ...
    param, ...
    save_poststepFigs, paths_NASTD_ECoG)

% set(0, 'DefaultFigureVisible', 'on');

%% 0.1) Specify vars, paths, and setup fieldtrip
%Add base dir and own script dir
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

path_fig = (['/isilon/LFMI/VMdrive/Lua/NASTD/Figs/Allsub_n' num2str(length(subs)) '_IdentvsSim_StimCorr_p' ...
    num2str(param.pval_plotting) '/']);
if (~exist(path_fig, 'dir')); mkdir(path_fig); end

%% 1) Load stimulus correlation data and aggregate relevant info across subjects
clear Param2plot

coords_sub                      = [];
coords_allsub                   = [];
SensorLabelsSubMarker_allsub    = [];
SensorLabels_allsub             = [];
AnatLabels_allsub               = [];
usedElecs_chanposIndex          = [];
Param2plot.all_subs.sub_index   = [];

for i_sub = 1:length(subs)
    usedElecs_chanposIndex = [];
    tic
    sub = subs{i_sub};
    disp(['-- Loading data for sub: ' sub ' --'])
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;
    
   for i_TD = 1:length(FuncInput_ToneDur_text)
        %Load stimullus correlation data and select current effect
        path_inputdata = [paths_NASTD_ECoG.ECoGdata_StimCorr '/' sub '/Data/'];
        load([path_inputdata sub '_StimCorr_' FuncInput_DataType '_' ...
            FuncInput_ToneDur_text{i_TD} 'sTD.mat'], ...
            'corrcoeff_identical', 'corrcoeff_similar', 'corr_ttest', ...
            'anova_F', 'anova_p', ...
            'SensorLabels');
        
        %Also load ECoG preproc data for channel labels and position
        loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];
        load(loadfile_ECoGpreprocdata);
        
        %Determine basic parameters
        SampleFreq  = DataClean_AllTrials.fsample;
        nSensors    = length(SensorLabels);
       
        %Read electrode labels, coordinates, and anatomical labels for all
        %analyzed electrodes and aggregate them across subjects
        for i_elec = 1:length(SensorLabels)
            usedElecs_chanposIndex(1,i_elec) = ...
                find(strcmp(SensorLabels{i_elec}, ...
                DataClean_AllTrials.elec.label));
        end
        
        coords_sub{i_sub, i_TD}   = ...
            DataClean_AllTrials.elec.chanpos(usedElecs_chanposIndex,:);
        coords_allsub       = [coords_allsub; coords_sub{i_sub, i_TD}];
        AnatLabels_allsub   = ...
            [AnatLabels_allsub; DataClean_AllTrials.elec.T1AnatLabel(usedElecs_chanposIndex)];
        
        for i_elec = 1:length(usedElecs_chanposIndex)
            SensorLabelsSubMarker_allsub{end+1,1}   = [SensorLabels{i_elec} ' ' sub];
            SensorLabels_allsub{end+1,1}            = [SensorLabels{i_elec}];
        end
        
        %Store results for each subject
        Param2plot.per_sub.corrcoeff_identical{i_sub,i_TD}   = nan(nSensors, 1);
        Param2plot.per_sub.corrcoeff_similar{i_sub,i_TD}     = nan(nSensors, 1);
        Param2plot.per_sub.corr_ttest_t{i_sub,i_TD}          = nan(nSensors, 1);
        Param2plot.per_sub.corr_ttest_p{i_sub,i_TD}          = nan(nSensors, 1);
        Param2plot.per_sub.anova_F{i_sub,i_TD}               = nan(nSensors, 1);
        Param2plot.per_sub.anova_p{i_sub,i_TD}               = nan(nSensors, 1);
        
        Param2plot.per_sub.corrcoeff_identical{i_sub,i_TD}(:,1) = corrcoeff_identical;
        Param2plot.per_sub.corrcoeff_similar{i_sub,i_TD}(:,1)   = corrcoeff_similar;
        Param2plot.per_sub.corr_ttest_t{i_sub,i_TD}(:,1)        = corr_ttest.t;
        Param2plot.per_sub.corr_ttest_p{i_sub,i_TD}(:,1)        = corr_ttest.p;
        Param2plot.per_sub.anova_F{i_sub,i_TD}(:,1)             = anova_F;
        Param2plot.per_sub.anova_p{i_sub,i_TD}(:,1)             = anova_p;
        
        %Array to differentiate subject entries
        Param2plot.all_subs.sub_index = ...
            [Param2plot.all_subs.sub_index; ones(length(coords_sub{i_sub, i_TD}),1)*i_sub];
        
        %data aggregation over subjects
        if i_sub == 1 && i_TD == 1
            Param2plot.all_subs.corrcoeff_identical     = ...
                Param2plot.per_sub.corrcoeff_identical{i_sub,i_TD};
            Param2plot.all_subs.corrcoeff_similar       = ...
                Param2plot.per_sub.corrcoeff_similar{i_sub,i_TD};
            Param2plot.all_subs.corr_ttest_t            = ...
                Param2plot.per_sub.corr_ttest_t{i_sub,i_TD};
            Param2plot.all_subs.corr_ttest_p            = ...
                Param2plot.per_sub.corr_ttest_p{i_sub,i_TD};
            Param2plot.all_subs.anova_F                 = ...
                Param2plot.per_sub.anova_F{i_sub,i_TD};
            Param2plot.all_subs.anova_p                 = ...
                Param2plot.per_sub.anova_p{i_sub,i_TD};
        else
            Param2plot.all_subs.corrcoeff_identical     = ...
                [Param2plot.all_subs.corrcoeff_identical; ...
                Param2plot.per_sub.corrcoeff_identical{i_sub,i_TD}];
            Param2plot.all_subs.corrcoeff_similar       = ...
                [Param2plot.all_subs.corrcoeff_similar; ...
                Param2plot.per_sub.corrcoeff_similar{i_sub,i_TD}];
            Param2plot.all_subs.corr_ttest_t            = ...
                [Param2plot.all_subs.corr_ttest_t; ...
                Param2plot.per_sub.corr_ttest_t{i_sub,i_TD}];
            Param2plot.all_subs.corr_ttest_p            = ...
                [Param2plot.all_subs.corr_ttest_p; ...
                Param2plot.per_sub.corr_ttest_p{i_sub,i_TD}];
            Param2plot.all_subs.anova_F                 = ...
                [Param2plot.all_subs.anova_F; ...
                Param2plot.per_sub.anova_F{i_sub,i_TD}];
            Param2plot.all_subs.anova_p                 = ...
                [Param2plot.all_subs.anova_p; ...
                Param2plot.per_sub.anova_p{i_sub,i_TD}];
        end
        
        %Cleanup
        usedElecs_chanposIndex = [];
        clear SensorLabels corr* anova* DataClean_AllTrials
    end
    disp(['-- done loading in ' num2str(toc) ' sec --'])
end

%% 2) Determine sign. electrodes

%% LK change Feb 18, 2025: Take as significant an electrode if it has significant stimulus correlation in at least 1 tone duration. Previous code 
%% computed number of significant electrodes across both tone durations so if one electrode was significant in both tone durations it was counted
%% twice in the number of significant electrodes. Similarly, the total number of electrodes was doubled (across both tone durs).
uniqueSensors = unique(SensorLabelsSubMarker_allsub);
temp_indexsignelec      = nan(length(uniqueSensors),1);
temp_indices = nan(length(uniqueSensors),1);

for i_elec = 1:length(uniqueSensors)
    elecName = uniqueSensors{i_elec};
    index_of_bothToneDurs = find(strcmp(SensorLabelsSubMarker_allsub, elecName));
    % Take the smallest p-value to evaluate if it is lower than
    % param.pval_plotting since only one p-val needs to be significant for
    % the electrode to count as significant (liberal criterion)
    [temp_indexsignelec(i_elec), idx] = min(Param2plot.all_subs.corr_ttest_p(index_of_bothToneDurs));
    temp_indices(i_elec) = index_of_bothToneDurs(idx);
end

SignElecs.array         = temp_indexsignelec < param.pval_plotting;
WhichSig  = find(SignElecs.array);
SignElecs.index = temp_indices(WhichSig);
SignElecs.array = logical(zeros(length(Param2plot.all_subs.corrcoeff_identical),1));
SignElecs.array(SignElecs.index) = 1;
SignElecs.num_elecs     = length(SignElecs.index);

disp(SignElecs.num_elecs);

%Determine labels of sign. elecs and add number to subtitle
SignElecs.labels        = SensorLabelsSubMarker_allsub(SignElecs.index);
SignElecs.anatlabels    = AnatLabels_allsub(SignElecs.index);
SignElecs.fulllabels    = [];
for i_elec = 1:length(SignElecs.labels)
    SignElecs.fulllabels{i_elec,1} = ...
        [SignElecs.labels{i_elec}, ' ', SignElecs.anatlabels{i_elec}];
end

if ~isempty(SignElecs.labels)
    sign_title = ...
        [num2str(SignElecs.num_elecs)...
        ' / ' num2str(length(SensorLabelsSubMarker_allsub)) ' sign. elecs']; %1 line
else
    sign_title = ...
        ['0 / ' num2str(length(SensorLabelsSubMarker_allsub)) ' sign. elecs']; %1 line
end

%Create numerical and name labels exclsuively for sign. elecs 
for i_elec = 1:length(SensorLabelsSubMarker_allsub)%No electrode labels
    EmptyLabels_plotting{i_elec}   = '';
end
counter_elecs           = 0;
NumberLabels_plotting   = EmptyLabels_plotting;
NameLabels_plotting     = EmptyLabels_plotting;
for i_elec = SignElecs.index'
    counter_elecs = counter_elecs +1;
    NumberLabels_plotting{i_elec}   = num2str(counter_elecs);
    NameLabels_plotting{i_elec}     = SensorLabels_allsub{i_elec};
end

%% 3. LK Edit June 9, 2025: For each input data type, check which electrodes are significant, and color code by subject ID. 
%  Goal is to see if significant correlations are driven by several
%  subjects or many
last5 = cellfun(@(x) x(end-4:end), SignElecs.labels, 'UniformOutput', false); % extract last five characters of labels which indicates subject ID
% Assign unique number for each subject ID
[uniqueTags, ~, tagIndices] = unique(last5, 'stable');
subjID_sigelecs      = nan(length(Param2plot.all_subs.corrcoeff_identical),1);
subjID_sigelecs(SignElecs.index) = tagIndices;

% Plot for different levels of significance
pvals = [0.05, 0.01, 0.001];
tvals = [1.9, 3, 4.79];

clear plot_struct

%7.0 Set up subplot structure
h = figure('visible','on'); %ensures figure doesn't pup during plotting
set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen

DimSubplot          = [2 2];
SubplotPosition     = [0 -0.05 0 0];
ColorbarPosition    = [0 0 0 0];
SizeFactor          = 2.5;
CounterSubplot      = 1;

sgtitle({...
    ['Across-trial correlation of neural data for identical vs. similar sequences'], ...
    ['n = ' num2str(length(subs)) ' - ' FuncInput_DataType ' - ' ...
    'All TD - (RH elecs projected on LH)']}, ...
    'Interpreter','none')

plot_struct.coords      = coords_allsub;

%Project all electrodes on left hemisphere,
plot_struct.coords(:,1) = abs(coords_allsub(:,1)) * -1;

plot_struct.dimord          = 'chan_time';
plot_struct.time            = 0;
plot_struct.sign_elecs      = ...
    logical(ones(length(Param2plot.all_subs.corrcoeff_identical),1));
plot_struct.clims           = [-0.25 0.25]; %free symmetric scaling
plot_struct.chanSize        = ...
    ones(1,length(plot_struct.sign_elecs))*SizeFactor; %electrode size (arbitrary)
plot_struct.cmap            = 'jet';
plot_struct.textcolor_rgb   = [0 0 0];

% %Remove electrode labels
for i_elec = 1:length(plot_struct.sign_elecs)
    plot_struct.label{i_elec} = '';
end

%7.3 Plot surface plot showing stats for identical vs. similar contrast
%Determine labels of sign. electrodes
% plot_struct.label           = NameLabels_plotting; %NumberLabels_plotting
plot_struct.avg             = subjID_sigelecs;
plot_struct.sign_elecs      = ...
    logical(zeros(length(Param2plot.all_subs.corrcoeff_identical),1));
plot_struct.sign_elecs(SignElecs.array)= 1;
plot_struct.clims           = [1 9];
plot_struct.textcolor_rgb   = [0 0 0];
plot_struct.chanSize        = ...
    ones(1,length(plot_struct.sign_elecs))*SizeFactor; %electrode size (arbitrary)

%Plot surface with highlighted sign electrodes
sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH...
    (plot_struct.coords, plot_struct.label, plot_struct.avg, plot_struct.sign_elecs,...
    plot_struct.chanSize, plot_struct.clims, plot_struct.cmap, plot_struct.textcolor_rgb, ...
    DimSubplot, CounterSubplot, SubplotPosition, ColorbarPosition, []);

sp_handle_surf{CounterSubplot} = sp_handle_surf_temp.L;
title({['Identical vs. similar sequences (p < ' num2str(param.pval_plotting) ' uncorr)'], ...
    ['# sign. elecs / all elecs: ' ...
    num2str(SignElecs.num_elecs) ' / ' num2str(length(uniqueSensors))]}, ...
    'FontSize',14)

% Plot only electrodes that are significant above a specific t-value
for tval_index = 2:3
    CounterSubplot = CounterSubplot + 1;
    
    data_to_plot = abs(Param2plot.all_subs.corr_ttest_t(SignElecs.index)) > tvals(tval_index);
    plot_struct.avg = subjID_sigelecs;
    plot_struct.sign_elecs      = logical(zeros(length(Param2plot.all_subs.corrcoeff_identical),1));
    plot_struct.sign_elecs(SignElecs.index(data_to_plot))= 1;
    num_elecs = sum(data_to_plot);

    sp_handle_surf_temp2 = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH...
    (plot_struct.coords, plot_struct.label, plot_struct.avg, plot_struct.sign_elecs,...
    plot_struct.chanSize, plot_struct.clims, plot_struct.cmap, plot_struct.textcolor_rgb, ...
    DimSubplot, CounterSubplot, SubplotPosition, ColorbarPosition, []);

    sp_handle_surf{CounterSubplot} = sp_handle_surf_temp2.L;
    title({['Identical vs. similar sequences (p < ' num2str(pvals(tval_index)) ' uncorr)'], ...
    ['# sign. elecs / all elecs: ' ...
    num2str(num_elecs) ' / ' num2str(length(uniqueSensors))]}, ...
    'FontSize',14)
end

%Add colorbar
ColorbarPosition = [sp_handle_surf{CounterSubplot-1}.Position(1)-0.05 ...
    sp_handle_surf{CounterSubplot-1}.Position(2)+0.02 0 0];
h = colorbar;
h.Position(1) = ColorbarPosition(1); %sets colorbar to the right
h.Position(2) = ColorbarPosition(2); %sets colorbar higher
h.Position(4) = h.Position(4)*0.8; %makes colorbar shorter
h.Label.String = ['subject ID'];
h.FontSize = 14;
caxis(plot_struct.clims)    

if save_poststepFigs == 1
    filename     = ['StimCorr_Surf1HSignElec_n' num2str(length(subs)) '_' ...
        FuncInput_DataType '_AllTD_BySubjColor.png'];
    figfile      = [path_fig filename];
    saveas(gcf, figfile, 'png'); %save png version
    close all;
end
end