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

path_fig = ([paths_NASTD_ECoG.ECoGdata_StimCorr ... 
    '/Allsub_n' num2str(length(subs)) '/Figs/IdentvsSim/p' ...
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
temp_indexsignelec      = nan(length(Param2plot.all_subs.corr_ttest_p),1);
for i_elec = 1:length(Param2plot.all_subs.corr_ttest_p)
    temp_indexsignelec(i_elec) = Param2plot.all_subs.corr_ttest_p(i_elec);
end
SignElecs.array         = temp_indexsignelec < param.pval_plotting;
SignElecs.index         = find(SignElecs.array);
SignElecs.num_elecs     = length(SignElecs.index);

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

%% 3) Plot summary figure containing:
%1) Correlation coefficient to identical sequences
%2) Correlation coefficient to similar sequences
%3) Statistic of indentical vs. similar comparison
%4) Table with anatomical labels for sign. electrodes

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

%7.1 Plot surface plot with correlation coefficient for identical sequences
plot_struct.dimord          = 'chan_time';
plot_struct.time            = 0;
plot_struct.sign_elecs      = ...
    logical(ones(length(Param2plot.all_subs.corrcoeff_identical),1));
plot_struct.clims           = [-0.25 0.25]; %free symmetric scaling
plot_struct.chanSize        = ...
    ones(1,length(plot_struct.sign_elecs))*SizeFactor; %electrode size (arbitrary)
plot_struct.cmap            = 'jet';
plot_struct.textcolor_rgb   = [0 0 0];

%Remove electrode labels
for i_elec = 1:length(plot_struct.sign_elecs)
    plot_struct.label{i_elec} = '';
end

plot_struct.avg             = Param2plot.all_subs.corrcoeff_identical;

%Plot surface with highlighted sign electrodes
sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH...
    (plot_struct.coords, plot_struct.label, plot_struct.avg, plot_struct.sign_elecs,...
    plot_struct.chanSize, plot_struct.clims, plot_struct.cmap, plot_struct.textcolor_rgb, ...
    DimSubplot, CounterSubplot, SubplotPosition, ColorbarPosition, []);

sp_handle_surf{CounterSubplot} = sp_handle_surf_temp.L;
CounterSubplot = CounterSubplot + 1;
title('Identical sequences','FontSize',14)

%7.2 Plot surface plot with correlation coefficient for identical sequences
plot_struct.avg         = Param2plot.all_subs.corrcoeff_similar;

%Plot surface with highlighted sign electrodes
sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH...
    (plot_struct.coords, plot_struct.label, plot_struct.avg, plot_struct.sign_elecs,...
    plot_struct.chanSize, plot_struct.clims, plot_struct.cmap, plot_struct.textcolor_rgb, ...
    DimSubplot, CounterSubplot, SubplotPosition, ColorbarPosition, []);

sp_handle_surf{CounterSubplot} = sp_handle_surf_temp.L;
CounterSubplot = CounterSubplot + 1;
title('Similar sequences','FontSize',14)

%Add colorbar
ColorbarPosition = [sp_handle_surf{CounterSubplot-1}.Position(1)-0.05 ...
    sp_handle_surf{CounterSubplot-1}.Position(2)+0.05 0 0];
h = colorbar;
h.Position(1) = ColorbarPosition(1); %sets colorbar to the right
h.Position(2) = ColorbarPosition(2); %sets colorbar higher
h.Position(4) = h.Position(4)*0.8; %makes colorbar shorter
h.Label.String = ['Pearsons r'];
h.FontSize = 14;
caxis(plot_struct.clims)

%7.3 Plot surface plot showing stats for identical vs. similar contrast
%Determine labels of sign. electrodes
% plot_struct.label           = NameLabels_plotting; %NumberLabels_plotting
plot_struct.avg             = Param2plot.all_subs.corr_ttest_t;
plot_struct.sign_elecs      = logical(SignElecs.array);
plot_struct.clims           = [0 5];
plot_struct.textcolor_rgb   = [0 0 0];
plot_struct.chanSize        = ...
    ones(1,length(plot_struct.sign_elecs))*SizeFactor; %electrode size (arbitrary)

%Plot surface with highlighted sign electrodes
sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH...
    (plot_struct.coords, plot_struct.label, plot_struct.avg, plot_struct.sign_elecs,...
    plot_struct.chanSize, plot_struct.clims, plot_struct.cmap, plot_struct.textcolor_rgb, ...
    DimSubplot, CounterSubplot, SubplotPosition, ColorbarPosition, []);

sp_handle_surf{CounterSubplot} = sp_handle_surf_temp.L;
CounterSubplot = CounterSubplot + 1;
title({['Identical vs. similar sequences (p < ' num2str(param.pval_plotting) ' uncorr)'], ...
    ['# sign. elecs / all elecs: ' ...
    num2str(SignElecs.num_elecs) ' / ' num2str(length(SensorLabelsSubMarker_allsub))]}, ...
    'FontSize',14)

%Add colorbar
ColorbarPosition = [sp_handle_surf{CounterSubplot-1}.Position(1)-0.05 ...
    sp_handle_surf{CounterSubplot-1}.Position(2)+0.02 0 0];
h = colorbar;
h.Position(1) = ColorbarPosition(1); %sets colorbar to the right
h.Position(2) = ColorbarPosition(2); %sets colorbar higher
h.Position(4) = h.Position(4)*0.8; %makes colorbar shorter
h.Label.String = ['t-value (identical vs. different)'];
h.FontSize = 14;
caxis(plot_struct.clims)

%7.4 Table with analysis & electrode information
textbox         = [];
textbox_info2   = {};
NumTextRows = ceil(SignElecs.num_elecs/25);
if NumTextRows < 3
    param.FontSize  = 10;
else
    param.FontSize  = 8;
end

% subplot(DimSubplot(1), DimSubplot(2), 4)
% textbox_info1 = {...
%     [sub] ...
%     ['# sign. elecs / all elecs: ' ...
%     num2str(length(SignElecs.index)) ' / ' num2str(length(SensorLabelsSubMarker_allsub))] ...
%     ''};
% if NumTextRows == 1
%     for i_elec = 1:SignElecs.num_elecs
%         sel_elec = SignElecs.index(i_elec);
%         textbox_info2{i_elec} = ...
%             [SignElecs.fulllabels{i_elec}];
%     end
%     textbox = [textbox_info1 textbox_info2];
%     t = text(0, 0.5, 0, textbox, ...
%         'FontSize',param.FontSize,'Interpreter','none');
% else
%     elec_counter = 0;
%     for i_NumTextRows = 1:NumTextRows
%         counter_placing = 1;
%         for i_elec = 1+elec_counter:25+elec_counter
%             while i_elec <= SignElecs.num_elecs
%                 sel_elec = SignElecs.index(i_elec);
%                 textbox_info2{i_NumTextRows}{counter_placing} = ...
%                     SignElecs.fulllabels{i_elec};
%                 counter_placing = counter_placing +1;
%                 break
%             end
%         end
%         elec_counter = elec_counter + 25;
%     end
%     
%     textbox = [textbox_info1 textbox_info2{1}];       
%     t = text(-0.3, 0.5, 0, textbox, ...
%         'FontSize', param.FontSize,'Interpreter','none');
%     for i_NumTextRows = 2:NumTextRows    
%         t2 = text(-0.525 + (0.33 * i_NumTextRows), 0.425, 0, textbox_info2{i_NumTextRows}, ...
%             'FontSize',param.FontSize,'Interpreter','none');
%     end    
% end
% % set(gca,'visible','on')

%7.5 Save Figure
if save_poststepFigs == 1
    filename     = ['StimCorr_Surf1HSignElec_n' num2str(length(subs)) '_' ...
        FuncInput_DataType '_AllTD.png'];
    figfile      = [path_fig filename];
    saveas(gcf, figfile, 'png'); %save png version
    close all;
end

end