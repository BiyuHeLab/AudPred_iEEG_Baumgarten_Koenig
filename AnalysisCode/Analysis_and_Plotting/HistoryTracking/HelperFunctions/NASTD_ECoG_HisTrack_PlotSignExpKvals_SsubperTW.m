function NASTD_ECoG_HisTrack_PlotSignExpKvals_SsubperTW...
    (sub, input_ToneDurLabel, input_DataLabel, ...
    SamplesTW, Label_TW,...
    param, ...
    FDR_correct, pval_plotting,...
    plot_SubplotperTW, save_poststepFigs, ...
    paths_NASTD_ECoG)

%Aim: Plot significant exp Kprime values for all electrodes on MNI volume for single
%subjects,per TW

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

path_inputdata = [paths_NASTD_ECoG.ECoGdata_HisTrack sub '/ExpvsShuff/' Label_TW '/'];
path_fig1 = [paths_NASTD_ECoG.ECoGdata_HisTrack sub '/ExpvsShuff/' Label_TW '/Figs/Kprime_Surf/'];
if (~exist(path_fig1, 'dir')); mkdir(path_fig1); end

%Tone duration condition for data load-in
if strcmp(input_ToneDurLabel,'0.2')
    tonedur_title = '200msTD';
elseif  strcmp(input_ToneDurLabel,'0.4')
    tonedur_title = '400msTD';
end

%% 1) Load preproc data (for electrode info) and combined expKprime data
NASTD_ECoG_subjectinfo %load subject info file (var: si)
loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];%path to indiv preprocessed ECoG data
tic
disp(['Loading preprocessed & combined exp Kprime data set for sub: ' sub])
load(loadfile_ECoGpreprocdata);

load([path_inputdata sub '_' input_DataLabel '_' tonedur_title ...
    '_HisTrack_ExpvsShuffKprimeCombFolds_unbalanced.mat']); %filename: Kprime_data, labels_loadedData
disp(['done loading in ' num2str(toc) ' sec'])

%% 2) Define analysis parameters
%2.1 Define analysis parameters
fsample         = DataClean_AllTrials.fsample;
toneDur_inSecs  = str2num(input_ToneDurLabel);
nSamplesPerTone = toneDur_inSecs * fsample;

win_size    = SamplesTW;  %25 data points with fs 512 = ~49ms as mentioned in manuscript
s_per_win   = (1/fsample)*win_size;
win_overlap = 0; %as mentioned in manuscript

windows     = [1 win_size];
while windows(end,end) < nSamplesPerTone
    windows = [windows; windows(end,:) + (win_size - win_overlap)];
end

if windows(end,end) > nSamplesPerTone
    windows(end,:) = [];
end

% windows_inms = (windows / fsample) * 1000;

%Optional: Load stimulus correlation data and select current effect
path_inputdata = [paths_NASTD_ECoG.ECoGdata_StimCorr  sub '/Data/'];

if strcmp(param.ElecSelect, 'StimCorr')
    load([path_inputdata sub '_StimCorr_' input_DataLabel '_' ...
        input_ToneDurLabel 'sTD.mat'], ...
        'corr_ttest', 'SensorLabels');
    filt_signelecs_StimCorr = corr_ttest.p < 0.05;%only sign. stim corr elecs
else
    filt_signelecs_StimCorr = true(length(labels_loadedData),1);  %all elecs
end
ind_selelecs  = find(filt_signelecs_StimCorr);
nSensors_sel  = sum(filt_signelecs_StimCorr);
%labels_loadedData(ind_selelecs)

%Determine selected electrode EDF-index to read out elec coords
ind_selelecs_MNI   = [];
for i_selelec = 1:length(ind_selelecs)
    ind_selelecs_MNI = [ind_selelecs_MNI; ...
        find(strcmp(DataClean_AllTrials.elec.label, ...
        labels_loadedData{ind_selelecs(i_selelec)}))];
end

%% 3) Plot exp Kprime values on surface

%3.1 find zlimits for constant plotting across time windows
dv_absmax = 0;
for i_win = 1:size(windows,1)
    dv_win_absmax = max(Kprime_data.Exp.avgRuns.Kprime_Avg{i_win});
    if dv_win_absmax > dv_absmax
        dv_absmax = dv_win_absmax;
    end
end

if plot_SubplotperTW == 1 %one subplot per TW
    figure('visible','off'); %ensures figure doesn't pop up during plotting
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    if strcmp(sub,'NY723')
        DimSubplot = [4,size(windows,1)/2];        
    else
        DimSubplot = [2,size(windows,1)/2];
    end
    CounterSubplot = 1;
end

%3.2 Plot data for current time window
for i_win = 1 : size(windows,1)
    
    %3.2.1 Determine significance (FDR or non-FDR thresholding)
    if FDR_correct == 1
        pval_FDRcorrected       = mafdr(...
            Kprime_data.Exp.avgRuns.pval_ExpvsShuff{i_win}(ind_selelecs),'BHFDR', true);
        Index_SuprathresElecs   = pval_FDRcorrected < pval_plotting;
        FDR_label               = 'FDRcorrected';
    else
        Index_SuprathresElecs   = ...
            Kprime_data.Exp.avgRuns.pval_ExpvsShuff{i_win}(ind_selelecs) < pval_plotting;
        FDR_label               = 'uncorrected';
    end
        sumSignELec             = sum(Index_SuprathresElecs);
    
    temp_Kprime_SuprathresElecs = [];        
    temp_Kprime_SuprathresElecs = Kprime_data.Exp.avgRuns.Kprime_Avg{i_win}(ind_selelecs);
    Kprime_SuprathresElecs      = temp_Kprime_SuprathresElecs(Index_SuprathresElecs);
    
    %3.2.1 Set up title information
    w1 = num2str( 1000 * windows(i_win, 1) / fsample , 3 );
    w2 = num2str( 1000 * windows(i_win, 2) / fsample , 3 );
    
    %     effect_title = 'Exp Kprime values (uncorrected; defined by minSumSquaredResiduals/ModelOrder for LinReg on tone pitch history))';
    effect_title = ['Thresholded Exp Kprime values (Exp vs. Shuff; p < ' ...
        num2str(pval_plotting) '; ' FDR_label ')'];
    if i_win == 1
        win_title = {['TW = ' w1 ' - ' w2 ' ms '],...
            ['min/max (mean +- SD) Kprime = ' ...
            num2str(min(Kprime_SuprathresElecs)) '/' num2str(max(Kprime_SuprathresElecs)) ' (' ...
            num2str(mean(Kprime_SuprathresElecs)) ' +- ' num2str(std(Kprime_SuprathresElecs)) ')']};
    else
        win_title = {['TW = ' w1 ' - ' w2 ' ms '],...
            [num2str(min(Kprime_SuprathresElecs)) '/' num2str(max(Kprime_SuprathresElecs)) ' (' ...
            num2str(mean(Kprime_SuprathresElecs)) ' +- ' num2str(std(Kprime_SuprathresElecs)) ')']};        
    end
    
    index_signelec_allelec  = ind_selelecs_MNI(Index_SuprathresElecs);
    label_signelec          = [];
    if ~isempty(index_signelec_allelec)
        for i_signelec = 1:length(index_signelec_allelec)
            label_signelec = [label_signelec ' ' ...
                DataClean_AllTrials.label{index_signelec_allelec(i_signelec)}];
        end
    end
    sign_title = [num2str(sumSignELec)...
        ' / ' num2str(length(ind_selelecs_MNI)) ' supra-thresh. elecs']; %1 line
    
    %3.2.2 Set up data struct for plotting    
    vals        = Kprime_data.Exp.avgRuns.Kprime_Avg{i_win}(ind_selelecs);%parameter of interest for resp. electrodes    
    sign_index  = Index_SuprathresElecs;    
    coords      = DataClean_AllTrials.elec.chanpos(ind_selelecs_MNI,:); %MNI coordinates for selected electrodes
    chanSize    = ones(1,length(vals))*3; %electrode size (arbitrary)    
% %Optional: plot higher Kprime values bigger
% for i_elec = 1:length(vals)
%     chanSize(i_elec,1) = ...
%         (vals(i_elec)./vals(i_elec)) * ...
%         (abs(log10(vals(i_elec)))*4);
% end

    clims       = [0 15]; %fixed scaling
    cmap        = 'jet';
    if strcmp(sub,'NY787')
        view_angle  = [90,0]; %Lat R        
    else
        view_angle  = [270,0]; %Lat L
    end    
    
    %3.2.1) Plot electrodes as spheres on MNI brain with each sphere color-coded
    %depending on Kprime-value    
    if plot_SubplotperTW == 0 %one separate figure per TW
        
        Figtitle = {[sub ' - ' input_DataLabel ' - ToneDur: ' ...
            tonedur_title ' - ' win_title] [effect_title] [sign_title]};
        
        if strcmp(sub,'NY723')
            NASTD_ECoG_Plot_PlotSignElecsSurf_SubplotHemis...
                (coords,vals,sign_index,chanSize,clims,cmap) %Both H
            suptitle(Figtitle)
        else
            NASTD_ECoG_Plot_PlotSignElecsSurf...
                (coords,vals,sign_index,chanSize,clims,cmap,view_angle) 
            title(Figtitle)
        end
        
        if save_poststepFigs == 1
            path_fig2 = [path_fig1 'signKprime3DSurf/'];
            mkdir([path_fig2]);
            filename     = ['3DSurf_SignExpKprime_' sub '_' ...
                input_DataLabel '_' tonedur_title '_' ...
                Label_TW num2str(i_win) '_p' num2str(pval_plotting) ...
                '_' param.ElecSelect 'elec.png'];
            figfile      = [path_fig2 filename];
            saveas(gcf, figfile, 'png'); %save png version
            close all;
            
        end
        
    elseif plot_SubplotperTW == 1 %one subplot per TW
        
        Figtitle = {[sub ' - ' input_DataLabel ' - ToneDur: ' ...
            tonedur_title] [effect_title]};
        
        if strcmp(sub,'NY723')
            NASTD_ECoG_Plot_SubplotSignElecsSurf_SubplotHemis(...
                coords,vals,sign_index,chanSize,clims,cmap,DimSubplot,...
                CounterSubplot,[win_title sign_title]);
            CounterSubplot = CounterSubplot+2;            
        else
            NASTD_ECoG_Plot_SubplotSignElecsSurf(...
                coords,vals,sign_index,chanSize,clims,cmap,view_angle,...
                DimSubplot,CounterSubplot);
            title([win_title sign_title],'Interpreter','none','FontSize',8);            
            CounterSubplot = CounterSubplot+1;           
                        
        end
    end
end

if plot_SubplotperTW == 1 %one subplot per TW
    sgtitle(Figtitle,'Interpreter', 'none');
    
    if save_poststepFigs == 1
        path_fig2 = [path_fig1 FDR_label '/'];
        if (~exist(path_fig2, 'dir')); mkdir(path_fig2); end
        
        filename     = ['3DSurf_SignExpKprime_' sub '_' input_DataLabel ...
            '_' tonedur_title '_' Label_TW 'all_p' num2str(pval_plotting) ...
            FDR_label '_' param.ElecSelect 'elec.png'];
        figfile      = [path_fig2 filename];
        saveas(gcf, figfile, 'png'); %save png version
        close all;
    end
end

%% 4) Plot exp Kprime values averaged across TW on surface

%Average exp and shuff Kprime values across all TW
for i_win = 1:size(windows,1)
    ExpKprime_data_perTW(:,i_win) = ...
        Kprime_data.Exp.avgRuns.Kprime_Avg{i_win};
    ShuffKprime_data_perTW(i_win, :, :) = ...
        Kprime_data.Shuff.perRepavgRuns.Kprime_NullDistribution{i_win};
end
ExpKprime_data_avgTW    = mean(ExpKprime_data_perTW,2);
ShuffKprime_data_avgTW  = squeeze(mean(ShuffKprime_data_perTW,1));

%Recompute p-values based on across-TW averages
for i_sensor = 1:length(ExpKprime_data_avgTW)
    Kprime_data.Exp.avgRuns.pval_ExpvsShuff_AvgTW(i_sensor,1) = ...
        sum(ShuffKprime_data_avgTW(i_sensor, :) >= ...
        ExpKprime_data_avgTW(i_sensor)) ...
        / length(ShuffKprime_data_avgTW(i_sensor, :));   
end

%Determine significance (FDR or non-FDR thresholding)
if FDR_correct == 1
    pval_FDRcorrected       = mafdr(...
        Kprime_data.Exp.avgRuns.pval_ExpvsShuff_AvgTW(ind_selelecs),'BHFDR', true);
    Index_SuprathresElecs   = pval_FDRcorrected < pval_plotting;
    FDR_label               = 'FDRcorrected';
else
    Index_SuprathresElecs   = ...
        Kprime_data.Exp.avgRuns.pval_ExpvsShuff_AvgTW(ind_selelecs) < pval_plotting;
    FDR_label               = 'uncorrected';
end
    sumSignELec             = sum(Index_SuprathresElecs);
  
    temp_Kprime_SuprathresElecs = [];
    temp_Kprime_SuprathresElecs = ExpKprime_data_avgTW(ind_selelecs);
    Kprime_SuprathresElecs      = temp_Kprime_SuprathresElecs(Index_SuprathresElecs);
   
%Determine title information
effect_title = ['Thresholded Exp Kprime values (Exp vs. Shuff; p < ' ...
    num2str(pval_plotting) '; ' FDR_label ')'];
win_title = ['Averaged across TW - max/mean(SD) Kprime = ' ...
    num2str(round(max(Kprime_SuprathresElecs),2)) ' / ' ...
    num2str(round(mean(Kprime_SuprathresElecs),2)) ...
    ' (' num2str(round(std(Kprime_SuprathresElecs),2)) ')'];
sign_title = [num2str(sumSignELec)...
    ' / ' num2str(length(ind_selelecs_MNI)) ' supra-thresh. elecs'];

%Set up data struct for plotting
%Plot electrodes as spheres on MNI brain with each sphere color-coded depending on Kprime-value
coords      = DataClean_AllTrials.elec.chanpos(ind_selelecs_MNI,:); %MNI coordinates for selected electrodes
vals        = ExpKprime_data_avgTW(ind_selelecs);
sign_index  = Index_SuprathresElecs;   
    %     clims       = [0 max([dv_absmax, 1.1])]; %free scaling
clims       = [0 15]; %fixed scaling
cmap        = 'jet';
chanSize    = (vals./vals)*3; %electrode size (arbitrary)
% %Optional: plot higher Kprime values bigger
% for i_elec = 1:length(vals)
%     chanSize(i_elec,1) = ...
%         (vals(i_elec)./vals(i_elec)) * ...
%         (abs(log10(vals(i_elec)))*4);
% end


if strcmp(sub,'NY787')
    view_angle  = [90,0]; %Lat R
else
    view_angle  = [270,0]; %Lat L
end

Figtitle = {[sub ' - ' input_DataLabel ' - ToneDur: ' tonedur_title] ...
    [win_title] [effect_title] [sign_title]};

if strcmp(sub,'NY723')
    NASTD_ECoG_Plot_PlotSignElecsSurf_SubplotHemis...
        (coords,vals,sign_index,chanSize,clims,cmap) %Both H
    sgtitle(Figtitle,'Interpreter','none')
else
    NASTD_ECoG_Plot_PlotSignElecsSurf...
        (coords,vals,sign_index,chanSize,clims,cmap,view_angle)
    title(Figtitle,'Interpreter','none')
end

if save_poststepFigs == 1
    filename     = ['3DSurf_SignExpKprime_' sub '_' input_DataLabel ...
        '_' tonedur_title '_' Label_TW 'Avgall_p' num2str(pval_plotting) ...
        FDR_label '_' param.ElecSelect 'elec.png'];
     figfile      = [path_fig2 filename];
    saveas(gcf, figfile, 'png'); %save png version
    close all;    
end

end