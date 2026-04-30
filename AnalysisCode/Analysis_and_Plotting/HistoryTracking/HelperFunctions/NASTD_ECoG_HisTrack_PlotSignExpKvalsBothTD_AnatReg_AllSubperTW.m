function NASTD_ECoG_HisTrack_PlotSignExpKvalsBothTD_AnatReg_AllSubperTW...
    (subs, input_ToneDurLabel, input_DataLabel, ...
    SamplesTW, Label_TW, ...
    FDR_correct, pval_threshold,...
    param, ...
    save_poststepFigs, ...
    paths_NASTD_ECoG)

%Aim: Contrast and plot group-level exp Kprime values across both TD for
%anatomical parcellations. Both per TW and averaged across TW.
%Both for all electrodes per parcellation or only sign. electrodes.

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

if save_poststepFigs == 1
    path_fig_elecparcel = [paths_NASTD_ECoG.ECoGdata_ElecParcellation ...
        'Allsub_n' num2str(length(subs)) '/Figs/'];
    if (~exist(path_fig_elecparcel, 'dir')); mkdir(path_fig_elecparcel); end
    path_fig_Kperelecparcel = [paths_NASTD_ECoG.ECoGdata_HisTrack ...
        'Allsub_n' num2str(length(subs)) '/ExpvsShuff/' Label_TW '/Figs/Kprime_AnatReg/'];
    if (~exist(path_fig_Kperelecparcel, 'dir')); mkdir(path_fig_Kperelecparcel); end
end

%Define colormap for later electrode-on-surface plotting
map = [1 1 0;%'FrontalLobe'
    0.93 0.69 0.12;%'PrecentralGyrus'
    0.85 0.33 0.1;%'InferiorFrontalGyrus'
    0 0 1;%'TemporalLobe'
    0.3 0.75 0.93;%'STG'
    0.49 0.18 0.56;%'MTG'
    0 1 0;%'ParietalLobe'
    0.47 0.67 0.19;%'PostcentralGyrus'
    0.25 0.6 0.5;%'SupramarginalGyrus'
    0 0 0];%'OccipitalLobe'

%% 1) Define analysis parameters and windows based on shorter TD
%Define time window (TW) parameters
for i_TD = 1:length(input_ToneDurLabel)
    fsample                 = 512;
    toneDur_inSecs{i_TD}    = str2num(input_ToneDurLabel{i_TD});
    nSamplesPerTone{i_TD}   = toneDur_inSecs{i_TD} * fsample;
    
    win_size                = SamplesTW;
    s_per_win               = (1/fsample)*win_size;
    win_overlap             = 0;
    
    windows{i_TD}           = [1 win_size];
    while windows{i_TD}(end,end) < nSamplesPerTone{i_TD}
        windows{i_TD}       = ...
            [windows{i_TD}; windows{i_TD}(end,:) + (win_size - win_overlap)];
    end
    
    if windows{i_TD}(end,end) > nSamplesPerTone{i_TD}
        windows{i_TD}(end,:) = [];
    end
    num_win(i_TD,1) = length(windows{i_TD});
    num_win_common = min(num_win);
end


%% 2) Load and aggregate Kprime data across TD and subjects
%Create empty proxy struct
data_selElec_allSubs                    = struct;
data_selElec_allSubs.label              = [];
data_selElec_allSubs.sub                = [];
data_selElec_allSubs.elecpos            = [];
data_selElec_allSubs.allelecs_persub    = [];
data_selElec_allSubs.selelecs_persub    = [];
data_selElec_allSubs.label_AnatCat      = [];
data_selElec_allSubs.index_AnatCat      = [];
data_selElec_allSubs.filt_signelecs_StimCorr_allsub = [];

for i_sub = 1:length(subs)
    
    sub = subs{i_sub};
    
    %Load preprocessed data for elec coords
    tic
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    tic
    disp([' -- Loading preprocessed data set for sub: ' sub])
    load([si.path_preprocdata_sub], 'DataClean_AllTrials');
    %     disp([' -- Done loading in ' num2str(toc) ' sec'])
    
    for i_TD = 1:length(input_ToneDurLabel)
        if strcmp(input_ToneDurLabel{i_TD},'0.2')
            tonedur_title = '200msTD';
        elseif  strcmp(input_ToneDurLabel{i_TD},'0.4')
            tonedur_title = '400msTD';
        end
        
        %Load Kprime data per TD
        path_inputdata = [paths_NASTD_ECoG.ECoGdata_HisTrack sub '/ExpvsShuff/' Label_TW '/'];
        load([path_inputdata sub '_' input_DataLabel '_' tonedur_title ...
            '_HisTrack_ExpvsShuffKprimeCombFolds_unbalanced.mat'], ...
            'Kprime_data', 'labels_loadedData'); %filename: Kprime_data
        
        %Place Kprime-values per subject, window, TD for every electrode in common struct
        for i_win = 1:length(Kprime_data.Exp.avgRuns.Kprime_Avg)
            data_selElec_allSubs.kprime_allelecs_perSubWinTD{i_sub}{i_win}{i_TD} = ...
                Kprime_data.Exp.avgRuns.Kprime_Avg{i_win};
            data_selElec_allSubs.pval_allelecs_perSubWinTD{i_sub}{i_win}{i_TD} = ...
                Kprime_data.Exp.avgRuns.pval_ExpvsShuff{i_win};
        end
        %Average exp and shuff Kprime values across windows to recompute
        %sign. when averaging across windows
        for i_win = 1:length(Kprime_data.Exp.avgRuns.Kprime_Avg)
            temp_ExpKprime_data_pTW(:,i_win) = ...
                Kprime_data.Exp.avgRuns.Kprime_Avg{i_win};
            temp_ShuffKprime_data_pTW(i_win, :, :) = ...
                Kprime_data.Shuff.perRepavgRuns.Kprime_NullDistribution{i_win};
        end
        data_selElec_allSubs.kprime_allelecs_perSubTD_avgWin{i_sub}{i_TD} = ...
            mean(temp_ExpKprime_data_pTW,2);
        temp_ShuffKprime_data_avgTW  = squeeze(mean(temp_ShuffKprime_data_pTW,1));
        %Recompute p-values based on across-TW averages
        for i_sensor = 1:length(data_selElec_allSubs.kprime_allelecs_perSubTD_avgWin{i_sub}{i_TD})
            data_selElec_allSubs.pval_allelecs_perSubTD_avgWin{i_sub}{i_TD}(i_sensor,1) = ...
                sum(temp_ShuffKprime_data_avgTW(i_sensor, :) >= ...
                data_selElec_allSubs.kprime_allelecs_perSubTD_avgWin{i_sub}{i_TD}(i_sensor)) ...
                / length(temp_ShuffKprime_data_avgTW(i_sensor, :));
        end
        clear temp*
        
        %Load sequence tracking data and select current respective sign. electrodes per TD
        if strcmp(param.ElecSelect, 'StimCorr')
            path_inputdata = [paths_NASTD_ECoG.ECoGdata_StimCorr  sub '/Data/'];
            load([path_inputdata sub '_StimCorr_' input_DataLabel '_' ...
                input_ToneDurLabel{i_TD} 'sTD.mat'], ...
                'corr_ttest', 'SensorLabels');
            filt_signelecs_StimCorr_perTD(:,i_TD) = corr_ttest.p < 0.05; %only sign. stim corr elecs
        else
            filt_signelecs_StimCorr_perTD(:,i_TD) = true(length(labels_loadedData),1);  %all elecs
        end
        clear ExpKprime_data corr_ttest SensorLabels
    end
    
    %Label electrodes used per subject (before Sequence tracking electrode selection)
    data_selElec_allSubs.allelecs_persub = ...
        [data_selElec_allSubs.allelecs_persub; ...
        ones(length(labels_loadedData),1)*i_sub];
    
    %Combine sequence tracking electrode-selection across TD,
    %Select electrodes with sign. sequence tracking effect in at least one TD (liberal criterion)
    filt_signelecs_StimCorr_perSub = ...
        filt_signelecs_StimCorr_perTD(:,1) | filt_signelecs_StimCorr_perTD(:,2);
    ind_selelecs  = find(filt_signelecs_StimCorr_perSub);
    %nSensors_sel  = sum(filt_signelecs_StimCorr);
    %labels_loadedData(ind_selelecs)
    clear filt_signelecs_StimCorr_perTD
    
    %Determine MNI index of selected electrodes to select MNI coordinates
    ind_selelecs_MNI   = [];
    for i_selelec = 1:length(ind_selelecs)
        ind_selelecs_MNI = [ind_selelecs_MNI; ...
            find(strcmp(DataClean_AllTrials.elec.label, ...
            labels_loadedData{ind_selelecs(i_selelec)}))];
    end
    data_selElec_allSubs.elecpos = ...
        [data_selElec_allSubs.elecpos; ...
        DataClean_AllTrials.elec.chanpos(ind_selelecs_MNI,:)];
    clear ind_selelecs_MNI
    
    %Read out selected elecs per sub, Add sub-identifier to each elec label,
    %and append selected electrodes to one common summary struct
    label_selelecs_currSub  = {};
    for i_selelec = 1:length(ind_selelecs)
        label_selelecs_currSub{i_selelec,1} = ...
            strcat(labels_loadedData{ind_selelecs(i_selelec)},'_', sub);
    end
    
    %Aggregate identifiers and labels across subs
    data_selElec_allSubs.label = ...
        [data_selElec_allSubs.label; label_selelecs_currSub];
    data_selElec_allSubs.sub = ...
        [data_selElec_allSubs.sub; ones(length(label_selelecs_currSub),1)*i_sub];
    data_selElec_allSubs.selelecs_persub{i_sub} = ...
        find(data_selElec_allSubs.sub == i_sub);
    
    clear ind_selelecs label_selelecs_currSub
    
    %Aggregate sequence tracking electrode selection filter over subjects
    data_selElec_allSubs.filt_signelecs_StimCorr_allsub = ...
        [data_selElec_allSubs.filt_signelecs_StimCorr_allsub; ...
        filt_signelecs_StimCorr_perSub];
    data_selElec_allSubs.filt_signelecs_StimCorr_allsub = ...
        logical(data_selElec_allSubs.filt_signelecs_StimCorr_allsub);
    
    %Aggregate sequence tracking selected electrodes over subjects for each TD
    %Per TW
    for i_TD = 1:length(input_ToneDurLabel)
        for i_win = 1:num_win(i_TD)
            if i_sub == 1
                data_selElec_allSubs.kprime_selelecs_perTD{i_win}{i_TD} = ...
                    data_selElec_allSubs.kprime_allelecs_perSubWinTD{i_sub}{i_win}{i_TD} ...
                    (filt_signelecs_StimCorr_perSub);
            else
                data_selElec_allSubs.kprime_selelecs_perTD{i_win}{i_TD} = ...
                    [data_selElec_allSubs.kprime_selelecs_perTD{i_win}{i_TD}; ...
                    data_selElec_allSubs.kprime_allelecs_perSubWinTD{i_sub}{i_win}{i_TD}(filt_signelecs_StimCorr_perSub)];
            end
        end
    end
    %Averaged across TW
    for i_TD = 1:length(input_ToneDurLabel)
        if i_sub == 1
            data_selElec_allSubs.kprime_selelecs_perTD_avgWin{i_TD} = ...
                data_selElec_allSubs.kprime_allelecs_perSubTD_avgWin{i_sub}{i_TD} ...
                (filt_signelecs_StimCorr_perSub);
        else
            data_selElec_allSubs.kprime_selelecs_perTD_avgWin{i_TD} = ...
                [data_selElec_allSubs.kprime_selelecs_perTD_avgWin{i_TD}; ...
                data_selElec_allSubs.kprime_allelecs_perSubTD_avgWin{i_sub}{i_TD}(filt_signelecs_StimCorr_perSub)];
        end
    end
    
    
    %% 3) Determine sign. SHI electrodes per subject
    for i_TD = 1:length(input_ToneDurLabel)
        %per TW
        for i_win = 1:num_win(i_TD)
            %Determine significance (FDR or non-FDR thresholding) for selected elecs
            if FDR_correct == 1
                FDR_label = 'FDRcorrected';
                pval.(FDR_label){i_TD}{i_win} = ...
                    mafdr(data_selElec_allSubs.pval_allelecs_perSubWinTD{i_sub}{i_win}{i_TD}...
                    (filt_signelecs_StimCorr_perSub), ...
                    'BHFDR', true);
                filt_SuprathresElecs{i_TD}{i_win} = pval.(FDR_label){i_TD}{i_win} < pval_threshold;
            else
                FDR_label = 'uncorrected';
                pval.(FDR_label){i_TD}{i_win} = ...
                    data_selElec_allSubs.pval_allelecs_perSubWinTD{i_sub}{i_win}{i_TD}...
                    (filt_signelecs_StimCorr_perSub);
                filt_SuprathresElecs{i_TD}{i_win} = pval.(FDR_label){i_TD}{i_win} < pval_threshold;
            end
            %Aggregate sign. filter and min p-values across subjects
            if i_sub == 1
                data_selElec_allSubs.filt_signElecs{i_TD}{i_win} = [];
                data_selElec_allSubs.pval_signElecs{i_TD}{i_win} = [];
            end
            data_selElec_allSubs.filt_signElecs{i_TD}{i_win} = ...
                [data_selElec_allSubs.filt_signElecs{i_TD}{i_win}; filt_SuprathresElecs{i_TD}{i_win}];
            data_selElec_allSubs.pval_signElecs{i_TD}{i_win} = ...
                [data_selElec_allSubs.pval_signElecs{i_TD}{i_win}; pval.(FDR_label){i_TD}{i_win}];
        end
        %Avg across TW
        %Determine significance (FDR or non-FDR thresholding) for selected elecs
        if FDR_correct == 1
            pval_avgWin.(FDR_label){i_TD} = ...
                mafdr(data_selElec_allSubs.pval_allelecs_perSubTD_avgWin{i_sub}{i_TD}...
                (filt_signelecs_StimCorr_perSub), ...
                'BHFDR', true);
            filt_SuprathresElecs_avgWin{i_TD} = ...
                pval_avgWin.(FDR_label){i_TD} < pval_threshold;
        else
            pval_avgWin.(FDR_label){i_TD} = ...
                data_selElec_allSubs.pval_allelecs_perSubTD_avgWin{i_sub}{i_TD}...
                (filt_signelecs_StimCorr_perSub);
            filt_SuprathresElecs_avgWin{i_TD} = pval_avgWin.(FDR_label){i_TD} < pval_threshold;
        end
        %Aggregate sign. filter and min p-values across subjects
        if i_sub == 1
            data_selElec_allSubs.filt_signElecs_avgWin{i_TD} = [];
            data_selElec_allSubs.pval_signElecs_avgWin{i_TD} = [];
        end
        data_selElec_allSubs.filt_signElecs_avgWin{i_TD} = ...
            [data_selElec_allSubs.filt_signElecs_avgWin{i_TD}; ...
            filt_SuprathresElecs_avgWin{i_TD}];
        data_selElec_allSubs.pval_signElecs_avgWin{i_TD} = ...
            [data_selElec_allSubs.pval_signElecs_avgWin{i_TD}; ...
            pval_avgWin.(FDR_label){i_TD}];
    end
    
    %Combine across TD and select all elecs that show sign. SHI in at least one TD (liberal criterion)
    %per win
    for i_win = 1:num_win_common
        temp_pval_bothTD = [pval.(FDR_label){1}{i_win}, pval.(FDR_label){2}{i_win}];
        pval.(FDR_label){3}{i_win} = min(temp_pval_bothTD,[],2);
        filt_SuprathresElecs{3}{i_win} = ...
            pval.(FDR_label){3}{i_win} < pval_threshold;
        %sumSignELec{i_win} = sum(filt_SuprathresElecs{i_win}{3});
        clear temp_pval_bothTD
        %Aggregate sign. filter and min p-values across subjects
        if i_sub == 1
            data_selElec_allSubs.filt_signElecs{3}{i_win} = [];
            data_selElec_allSubs.pval_signElecs{3}{i_win} = [];
        end
        data_selElec_allSubs.filt_signElecs{3}{i_win} = ...
            [data_selElec_allSubs.filt_signElecs{3}{i_win}; ...
            filt_SuprathresElecs{3}{i_win}];
        data_selElec_allSubs.pval_signElecs{3}{i_win} = ...
            [data_selElec_allSubs.pval_signElecs{3}{i_win}; ...
            pval.(FDR_label){3}{i_win}];
    end
    clear filt_SuprathresElecs
    
    %Avg across win and select all elecs that show sign. SHI in at least one TW (liberal criterion)
    temp_pval_bothTD = [pval_avgWin.(FDR_label){1}, pval_avgWin.(FDR_label){2}];
    pval_avgWin.(FDR_label){3} = min(temp_pval_bothTD,[],2);
    filt_SuprathresElecs_avgWin{3} = ...
        pval_avgWin.(FDR_label){3} < pval_threshold;
    clear temp_pval_bothTD
    %Aggregate sign. filter and min p-values across subjects
    if i_sub == 1
        data_selElec_allSubs.filt_signElecs_avgWin{3} = [];
        data_selElec_allSubs.pval_signElecs_avgWin{3} = [];
    end
    data_selElec_allSubs.filt_signElecs_avgWin{3} = ...
        [data_selElec_allSubs.filt_signElecs_avgWin{3}; ...
        filt_SuprathresElecs_avgWin{3}];
    data_selElec_allSubs.pval_signElecs_avgWin{3} = ...
        [data_selElec_allSubs.pval_signElecs_avgWin{3}; ...
        pval_avgWin.(FDR_label){3}];
    
    clear filt_SuprathresElecs_avgWin
    
    
    %% 4) Categorize electrodes according to anatomical regions
    AnatReg_allSubs{i_sub} = ...
        NASTD_ECoG_AssignAnatRegions(DataClean_AllTrials, labels_loadedData);
    clear DataClean_AllTrials labels_loadedData
    
    %Restrict anatomical parcellation output to selected electrodes
    %(e.g., sequence tracking selection)
    temp_filt_allelecs_currsub = ...
        find(data_selElec_allSubs.allelecs_persub == i_sub);
    temp_filt_selelecs_currsub = ...
        data_selElec_allSubs.filt_signelecs_StimCorr_allsub(temp_filt_allelecs_currsub);
    
    AnatReg_allSubs{i_sub}.CatIndex = ...
        AnatReg_allSubs{i_sub}.CatIndex(temp_filt_selelecs_currsub,:);
    AnatReg_allSubs{i_sub}.Info_perelec = ...
        AnatReg_allSubs{i_sub}.Info_perelec(temp_filt_selelecs_currsub,:);
    clear temp_filt*
    
    %Aggregate category labels and indices across subjects
    data_selElec_allSubs.label_AnatCat  = ...
        [data_selElec_allSubs.label_AnatCat; ...
        AnatReg_allSubs{i_sub}.Info_perelec(:,3)];
    data_selElec_allSubs.index_AnatCat  = ...
        [data_selElec_allSubs.index_AnatCat; ...
        AnatReg_allSubs{i_sub}.CatIndex];
    
    
    clear pval filt_SuprathresElecs filt_signelecs_StimCorr_perSub
    disp(['  -- ' sub ' Kprime processing finished in ' num2str(round(toc)) ' sec --'])
    
    %% 5) Plot surface with electrode names and categories per subject
    %     %Create a colormap that nicely distinguishes categories but keeps
    %     %related categories similar
    %
    %     %set up plotting struct
    %     clear dat
    %     dat.dimord  = 'chan_time';
    %     dat.time    = 0;
    %     coords      = data_selElec_allSubs.elecpos(data_selElec_allSubs.selelecs_persub{i_sub},:); %MNI coordinates for selected electrodes
    %     vals        = max(data_selElec_allSubs.index_AnatCat(data_selElec_allSubs.selelecs_persub{i_sub},:),[],2); %Finest category distinction
    % %     vals        = data_selElec_allSubs.index_AnatCat(data_selElec_allSubs.selelecs_persub{i_sub},:);%coarsest category distinction
    %     sign_index  = logical(ones(length(data_selElec_allSubs.selelecs_persub{i_sub}),1));
    %     label       = data_selElec_allSubs.label_AnatCat(data_selElec_allSubs.selelecs_persub{i_sub},:);
    %     chanSize    = (vals./vals)*2.5; %electrode size (arbitrary)
    %     clims       = [1 10]; %fixed scaling
    %     cmap        = map;
    %     if strcmp(sub,'NY787')
    %         view_angle  = [90,0]; %Lat R
    %     else
    %         view_angle  = [270,0]; %Lat L
    %     end
    %
    %     %plot surface brain
    %     f = figure;
    %     set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    %     NASTD_ECoG_Plot_SubplotSignElecsSurf_Label...
    %         (coords, label, vals, sign_index,...
    %         chanSize, clims, cmap, view_angle, ...
    %         [1 1], 1, [0 0 0 0], [0 0 0 0]);
    %
    %     %Add colorbar
    %     c = colorbar;
    %     c.Limits = [1 10];
    %     colorbar_labels = {};
    %     for i_cat = 1:length(AnatReg_allSubs{i_sub}.CatLabels)
    %         colorbar_labels{end +1,1} = [AnatReg_allSubs{i_sub}.CatLabels{i_cat} ...
    %             ' - (' num2str(sum(max(...
    %             data_selElec_allSubs.index_AnatCat(...
    %             data_selElec_allSubs.selelecs_persub{i_sub},:),[],2) == i_cat)) ...
    %             ') elecs'];
    %     end
    %     c.TickLabels = colorbar_labels;
    %     c.Position(1) = c.Position(1)-0.1;
    %
    %     %Add title
    %     title({[sub ' - Anatomical parcellation'] ...
    %         [num2str(length(data_selElec_allSubs.selelecs_persub{i_sub})) ' electrodes total']})
    %
    %     if save_poststepFigs == 1
    %         filename     = ['3DSurf_' sub '_ElecParcel_' param.ElecSelect 'elec.png'];
    %         figfile      = [path_fig filename];
    %         saveas(gcf, figfile, 'png'); %save png version
    %         close all;
    %     end
    
    
end

%Convert Filter to logical and count sign. elecs
for i_TD = 1:3
    %Per TW
    for i_win = 1:length(data_selElec_allSubs.filt_signElecs{i_TD})
        data_selElec_allSubs.filt_signElecs{i_TD}{i_win} = ...
            logical(data_selElec_allSubs.filt_signElecs{i_TD}{i_win});
        sumSignELec{i_TD}{i_win} = sum(data_selElec_allSubs.filt_signElecs{i_TD}{i_win});
    end
    %Averaged across TW
    data_selElec_allSubs.filt_signElecs_avgWin{i_TD} = ...
        logical(data_selElec_allSubs.filt_signElecs_avgWin{i_TD});
    sumSignELec_avgWin{i_TD} = sum(data_selElec_allSubs.filt_signElecs_avgWin{i_TD});
end

% %% 7) Plot anatomical electrode parcellation across all subjects
% %set up plotting struct
% clear dat
% dat.dimord  = 'chan_time';
% dat.time    = 0;
% coords      = data_selElec_allSubs.elecpos; %MNI coordinates for selected electrodes
% vals        = max(data_selElec_allSubs.index_AnatCat,[],2); %Finest category distinction
% %     vals        = data_selElec_allSubs.index_AnatCat;%coarsest category distinction
% sign_index  = logical(ones(length(data_selElec_allSubs.elecpos),1));
% label       = data_selElec_allSubs.label_AnatCat;
% chanSize    = (vals./vals)*2; %electrode size (arbitrary)
% clims       = [1 10]; %fixed scaling
% cmap        = map;
% 
% %plot surface brain
% f = figure; hold on;
% 
% view_angle  = [270,0]; %Lat L
% set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
% NASTD_ECoG_Plot_SubplotSignElecsSurf_Label...
%     (coords, label, vals, sign_index,...
%     chanSize, clims, cmap, view_angle, ...
%     [1 2], 1, [0 0 0 0], [0 0 0 0]);
% 
% view_angle  = [90,0]; %Lat R
% NASTD_ECoG_Plot_SubplotSignElecsSurf_Label...
%     (coords, label, vals, sign_index,...
%     chanSize, clims, cmap, view_angle, ...
%     [1 2], 2, [0 0 0 0], [0 0 0 0]);
% 
% %Add colorbar
% c               = colorbar;
% c.Limits        = [1 10];
% colorbar_labels = {};
% for i_cat = 1:length(AnatReg_allSubs{i_sub}.CatLabels)
%     colorbar_labels{end +1,1} = [AnatReg_allSubs{i_sub}.CatLabels{i_cat} ...
%         ' - (' num2str(sum(max(data_selElec_allSubs.index_AnatCat,[],2) == i_cat)) ') elecs'];
% end
% c.TickLabels    = colorbar_labels;
% c.Position(1)   = c.Position(1)-0.4;
% c.Position(4)   = c.Position(1)*1.6;
% c.AxisLocation  = 'in';
% c.FontSize      = 10;
% 
% %Add title
% sgtitle({['All subs (n = ' num2str(length(subs)) ', ' FDR_label ') - Anatomical parcellation'] ...
%     [num2str(length(data_selElec_allSubs.index_AnatCat)) ' electrodes total']})
% 
% if save_poststepFigs == 1
%     filename     = ['3DSurf_Allsubn' num2str(length(subs)) '_ElecParcel_' param.ElecSelect 'elec.png'];
%     figfile      = [path_fig_elecparcel filename];
%     saveas(gcf, figfile, 'png'); %save png version
%     close all;
% end


%% 8) Plot Kprime values per anatomical parcellation for each TW and TD
%Plot barcharts of group-average and single-subject Kprime per time window across parcellations

%Read out Kprime data per win and per electrode parcellation
%8.1 Per subject
plotstruct = [];
clear sum_signelec_percat_allsub
for i_TD = 1:length(input_ToneDurLabel)
    
    sum_signelec_percat_allsub{i_TD} = ...
        zeros(AnatReg_allSubs{i_sub}.Num_Cat, num_win(i_TD));
    
    for i_sub = 1:length(subs)
        for i_win = 1:num_win(i_TD)
            for i_cat = 1:AnatReg_allSubs{i_sub}.Num_Cat
                
                %All elecs per cat
                filt_elecpercat = max(AnatReg_allSubs{i_sub}.CatIndex,[],2) == i_cat; %Parcellation filter
                
                plotstruct.avgKprime_percatsub.allelecs{i_TD}{i_sub}(i_cat, i_win) = ...
                    nanmean(data_selElec_allSubs.kprime_selelecs_perTD{i_win}{i_TD} ...
                    (data_selElec_allSubs.selelecs_persub{i_sub}(filt_elecpercat)));
                plotstruct.STDKprime_percatsub.allelecs{i_TD}{i_sub}(i_cat, i_win) = ...
                    nanstd(data_selElec_allSubs.kprime_selelecs_perTD{i_win}{i_TD} ...
                    (data_selElec_allSubs.selelecs_persub{i_sub}(filt_elecpercat)));
                
                %Sign. elecs per cat
                filt_elecpercat = max(AnatReg_allSubs{i_sub}.CatIndex,[],2) == i_cat; %Parcellation filter
                filt_signelec   = data_selElec_allSubs.filt_signElecs{i_TD}{i_win}(data_selElec_allSubs.selelecs_persub{i_sub}); %Sign. elec filter
                filt_signcat    = and(filt_elecpercat, filt_signelec); %Combined parcellation + sign. filter
                
                sum_signelec_percat_allsub{i_TD}(i_cat, i_win) = ...
                    sum_signelec_percat_allsub{i_TD}(i_cat, i_win) + sum(filt_signcat);
                
                plotstruct.avgKprime_percatsub.signelecs{i_TD}{i_sub}(i_cat, i_win) = ...
                    nanmean(data_selElec_allSubs.kprime_selelecs_perTD{i_win}{i_TD} ...
                    (data_selElec_allSubs.selelecs_persub{i_sub}(filt_signcat)));
                plotstruct.STDKprime_percatsub.signelecs{i_TD}{i_sub}(i_cat, i_win) = ...
                    nanstd(data_selElec_allSubs.kprime_selelecs_perTD{i_win}{i_TD} ...
                    (data_selElec_allSubs.selelecs_persub{i_sub}(filt_elecpercat)));
                
            end
        end
    end
end

%8.2 Group average across all subs
for i_TD = 1:length(input_ToneDurLabel)
    for i_cat = 1:AnatReg_allSubs{i_sub}.Num_Cat
        
        %All elecs per cat
        temp_Kprime_allsub = [];
        for i_sub = 1:length(subs)
            %Aggregate avg Kprime across elecs per cat for all subs
            temp_Kprime_allsub(i_sub,:) = ...
                plotstruct.avgKprime_percatsub.allelecs{i_TD}{i_sub}(i_cat,:);
        end
        %Average across subs
        plotstruct.avgKprime_percatallsub.allelecs{i_TD}(i_cat,:) = ...
            nanmean(temp_Kprime_allsub,1);
        %Compute SEM across subs. If only 1 sub, use STD across elecs.
        if length(subs) > 1
            plotstruct.SEMKprime_percatallsub.allelecs{i_TD}(i_cat,:) = ...
                nanstd(temp_Kprime_allsub,[],1) / sqrt(length(subs));
        else
            plotstruct.SEMKprime_percatallsub.allelecs{i_TD}(i_cat,:) = ...
                plotstruct.STDKprime_percatsub.allelecs{i_TD}{i_sub}(i_cat, :);
        end
        
        %Sign. elecs per cat
        temp_Kprime_allsub = [];
        for i_sub = 1:length(subs)
            %Aggregate avg Kprime across elecs per cat for all sub
            temp_Kprime_allsub(i_sub,:) = ...
                plotstruct.avgKprime_percatsub.signelecs{i_TD}{i_sub}(i_cat,:);
        end
        %Average across subs
        plotstruct.avgKprime_percatallsub.signelecs{i_TD}(i_cat,:) = ...
            nanmean(temp_Kprime_allsub,1);
        %Compute SEM across subs. If only 1 sub, use STD across elecs.
        if length(subs) > 1
            plotstruct.SEMKprime_percatallsub.signelecs{i_TD}(i_cat,:) = ...
                nanstd(temp_Kprime_allsub,[],1) / sqrt(length(subs));
        else
            plotstruct.SEMKprime_percatallsub.signelecs{i_TD}(i_cat,:) = ...
                plotstruct.STDKprime_percatsub.signelecs{i_TD}{i_sub}(i_cat, :);
        end
    end
end
clear temp_Kprime_allsub

%8.3 Plot barchart showing Kprime per parcellation and win, per TD
for i_TD = 1:length(input_ToneDurLabel)

    FigName = ['Group_n', num2str(length(subs)), ' AvgKprime per anat parcel - TD: ' input_ToneDurLabel{i_TD} ' s'];
    figure('NumberTitle','off','Name', FigName, ...
        'units','normalized','outerposition',[0 0 1 1]);
    fig_counter = 0;
    temp_num_win = num_win(i_TD);

    for i_win = 1:temp_num_win

        x_labels = {};
        for i_cat = 1:length(AnatReg_allSubs{i_sub}.CatLabels)
            x_labels{end +1,1} = [AnatReg_allSubs{i_sub}.CatLabels{i_cat} ...
                ' (' num2str(sum(max(data_selElec_allSubs.index_AnatCat,[],2) == i_cat)) '/' ...
                num2str(sum_signelec_percat_allsub{i_TD}(i_cat,i_win)) ')'];
        end

        fig_counter = fig_counter +1;

        subplot(round(sqrt(temp_num_win)), round(sqrt(temp_num_win)), fig_counter);
        h = barwitherr( ... %Plots Kprime values for both all and sign. elecs
            [plotstruct.SEMKprime_percatallsub.allelecs{i_TD}(:,i_win), ...
            plotstruct.SEMKprime_percatallsub.signelecs{i_TD}(:,i_win)], ...
            [plotstruct.avgKprime_percatallsub.allelecs{i_TD}(:,i_win), ...
            plotstruct.avgKprime_percatallsub.signelecs{i_TD}(:,i_win)]);
        title({'Kprime (Group Avg + SEM) per anatomical parcellation'},'FontSize',10)

        set(gca,'xticklabel',x_labels,'fontsize',8)
        xtickangle(45)
        ylim([0 10])
        set(gca, 'YGrid', 'on');

        hold on; %Plot single subject data
        for i_sub = 1:length(subs)
            scatter((h(1).XEndPoints(1):AnatReg_allSubs{1}.Num_Cat),...
                plotstruct.avgKprime_percatsub.allelecs{i_TD}{i_sub}(:,i_win)',...
                15,[0 0 0.2],'s', 'filled','MarkerFaceAlpha', 0.5)
            scatter((h(1).XEndPoints(1):AnatReg_allSubs{1}.Num_Cat) + 0.3,...
                plotstruct.avgKprime_percatsub.signelecs{i_TD}{i_sub}(:,i_win)',...
                15,[0.2 0 0],'filled','MarkerFaceAlpha', 0.5)
        end

        title(['TW: ' num2str(i_win)])

        if fig_counter == 1
            f=get(gca,'Children');
            legend([f(end-1),f(end)], 'Sign. elecs', 'All elecs')
        end

    end

    sgtitle({['Kprime (Group Avg (n = ' num2str(length(subs)) ...
        ') + SEM/SD) per anatomical parcellation'], ...
        [input_DataLabel ', TD: ' input_ToneDurLabel{i_TD} ' s, ' ...
        FDR_label ', per TW']}, ...
        'FontSize',14, 'Interpreter', 'none')

    if save_poststepFigs == 1
        filename     = ['ErrBar_Allsubn' num2str(length(subs)) '_' input_DataLabel...
            '_AvgKperElecParcelperTW_TD' num2str(i_TD) '_' ...
            FDR_label '_' param.ElecSelect 'elec.png'];
        path_fig_curr = [path_fig_Kperelecparcel '/Kprime_perTD/'];
        if (~exist(path_fig_curr, 'dir')); mkdir(path_fig_curr); end        
        figfile      = [path_fig_curr filename];
        saveas(gcf, figfile, 'png');
    end
end


%% 9) Plot Kprime values per anatomical parcellation averaged across TW for each TD
%Read out Kprime data per win and per electrode parcellation
%9.1 Per subject
plotstruct = [];
clear sum_signelec_percat_allsub

for i_TD = 1:length(input_ToneDurLabel)
    sum_signelec_percat_allsub{i_TD} = ...
        zeros(AnatReg_allSubs{i_sub}.Num_Cat, 1);
    for i_sub = 1:length(subs)
        for i_cat = 1:AnatReg_allSubs{i_sub}.Num_Cat
            %All elecs per cat
            filt_elecpercat = max(AnatReg_allSubs{i_sub}.CatIndex,[],2) == i_cat; %Parcellation filter
            plotstruct.avgKprime_percatsub_avgWin.allelecs{i_TD}(i_cat, i_sub) = ...
                nanmean(data_selElec_allSubs.kprime_selelecs_perTD_avgWin{i_TD} ...
                (data_selElec_allSubs.selelecs_persub{i_sub}(filt_elecpercat)));
            plotstruct.STDKprime_percatsub_avgWin.allelecs{i_TD}(i_cat, i_sub) = ...
                nanstd(data_selElec_allSubs.kprime_selelecs_perTD_avgWin{i_TD} ...
                (data_selElec_allSubs.selelecs_persub{i_sub}(filt_elecpercat)));
            
            %Sign. elecs per cat
            filt_elecpercat = max(AnatReg_allSubs{i_sub}.CatIndex,[],2) == i_cat; %Parcellation filter
            filt_signelec   = data_selElec_allSubs.filt_signElecs_avgWin{i_TD}(data_selElec_allSubs.selelecs_persub{i_sub}); %Sign. elec filter
            filt_signcat    = and(filt_elecpercat, filt_signelec); %Combined parcellation + sign. filter
            %Count sign. electrodes per cat
            sum_signelec_percat_allsub{i_TD}(i_cat) = ...
                sum_signelec_percat_allsub{i_TD}(i_cat) + sum(filt_signcat);
            plotstruct.avgKprime_percatsub_avgWin.signelecs{i_TD}(i_cat, i_sub) = ...
                nanmean(data_selElec_allSubs.kprime_selelecs_perTD_avgWin{i_TD} ...
                (data_selElec_allSubs.selelecs_persub{i_sub}(filt_signcat)));
            plotstruct.STDKprime_percatsub_avgWin.signelecs{i_TD}(i_cat, i_sub) = ...
                nanstd(data_selElec_allSubs.kprime_selelecs_perTD_avgWin{i_TD} ...
                (data_selElec_allSubs.selelecs_persub{i_sub}(filt_elecpercat)));
        end
    end
end

%9.2 Group average across all subs
for i_TD = 1:length(input_ToneDurLabel)
    for i_cat = 1:AnatReg_allSubs{i_sub}.Num_Cat
        
        %All elecs per cat
        %Average across subs
        plotstruct.avgKprime_percatallsub_avgWin.allelecs{i_TD}(i_cat,1) = ...
            nanmean(plotstruct.avgKprime_percatsub_avgWin.allelecs{i_TD}(i_cat,:),2);
        %Compute SEM across subs. If only 1 sub, use STD across elecs.
        if length(subs) > 1
            plotstruct.SEMKprime_percatallsub_avgWin.allelecs{i_TD}(i_cat,1) = ...
                nanstd(plotstruct.avgKprime_percatsub_avgWin.allelecs{i_TD}(i_cat,:),[],2) / sqrt(length(subs));
        else
            plotstruct.SEMKprime_percatallsub_avgWin.allelecs{i_TD}(i_cat,1) = ...
                plotstruct.STDKprime_percatsub_avgWin.allelecs{i_TD}(i_cat,i_sub);
        end
        
        %Sign. elecs per cat
        %Average across subs
        plotstruct.avgKprime_percatallsub_avgWin.signelecs{i_TD}(i_cat,1) = ...
            nanmean(plotstruct.avgKprime_percatsub_avgWin.signelecs{i_TD}(i_cat,:),2);
        %Compute SEM across subs. If only 1 sub, use STD across elecs.
        if length(subs) > 1
            plotstruct.SEMKprime_percatallsub_avgWin.signelecs{i_TD}(i_cat,1) = ...
                nanstd(plotstruct.avgKprime_percatsub_avgWin.signelecs{i_TD}(i_cat,:),[],2) / sqrt(length(subs));
        else
            plotstruct.SEMKprime_percatallsub_avgWin.signelecs{i_TD}(i_cat,1) = ...
                plotstruct.STDKprime_percatsub_avgWin.signelecs{i_TD}(i_cat,i_sub);
        end
    end
end
clear temp_Kprime_allsub

%9.3 Plot barchart showing Kprime per parcellation per TD, averaged across TW
for i_TD = 1:length(input_ToneDurLabel)

    FigName = ['Group_n', num2str(length(subs)), ' AvgKprime per anat parcel - TD: ' input_ToneDurLabel{i_TD} ' s'];
    figure('NumberTitle','off','Name', FigName, ...
        'units','normalized','outerposition',[0 0 1 1]);

    x_labels = {};
    for i_cat = 1:length(AnatReg_allSubs{i_sub}.CatLabels)
        x_labels{end +1,1} = [AnatReg_allSubs{i_sub}.CatLabels{i_cat} ...
            ' (' num2str(sum(max(data_selElec_allSubs.index_AnatCat,[],2) == i_cat)) '/' ...
            num2str(round(nanmean(sum_signelec_percat_allsub{i_TD}(i_cat,:)))) ')'];
    end

    h = barwitherr( ... %Plots Kprime values for both all and sign. elecs
        [plotstruct.SEMKprime_percatallsub_avgWin.allelecs{i_TD}, ...
        plotstruct.SEMKprime_percatallsub_avgWin.signelecs{i_TD}], ...
        [plotstruct.avgKprime_percatallsub_avgWin.allelecs{i_TD}, ...
        plotstruct.avgKprime_percatallsub_avgWin.signelecs{i_TD}]);

    set(gca,'xticklabel',x_labels,'fontsize',8)
    xtickangle(45)
    ylim([0 10])
    set(gca, 'YGrid', 'on');

    hold on; %Plot single subject data
    for i_sub = 1:length(subs)
        scatter((h(1).XEndPoints(1):AnatReg_allSubs{1}.Num_Cat),...
            plotstruct.avgKprime_percatsub_avgWin.allelecs{i_TD}(:,i_sub)',...
            15,[0 0 0.2],'s', 'filled','MarkerFaceAlpha', 0.5)
        scatter((h(1).XEndPoints(1):AnatReg_allSubs{1}.Num_Cat) + 0.3,...
            plotstruct.avgKprime_percatsub_avgWin.signelecs{i_TD}(:,i_sub)',...
            15,[0.2 0 0],'filled','MarkerFaceAlpha', 0.5)
    end

%     hold on; %Connect single subject dots
%     for i_sub = 1:length(subs)
%         line((h(1).XEndPoints(1):AnatReg_allSubs{2}.Num_Cat),...
%             plotstruct.avgKprime_percatsub_avgWin.allelecs{i_TD}(:,i_sub)','Color',[0 0 0], 'LineWidth', 0.2)
%         line((h(1).XEndPoints(1):AnatReg_allSubs{2}.Num_Cat) + 0.3,...
%             plotstruct.avgKprime_percatsub_avgWin.signelecs{i_TD}(:,i_sub)','Color',[0.3 0.3 0.3], 'LineWidth', 0.2, 'LineStyle', '--')
%     end

    f = get(gca,'Children');
    legend([f(end-1),f(end)], 'Sign. elecs', 'All elecs')

    title({['Kprime (Group Avg (n = ' num2str(length(subs)) ...
        ') + SEM/SD) per anatomical parcellation'], ...
        [input_DataLabel ', TD: ' input_ToneDurLabel{i_TD} ' s, ' ...
        FDR_label ', averaged across TW']}, ...
        'FontSize',14, 'Interpreter', 'none')

    if save_poststepFigs == 1
        filename     = ['ErrBar_Allsubn' num2str(length(subs)) '_' input_DataLabel...
            '_AvgKperElecParcelallTW_TD' num2str(i_TD) '_' ...
            FDR_label '_' param.ElecSelect 'elec.png'];
        path_fig_curr = [path_fig_Kperelecparcel '/Kprime_perTD/'];
        if (~exist(path_fig_curr, 'dir')); mkdir(path_fig_curr); end        
        figfile      = [path_fig_curr filename];
        saveas(gcf, figfile, 'png');       
    end
end

%% 10) Plot Kprime ratios between TD per anatomical parcellation
%Plot barcharts of group-average and single-subject Kprime ratios across parcellations

%10.1 Read out ratios per parcellation for each time window
%Per subject
plotstruct = [];
clear sum_signelec_percat_allsub

for i_TD = 1:length(input_ToneDurLabel)
    
    sum_signelec_percatwin_allsub = ...
        zeros(AnatReg_allSubs{i_sub}.Num_Cat, num_win_common);
    
    for i_sub = 1:length(subs)
        for i_win = 1:num_win_common
            for i_cat = 1:AnatReg_allSubs{i_sub}.Num_Cat
                
                %All elecs per cat
                filt_elecpercat = max(AnatReg_allSubs{i_sub}.CatIndex,[],2) == i_cat; %Parcellation filter
                
                %Read out Kprime values per TD for each selected electrode (for later use)
                plotstruct.Kprime_perTDsubcatwin.allelecs{i_sub}{i_cat}{i_win}(:, i_TD) = ...
                    data_selElec_allSubs.kprime_selelecs_perTD{i_win}{i_TD} ...
                    (data_selElec_allSubs.selelecs_persub{i_sub}(filt_elecpercat));
                
                %Compute ratio of across-electrode Kprime averages (not average of ratio!)
                plotstruct.RatioKprime_avgelecs_persubcatwin.allelecs{i_sub}(i_cat, i_win) = ...
                    nanmean(data_selElec_allSubs.kprime_selelecs_perTD{i_win}{1} ...
                    (data_selElec_allSubs.selelecs_persub{i_sub}(filt_elecpercat))) ...
                    / ...
                    nanmean(data_selElec_allSubs.kprime_selelecs_perTD{i_win}{2} ...
                    (data_selElec_allSubs.selelecs_persub{i_sub}(filt_elecpercat)));
                
                %Compute Kprime ratio for each electrode seperately
                plotstruct.RatioKprime_perelec_subcatwin.allelecs{i_sub}{i_cat, i_win} = ...
                    data_selElec_allSubs.kprime_selelecs_perTD{i_win}{1} ...
                    (data_selElec_allSubs.selelecs_persub{i_sub}(filt_elecpercat)) ...
                    ./ ...
                    data_selElec_allSubs.kprime_selelecs_perTD{i_win}{2} ...
                    (data_selElec_allSubs.selelecs_persub{i_sub}(filt_elecpercat));
                
                %Sign. elecs per cat
                filt_elecpercat = max(AnatReg_allSubs{i_sub}.CatIndex,[],2) == i_cat; %Parcellation filter
                filt_signelec   = data_selElec_allSubs.filt_signElecs{3}{i_win}(data_selElec_allSubs.selelecs_persub{i_sub}); %Common, not TD-wise sign.
                filt_signcat    = and(filt_elecpercat, filt_signelec); %Combined parcellation + sign. filter
                
                sum_signelec_percatwin_allsub(i_cat, i_win) = ...
                    sum_signelec_percatwin_allsub(i_cat, i_win) + sum(filt_signcat);
                
                %Read out Kprime values per TD for each selected electrode (for later use)
                plotstruct.Kprime_perTDsubcatwin.signelecs{i_sub}{i_cat}{i_win}(:, i_TD) = ...
                    data_selElec_allSubs.kprime_selelecs_perTD{i_win}{i_TD} ...
                    (data_selElec_allSubs.selelecs_persub{i_sub}(filt_signcat));
                
                %Compute ratio of across-electrode Kprime averages (not average of ratio!)
                plotstruct.RatioKprime_avgelecs_persubcatwin.signelecs{i_sub}(i_cat, i_win) = ...
                    nanmean(data_selElec_allSubs.kprime_selelecs_perTD{i_win}{1} ...
                    (data_selElec_allSubs.selelecs_persub{i_sub}(filt_signcat))) / ...
                    nanmean(data_selElec_allSubs.kprime_selelecs_perTD{i_win}{2} ...
                    (data_selElec_allSubs.selelecs_persub{i_sub}(filt_signcat)));
                
                %Compute Kprime ratio for each electrode seperately
                plotstruct.RatioKprime_perelec_subcatwin.signelecs{i_sub}{i_cat, i_win} = ...
                    data_selElec_allSubs.kprime_selelecs_perTD{i_win}{1} ...
                    (data_selElec_allSubs.selelecs_persub{i_sub}(filt_signcat)) ./ ...
                    data_selElec_allSubs.kprime_selelecs_perTD{i_win}{2} ...
                    (data_selElec_allSubs.selelecs_persub{i_sub}(filt_signcat));
                
            end
        end
    end
end

%Compute Kprime Ratio for across-electrode averages per parcellation and TW
for i_win = 1:num_win_common
    for i_cat = 1:AnatReg_allSubs{i_sub}.Num_Cat
        
        %All elecs per cat
        %Read out Kprime values for both TD for all elecs, per cat and win, and aggregate across subs
        temp_RatioKprime_allsub = [];
        for i_sub = 1:length(subs)
            temp_RatioKprime_allsub = ...
                [temp_RatioKprime_allsub; ...
                plotstruct.Kprime_perTDsubcatwin.allelecs{i_sub}{i_cat}{i_win}];
        end
        %Compute ratio of across-electrode Kprime averages (not average of ratio!)
        plotstruct.RatioKprime_avgelecs_percatwin.allelecs(i_cat,i_win) = ...
            nanmean(temp_RatioKprime_allsub(:,1)) / nanmean(temp_RatioKprime_allsub(:,2));        
        
        %Sign. elecs per cat
        temp_RatioKprime_allsub = [];
        for i_sub = 1:length(subs)
            temp_RatioKprime_allsub = ...
                [temp_RatioKprime_allsub; ...
                plotstruct.Kprime_perTDsubcatwin.signelecs{i_sub}{i_cat}{i_win}];
        end
        %Compute ratio of across-electrode Kprime averages (not average of ratio!
        plotstruct.RatioKprime_avgelecs_percatwin.signelecs(i_cat,i_win) = ...
            nanmean(temp_RatioKprime_allsub(:,1)) / nanmean(temp_RatioKprime_allsub(:,2));
        
    end
end
clear temp_Kprime_allsub

%Aggregate single-electrode Kprime ratios for each single electrode per TW
for i_win = 1:num_win_common
    for i_cat = 1:AnatReg_allSubs{i_sub}.Num_Cat
        
        %All elecs per cat for specific win
        temp_RatioKprime_pereleccatwin = [];
        for i_sub = 1:length(subs)
            temp_RatioKprime_pereleccatwin = ...
                [temp_RatioKprime_pereleccatwin;
                plotstruct.RatioKprime_perelec_subcatwin.allelecs{i_sub}{i_cat, i_win}];
        end
        plotstruct.RatioKprime_perelec_catwin.allelecs{i_cat, i_win} = ...
            temp_RatioKprime_pereleccatwin;
        %Compute SEM across all electrodes per parcel
        plotstruct.SEM_RatioKprime_perelec_wincat.allelecs{i_win}(1,i_cat) = ...
            nanstd(temp_RatioKprime_pereleccatwin) / sqrt(length(temp_RatioKprime_pereleccatwin));
        
        %Sign elecs per cat
        temp_RatioKprime_pereleccatwin = [];
        for i_sub = 1:length(subs)
            temp_RatioKprime_pereleccatwin = ...
                [temp_RatioKprime_pereleccatwin;
                plotstruct.RatioKprime_perelec_subcatwin.signelecs{i_sub}{i_cat, i_win}];
        end
        plotstruct.RatioKprime_perelec_catwin.signelecs{i_cat, i_win} = ...
            temp_RatioKprime_pereleccatwin;
        %Compute SEM across electrodes
        plotstruct.SEM_RatioKprime_perelec_wincat.signelecs{i_win}(1,i_cat) = ...
            nanstd(temp_RatioKprime_pereleccatwin) / sqrt(length(temp_RatioKprime_pereleccatwin));
        
    end
end

%10.2 Compute ratios per parcellation across all subjects, averaged across all common TW
%Per Subject
sum_signelec_percat_avgwin_allsub = ...
    zeros(AnatReg_allSubs{i_sub}.Num_Cat, 1);

for i_sub = 1:length(subs)
    for i_cat = 1:AnatReg_allSubs{i_sub}.Num_Cat
        
        %All elecs per cat
        filt_elecpercat = max(AnatReg_allSubs{i_sub}.CatIndex,[],2) == i_cat; %Parcellation filter
        %Read out Kprime values per TD for each selected electrode (for later use)
        for i_TD = 1:length(input_ToneDurLabel)
            plotstruct.Kprime_persubcatTD_avgWin.allelecs{i_sub}{i_cat}(:, i_TD) = ...
                data_selElec_allSubs.kprime_selelecs_perTD_avgWin{i_TD} ...
                (data_selElec_allSubs.selelecs_persub{i_sub}(filt_elecpercat));
        end        
        %Compute ratio of across-electrode Kprime averages (not average of ratio!)
        plotstruct.RatioKprime_avgelecs_percatsub_avgWin.allelecs(i_cat, i_sub) = ...
            nanmean(data_selElec_allSubs.kprime_selelecs_perTD_avgWin{1} ...
            (data_selElec_allSubs.selelecs_persub{i_sub}(filt_elecpercat))) ...
            / ...
            nanmean(data_selElec_allSubs.kprime_selelecs_perTD_avgWin{2} ...
            (data_selElec_allSubs.selelecs_persub{i_sub}(filt_elecpercat)));
        %Compute ratio for each single electrodes
        plotstruct.RatioKprime_perelec_catsub_avgWin.allelecs{i_cat, i_sub} = ...
            data_selElec_allSubs.kprime_selelecs_perTD_avgWin{1} ...
            (data_selElec_allSubs.selelecs_persub{i_sub}(filt_elecpercat)) ...
            ./ ...
            data_selElec_allSubs.kprime_selelecs_perTD_avgWin{2} ...
            (data_selElec_allSubs.selelecs_persub{i_sub}(filt_elecpercat));      
        
        %Sign. elecs per cat
        filt_signelec   = data_selElec_allSubs.filt_signElecs_avgWin{3}(data_selElec_allSubs.selelecs_persub{i_sub});
        filt_signcat    = and(filt_elecpercat, filt_signelec); %Combined parcellation + sign. filter
        sum_signelec_percat_avgwin_allsub(i_cat, 1) = ...
            sum_signelec_percat_avgwin_allsub(i_cat, 1) + sum(filt_signcat);
        for i_TD = 1:length(input_ToneDurLabel)
            %Read out Kprime values per TD for each selected electrode (for later use)
            plotstruct.Kprime_persubcatTD_avgWin.signelecs{i_sub}{i_cat}(:, i_TD) = ...
                data_selElec_allSubs.kprime_selelecs_perTD_avgWin{i_TD} ...
                (data_selElec_allSubs.selelecs_persub{i_sub}(filt_signcat));
        end
        %Compute ratio of across-electrode Kprime averages (not average of ratio!)
        plotstruct.RatioKprime_avgelecs_percatsub_avgWin.signelecs(i_cat, i_sub) = ...
            nanmean(data_selElec_allSubs.kprime_selelecs_perTD_avgWin{1} ...
            (data_selElec_allSubs.selelecs_persub{i_sub}(filt_signcat))) ...
            / ...
            nanmean(data_selElec_allSubs.kprime_selelecs_perTD_avgWin{2} ...
            (data_selElec_allSubs.selelecs_persub{i_sub}(filt_signcat)));
        %Compute ratio for each single electrodes
        plotstruct.RatioKprime_perelec_catsub_avgWin.signelecs{i_cat, i_sub} = ...
            data_selElec_allSubs.kprime_selelecs_perTD_avgWin{1} ...
            (data_selElec_allSubs.selelecs_persub{i_sub}(filt_signcat)) ...
            ./ ...
            data_selElec_allSubs.kprime_selelecs_perTD_avgWin{2} ...
            (data_selElec_allSubs.selelecs_persub{i_sub}(filt_signcat));        
        
    end
end

%Aggregate single electrode KprimeRatios across subjects for each anatomical parcel
for i_cat = 1:AnatReg_allSubs{i_sub}.Num_Cat
    
    %All elecs per cat
    temp_RatioKprime_allsub = [];
    for i_sub = 1:length(subs)
        temp_RatioKprime_allsub = ...
            [temp_RatioKprime_allsub; ...
            plotstruct.RatioKprime_perelec_catsub_avgWin.allelecs{i_cat, i_sub}];
    end
    plotstruct.RatioKprime_perelec_cat_avgWin.allelecs{i_cat} = ...
        temp_RatioKprime_allsub;
    
    %Sign. elecs per cat
    temp_RatioKprime_allsub = [];
    for i_sub = 1:length(subs)
        temp_RatioKprime_allsub = ...
            [temp_RatioKprime_allsub; ...
            plotstruct.RatioKprime_perelec_catsub_avgWin.signelecs{i_cat, i_sub}];
    end
    plotstruct.RatioKprime_perelec_cat_avgWin.signelecs{i_cat} = ...
        temp_RatioKprime_allsub;    
    
end

%Average Kprime Ratio per anatomical parcel across all subjects (for plotting)
for i_cat = 1:AnatReg_allSubs{i_sub}.Num_Cat
    
    %All elecs per cat
    %Read out Kprime values for both TD for all elecs, per cat and win, and aggregate across subs
    temp_RatioKprime_allsub = [];
    for i_sub = 1:length(subs)
        temp_RatioKprime_allsub = ...
            [temp_RatioKprime_allsub; ...
            plotstruct.Kprime_persubcatTD_avgWin.allelecs{i_sub}{i_cat}];
    end
    %Compute ratio of across-electrode Kprime averages (not average of ratio!)
    plotstruct.RatioKprime_avgelecs_percat_avgSubWin.allelecs(i_cat,1) = ...
        nanmean(temp_RatioKprime_allsub(:,1)) / nanmean(temp_RatioKprime_allsub(:,2));
    
    %Sign. elecs per cat
    temp_RatioKprime_allsub = [];
    for i_sub = 1:length(subs)
        temp_RatioKprime_allsub = ...
            [temp_RatioKprime_allsub; ...
            plotstruct.Kprime_persubcatTD_avgWin.signelecs{i_sub}{i_cat}];
    end
    %Compute ratio of across-electrode Kprime averages (not average of ratio!
    plotstruct.RatioKprime_avgelecs_percat_avgSubWin.signelecs(i_cat,1) = ...
        nanmean(temp_RatioKprime_allsub(:,1)) / nanmean(temp_RatioKprime_allsub(:,2));
    
end
clear temp_Kprime_allsub


%Statistically test if median of Kprime ratio for all elecs within a parcel is sign. different from 1
stats_signrank = struct;
for i_cat = 1:AnatReg_allSubs{i_sub}.Num_Cat
    %WSR-test
    if any(~isnan(plotstruct.RatioKprime_perelec_cat_avgWin.allelecs{i_cat})) ...
            && length(plotstruct.RatioKprime_perelec_cat_avgWin.allelecs{i_cat}) > 1
        [stats_signrank.allelecs.pval_uncorrected(i_cat,1), ~, stats_signrank.allelecs.stats{i_cat}] = ...
            signrank(plotstruct.RatioKprime_perelec_cat_avgWin.allelecs{i_cat}, 1,'tail', 'both');
    %     [~, stats_signrank.allelecs.pval_uncorrected(i_cat,1)] = ...
    %         ttest(plotstruct.RatioKprime_perelec_cat_avgWin.allelecs{i_cat}, 1,'tail', 'both');
        
    end
    if any(~isnan(plotstruct.RatioKprime_perelec_cat_avgWin.signelecs{i_cat})) ...
            && length(plotstruct.RatioKprime_perelec_cat_avgWin.allelecs{i_cat}) > 1            
        [stats_signrank.signelecs.pval_uncorrected(i_cat,1), ~, stats_signrank.signelecs.stats{i_cat}] = ...
            signrank(plotstruct.RatioKprime_perelec_cat_avgWin.signelecs{i_cat}, 1,'tail', 'both');
    %     [~, stats_signrank.signelecs.pval_uncorrected(i_cat,1)] = ...
    %         ttest(plotstruct.RatioKprime_perelec_cat_avgWin.signelecs{i_cat}, 1,'tail', 'both');
        
    end
end

%FDR-correct across parcels
stats_signrank.allelecs.pval_FDR = ...
    mafdr(stats_signrank.allelecs.pval_uncorrected,'BHFDR', true);
stats_signrank.signelecs.pval_FDR = ...
    mafdr(stats_signrank.signelecs.pval_uncorrected,'BHFDR', true);


%10.3 Plot barchart showing Kprime ratio per parcellation avg across TW
x_labels = {};
for i_cat = 1:length(AnatReg_allSubs{i_sub}.CatLabels)
    x_labels{end +1,1} = [AnatReg_allSubs{i_sub}.CatLabels{i_cat} ...
        ' (' num2str(sum(max(data_selElec_allSubs.index_AnatCat,[],2) == i_cat)) '/' ...
        num2str(round(nanmean(sum_signelec_percat_avgwin_allsub(i_cat)))) ')'];
end

FigName = ['Group_n', num2str(length(subs)), ' Kprime ratio across TD per anat parcel'];
figure('NumberTitle','off','Name', FigName, ...
    'units','normalized','outerposition',[0 0 1 1]);

yline(1, 'k-', 'LineWidth', 2)

hold on;
h = barwitherr( ... %Plots Kprime values for both all and sign. elecs
    [nanstd(plotstruct.RatioKprime_avgelecs_percatsub_avgWin.allelecs,[],2) / sqrt(length(subs)), ...
    nanstd(plotstruct.RatioKprime_avgelecs_percatsub_avgWin.signelecs,[],2) / sqrt(length(subs))], ...
    [plotstruct.RatioKprime_avgelecs_percat_avgSubWin.allelecs, ...
    plotstruct.RatioKprime_avgelecs_percat_avgSubWin.signelecs]);

set(gca,'xticklabel',x_labels,'fontsize',8)
xtickangle(45)
ylim([0 2])
set(gca, 'YGrid', 'on');

hold on; %Plot single electrode data points
for i_cat = 1:AnatReg_allSubs{i_sub}.Num_Cat
    scatter((ones(length(plotstruct.RatioKprime_perelec_cat_avgWin.allelecs{i_cat}),1) * h(1).XEndPoints(i_cat)),...
        plotstruct.RatioKprime_perelec_cat_avgWin.allelecs{i_cat},...
        5,[0 0 0.2],'o', 'filled','MarkerFaceAlpha', 0.5)
    scatter((ones(length(plotstruct.RatioKprime_perelec_cat_avgWin.signelecs{i_cat}),1) * h(1).XEndPoints(i_cat)) + 0.3,...
        plotstruct.RatioKprime_perelec_cat_avgWin.signelecs{i_cat},...
        5,[0 0 0.2],'o', 'filled','MarkerFaceAlpha', 0.5)
end

hold on; %Plot p-values in case of significance
for i_cat = 1:length(AnatReg_allSubs{i_sub}.CatLabels)
    if stats_signrank.allelecs.pval_uncorrected(i_cat) < 0.05 || stats_signrank.allelecs.pval_FDR(i_cat) < 0.05
        temp_text = {'WSR-test (elec) vs. 1:', ...
            ['p = ' num2str(round(stats_signrank.allelecs.pval_uncorrected(i_cat),2)) ...
            ' (uncorr)'],['p = ' num2str(round(stats_signrank.allelecs.pval_FDR(i_cat),2)) ' (FDR)']};
        text(h(1).XEndPoints(i_cat), 1.9, ...
            temp_text, 'FontSize',10)
        temp_text = {['p = ' num2str(round(stats_signrank.signelecs.pval_uncorrected(i_cat),2)) ...
            ' (uncorr)'],['p = ' num2str(round(stats_signrank.signelecs.pval_FDR(i_cat),2)) ' (FDR)']};
        text(h(1).XEndPoints(i_cat) + 0.3, 1.7, ...
            temp_text, 'FontSize',10)
    end
end

f=get(gca,'Children');
legend([f(end-2),f(end-1)], 'Sign. elecs', 'All elecs')

title({['Kprime-Ratio across TD (Group Avg (n = ' num2str(length(subs)) ...
    ') + SEM across subjects) per anatomical parcellation'], ...
    [input_DataLabel ', ' ...
    FDR_label ', averaged across TW']}, ...
    'FontSize',14, 'Interpreter', 'none')

if save_poststepFigs == 1
    filename     = ['ErrBar_Allsubn' num2str(length(subs)) '_' input_DataLabel...
        '_KRatioperElecParcel_allTW_' ...
        FDR_label '_' param.ElecSelect 'elec.png'];
    path_fig_curr = [path_fig_Kperelecparcel '/KprimeRatio/'];
    if (~exist(path_fig_curr, 'dir')); mkdir(path_fig_curr); end
    figfile      = [path_fig_curr filename];
    saveas(gcf, figfile, 'png');    
end

close all


%10.4 Plot barchart showing Kprime ratio per parcellation for each TW, along
%with stat results for across-electrode-per-parcel testing

%Statistically compare KprimeRatio against 1 across all elecs of a parcel
stats_signrank = struct;
for i_win = 1:num_win_common
    for i_cat = 1:AnatReg_allSubs{i_sub}.Num_Cat
        %WSR-test
        if any(~isempty(plotstruct.RatioKprime_perelec_catwin.allelecs{i_cat, i_win}))
            [stats_signrank.allelecs.pval_uncorrected(i_cat, i_win), ~, stats_signrank.allelecs.stats{i_cat, i_win}] = ...
                signrank(plotstruct.RatioKprime_perelec_catwin.allelecs{i_cat, i_win}, 1,'tail', 'both');
        %             [~, stats_signrank.allelecs.pval_uncorrected(i_cat,i_win)] = ...
        %                 ttest(plotstruct.RatioKprime_pereleccat_allsub.allelecs{i_win}{i_cat}, 1,'tail', 'both');
            
        end
        if any(~isempty(plotstruct.RatioKprime_perelec_catwin.signelecs{i_cat, i_win}))
            [stats_signrank.signelecs.pval_uncorrected(i_cat, i_win), ~, stats_signrank.signelecs.stats{i_cat, i_win}] = ...
                signrank(plotstruct.RatioKprime_perelec_catwin.signelecs{i_cat, i_win}, 1,'tail', 'both');
        %             [~, stats_signrank.signelecs.pval_uncorrected(i_cat,i_win)] = ...
        %                 ttest(plotstruct.RatioKprime_pereleccat_allsub.allelecs{i_win}{i_cat}, 1,'tail', 'both');
            
        end
    end
    %FDR-correct across parcels
    stats_signrank.allelecs.pval_FDR(:,i_win) = ...
        mafdr(stats_signrank.allelecs.pval_uncorrected(:,i_win),'BHFDR', true);
    stats_signrank.signelecs.pval_FDR(:,i_win) = ...
        mafdr(stats_signrank.signelecs.pval_uncorrected(:,i_win),'BHFDR', true);
end

%Plot barchart showing Kprime per parcellation and win, per TD

FigName = ['Group_n', num2str(length(subs)), ' Avg KprimeRatio per anat parcel'];
figure('NumberTitle','off','Name', FigName, ...
    'units','normalized','outerposition',[0 0 1 1]);
fig_counter = 0;
temp_num_win = num_win(1);

for i_win = 1:temp_num_win
    
    x_labels = {};
    for i_cat = 1:length(AnatReg_allSubs{i_sub}.CatLabels)
        x_labels{end +1,1} = [AnatReg_allSubs{i_sub}.CatLabels{i_cat} ...
            ' (' num2str(sum(max(data_selElec_allSubs.index_AnatCat,[],2) == i_cat)) '/' ...
            num2str(sum_signelec_percatwin_allsub(i_cat,i_win)) ')'];
    end
    
    fig_counter = fig_counter +1;
    
    subplot(round(sqrt(temp_num_win)), round(sqrt(temp_num_win)), fig_counter);
    h = barwitherr( ... %Plots Kprime values for both all and sign. elecs
        [plotstruct.SEM_RatioKprime_perelec_wincat.allelecs{i_win}', ...
        plotstruct.SEM_RatioKprime_perelec_wincat.signelecs{i_win}'], ...
        [plotstruct.RatioKprime_avgelecs_percatwin.allelecs(:,i_win), ...
        plotstruct.RatioKprime_avgelecs_percatwin.signelecs(:,i_win)]);
    title({'Kprime Ratio (Avg + SEM across electrodes) per anatomical parcellation'},'FontSize',10)
    
    yline(1, 'k-', 'LineWidth', 2)    
    
    set(gca,'xticklabel',x_labels,'fontsize',8)
    xtickangle(45)
    ylim([0 2.5])
    set(gca, 'YGrid', 'on');
    
    hold on; %Plot single electrode
    for i_cat = 1:AnatReg_allSubs{i_sub}.Num_Cat
        scatter((ones(length(plotstruct.RatioKprime_perelec_catwin.allelecs{i_cat, i_win}),1) * h(1).XEndPoints(i_cat)),...
            plotstruct.RatioKprime_perelec_catwin.allelecs{i_cat, i_win},...
            5,[0 0 0.2],'o', 'filled','MarkerFaceAlpha', 0.5)
        scatter((ones(length(plotstruct.RatioKprime_perelec_catwin.signelecs{i_cat, i_win}),1) * h(1).XEndPoints(i_cat) + 0.3),...
            plotstruct.RatioKprime_perelec_catwin.signelecs{i_cat, i_win},...
            5,[0 0 0.2],'o', 'filled','MarkerFaceAlpha', 0.5)
    end
    
    hold on; %Plot p-values in case of significance
    for i_cat = 1:length(AnatReg_allSubs{i_sub}.CatLabels)
        text_height = (-1)^i_cat * 0.3;
        if stats_signrank.allelecs.pval_uncorrected(i_cat) < 0.05 || stats_signrank.allelecs.pval_FDR(i_cat) < 0.05
            temp_text = {'WSR vs. 1:', ...
                ['p = ' num2str(round(stats_signrank.allelecs.pval_uncorrected(i_cat, i_win),2)) ...
                ' (uncorr)'],['p = ' num2str(round(stats_signrank.allelecs.pval_FDR(i_cat, i_win),2)) ' (FDR)']};
            text(h(1).XEndPoints(i_cat), 1.9 + text_height, ...
                temp_text, 'FontSize',8)
            temp_text = {['p = ' num2str(round(stats_signrank.signelecs.pval_uncorrected(i_cat, i_win),2)) ...
                ' (uncorr)'],['p = ' num2str(round(stats_signrank.signelecs.pval_FDR(i_cat, i_win),2)) ' (FDR)']};
            text(h(1).XEndPoints(i_cat) + 0.3, 1.6 + text_height, ...
                temp_text, 'FontSize',8)
        end
    end
    
    title(['TW: ' num2str(i_win)])
    
    if fig_counter == 1
        f=get(gca,'Children');
        legend([f(end-1),f(end)], 'Sign. elecs', 'All elecs')
    end
    
end

sgtitle({['Kprime Ratio (Avg + SEM across all elecs) per anatomical parcellation'], ...
    [input_DataLabel ', ' ...
    FDR_label ', per TW']}, ...
    'FontSize',14, 'Interpreter', 'none')

if save_poststepFigs == 1
    filename     = ['ErrBar_Allsubn' num2str(length(subs)) '_' input_DataLabel...
        '_AvgKRatioperElecParcelperTW_' ...
        FDR_label '_' param.ElecSelect 'elec.png'];
    path_fig_curr = [path_fig_Kperelecparcel '/KprimeRatio/'];
    if (~exist(path_fig_curr, 'dir')); mkdir(path_fig_curr); end
    figfile      = [path_fig_curr filename];
    saveas(gcf, figfile, 'png'); 
    close
end


end