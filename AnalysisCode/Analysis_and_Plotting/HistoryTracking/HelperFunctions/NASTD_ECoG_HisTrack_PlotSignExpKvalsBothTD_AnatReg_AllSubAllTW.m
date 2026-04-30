function NASTD_ECoG_HisTrack_PlotSignExpKvalsBothTD_AnatReg_AllSubAllTW...
    (subs, input_ToneDurLabel, input_DataLabel, ...
    SamplesTW, Label_TW, NumRuns, ...
    FDR_correct, pval_threshold,...
    param, ...
    plot_SubplotperTW, save_poststepFigs, ...
    paths_NASTD_ECoG)

%Aim: Contrast and plot exp Kprime values across both TD, averaged across
%all common TW, for anatomical regions defined by T1 electrode labels. 
%For all subjects, optionally restricted to sign. electrodes

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

path_fig = [paths_NASTD_ECoG.ECoGdata_HisTrack ...
    'Allsub_n' num2str(length(subs)) '/ExpvsShuff/' Label_TW '/Figs/Kprime_AnatReg/'];
if (~exist(path_fig, 'dir')); mkdir(path_fig); end


%% 1) Define analysis parameters and windows based on shorter TD
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


%% 2) Load and aggregate data across TD and subjects
data_selElec_allSubs                = struct;
data_selElec_allSubs.label          = [];
data_selElec_allSubs.sub            = [];
data_selElec_allSubs.elecpos        = [];
data_selElec_allSubs.elecs_persub   = [];
data_selElec_allSubs.label_AnatCat  = [];
data_selElec_allSubs.index_AnatCat  = [];

data_selElec_allSubs.filt_signelecs_StimCorr_allsub = [];

for i_sub = 1:length(subs)
    
    sub = subs{i_sub};
    
    tic
    %Load preprocessed data for elec coords
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    tic
    disp([' -- Loading preprocessed data set for sub: ' sub])
    load([si.path_preprocdata_sub], 'DataClean_AllTrials');
    disp([' -- Done loading in ' num2str(toc) ' sec'])
    
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
        
        %Load sequence tracking data and select current respective sign.
        %electrodes per TD
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
    
    %Combine sequence electrode electrode-selection across TD,
    %Select electrodes with sign. sign. sequence tracking effect in at least one TD
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
    data_selElec_allSubs.elecs_persub{i_sub} = ...
        find(data_selElec_allSubs.sub == i_sub);
    
    clear ind_selelecs label_selelecs_currSub
    
    %Aggregate electrode selection filter over subjects
    data_selElec_allSubs.filt_signelecs_StimCorr_allsub = ...
        [data_selElec_allSubs.filt_signelecs_StimCorr_allsub; ...
        filt_signelecs_StimCorr_perSub];
    data_selElec_allSubs.filt_signelecs_StimCorr_allsub = ...
        logical(data_selElec_allSubs.filt_signelecs_StimCorr_allsub);
    
    %Aggregate selected electrodes over subjects for each TD
    for i_TD = 1:length(input_ToneDurLabel)
        for i_win = 1:num_win_common
            if i_sub == 1
                data_selElec_allSubs.kprime_selelecs_perTD{i_win}{i_TD} = [];
            end
            data_selElec_allSubs.kprime_selelecs_perTD{i_win}{i_TD} = ...
                [data_selElec_allSubs.kprime_selelecs_perTD{i_win}{i_TD}; ...
                data_selElec_allSubs.kprime_allelecs_perSubWinTD{i_sub}{i_win}{i_TD}(filt_signelecs_StimCorr_perSub)];
        end
    end
    
    
    %% 3) Determine sign. SHI electrodes per subject
    for i_win = 1:num_win_common
        %Determine significance (FDR or non-FDR thresholding) for selected elecs
        for i_TD = 1:length(input_ToneDurLabel)
            if FDR_correct == 1
                FDR_label = 'FDRcorrected';
                pval.(FDR_label){i_win}{i_TD} = ...
                    mafdr(data_selElec_allSubs.pval_allelecs_perSubWinTD{i_sub}{i_win}{i_TD}...
                    (filt_signelecs_StimCorr_perSub), ...
                    'BHFDR', true);
                filt_SuprathresElecs{i_win}{i_TD} = pval.(FDR_label){i_win}{i_TD} < pval_threshold;
            else
                FDR_label = 'uncorrected';
                pval.(FDR_label){i_win}{i_TD} = ...
                    data_selElec_allSubs.pval_allelecs_perSubWinTD{i_sub}{i_win}{i_TD}...
                    (filt_signelecs_StimCorr_perSub);
                filt_SuprathresElecs{i_win}{i_TD} = pval.(FDR_label){i_win}{i_TD} < pval_threshold;
            end
        end
        
        %Combine across TD and select all elecs that show sign. SHI in at least
        %one TD (not in both)
        temp_pval_bothTD = [pval.(FDR_label){i_win}{1}, pval.(FDR_label){i_win}{2}];
        pval.(FDR_label){i_win}{3} = min(temp_pval_bothTD,[],2);
        filt_SuprathresElecs{i_win}{3} = ...
            pval.(FDR_label){i_win}{3} < pval_threshold;
        %sumSignELec{i_win} = sum(filt_SuprathresElecs{i_win}{3});
        clear temp_pval_bothTD
        
        %Aggregate sign. filter and min p-values across subjects
        if i_sub == 1
            data_selElec_allSubs.filt_signElecs{i_win} = [];
            data_selElec_allSubs.pval_signElecs{i_win} = [];
        end
        data_selElec_allSubs.filt_signElecs{i_win} = ...
            [data_selElec_allSubs.filt_signElecs{i_win}; filt_SuprathresElecs{i_win}{3}];
        data_selElec_allSubs.pval_signElecs{i_win} = ...
            [data_selElec_allSubs.pval_signElecs{i_win}; pval.(FDR_label){i_win}{3}];
    end
    
    
    %% 4) Categorize electrodes according to anatomical regions
    AnatReg_allSubs{i_sub} = ...
        NASTD_ECoG_AssignAnatRegions(DataClean_AllTrials, labels_loadedData);
    clear DataClean_AllTrials labels_loadedData
    
    data_selElec_allSubs.label_AnatCat  = ...
        [data_selElec_allSubs.label_AnatCat; ...
        AnatReg_allSubs{i_sub}.Info_perelec(:,3)];
    data_selElec_allSubs.index_AnatCat  = ...
        [data_selElec_allSubs.index_AnatCat; ...
        AnatReg_allSubs{i_sub}.CatIndex];
    
    clear pval filt_SuprathresElecs filt_signelecs_StimCorr_perSub
    disp(['-- ' sub ' finished in ' num2str(round(toc)) ' sec --'])
    
    
    %% 5) Plot surface with electrode names and categories per subject
    clear dat
    dat.dimord  = 'chan_time';
    dat.time    = 0;
    coords      = data_selElec_allSubs.elecpos(data_selElec_allSubs.elecs_persub{i_sub},:); %MNI coordinates for selected electrodes
    vals        = data_selElec_allSubs.index_AnatCat(data_selElec_allSubs.elecs_persub{i_sub});
    sign_index  = logical(ones(length(data_selElec_allSubs.elecs_persub{i_sub}),1));
    label       = data_selElec_allSubs.label_AnatCat(data_selElec_allSubs.elecs_persub{i_sub},:);
    chanSize    = (vals./vals)*3; %electrode size (arbitrary)
    clims       = [1 12]; %fixed scaling
    cmap        = 'jet';
    view_angle  = [270,0];
    
    figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    NASTD_ECoG_Plot_SubplotSignElecsSurf_Label...
        (coords, label, vals, sign_index,...
        chanSize, clims, cmap, view_angle, ...
        [1 1], 1, [0 0 0 0], [0 0 0 0]);
    c = colorbar('Location', 'West') ;
    c.TickLabels = AnatReg_allSubs{i_sub}.CatLabels;
    c.FontSize = 12;

%% 6) Average Kprime across all common time windows and combine sign. filter across all time windows
for i_win = 1:num_win_common
    temp_kprime_selelecs_perTDwin{1}(:,i_win) = ...
        data_selElec_allSubs.kprime_selelecs_perTD{i_win}{1};
    temp_kprime_selelecs_perTDwin{2}(:,i_win) = ...
        data_selElec_allSubs.kprime_selelecs_perTD{i_win}{2};   
    temp_filt_signelecs_allwin(:,i_win) = ...
        data_selElec_allSubs.filt_signElecs{i_win};    
end
data_selElec_allSubs.kprime_selelecs_avgwinperTD(:,1) = ...
    nanmean(temp_kprime_selelecs_perTDwin{1},2);
data_selElec_allSubs.kprime_selelecs_avgwinperTD(:,2) = ...
    nanmean(temp_kprime_selelecs_perTDwin{2},2);
data_selElec_allSubs.filt_signElecs_allwin = ...
    logical(sum(temp_filt_signelecs_allwin,2));
    
   
%% 7) Compute Kprime ratio across TD for sign. electrodes
    data_selElec_allSubs.RatioTD_ExpKprime_selElecs{i_win} = [];
    for i_win = 1:num_win_common
        data_selElec_allSubs.RatioTD_ExpKprime_selElecs{i_win} = ...
            data_selElec_allSubs.kprime_selelecs_perTD{i_win}{1} ./ ...
            data_selElec_allSubs.kprime_selelecs_perTD{i_win}{2};
    end   



end

%Convert Filter to logical and count sign. elecs
for i_win = 1:num_win_common
    data_selElec_allSubs.filt_signElecs{i_win} = ...
        logical(data_selElec_allSubs.filt_signElecs{i_win});
    sumSignELec{i_win} = sum(data_selElec_allSubs.filt_signElecs{i_win});
end

%% 7) Plot Kprime values per anatomical category

figure
set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen

%7.1 kprime values per TD - single-subjects and all subjects per category - all elecs
for i_sub = 1:length(subs)
        for i_cat = 1:AnatReg_allSubs{i_sub}.Num_IndivCategories
            temp_filt = AnatReg_allSubs{i_sub}.CatIndex == i_cat;
            temp_avgKprime_percat{i_sub}{1}(i_cat, :) = ...
                nanmean(data_selElec_allSubs.kprime_selelecs_perTD{i_win}{1}(temp_filt));
            temp_avgKprime_percat{i_sub}{2}(i_cat, :) = ...
                nanmean(data_selElec_allSubs.kprime_selelecs_perTD{i_win}{2}(temp_filt));
%             temp_stdKprime_percat{i_sub}{i_win}(i_cat,:) = ...
%                [nanstd(data_selElec_allSubs.kprime_selelecs_perTD{i_win}{1}(temp_filt)), ...
%                 nanstd(data_selElec_allSubs.kprime_selelecs_perTD{i_win}{2}(temp_filt))];
            
        end
end

figure
set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen

for i_sub = 1:length(subs)
    subplot(2,5,i_sub)
    bar([nanmean(temp_avgKprime_percat{i_sub}{1},1)', nanmean(temp_avgKprime_percat{i_sub}{2},1)'])
    xticklabels(AnatReg_allSubs{i_sub}.CatLabels)
    xtickangle(45)
end

%Note: Problem - not all subjects have identical category names - problem
%with averaging - ensure identical names, even if ther eis no entry in said
%field
