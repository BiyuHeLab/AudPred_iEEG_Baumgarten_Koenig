function NASTD_ECoG_HisTrack_PlotSignExpKvalsBothTD_AllSubAvgTW...
    (subs, input_ToneDurLabel, input_DataLabel, ...
    SamplesTW, Label_TW, ...
    FDR_correct, pval_threshold,...
    param, ...
    save_poststepFigs, ...
    paths_NASTD_ECoG)

%Aim: Contrast and plot exp Kprime ratios across TD for all subjects.
%Averaged across common TW. Optionally restricted to sign. electrodes
%(sign. electrodes determine by across-TW-averaged exp vs. shuff Kprime values)

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

path_fig = [paths_NASTD_ECoG.ECoGdata_HisTrack ...
    'Allsub_n' num2str(length(subs)) '/ExpvsShuff/' Label_TW '/Figs/KprimeTDRatio_Surf/'];
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
        
        %Average exp and shuff Kprime values across TW
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
    
    %Combine sequence electrode electrode-selection across TD,
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
    %Determine significance (FDR or non-FDR thresholding) for selected elecs
    for i_TD = 1:length(input_ToneDurLabel)
        %Determine significance (FDR or non-FDR thresholding) for selected elecs
        if FDR_correct == 1
            FDR_label = 'FDRcorrected';
            pval_avgWin.(FDR_label){i_TD} = ...
                mafdr(data_selElec_allSubs.pval_allelecs_perSubTD_avgWin{i_sub}{i_TD}...
                (filt_signelecs_StimCorr_perSub), ...
                'BHFDR', true);
            filt_SuprathresElecs_avgWin{i_TD} = ...
                pval_avgWin.(FDR_label){i_TD} < pval_threshold;
        else
            FDR_label = 'uncorrected';
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
    
    %Select all elecs that show sign. SHI in at least one TW (liberal criterion)
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
    
end

for i_filter = 1:length(data_selElec_allSubs.filt_signElecs_avgWin)
    %Change sign. filter to logical
    data_selElec_allSubs.filt_signElecs_avgWin{i_filter} = ...
        logical(data_selElec_allSubs.filt_signElecs_avgWin{i_filter});
    %Compute overall number of sign. electrodes
    sumSignELec(i_filter) = sum(data_selElec_allSubs.filt_signElecs_avgWin{i_filter});
end


%% 4) Compute Kprime ratio across TD
%For each electrode
data_selElec_allSubs.RatioTD_ExpKprime_selElecs_avgWin = ...
    data_selElec_allSubs.kprime_selelecs_perTD_avgWin{1} ./ ...
    data_selElec_allSubs.kprime_selelecs_perTD_avgWin{2};

%Optional: Remove entries with kprime = 0 or inf, since these mess up
%min/max estimates and averages
ind_infentries = [];
ind_infentries = ...
    find(data_selElec_allSubs.RatioTD_ExpKprime_selElecs_avgWin == inf | ...
    data_selElec_allSubs.RatioTD_ExpKprime_selElecs_avgWin == 0);
data_selElec_allSubs.RatioTD_ExpKprime_selElecs_avgWin(ind_infentries) = NaN;

%Compute average ratio across sign. electrodes (using ratio of averages)
mean_KprimeRatio = ...
    nanmean(data_selElec_allSubs.kprime_selelecs_perTD_avgWin{1}...
    (data_selElec_allSubs.filt_signElecs_avgWin{3})) ...
    / ...
    nanmean(data_selElec_allSubs.kprime_selelecs_perTD_avgWin{2}...
    (data_selElec_allSubs.filt_signElecs_avgWin{3}));

%Compute std ratio across sign. electrodes (using std across electrode-specific ratios)
std_KprimeRatio = ...
    std(data_selElec_allSubs.kprime_selelecs_perTD_avgWin{1}...
    (data_selElec_allSubs.filt_signElecs_avgWin{3}) ...
    ./ ...
    data_selElec_allSubs.kprime_selelecs_perTD_avgWin{2}...
    (data_selElec_allSubs.filt_signElecs_avgWin{3}));


%Compute number of sign. elecs above and below 1
num_elecs_below1 = sum(data_selElec_allSubs.RatioTD_ExpKprime_selElecs_avgWin...
    (data_selElec_allSubs.filt_signElecs_avgWin{3}) < 1);
num_elecs_above1 = sum(data_selElec_allSubs.RatioTD_ExpKprime_selElecs_avgWin...
    (data_selElec_allSubs.filt_signElecs_avgWin{3}) > 1);


%% 5) Plot figure (surface + scatter plot)
% %Optional: Restrict plotting only to positive or negative ratios
% for i_win = 1:num_win_common
%     ind_selentries = [];
%     ind_selentries = ...
%         find(data_selElec_allSubs.RatioTD_ExpKprime_selElecs_avgWin < 1);
%     data_selElec_allSubs.RatioTD_ExpKprime_selElecs_avgWin(ind_selentries) = NaN;
% end

%Find zlimits (for difference in K-prime between TD) for constant plotting across time windows
dv_absmax = 0;
dv_absmin = 0;
dv_win_absmax = max(data_selElec_allSubs.RatioTD_ExpKprime_selElecs_avgWin...
    (data_selElec_allSubs.filt_signElecs_avgWin{3}));
if dv_win_absmax > dv_absmax
    dv_absmax = dv_win_absmax;
end
dv_win_absmin = min(data_selElec_allSubs.RatioTD_ExpKprime_selElecs_avgWin...
    (data_selElec_allSubs.filt_signElecs_avgWin{3}));
if dv_win_absmin < dv_absmin
    dv_absmin = dv_win_absmin;
end
dv_abs = max([abs(dv_absmax) abs(dv_absmin)]);

%Set up figure
fig = figure('visible','on');
set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen

%Set up data struct for plotting
clear dat
dat.dimord  = 'chan_time';
dat.time    = 0;
coords      = data_selElec_allSubs.elecpos; %MNI coordinates for selected electrodes
vals        = data_selElec_allSubs.RatioTD_ExpKprime_selElecs_avgWin; %parameter of interest for resp. electrodes
sign_index  = data_selElec_allSubs.filt_signElecs_avgWin{3};
chanSize    = (vals./vals)*2.5; %electrode size (arbitrary)
%clims       = [0 dv_abs]; %free scaling
clims       = [0 2]; %fixed scaling
cmap        = 'jet';

% %Both Hemispheres separately
% DimSubplot = [1, 3]; %Surface left, 2D plot right
% num_subplot = 2;
% NASTD_ECoG_Plot_SubplotSignElecsSurf_SubplotHemis_avgTW...
%     (coords, vals, sign_index, ...
%     chanSize, clims, cmap, DimSubplot, 1,...
%     []);

%Projected on left hemisphere
DimSubplot = [1, 1]; %Surface left, 2D plot right
num_subplot = 1;
%Project all electrodes on one hemisphere,
coords_projectedL = coords;
coords_projectedL(:,1) = abs(coords(:,1)) * -1;

NASTD_ECoG_Plot_SubplotSignElecsSurf_SubplotLH_avgTW...
    (coords_projectedL, vals, sign_index, ...
    chanSize, clims, cmap, DimSubplot, 1,...
    []);
fig.Children(2).Ticks = [0 0.5 1 1.5 2];
fig.Children(2).FontSize = 16;

%Add scatterplot showing k' values per TD
subplot(DimSubplot(1),DimSubplot(2),num_subplot+1);
plot(0:15, 0:15, 'g-','LineWidth',2) %unity line/information
hold on;
plot(0:1:15, 0:0.5:7.5, 'r-','LineWidth',2) %duration line

hold on;
scatter(...
    data_selElec_allSubs.kprime_selelecs_perTD_avgWin{1}...
    (data_selElec_allSubs.filt_signElecs_avgWin{3}), ...
    data_selElec_allSubs.kprime_selelecs_perTD_avgWin{2}...
    (data_selElec_allSubs.filt_signElecs_avgWin{3}), ...
    100, ...
    data_selElec_allSubs.RatioTD_ExpKprime_selElecs_avgWin...
    (data_selElec_allSubs.filt_signElecs_avgWin{3}), ...
    '.')
%Avg based on computing ratio from across-electrode averaged Kprime values
scatter(...
    nanmean(data_selElec_allSubs.kprime_selelecs_perTD_avgWin{1}...
    (data_selElec_allSubs.filt_signElecs_avgWin{3})), ...
    nanmean(data_selElec_allSubs.kprime_selelecs_perTD_avgWin{2}...
    (data_selElec_allSubs.filt_signElecs_avgWin{3})), ...
    200, ...
    nanmean(data_selElec_allSubs.kprime_selelecs_perTD_avgWin{1}...
    (data_selElec_allSubs.filt_signElecs_avgWin{3})) ./ ...
    nanmean(data_selElec_allSubs.kprime_selelecs_perTD_avgWin{2}...
    (data_selElec_allSubs.filt_signElecs_avgWin{3})), ...
    'd', 'filled', 'MarkerEdgeColor',[0.2 0.2 0.2], 'LineWidth', 2)
xlim([0 15])
ylim([0 15])
caxis([0 2])
c = colorbar;
c.Label.String = 'Kprime Ratio across TD';
c.Ticks = [0 0.5 1 1.5 2];
c.FontSize = 16;
colormap('jet');
xlabel('kPrime 200ms TD','FontSize',16)
ylabel('kPrime 400ms TD','FontSize',16)
%pbaspect([1 1 1])
axis('square')
grid on

%Set up and add title
title1      = ['Kprime Ratio across TD (200ms TD / 400 ms TD) for superthresh elecs (p < ' ...
    num2str(pval_threshold) ', ' FDR_label ')'];
w1          = num2str( 1000 * windows{1}(1, 1) / fsample , 3 );
w2          = num2str( 1000 * windows{1}(num_win_common, 2) / fsample , 3 );
title2      = ['Averaged across all common windows (TW = ' w1 ' - ' w2 ' ms)'];
title3      = [num2str(sumSignELec(3))...
    ' / ' num2str(length(data_selElec_allSubs.elecpos)) ' supra-thresh. elecs; '];
title4      = [num2str(num_elecs_above1) ' elecs > 1; ' num2str(num_elecs_below1) ' elecs < 1'];
title5      = ['max/min: ' num2str(round(dv_absmax,2)) '/' num2str(round(dv_absmin,2)) ...
    '; mean across elecs: ' num2str(round(mean_KprimeRatio,2))];
title6      = [num2str(length(ind_infentries)) ' elecs removed due to ratio values of 0 or Inf'];

Figtitle = {['Allsub (n = ' num2str(length(subs)) ') - '  input_DataLabel] ...
    title1 title2 [] title3 title4 title5 title6};
sgtitle(Figtitle,'Interpreter','none')

%Save
if save_poststepFigs == 1
    filename     = ['3DSurf_KprimeRatio_Allsubn' num2str(length(subs)) ...
        '_' input_DataLabel '_p' num2str(pval_threshold) FDR_label ...
        '_' param.ElecSelect 'elec_AvgWin.png'];
    figfile      = [path_fig filename];
    saveas(gcf, figfile, 'png'); %save png version
    close all;
end

end