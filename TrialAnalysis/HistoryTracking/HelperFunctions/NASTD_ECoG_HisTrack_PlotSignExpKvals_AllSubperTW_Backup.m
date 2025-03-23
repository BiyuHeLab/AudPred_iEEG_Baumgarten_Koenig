function NASTD_ECoG_HisTrack_PlotSignExpKvals_AllSubperTW...
    (subs, input_ToneDurLabel, input_DataLabel, ...
    SamplesTW, Label_TW,...
    param, ...
    FDR_correct, pval_plotting,...
    plot_SubplotperTW, save_poststepFigs, ...
    paths_NASTD_ECoG)

%Aim: Plot exp Kprime values for all electrodes and all subjects on MNI volume,per TW

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

path_fig = [paths_NASTD_ECoG.ECoGdata_HisTrack ...
    'Allsub_n' num2str(length(subs)) '/ExpvsShuff/' Label_TW '/Figs/Kprime_Surf/'];
if (~exist(path_fig, 'dir')); mkdir(path_fig); end

%Tone duration condition for data load-in
if strcmp(input_ToneDurLabel,'0.2')
    tonedur_title = '200msTD';
elseif  strcmp(input_ToneDurLabel,'0.4')
    tonedur_title = '400msTD';
end

%% 1) Load in data
%1.1 Load in ECoG preproc data for channel labels and position
data_selElec_allSubs                = struct;
data_selElec_allSubs.label          = [];
data_selElec_allSubs.sub            = [];
data_selElec_allSubs.elecpos        = [];
data_selElec_allSubs.elecs_persub   = [];
    
for i_sub = 1:length(subs)
    
    tic    
    sub = subs{i_sub};
    NASTD_ECoG_subjectinfo
    %Load preprocessed data    
    load([si.path_preprocdata_sub]);
    fsample = DataClean_AllTrials.fsample; %read this out while preproc data in workspace
    
    %Load Kprime data
    path_inputdata = [paths_NASTD_ECoG.ECoGdata_HisTrack sub '/ExpvsShuff/' Label_TW '/'];
    load([path_inputdata sub '_' input_DataLabel '_' tonedur_title ...
        '_HisTrack_ExpvsShuffKprimeCombFolds_unbalanced.mat']); %filename: Kprime_data, labels_loadedData
    
    %Optional: Load stimulus correlation data and select current effect
    path_inputdata = [paths_NASTD_ECoG.ECoGdata_StimCorr  sub '/Data/'];    
    if strcmp(param.ElecSelect, 'StimCorr')
        load([path_inputdata sub '_StimCorr_' input_DataLabel '_' ...
            input_ToneDurLabel 'sTD.mat'], ...
            'corr_ttest');
        filt_signelecs_StimCorr = corr_ttest.p < 0.05;%only sign. stim corr elecs
    else
        filt_signelecs_StimCorr = true(length(labels_loadedData),1);  %all elecs
    end
    ind_selelecs  = find(filt_signelecs_StimCorr);
    nSensors_sel  = sum(filt_signelecs_StimCorr);
    %labels_loadedData(ind_selelecs)
    
    %Determine selected electrode index to read out elec coords in MNI space
    ind_selelecs_MNI   = [];
    for i_selelec = 1:length(ind_selelecs)
        ind_selelecs_MNI = [ind_selelecs_MNI; ...
            find(strcmp(DataClean_AllTrials.elec.label, ...
            labels_loadedData{ind_selelecs(i_selelec)}))];
    end
    
     data_selElec_allSubs.elecpos = ...
        [data_selElec_allSubs.elecpos; ...
        DataClean_AllTrials.elec.chanpos(ind_selelecs_MNI,:)];    
   
    %Read out selected elecs per sub,and append selected electrodes 
    %to one common summary struct
    indiv_label = {};
    sub = subs{i_sub};
    
    for i_elec = 1:length(ind_selelecs)
        indiv_label(i_elec,1) = ...
            strcat(labels_loadedData(ind_selelecs(i_elec)),'_', sub); 
        %Add sub-identifier to each elec label
    end
    
    data_selElec_allSubs.label                  = ...
        [data_selElec_allSubs.label; indiv_label];
    data_selElec_allSubs.sub                    = ...
        [data_selElec_allSubs.sub; ones(length(indiv_label),1)*i_sub];
    data_selElec_allSubs.elecs_persub{i_sub}    = ...
        find(data_selElec_allSubs.sub == i_sub);    
    
    %Read out exp Kprime values for selected electrodes
    for i_win = 1:length(Kprime_data.Exp.avgRuns.Kprime_Avg)
        if i_sub == 1
            data_selElec_allSubs.kprime{i_win} = [];
            data_selElec_allSubs.index_signkprime{i_win} = [];
        end        
        data_selElec_allSubs.kprime{i_win} = ...
            [data_selElec_allSubs.kprime{i_win}; ...
            Kprime_data.Exp.avgRuns.Kprime_Avg{i_win}(ind_selelecs)];
        
    %Read out shuff Kprime values for selected electrodes to later re-compute pvals for across TW averaged data       
        if i_sub == 1
            data_selElec_allSubs.kprime_shuff{i_win} = [];
        end        
        data_selElec_allSubs.kprime_shuff{i_win} = ...
            [data_selElec_allSubs.kprime_shuff{i_win}; ...
            Kprime_data.Shuff.perRepavgRuns.Kprime_NullDistribution{i_win}(ind_selelecs,:)];        
    
    %Thresholding and FDR
         if FDR_correct == 1
            pval_FDRcorrected = mafdr(...
                Kprime_data.Exp.avgRuns.pval_ExpvsShuff{i_win}(ind_selelecs), ...
                'BHFDR', true);
            data_selElec_allSubs.index_signkprime{i_win} = ...
                [data_selElec_allSubs.index_signkprime{i_win}; ...
                (pval_FDRcorrected < pval_plotting)];
            FDR_label = 'FDRcorrected';
        else
            data_selElec_allSubs.index_signkprime{i_win} = ...
                [data_selElec_allSubs.index_signkprime{i_win}; ...
                (Kprime_data.Exp.avgRuns.pval_ExpvsShuff{i_win}(ind_selelecs)...
                < pval_plotting)];
            FDR_label = 'uncorrected';
         end 
         sumSignELec = sum(data_selElec_allSubs.index_signkprime{i_win});
    end
    
    clear DataClean_AllTrials Kprime_data labels_loadedData 'corr_ttest'
    disp(['-- ' sub ' finished in ' num2str(round(toc)) ' sec --']) 
end


%% 2) Define analysis parameters
%2.1 Define analysis parameters
toneDur_inSecs  = str2num(input_ToneDurLabel);
nSamplesPerTone = toneDur_inSecs * fsample;

win_size    = SamplesTW;  %25 data points with fs 512 = ~49ms as mentioned in manuscript
s_per_win   = (1/fsample)*win_size;
win_overlap = 0; %as mentioned in manuscript

windows = [1 win_size];
while windows(end,end) < nSamplesPerTone
    windows = [windows; windows(end,:) + (win_size - win_overlap)];
end

if windows(end,end) > nSamplesPerTone
    windows(end,:) = [];
end


%% 3) Plot exp Kprime values on surface
%3.1 find zlimits for constant plotting across time windows
dv_absmax = 0;
dv_absmin = 16;

for i_win = 1:size(windows,1)
    dv_win_absmax = ceil(max(data_selElec_allSubs.kprime{i_win}...
        (logical(data_selElec_allSubs.index_signkprime{i_win}))));
    dv_win_absmin = floor(min(data_selElec_allSubs.kprime{i_win}...
        (logical(data_selElec_allSubs.index_signkprime{i_win}))));
    
    if dv_win_absmax > dv_absmax
        dv_absmax = dv_win_absmax;
    end
    if dv_win_absmin < dv_absmin
        dv_absmin = dv_win_absmin;
    end
end

if plot_SubplotperTW == 1 %one subplot per TW
    figure('visible','on'); %ensures figure doesn't pop up during plotting
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    if size(windows,1) == 4
        DimSubplot = [size(windows,1)/2, size(windows,1)];
    else
        DimSubplot = [size(windows,1)/2, size(windows,1)/2];        
    end
    CounterSubplot = 1;
end

%3.2 Plot data for current time window
for i_win = 1:size(windows,1)        
    kprime_signelecs    = ...
        data_selElec_allSubs.kprime{i_win}...
        (logical(data_selElec_allSubs.index_signkprime{i_win}));
    index_signelecs     = find(logical(...
        data_selElec_allSubs.index_signkprime{i_win}));
    
    %Set up title information
    w1 = num2str( 1000 * windows(i_win, 1) / fsample , 3 );
    w2 = num2str( 1000 * windows(i_win, 2) / fsample , 3 );
    
    effect_title1   = 'Thresholded Exp Kprime values';
    effect_title2   = ['(Exp vs. Shuff; p < ' num2str(pval_plotting) '; ' FDR_label ')'];
    win_title       = ['TW = ' w1 ' - ' w2 ' ms'];
    
    if i_win == 1    
        stat_title = ['min/max(mean+-SD) Kprime = ' ...
            num2str(round(min(kprime_signelecs),2)) '/' num2str(round(max(kprime_signelecs),2)) ' (' ...
            num2str(round(mean(kprime_signelecs),2)) ' +- ' ...
            num2str(round(std(kprime_signelecs),2)) ')'];
    else
        stat_title = [num2str(min(kprime_signelecs)) '/' ...
            num2str(round(max(kprime_signelecs),2)) ' (' ...
            num2str(round(mean(kprime_signelecs),2)) ' +- ' ...
            num2str(round(std(kprime_signelecs),2)) ')'];
    end
    
%     label_signelecs = [];
%     if ~isempty(index_signelecs)
%         for i_signelec = index_signelecs'
%             label_signelecs = [label_signelecs, ...
%                 strcat(data_selElec_allSubs.label{i_signelec}, ',-')];
%         end
%     end
    sign_title = [num2str(length(index_signelecs))...
        ' / ' num2str(length(data_selElec_allSubs.kprime{i_win})) ...
        ' supra-thresh. elecs'];    
    
    %Set up data struct for plotting
    vals        = data_selElec_allSubs.kprime{i_win}; %parameter of interest for resp. electrodes
    sign_index  = data_selElec_allSubs.index_signkprime{i_win};    
    coords      = data_selElec_allSubs.elecpos; %MNI coordinates for selected electrodes
    chanSize    = (vals./vals)*3; %electrode size (arbitrary, 3 is usually good)
    %Optional: plot higher Kprime values bigger
    for i_elec = 1:length(vals)
        chanSize(i_elec,1) = ...
            (vals(i_elec)./vals(i_elec)) * ...
            (abs(log10(vals(i_elec)))*4);
    end    
    %clims = [0 max([dv_absmax, 1.1])]; %free scaling
%     clims       = [0 15]; %fixed scaling
    clims       = [0 10]; %fixed scaling
    cmap        = 'jet';
    
    if plot_SubplotperTW == 0 %one separate figure per TW
        
        Figtitle = {['Allsub (n = ' num2str(length(subs)) ') - ' ...
            input_DataLabel ' - ToneDur: ' tonedur_title ' - ' win_title]...
            [effect_title1] [effect_title2] [sign_title ' - ' stat_title]};      

            NASTD_ECoG_Plot_PlotSignElecsSurf_SubplotHemis...
                (coords, vals, sign_index , ...
                chanSize, clims, cmap) %Both H
            sgtitle(Figtitle, 'Interpreter','none')
            
        if save_poststepFigs == 1
            path_fig1 = [path_fig FDR_label '/'];
            if (~exist(path_fig1, 'dir')); mkdir(path_fig1); end
            filename     = ['3DSurf_SignExpKprime_Allsub__' input_DataLabel ...
                '_' tonedur_title '_' Label_TW num2str(i_win) ...
                '_p' num2str(pval_plotting) FDR_label ...
                '_' param.ElecSelect 'elec.png'];
            figfile      = [path_fig1 filename];
            saveas(gcf, [figfile], 'png'); %save png version
            close all;
        end
        
    elseif plot_SubplotperTW == 1 %one subplot per TW
        
        Figtitle = {['Allsub (n = ' num2str(length(subs)) ') - '  effect_title1 ' - ' ...
            input_DataLabel ' - ToneDur: ' tonedur_title ' - ' effect_title2], [' '], [' ']};
        
            NASTD_ECoG_Plot_SubplotSignElecsSurf_SubplotHemis...
                (coords, vals, sign_index, ...
                chanSize, clims, cmap, ...
                DimSubplot, CounterSubplot,...
                {[win_title], [sign_title], [stat_title]});        
            CounterSubplot = CounterSubplot + 2;            
    end
end

%Increase and label colorbar
C = findall(gcf,'type','ColorBar');
C.Label.String = 'Kprime value';
C.Label.FontSize = 12;
if size(windows,1) > 4
    C.Position(4) = C.Position(4)*3;
    C.Position(2) = C.Position(2)*0.9; %lower Colorbar
end

if plot_SubplotperTW == 1 %one subplot per TW
    sgtitle(Figtitle, 'Interpreter','none')
    
    if save_poststepFigs == 1
        path_fig1 = [path_fig FDR_label '/'];
        if (~exist(path_fig1, 'dir')); mkdir(path_fig1); end
        filename     = ['3DSurf_SignExpKprime_Allsub_' input_DataLabel '_' ...
            tonedur_title '_' Label_TW 'all_p' num2str(pval_plotting) FDR_label ...
            '_' param.ElecSelect 'elec.png'];
        figfile      = [path_fig1 filename];
        saveas(gcf, [figfile], 'png'); %save png version
        close all;
    end
end


%% 4) Plot exp Kprime values averaged across TW on surface
%Average exp and shuff Kprime values across all TW
for i_win = 1:size(windows,1)
    ExpKprime_data_perTW(:,i_win) = ...
        data_selElec_allSubs.kprime{i_win};
    ShuffKprime_data_perTW(i_win, :, :) = ...
        data_selElec_allSubs.kprime_shuff{i_win};
end
ExpKprime_data_avgTW    = mean(ExpKprime_data_perTW,2);
ShuffKprime_data_avgTW  = squeeze(mean(ShuffKprime_data_perTW,1));

%Recompute p-values based on across-TW averages
for i_sensor = 1:length(ExpKprime_data_avgTW)
    data_selElec_allSubs.pval_avgTW(i_sensor,1) = ...
        sum(ShuffKprime_data_avgTW(i_sensor, :) >= ...
        ExpKprime_data_avgTW(i_sensor)) ...
        / length(ShuffKprime_data_avgTW(i_sensor, :));
end

%Determine significance (FDR or non-FDR thresholding)
if FDR_correct == 1
    pval_FDRcorrected       = mafdr(...
        data_selElec_allSubs.pval_avgTW, 'BHFDR', true);
    data_selElec_allSubs.index_signkprime_avgTW   = pval_FDRcorrected < pval_plotting;
    sumSignELec             = sum(data_selElec_allSubs.index_signkprime_avgTW);
    FDR_label               = 'FDRcorrected';
else
    data_selElec_allSubs.index_signkprime_avgTW   = ...
        data_selElec_allSubs.pval_avgTW < pval_plotting;
    sumSignELec             = sum(data_selElec_allSubs.index_signkprime_avgTW);
    FDR_label               = 'uncorrected';
end
index_signelecs = find(logical(...
    data_selElec_allSubs.index_signkprime_avgTW));

%Plot avg Kprime
effect_title    = ['Thresholded Exp Kprime values averaged across TW (Exp vs. Shuff; p < ' ...
    num2str(pval_plotting) '; ' FDR_label ')'];

stat_title = ['min/max(mean+-SD) Kprime = ' ...
    num2str(round(min(ExpKprime_data_avgTW(index_signelecs)),2)) '/' ...
    num2str(round(max(ExpKprime_data_avgTW(index_signelecs)),2)) ' (' ...
    num2str(round(mean(ExpKprime_data_avgTW(index_signelecs)),2)) ' +- ' ...
    num2str(round(std(ExpKprime_data_avgTW(index_signelecs)),2)) ')'];

% label_signelecs = [];
% if ~isempty(index_signelecs)
%     for i_signelec = index_signelecs'
%         label_signelecs = [label_signelecs, ...
%             strcat(data_selElec_allSubs.label{i_signelec}, ',-')];
%     end
% end
sign_title = [num2str(length(index_signelecs))...
    ' / ' num2str(length(data_selElec_allSubs.kprime{1})) ...
    ' supra-thresh. elecs'];

Figtitle = {['Allsub (n = ' num2str(length(subs)) ') - ' ...
    input_DataLabel ' - ToneDur: ' tonedur_title]...
    [effect_title]  [sign_title ' - ' stat_title]};

%Set up data struct for plotting
%Plot electrodes as spheres on MNI brain with each sphere color-coded depending on Kprime-value
coords      = data_selElec_allSubs.elecpos; %MNI coordinates for selected electrodes
vals        = ExpKprime_data_avgTW;
sign_index  = data_selElec_allSubs.index_signkprime_avgTW;
%     clims       = [0 max([dv_absmax, 1.1])]; %free scaling
% clims       = [0 15]; %fixed scaling
clims       = [0 10]; %fixed scaling
cmap        = 'jet';
chanSize    = (vals./vals)*3; %electrode size (arbitrary)
%Optional: plot higher Kprime values bigger
for i_elec = 1:length(vals)
    chanSize(i_elec,1) = ...
        (vals(i_elec)./vals(i_elec)) * ...
        (abs(log10(vals(i_elec)))*4);
end

NASTD_ECoG_Plot_PlotSignElecsSurf_SubplotHemis...
    (coords, vals, sign_index , ...
    chanSize, clims, cmap) %Both H
sgtitle(Figtitle, 'Interpreter','none')

%Increase and label colorbar
C = findall(gcf,'type','ColorBar');
C.Label.String = 'Kprime value';
C.Label.FontSize = 12;
C.Position(4) = C.Position(4)*0.8;
C.Position(2) = C.Position(2)*0.7; %lower Colorbar

if save_poststepFigs == 1
    filename     = ['3DSurf_SignExpKprime_Allsub__' input_DataLabel ...
        '_' tonedur_title '_' Label_TW 'Avgall_p' ...
        num2str(pval_plotting) FDR_label ...
        '_' param.ElecSelect 'elec.png'];
    figfile      = [path_fig1 filename];
    saveas(gcf, [figfile], 'png'); %save png version
    close all;
end

end