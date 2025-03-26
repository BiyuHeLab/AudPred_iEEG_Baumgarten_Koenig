function NASTD_ECoG_HisTrack_PlotSignExpKvals_StimCorrElec_AllSubperTW...
    (subs, FuncInput_ToneDur_text, FuncInput_DataType, ...
    SamplesTW, Label_TW,...
    FDR_correct, pval_plotting,...
    plot_SubplotperTW, save_poststepFigs, ...
    paths_NASTD_ECoG)

%Aim: Plot exp Kprime values for all electrodes and all subjects on MNI volume,per TW

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

path_fig_KprimeSurf = [paths_NASTD_ECoG.ECoGdata_HisTrack ...
    'Allsub_n' num2str(length(subs)) '/Exp/' Label_TW '/Figs/Kprime3DSurf_signElecs_StimCorrROI/'];

%Tone duration condition for data load-in
if strcmp(FuncInput_ToneDur_text,'0.2')
    tonedur_title = '200msTD';
elseif  strcmp(FuncInput_ToneDur_text,'0.4')
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
    load([path_inputdata sub '_' FuncInput_DataType '_' tonedur_title ...
        '_HisTrack_ExpvsShuffKprimeCombFolds_unbalanced.mat']); %filename: Kprime_data

    %Load stimulus correlation data and select current effect
    path_inputdata = [paths_NASTD_ECoG.ECoGdata_StimCorr '/' sub '/Data/'];
    load([path_inputdata sub '_StimCorr_' FuncInput_DataType '_' ...
        FuncInput_ToneDur_text 'sTD.mat'], ...
        'corr_ttest', 'SensorLabels');
    
    filt_signelecs_StimCorr = corr_ttest.p < 0.05;
    ind_signelecs_StimCorr  = find(filt_signelecs_StimCorr);
    nSensors_sel            = sum(filt_signelecs_StimCorr);
        
    %Read out selected elecs per sub,and append selected electrodes 
    %to one common summary struct
    selElec_currSub = [];
    indiv_label = {};
    sub = subs{i_sub};
    
    for i_elec = 1:nSensors_sel
        indiv_label(i_elec,1) = ...
            strcat(SensorLabels(ind_signelecs_StimCorr(i_elec)),'_', sub); 
        %Add sub-identifier to each elec label
    end
    
    data_selElec_allSubs.label                  = ...
        [data_selElec_allSubs.label; indiv_label];
    data_selElec_allSubs.sub                    = ...
        [data_selElec_allSubs.sub; ...
        ones(length(indiv_label),1)*i_sub];
    data_selElec_allSubs.elecs_persub{i_sub}    = ...
        find(data_selElec_allSubs.sub == i_sub);
    
    SelElec_MNIindex   = [];
    for i_elec = 1:nSensors_sel
        SelElec_MNIindex = [SelElec_MNIindex, ...
            find(strcmp(DataClean_AllTrials.elec.label, ...
            SensorLabels{ind_signelecs_StimCorr(i_elec)}))];
    end
    
    data_selElec_allSubs.elecpos = ...
        [data_selElec_allSubs.elecpos; ...
        DataClean_AllTrials.elec.chanpos(SelElec_MNIindex,:)];    
    
    for i_win = 1:length(Kprime_data.Exp.avgFolds.Kprime_Avg)
        if i_sub == 1
            data_selElec_allSubs.kprime{i_win} = [];
            data_selElec_allSubs.index_signkprime{i_win} = [];
        end        
        data_selElec_allSubs.kprime{i_win} = ...
            [data_selElec_allSubs.kprime{i_win}; ...
            Kprime_data.Exp.avgFolds.Kprime_Avg{i_win}(ind_signelecs_StimCorr)'];
    
    %Thresholding and FDR
         if FDR_correct == 1
            pval_FDRcorrected = mafdr(...
                Kprime_data.Exp.avgFolds.pval_ExpvsShuff{i_win}(ind_signelecs_StimCorr), ...
                'BHFDR', true);
          
            data_selElec_allSubs.index_signkprime{i_win} = ...
                [data_selElec_allSubs.index_signkprime{i_win}; ...
                (pval_FDRcorrected < pval_plotting)'];
            sumSignELec = sum(pval_FDRcorrected < pval_plotting);
            FDR_label = 'FDRcorrected';
        else
            data_selElec_allSubs.index_signkprime{i_win} = ...
                [data_selElec_allSubs.index_signkprime{i_win}; ...
                (Kprime_data.Exp.avgFolds.pval_ExpvsShuff{i_win}(ind_signelecs_StimCorr)...
                < pval_plotting)'];
            sumSignELec = sum(Kprime_data.Exp.avgFolds.pval_ExpvsShuff{i_win}(ind_signelecs_StimCorr)...
                < pval_plotting);
            FDR_label = 'uncorrected';
         end         
    end
    
    clear DataClean_AllTrials Kprime_data labels_loadedData corr_ttest SensorLabels
    disp(['-- ' sub ' finished in ' num2str(round(toc)) ' sec --']) 
end

%% 2) Define analysis parameters
%2.1 Define analysis parameters
toneDur_inSecs  = str2num(FuncInput_ToneDur_text);
nSamplesPerTone = toneDur_inSecs * fsample;

win_size    = SamplesTW;  %25 data points with fs 512 = ~49ms as mentioned in manuscript
s_per_win = (1/fsample)*win_size;
win_overlap = 0; %as mentioned in manuscript

windows = [1 win_size];
while windows(end,end) < nSamplesPerTone
    windows = [windows; windows(end,:) + (win_size - win_overlap)];
end

if windows(end,end) > nSamplesPerTone
    windows(end,:) = [];
end

windows_inms = (windows / fsample) * 1000;

%% 3) Plot exp Kprime values on surface

%3.1 find zlimits for constant plotting across time windows
% dv_absmax = 0;
% dv_absmin = 15;
% 
% for i_win = 1:size(windows,1)
%     dv_win_absmax = ceil(max(data_selElec_allSubs.kprime{i_win}...
%         (logical(data_selElec_allSubs.index_signkprime{i_win}))));
%     dv_win_absmin = floor(min(data_selElec_allSubs.kprime{i_win}...
%         (logical(data_selElec_allSubs.index_signkprime{i_win}))));
%     
%     if dv_win_absmax > dv_absmax
%         dv_absmax = dv_win_absmax;
%     end
%     if dv_win_absmin < dv_absmin
%         dv_absmin = dv_win_absmin;
%     end
% end

if plot_SubplotperTW == 1 %one subplot per TW
    h = figure('visible','on'); %ensures figure doesn't pop up during plotting
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    
    DimSubplot = [2,size(windows,1)/2];
    CounterSubplot = 1;
    SubplotPosition = [0 -0.05 0 0];    
end

%3.2 Plot data for current time window
for i_win = 1:size(windows,1)    
    
    kprime_signelecs = ...
        data_selElec_allSubs.kprime{i_win}...
        (logical(data_selElec_allSubs.index_signkprime{i_win}));
    index_signelecs = ...
        find(logical(data_selElec_allSubs.index_signkprime{i_win}));
    
    %3.2.1 Set up title information
    w1 = num2str( 1000 * windows(i_win, 1) / fsample , 3 );
    w2 = num2str( 1000 * windows(i_win, 2) / fsample , 3 );
    
    effect_title = ['Thresholded Exp Kprime values for all StimCorr sign. elecs (Exp vs. Shuff; p < ' num2str(pval_plotting) '; ' FDR_label ')'];
    win_title = ['TW = ' w1 ' - ' w2 ' ms'];
    
    if i_win == 1    
        stat_title = ['min/max(mean+-SD) Kprime = ' num2str(min(kprime_signelecs)) '/' num2str(max(kprime_signelecs)) ' (' ...
            num2str(round(mean(kprime_signelecs),2)) ' +- ' num2str(round(std(kprime_signelecs),2)) ')'];
    else
        stat_title = [num2str(min(kprime_signelecs)) '/' num2str(max(kprime_signelecs)) ' (' ...
            num2str(round(mean(kprime_signelecs),2)) ' +- ' num2str(round(std(kprime_signelecs),2)) ')'];
    end
    
    if ~isempty(index_signelecs)
        for i_elec = 1:length(data_selElec_allSubs.label)
            label_elecs{i_elec,1} = '';
%             label_elecs{i_elec,1} = data_selElec_allSubs.label{i_elec};
        end
    end
    
    sign_title = [num2str(length(index_signelecs))...
        ' / ' num2str(length(data_selElec_allSubs.kprime{i_win})) ...
        ' supra-thresh. elecs'];
    
    %3.2.2 Set up data struct for plotting
    coords      = data_selElec_allSubs.elecpos; %MNI coordinates for selected electrodes
    %Project all electrodes on one hemisphere, 
    coords(:,1) = abs(coords(:,1)) * -1; 
    
    vals        = data_selElec_allSubs.kprime{i_win}; %parameter of interest for resp. electrodes
    sign_index  = logical(data_selElec_allSubs.index_signkprime{i_win});    
    chanSize    = (vals./vals)*3; %electrode size (arbitrary, 3 is usually good)
    %clims = [0 max([dv_absmax, 1.1])]; %free scaling
    clims       = [0 15]; %fixed scaling
    cmap        = 'jet';
    view_angle  = [270,0];    
    textcolor_rgb   = [1 0 0];
    
    if plot_SubplotperTW == 0 %one separate figure per TW
        
        Figtitle = {['Allsub (n = ' num2str(length(subs)) ') - ' ...
            FuncInput_DataType ' - ToneDur: ' tonedur_title ' - ' win_title]...
            [effect_title]  [sign_title ' - ' stat_title]};      


        NASTD_ECoG_Plot_PlotSignElecsSurf_SubplotHemis...
            (coords,vals,sign_index,chanSize,clims,cmap) %Both H
        sgtitle(Figtitle)
            
        if save_poststepFigs == 1
            path_fig1 = [path_fig_KprimeSurf FDR_label '/'];
            if (~exist(path_fig1, 'dir')); mkdir(path_fig1); end
            filename     = ['3DSurf_SignExpKprime_Allsub__' FuncInput_DataType ...
                '_' tonedur_title '_' Label_TW num2str(i_win) ...
                '_p' num2str(pval_plotting) FDR_label '.png'];
            figfile      = [path_fig1 filename];
            saveas(gcf, [figfile], 'png'); %save png version
            close all;
        end
        
    elseif plot_SubplotperTW == 1 %one subplot per TW
        
        Figtitle = {['Allsub (n = ' num2str(length(subs)) ') - '  ...
            FuncInput_DataType ' - ToneDur: ' tonedur_title] [effect_title]};
        
        sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH...
            (coords, label_elecs, vals, sign_index,...
            chanSize, clims, cmap, textcolor_rgb, ...
            DimSubplot, CounterSubplot, SubplotPosition, [],...
            {[win_title],[sign_title],[stat_title]});

%         NASTD_ECoG_Plot_SubplotSignElecsSurf_SubplotHemis...
%             (coords,vals,sign_index,chanSize,clims,cmap,DimSubplot,CounterSubplot,...
%             {[win_title],[sign_title],[stat_title]});
%             CounterSubplot = CounterSubplot+2;     

        %Colorbar
        if CounterSubplot == 1
            ColorbarPosition    = [0.1 0.35 0 0]; %Left of surface
            c                   = colorbar;
            c.Position(1)       = ColorbarPosition(1); %sets colorbar to the right
            c.Position(2)       = ColorbarPosition(2); %sets colorbar higher
            c.Label.String      = 'Kprime';
            c.FontSize          = 12;
            caxis(clims)
        end
        
        CounterSubplot = CounterSubplot + 1;            

    end
end

if plot_SubplotperTW == 1 %one subplot per TW
    sgtitle(Figtitle)
    
    if save_poststepFigs == 1
        path_fig1 = [path_fig_KprimeSurf FDR_label '/'];
        if (~exist(path_fig1, 'dir')); mkdir(path_fig1); end
        filename     = ['3DSurf_SignExpKprime_Allsub_' FuncInput_DataType '_' ...
            tonedur_title '_' Label_TW 'all_p' num2str(pval_plotting) FDR_label '.png'];
        figfile      = [path_fig1 filename];
        saveas(gcf, [figfile], 'png'); %save png version
        close all;
    end
end

end