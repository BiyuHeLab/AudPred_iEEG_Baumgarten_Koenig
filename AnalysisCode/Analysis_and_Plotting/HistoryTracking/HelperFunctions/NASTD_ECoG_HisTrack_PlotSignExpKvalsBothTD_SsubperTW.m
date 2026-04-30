function NASTD_ECoG_HisTrack_PlotSignExpKvalsBothTD_SsubperTW...
    (sub, input_ToneDurLabel, input_DataLabel, ...
    SamplesTW, Label_TW, NumRuns, ...
    FDR_correct, pval_threshold,...
    param, ...
    plot_SubplotperTW, save_poststepFigs, ...
    paths_NASTD_ECoG)

%Aim: Contrast and plot exp Kprime values across both TD for electrodes
%showing a sign. SHI effect (exp Kprime vs. shuffled Kprime)
%on 3D surface plot for single subjects, per TW

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

path_inputdata = [paths_NASTD_ECoG.ECoGdata_HisTrack sub '/ExpvsShuff/' Label_TW '/'];
path_fig = [paths_NASTD_ECoG.ECoGdata_HisTrack sub '/ExpvsShuff/' Label_TW '/Figs/KprimeTDRatio_Surf/'];
if (~exist(path_fig, 'dir')); mkdir(path_fig); end


%% 1) Load preprocessed data for elec coords
NASTD_ECoG_subjectinfo %load subject info file (var: si)
tic
disp([' -- Loading preprocessed data set for sub: ' sub])
load([si.path_preprocdata_sub]);
disp([' -- Done loading in ' num2str(toc) ' sec'])

%% 2) Load and aggregate Kprime values across TD 
Kprime_data_BothTD              = struct;
filt_signelecs_StimCorr_perTD   = [];

for i_TD = 1:length(input_ToneDurLabel) %tone duration condition
    %Tone duration condition for data load-in
    if strcmp(input_ToneDurLabel{i_TD},'0.2')
        tonedur_title = '200msTD';
        tonedur_subfield = 'TD200ms';
    elseif  strcmp(input_ToneDurLabel{i_TD},'0.4')
        tonedur_title = '400msTD';
        tonedur_subfield = 'TD400ms';
    end
    
    %Define analysis parameters
    fsample                 = DataClean_AllTrials.fsample;
    toneDur_inSecs{i_TD}    = str2num(input_ToneDurLabel{i_TD});
    nSamplesPerTone{i_TD}   = toneDur_inSecs{i_TD} * fsample;
    
    win_size                = SamplesTW;  %25 data points with fs 512 = ~49ms
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
    
    load([path_inputdata sub '_' input_DataLabel '_' tonedur_title ...
        '_HisTrack_ExpvsShuffKprimeCombFolds_unbalanced.mat'], ...
        'Kprime_data', 'labels_loadedData'); %filename: Kprime_data
    
    Kprime_data_BothTD.(tonedur_subfield) = Kprime_data.Exp.avgRuns;
    clear Kprime_data
    
%% 3) Determine sequence-tracking selected electrodes and their EDF-index
    if strcmp(param.ElecSelect, 'StimCorr')
        path_inputdata_StimCorr = [paths_NASTD_ECoG.ECoGdata_StimCorr  sub '/Data/'];
        load([path_inputdata_StimCorr sub '_StimCorr_' input_DataLabel '_' ...
            input_ToneDurLabel{i_TD} 'sTD.mat'], ...
            'corr_ttest', 'SensorLabels');
        filt_signelecs_StimCorr_perTD(:,i_TD) = corr_ttest.p < 0.05;%only sign. stim corr elecs
    else
        filt_signelecs_StimCorr_perTD(:,i_TD) = true(length(labels_loadedData),1);  %all elecs
    end
end

num_win_common = min(num_win);

%Select electrodes that show a sign. sequence tracking effect in at least one TD
filt_signelecs_StimCorr = ...
    filt_signelecs_StimCorr_perTD(:,1) | filt_signelecs_StimCorr_perTD(:,2);
ind_selelecs  = find(filt_signelecs_StimCorr);
%nSensors_sel  = sum(filt_signelecs_StimCorr);
%labels_loadedData(ind_selelecs)

%Determine MNI indices for selected electrodes
ind_selelecs_MNI   = [];
for i_selelec = 1:length(ind_selelecs)
    ind_selelecs_MNI = [ind_selelecs_MNI; ...
        find(strcmp(DataClean_AllTrials.elec.label, ...
        labels_loadedData{ind_selelecs(i_selelec)}))];
end    
    
disp(['done loading in ' num2str(toc) ' sec'])


%% 4) Determine superthreshold electrodes from selected electrodes
%(as those that show (FDR-corrected) sign. for exp vs. shuff p < 0.05))
for i_win = 1:num_win_common
    %Determine significance (FDR or non-FDR thresholding) for selected elecs
    for i_TD = 1:length(input_ToneDurLabel) 
        if strcmp(input_ToneDurLabel{i_TD},'0.2')
            tonedur_subfield = 'TD200ms';
        elseif  strcmp(input_ToneDurLabel{i_TD},'0.4')
            tonedur_subfield = 'TD400ms';
        end
        
        if FDR_correct == 1
            FDR_label = 'FDRcorrected';
            pval.(FDR_label).(tonedur_subfield){i_win} = ...
                mafdr(Kprime_data_BothTD.(tonedur_subfield).pval_ExpvsShuff{i_win}(filt_signelecs_StimCorr),'BHFDR', true);
            filt_SuprathresElecs.(tonedur_subfield){i_win} = pval.(FDR_label).(tonedur_subfield){i_win} < pval_threshold;
            %sumSignELec{i_win} = sum(pval.(FDR_label).(tonedur_subfield){i_win} < pval_threshold);
        else
            FDR_label = 'uncorrected';
            pval.(FDR_label).(tonedur_subfield){i_win} = ...
                Kprime_data_BothTD.(tonedur_subfield).pval_ExpvsShuff{i_win}(filt_signelecs_StimCorr);
            filt_SuprathresElecs.(tonedur_subfield){i_win} = pval.(FDR_label).(tonedur_subfield){i_win} < pval_threshold;
            %sumSignELec{i_win} = sum(pval.(FDR_label).(tonedur_subfield){i_win} < pval_threshold);
        end
    end
    
    %Combine across TD and select all elecs that show sign. SHI in at least
    %one TD (not in both)
    temp_bothpval = [pval.(FDR_label).TD200ms{i_win}, pval.(FDR_label).TD400ms{i_win}];
    pval.(FDR_label).BothTD{i_win} = min(temp_bothpval,[],2);
    filt_SuprathresElecs.BothTD{i_win} = ...
        pval.(FDR_label).BothTD{i_win} < pval_threshold;
    sumSignELec{i_win} = sum(filt_SuprathresElecs.BothTD{i_win});
    clear temp_bothpval
    
    %Read out corresponding Kprime values for selected elecs
    for i_TD = 1:length(input_ToneDurLabel) %tone duration condition
        if strcmp(input_ToneDurLabel{i_TD},'0.2')
            tonedur_subfield = 'TD200ms';
        elseif  strcmp(input_ToneDurLabel{i_TD},'0.4')
            tonedur_subfield = 'TD400ms';
        end
                
        Kprime_SelectedElecs{i_win}(:,i_TD) = ...
            Kprime_data_BothTD.(tonedur_subfield).Kprime_Avg{i_win}(filt_signelecs_StimCorr);
    end
    
%     %Optional: Remove entries with kprime = 0 in any TD, since these mess up min/max estimates
%     for i_win = 1:num_win_common
%         ind_infentries = [];
%         ind_infentries = ...
%             find(Kprime_SelectedElecs{i_win}(:,1) == inf | ...
%             Kprime_SelectedElecs{i_win}(:,1) == 0);
%         Kprime_SelectedElecs{i_win}(ind_infentries,1) = NaN;
%         ind_infentries = [];
%         ind_infentries = ...
%             find(Kprime_SelectedElecs{i_win}(:,2) == inf | ...
%             Kprime_SelectedElecs{i_win}(:,2) == 0);
%         Kprime_SelectedElecs{i_win}(ind_infentries,2) = NaN;
%     end
    
    
%% 5) Compute ratio between TD for plotting
    RatioTD_ExpKprime{i_win} = ...
        Kprime_SelectedElecs{i_win}(:,1) ./ ...
        Kprime_SelectedElecs{i_win}(:,2);
end


%% 6) Plot figure (surface + scatter plot)
%Find color-limits for constant plotting across time windows
dv_absmax = 0;
dv_absmin = 0;
for i_win = 1:num_win_common
    dv_win_absmax = max(RatioTD_ExpKprime{i_win}(filt_SuprathresElecs.BothTD{i_win}));
    if dv_win_absmax > dv_absmax
        dv_absmax = dv_win_absmax;
    end
    dv_win_absmin = min(RatioTD_ExpKprime{i_win}(filt_SuprathresElecs.BothTD{i_win}));
    if dv_win_absmin < dv_absmin
        dv_absmin = dv_win_absmin;
    end   
end
dv_abs = max([abs(dv_absmax) abs(dv_absmin)]);
mean_DiffKprime_acrosswin = ...
    nanmean([nanmean(RatioTD_ExpKprime{1}(filt_SuprathresElecs.BothTD{1},:)) ...
    nanmean(RatioTD_ExpKprime{2}(filt_SuprathresElecs.BothTD{2},:)) ...
    nanmean(RatioTD_ExpKprime{3}(filt_SuprathresElecs.BothTD{3},:)) ...
    nanmean(RatioTD_ExpKprime{4}(filt_SuprathresElecs.BothTD{4},:))]);
std_DiffKprime_acrosswin = ...
    std([nanmean(RatioTD_ExpKprime{1}(filt_SuprathresElecs.BothTD{1},:)) ...
    nanmean(RatioTD_ExpKprime{2}(filt_SuprathresElecs.BothTD{2},:)) ...
    nanmean(RatioTD_ExpKprime{3}(filt_SuprathresElecs.BothTD{3},:)) ...
    nanmean(RatioTD_ExpKprime{4}(filt_SuprathresElecs.BothTD{4},:))]);    

%Set up figure
if plot_SubplotperTW == 1 %one subplot per TW
    fig = figure('visible','off');
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen    
    if strcmp(sub,'NY723') %sub requires both hemispheres
        DimSubplot = [2 6];
    else
        DimSubplot = [2 4];
    end
    CounterSubplot = 1;
end

%Set up title
effect_title1 = ['Exp Kprime ratio (200ms TD / 400 ms TD) for superthresh elecs (p < ' ...
    num2str(pval_threshold) ', ' FDR_label ')'];
effect_title2 = ['max/min & mean +- STD across win = ' ...
    num2str(round(dv_absmax,2)) '/' num2str(round(dv_absmin,2)) ' & ' ...
    num2str(round(mean_DiffKprime_acrosswin,2)) ' +- ' ...
    num2str(round(std_DiffKprime_acrosswin,2))];

for i_win = 1:num_win_common
    
    %Set up title information
    w1 = num2str( 1000 * windows{1}(i_win, 1) / fsample , 3 );
    w2 = num2str( 1000 * windows{1}(i_win, 2) / fsample , 3 );  
    win_title = {['TW = ' w1 ' - ' w2 ' ms '],...
        ['max/min & mean +- SD = ' ...
        num2str(max(round(RatioTD_ExpKprime{i_win}(filt_SuprathresElecs.BothTD{i_win},:),2))) ' / ' ...
        num2str(min(round(RatioTD_ExpKprime{i_win}(filt_SuprathresElecs.BothTD{i_win},:),2))) ' & ' ...
        num2str(round(mean(RatioTD_ExpKprime{i_win}(filt_SuprathresElecs.BothTD{i_win},:)),2)) ' +- '...
        num2str(round(std(RatioTD_ExpKprime{i_win}(filt_SuprathresElecs.BothTD{i_win},:)),2))]};    
    sign_title = [num2str(sumSignELec{i_win})...
        ' / ' num2str(length(ind_selelecs_MNI)) ' supra-thresh. elecs']; %1 line
    
    %Set up data struct for plotting
    %Note: Input is restricted selected electrodes and sign. is provided additionally
    clear dat
    dat.dimord  = 'chan_time';
    dat.time    = 0;
    dat.label   = labels_loadedData;
    coords      = DataClean_AllTrials.elec.chanpos(ind_selelecs_MNI,:); 
    vals        = RatioTD_ExpKprime{i_win};    
    sign_index  = filt_SuprathresElecs.BothTD{i_win};    
    chanSize    = (vals./vals)*3; 
    %clims       = [0 dv_abs]; %free scaling
    clims       = [0 2]; %fixed scaling
    cmap        = 'jet';
    if strcmp(sub,'NY787')
        view_angle  = [90,0]; %Lat R
    else
        view_angle  = [270,0]; %Lat L
    end
    
    if plot_SubplotperTW == 0 %one separate figure per TW
        
        Figtitle = {[sub ' - ' input_DataLabel ...
            ' - ' win_title ' - ' sign_title] ...
            [effect_title1] [effect_title2]};
        
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
            filename     = ['3DSurf_DiffTD_SignExpKprime' num2str(NumRuns) 'runs_' ...
                sub '_' input_DataLabel '_' num2str(i_win) 'TW_' ...
                param.ElecSelect 'elec.png'];
            figfile      = [path_fig filename];
            saveas(gcf, [figfile], 'png'); %save png version
            close all;            
        end
        
    elseif plot_SubplotperTW == 1 %one subplot per TW
        
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
        
    %Add 2D-scatterplot showing sign. k' values per TD and average sign. k' across TD
    subplot(DimSubplot(1),DimSubplot(2),CounterSubplot);      
    plot(0:15, 0:15, 'k-', 'LineWidth', 2)
    hold on;
    scatter(Kprime_SelectedElecs{i_win}(filt_SuprathresElecs.BothTD{i_win},1), ...
        Kprime_SelectedElecs{i_win}(filt_SuprathresElecs.BothTD{i_win},2), ...
        200, ...
        RatioTD_ExpKprime{i_win}(filt_SuprathresElecs.BothTD{i_win},:),'.')
    scatter(nanmean(Kprime_SelectedElecs{i_win}(filt_SuprathresElecs.BothTD{i_win},1)), ...
        nanmean(Kprime_SelectedElecs{i_win}(filt_SuprathresElecs.BothTD{i_win},2)), ...
        200, nanmean(RatioTD_ExpKprime{i_win}(filt_SuprathresElecs.BothTD{i_win},:)), ...
        'd', 'filled', 'MarkerEdgeColor',[0.2 0.2 0.2], 'LineWidth', 2)    
    xlim([0 15])
    ylim([0 15])
    caxis([0 2])
    c = colorbar;
    c.Label.String = 'Ratio Kprime';
    c.Ticks = [0 0.5 1 1.5 2];
    colormap('jet');
    xlabel('kPrime 200ms TD')
    ylabel('kPrime 400ms TD')
    %pbaspect([1 1 1])   
    axis('square')
    grid minor   
     
    CounterSubplot = CounterSubplot+1;

    end
end

if plot_SubplotperTW == 1 %one subplot per TW
    Figtitle = {[sub ' - ' input_DataLabel] [effect_title1] [effect_title2]};    
    sgtitle(Figtitle,'Interpreter','none')

    if save_poststepFigs == 1
        filename     = ['3DSurf_DiffTD_ExpKprime' num2str(NumRuns) 'runs_' ...
            sub '_' input_DataLabel '_p' ...
            num2str(pval_threshold) FDR_label '_' param.ElecSelect 'elec.png'];
        figfile      = [path_fig filename];
        saveas(gcf, figfile, 'png'); %save png version
        close all;
    end
end

end