function NASTD_ECoG_HisTrack_PlotExpKvalsBothTD_AllSubperTW...
    (subs, input_ToneDurLabel, input_DataLabel, ...
    SamplesTW, Label_TW, NumRuns, ...
    param, ...
    plot_SubplotperTW, save_poststepFigs, paths_NASTD_ECoG)

%Aim: Contrast and plot exp Kprime values across both TD on 3D surface plot
%combined for all subjects, per TW

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

path_fig = [paths_NASTD_ECoG.ECoGdata_HisTrack ...
    'Allsub_n' num2str(length(subs)) '/Exp/' Label_TW '/Figs/Kprime_DiffTD_Surf/'];
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

data_selElec_allSubs.filt_signelecs_StimCorr_allsub = [];

for i_sub = 1:length(subs)
    
    sub = subs{i_sub};
    
    tic
    %Load preprocessed data for elec coords
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    tic
    disp([' -- Loading preprocessed data set for sub: ' sub])
    load([si.path_preprocdata_sub]);
    disp([' -- Done loading in ' num2str(toc) ' sec'])
    
    for i_TD = 1:length(input_ToneDurLabel)
        if strcmp(input_ToneDurLabel{i_TD},'0.2')
            tonedur_title = '200msTD';
        elseif  strcmp(input_ToneDurLabel{i_TD},'0.4')
            tonedur_title = '400msTD';
        end
        
        %Load Kprime data per TD
        path_inputdata = [paths_NASTD_ECoG.ECoGdata_HisTrack sub '/Exp/' Label_TW '/'];
        load([path_inputdata sub '_' input_DataLabel '_' tonedur_title ...
            '_HisTrack_ExpKprimeCombFolds_unbalanced_' num2str(NumRuns) 'runs.mat'], ...
            'ExpKprime_data', 'labels_loadedData');
        
        %Place Kprime-values per TD for every electrode in common struct
        for i_win = 1:num_win_common
            if i_sub == 1
                data_selElec_allSubs.kprime_perTD{i_win}{i_TD} = [];
            end
            data_selElec_allSubs.kprime_perTD{i_win}{i_TD} = ...
                [data_selElec_allSubs.kprime_perTD{i_win}{i_TD}; ...
                ExpKprime_data.avgRuns.Kprime_Avg{i_win}];
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
        clear ExpKprime_data
    end
    
    %Combine sequence electrode electrode-selection across TD,
    %Select electrodes with sign. sign. sequence tracking effect in at least one TD
    filt_signelecs_StimCorr = ...
        filt_signelecs_StimCorr_perTD(:,1) | filt_signelecs_StimCorr_perTD(:,2);
    ind_selelecs  = find(filt_signelecs_StimCorr);
    %nSensors_sel  = sum(filt_signelecs_StimCorr);
    %labels_loadedData(ind_selelecs)
    
    
    %Aggregate electrode selection over subjects    
    data_selElec_allSubs.filt_signelecs_StimCorr_allsub = ...
        [data_selElec_allSubs.filt_signelecs_StimCorr_allsub; filt_signelecs_StimCorr];
    data_selElec_allSubs.filt_signelecs_StimCorr_allsub = ...
        logical(data_selElec_allSubs.filt_signelecs_StimCorr_allsub);
    
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
    
    disp(['-- ' sub ' finished in ' num2str(round(toc)) ' sec --'])
    
    clear DataClean_AllTrials labels_loadedData SensorLabels ...
        ind_selelecs ind_selelecs_MNI ...
        filt_signelecs_StimCorr_perTD filt_signelecs_StimCorr
    
end


%% 3) Compute Kprime ratio across TD for selected electrodes
data_selElec_allSubs.RatioTD_ExpKprime{i_win} = [];
for i_win = 1:num_win_common
    data_selElec_allSubs.RatioTD_ExpKprime{i_win} = ...
        data_selElec_allSubs.kprime_perTD{i_win}{1}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub) ./ ...
        data_selElec_allSubs.kprime_perTD{i_win}{2}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub);
end

% [data_selElec_allSubs.RatioTD_ExpKprime{i_win}, ...
%     data_selElec_allSubs.kprime_perTD{i_win}{1}, ...
%     data_selElec_allSubs.kprime_perTD{i_win}{2}]


%% 4) Plot figure (surface + scatter plot)
%Optional: Remove entries with kprime = 0 in any TD, since these mess up min/max
%estimates
for i_win = 1:num_win_common
    ind_infentries = [];
    ind_infentries = ...
        find(data_selElec_allSubs.RatioTD_ExpKprime{i_win} == inf | ...
        data_selElec_allSubs.RatioTD_ExpKprime{i_win} == 0);
    data_selElec_allSubs.RatioTD_ExpKprime{i_win}(ind_infentries) = NaN;
end

% %Optional: Restrict plotting only to positive or negative ratios
% for i_win = 1:num_win_common
%     ind_selentries = [];
%     ind_selentries = ...
%         find(data_selElec_allSubs.RatioTD_ExpKprime{i_win} < 1);
%     data_selElec_allSubs.RatioTD_ExpKprime{i_win}(ind_selentries) = NaN;
% end

%Find zlimits (for difference in K-prime between TD) for constant plotting across time windows
dv_absmax = 0;
dv_absmin = 0;
for i_win = 1:num_win_common
    dv_win_absmax = max(data_selElec_allSubs.RatioTD_ExpKprime{i_win});
    if dv_win_absmax > dv_absmax
        dv_absmax = dv_win_absmax;
    end
    dv_win_absmin = min(data_selElec_allSubs.RatioTD_ExpKprime{i_win});
    if dv_win_absmin < dv_absmin
        dv_absmin = dv_win_absmin;
    end
end
dv_abs = max([abs(dv_absmax) abs(dv_absmin)]);
mean_DiffKprime_acrosswin = ...
    nanmean(nanmean([...
    data_selElec_allSubs.kprime_perTD{1}{1}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub), ...
    data_selElec_allSubs.kprime_perTD{2}{1}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub), ...
    data_selElec_allSubs.kprime_perTD{3}{1}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub), ...
    data_selElec_allSubs.kprime_perTD{4}{1}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub)])) / ...
    nanmean(nanmean([...
    data_selElec_allSubs.kprime_perTD{1}{2}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub), ...
    data_selElec_allSubs.kprime_perTD{2}{2}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub), ...
    data_selElec_allSubs.kprime_perTD{3}{2}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub), ...
    data_selElec_allSubs.kprime_perTD{4}{2}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub)]));
std_DiffKprime_acrosswin = ...
    nanstd([nanmean(nanmean([...
    data_selElec_allSubs.kprime_perTD{1}{1}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub), ...
    data_selElec_allSubs.kprime_perTD{2}{1}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub), ...
    data_selElec_allSubs.kprime_perTD{3}{1}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub), ...
    data_selElec_allSubs.kprime_perTD{4}{1}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub)])), ...
    nanmean(nanmean([...
    data_selElec_allSubs.kprime_perTD{1}{2}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub), ...
    data_selElec_allSubs.kprime_perTD{2}{2}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub), ...
    data_selElec_allSubs.kprime_perTD{3}{2}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub), ...
    data_selElec_allSubs.kprime_perTD{4}{2}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub)]))]);

%Set up figure
if plot_SubplotperTW == 1 %All TW in 1 plot (2 subplot (surface + 2D plot) per TW)
    fig = figure('visible','off');
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    DimSubplot = [2, 6]; %Surface left, 2D plot right
    CounterSubplot = 1;
end
    
%Plot data for each time window
for i_win = 1:num_win_common
    
    %Set up title information
    w1 = num2str( 1000 * windows{1}(i_win, 1) / fsample , 3 );
    w2 = num2str( 1000 * windows{1}(i_win, 2) / fsample , 3 );
    win_title = {['TW = ' w1 ' - ' w2 ' ms '],...
        ['max/min: ' ...
        num2str(max(round(data_selElec_allSubs.RatioTD_ExpKprime{i_win},2))) ' / ' ...
        num2str(min(round(data_selElec_allSubs.RatioTD_ExpKprime{i_win},2))) ...
        '; mean +- SD across TD: ' ...
        num2str(round(...
        nanmean(data_selElec_allSubs.kprime_perTD{i_win}{1}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub)) / ...
        nanmean(data_selElec_allSubs.kprime_perTD{i_win}{2}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub)),2)) ' +- ' ...
        num2str(round(nanstd([...
        nanmean(data_selElec_allSubs.kprime_perTD{i_win}{1}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub)), ...
        nanmean(data_selElec_allSubs.kprime_perTD{i_win}{2}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub))]),2))], ...
        [num2str(sum(data_selElec_allSubs.RatioTD_ExpKprime{i_win} > 1)) ' elecs > 1; '...
        num2str(sum(data_selElec_allSubs.RatioTD_ExpKprime{i_win} < 1)) ' elecs < 1'...
        ]};
           
    %Set up data struct for plotting
    clear dat
    dat.dimord  = 'chan_time';
    dat.time    = 0;
    coords      = data_selElec_allSubs.elecpos; %MNI coordinates for selected electrodes
    vals        = data_selElec_allSubs.RatioTD_ExpKprime{i_win}; %parameter of interest for resp. electrodes
    chanSize    = (vals./vals)*2; %electrode size (arbitrary)
    %clims       = [0 dv_abs]; %free scaling
    clims       = [0 2]; %fixed scaling
    cmap        = 'jet';
    
    if plot_SubplotperTW == 0 %one separate figure per TW
        
        Figtitle = {['Allsub (n = ' num2str(length(subs)) ') - ' ...
            input_DataLabel ' - ' win_title]};
        
        NASTD_ECoG_Plot_PlotElecsSurf_SubplotHemis(coords,vals,chanSize,clims,cmap);
        suptitle(Figtitle)
        
        if save_poststepFigs == 1
            filename     = ['3DSurf_DiffTD_ExpKprime' num2str(NumRuns) 'runs_Allsubn' ...
                num2str(NumRuns) '_' input_DataLabel '_' num2str(i_win) 'TW_' ...
                param.ElecSelect 'elec.png'];
            figfile      = [path_fig filename];
            saveas(gcf, [figfile], 'png'); %save png version
            close all;
        end
        
    elseif plot_SubplotperTW == 1 %one subplot per TW
                
        NASTD_ECoG_Plot_SubplotElecsSurf_SubplotHemis...
            (coords,vals,chanSize,clims,cmap,DimSubplot,CounterSubplot,...
            [win_title]);
        CounterSubplot = CounterSubplot+2;
        
        %Add scatterplot showing k' values per TD
        subplot(DimSubplot(1),DimSubplot(2),CounterSubplot);
        plot(0:length(data_selElec_allSubs.kprime_perTD{i_win}{1}), 0:length(data_selElec_allSubs.kprime_perTD{i_win}{1}),'k-','LineWidth',2)
        hold on;
        scatter(data_selElec_allSubs.kprime_perTD{i_win}{1}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub), ...
            data_selElec_allSubs.kprime_perTD{i_win}{2}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub), ...
            20, ...
            data_selElec_allSubs.RatioTD_ExpKprime{i_win},'.')
%         %Avg based on averaging ratios across sign elecs
%         scatter(nanmean(data_selElec_allSubs.kprime_perTD{i_win}{1}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub)), ...
%             nanmean(data_selElec_allSubs.kprime_perTD{i_win}{2}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub)), ...
%             200, nanmean(data_selElec_allSubs.RatioTD_ExpKprime{i_win}), ...
%             'd', 'filled', 'MarkerEdgeColor',[0.2 0.2 0.2], 'LineWidth', 2)
        %Avg based on computing ratio from across-electrode averaged Kprime values 
        scatter(nanmean(data_selElec_allSubs.kprime_perTD{i_win}{1}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub)), ...
            nanmean(data_selElec_allSubs.kprime_perTD{i_win}{2}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub)), ...
            200, ...
            nanmean(data_selElec_allSubs.kprime_perTD{i_win}{1}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub)) ./ ...
            nanmean(data_selElec_allSubs.kprime_perTD{i_win}{2}(data_selElec_allSubs.filt_signelecs_StimCorr_allsub)), ...
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
        grid on   
        
        CounterSubplot = CounterSubplot+1;
        
    end
end

if plot_SubplotperTW == 1 %one subplot per TW
    %Set up and add title
    effect_title1 = 'Exp Kprime ratio (200ms TD / 400 ms TD) for all elecs (no threshold)';
    effect_title2 = ['max/min & mean +- STD across win = ' ...
        num2str(round(dv_absmax,2)) '/' num2str(round(dv_absmin,2)) ' & ' ...
        num2str(round(mean_DiffKprime_acrosswin,2)) ' +- ' ...
        num2str(round(std_DiffKprime_acrosswin,2))];    
    Figtitle = {['Allsub (n = ' num2str(length(subs)) ') - '  ...
        input_DataLabel] [effect_title1] [effect_title2]};
    sgtitle(Figtitle,'Interpreter','none')
    
    if save_poststepFigs == 1
        filename     = ['3DSurf_DiffTD_ExpKprime' num2str(NumRuns) 'runs_Allsubn' ...
            num2str(length(subs)) '_' input_DataLabel '_' param.ElecSelect 'elec.png'];
        figfile      = [path_fig filename];
        saveas(gcf, [figfile], 'png'); %save png version
        close all;
    end
end

end