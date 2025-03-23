function NASTD_ECoG_HisTrack_PlotExpKvalsBothTD_SsubperTW...
    (sub, input_ToneDurLabel, input_DataLabel, ...
    SamplesTW, Label_TW, NumRuns, ...
    param, ...
    plot_SubplotperTW, save_poststepFigs, ...
    paths_NASTD_ECoG)

%Aim: Contrast and plot exp Kprime values across both TD for all electrodes 
%on 3D surface plot for single subjects, per TW

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

path_inputdata = [paths_NASTD_ECoG.ECoGdata_HisTrack sub '/Exp/' Label_TW '/'];
path_fig = [paths_NASTD_ECoG.ECoGdata_HisTrack sub '/Exp/' Label_TW '/Figs/Kprime_DiffTD_Surf/'];
if (~exist(path_fig, 'dir')); mkdir(path_fig); end


%% 1) Load and contrast Kprime values across TD 
%Load preprocessed data for elec coords
NASTD_ECoG_subjectinfo %load subject info file (var: si)
tic
disp([' -- Loading preprocessed data set for sub: ' sub])
load([si.path_preprocdata_sub]);
disp([' -- Done loading in ' num2str(toc) ' sec'])

for i_TD = 1:length(input_ToneDurLabel)
    %Tone duration condition for data load-in
    if strcmp(input_ToneDurLabel{i_TD},'0.2')
        tonedur_title = '200msTD';
    elseif  strcmp(input_ToneDurLabel{i_TD},'0.4')
        tonedur_title = '400msTD';
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
    
    % Load expKprime data per TD
    load([path_inputdata sub '_' input_DataLabel '_' tonedur_title ...
        '_HisTrack_ExpKprimeCombFolds_unbalanced_' num2str(NumRuns) 'runs.mat'], ...
        'ExpKprime_data', 'labels_loadedData');
    
    %Aggregate Kprime across TD
    for i_win = 1:length(ExpKprime_data.avgRuns.Kprime_Avg)
        ExpKprime_bothTD{i_win}(:,i_TD) = ExpKprime_data.avgRuns.Kprime_Avg{i_win};
    end
    clear ExpKprime_data
end

num_win_common = min(num_win);

%Optional: Remove entries with kprime = 0 in any TD, since these mess up min/max
%estimates
for i_win = 1:num_win_common
    ind_infentries = [];
    ind_infentries = ...
        find(ExpKprime_bothTD{i_win}(:,1) == inf | ...
        ExpKprime_bothTD{i_win}(:,1) == 0);
    ExpKprime_bothTD{i_win}(ind_infentries,1) = NaN;
    ind_infentries = [];
    ind_infentries = ...
        find(ExpKprime_bothTD{i_win}(:,2) == inf | ...
        ExpKprime_bothTD{i_win}(:,2) == 0);
    ExpKprime_bothTD{i_win}(ind_infentries,2) = NaN;
end


%% 3) Compute Kprime ratio across tone durations
for i_win = 1:num_win_common
    RatioTD_ExpKprime{i_win} = ...
        ExpKprime_bothTD{i_win}(:,1) ./ ExpKprime_bothTD{i_win}(:,2);
end
    

%% 4) Determine sequence-tracking selected electrodes and their EDF-index
for i_TD = 1:length(input_ToneDurLabel)
    if strcmp(param.ElecSelect, 'StimCorr')
        path_inputdata = [paths_NASTD_ECoG.ECoGdata_StimCorr  sub '/Data/'];
        load([path_inputdata sub '_StimCorr_' input_DataLabel '_' ...
            input_ToneDurLabel{i_TD} 'sTD.mat'], ...
            'corr_ttest', 'SensorLabels');
        filt_signelecs_StimCorr_perTD(:,i_TD) = corr_ttest.p < 0.05; %only sign. sequence tracking elecs
    else
        filt_signelecs_StimCorr_perTD(:,i_TD) = true(length(labels_loadedData),1);  %all elecs
    end
end

%Select electrodes with sign. sign. sequence tracking effect in at least one TD
filt_signelecs_StimCorr = ...
    filt_signelecs_StimCorr_perTD(:,1) | filt_signelecs_StimCorr_perTD(:,2);
ind_selelecs  = find(filt_signelecs_StimCorr);
%nSensors_sel  = sum(filt_signelecs_StimCorr);
%labels_loadedData(ind_selelecs)

%Determine MNI index to select MNI coordinates
ind_selelecs_MNI   = [];
for i_selelec = 1:length(ind_selelecs)
    ind_selelecs_MNI = [ind_selelecs_MNI; ...
        find(strcmp(DataClean_AllTrials.elec.label, ...
        labels_loadedData{ind_selelecs(i_selelec)}))];
end


%% 5) Plot exp Kprime differences on surface
%Find color-limits for constant plotting across time windows
dv_absmax = 0;
dv_absmin = 0;
for i_win = 1:num_win_common
    dv_win_absmax = max(RatioTD_ExpKprime{i_win}(ind_selelecs,:));
    if dv_win_absmax > dv_absmax
        dv_absmax = dv_win_absmax;
    end
    dv_win_absmin = min(RatioTD_ExpKprime{i_win}(ind_selelecs,:));
    if dv_win_absmin < dv_absmin
        dv_absmin = dv_win_absmin;
    end   
end
dv_abs = max([abs(dv_absmax) abs(dv_absmin)]);
mean_DiffKprime_acrosswin = ...
    nanmean([nanmean(RatioTD_ExpKprime{1}(ind_selelecs,:)) ...
    nanmean(RatioTD_ExpKprime{2}(ind_selelecs,:)) ...
    nanmean(RatioTD_ExpKprime{3}(ind_selelecs,:)) ...
    nanmean(RatioTD_ExpKprime{4}(ind_selelecs,:))]);
std_DiffKprime_acrosswin = ...
    nanstd([nanmean(RatioTD_ExpKprime{1}(ind_selelecs,:)) ...
    nanmean(RatioTD_ExpKprime{2}(ind_selelecs,:)) ...
    nanmean(RatioTD_ExpKprime{3}(ind_selelecs,:)) ...
    nanmean(RatioTD_ExpKprime{4}(ind_selelecs,:))]);    

%Set up figure
if plot_SubplotperTW == 1 %one subplot per TW
    figure('visible','off');
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen    
    if strcmp(sub,'NY723') %sub requires both hemispheres
        DimSubplot = [2 6];
    else
        DimSubplot = [2 4];
    end
    CounterSubplot = 1;
end

%Set up title
effect_title1 = 'Exp Kprime ratio (200ms TD / 400 ms TD) for all elecs (no threshold)';
effect_title2 = ['max/min & mean +- STD across win = ' ...
    num2str(round(dv_absmax,2)) '/' num2str(round(dv_absmin,2)) ' & ' ...
    num2str(round(mean_DiffKprime_acrosswin,2)) ' +- ' ...
    num2str(round(std_DiffKprime_acrosswin,2))];
    
%Plot data for each time window
for i_win = 1:num_win_common
    
    %Set up title information    
    w1 = num2str( 1000 * windows{1}(i_win, 1) / fsample , 3 );
    w2 = num2str( 1000 * windows{1}(i_win, 2) / fsample , 3 );    
    win_title = {['TW = ' w1 ' - ' w2 ' ms '],...
        ['max/min & mean +- SD = ' ...
        num2str(max(round(RatioTD_ExpKprime{i_win}(ind_selelecs,:),2))) ' / ' ...
        num2str(min(round(RatioTD_ExpKprime{i_win}(ind_selelecs,:),2))) ' & ' ...
        num2str(round(nanmean(RatioTD_ExpKprime{i_win}(ind_selelecs,:)),2)) ' +- '...
        num2str(round(nanstd(RatioTD_ExpKprime{i_win}(ind_selelecs,:)),2))]};
    
    %Set up data struct for plotting
    clear dat
    dat.dimord  = 'chan_time';
    dat.time    = 0;
    dat.label   = labels_loadedData;
    dat.avg     = RatioTD_ExpKprime{i_win}(ind_selelecs,:);    
    coords      = DataClean_AllTrials.elec.chanpos(ind_selelecs_MNI,:); %MNI coordinates for selected electrodes
    vals        = dat.avg; %parameter of interest for resp. electrodes
    chanSize    = (vals./vals)*3; %electrode size (arbitrary)
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
            ' - ' win_title] [effect_title1] [effect_title2]};        
        if strcmp(sub,'NY723') %sub requires both hemispheres
            NASTD_ECoG_Plot_PlotElecsSurf_SubplotHemis...
                (coords,vals,chanSize,clims,cmap);
            sgtitle(Figtitle)
        else %1 hemisphere
            NASTD_ECoG_Plot_PlotElecsSurf...
                (coords,vals,chanSize,clims,cmap,view_angle);
            title(Figtitle, 'FontSize', 8)
        end
        
        if save_poststepFigs == 1
            filename     = ['3DSurf_DiffTD_ExpKprime' num2str(NumRuns) 'runs_' ...
                sub '_' input_DataLabel '_' num2str(i_win) 'TW_' ...
                param.ElecSelect 'elec.png'];
            figfile      = [path_fig filename];
            saveas(gcf, [figfile], 'png'); %save png version
            close all;            
        end
        
    elseif plot_SubplotperTW == 1 %one subplot per TW
        
        if strcmp(sub,'NY723')
            NASTD_ECoG_Plot_SubplotElecsSurf_SubplotHemis...
                (coords,vals,chanSize,clims,cmap,DimSubplot,CounterSubplot,win_title);
            CounterSubplot = CounterSubplot + 2;
        else
            NASTD_ECoG_Plot_SubplotElecsSurf...
                (coords,vals,chanSize,clims,cmap,view_angle,DimSubplot,CounterSubplot);
            title(win_title,'FontSize',8)
            CounterSubplot = CounterSubplot+1;
        end
        
    %Add 2D-scatterplot showing k' values per TD and average k' across TD
    subplot(DimSubplot(1),DimSubplot(2),CounterSubplot);
    plot(0:length(ExpKprime_bothTD{i_win}), 0:length(ExpKprime_bothTD{i_win}), ...
        'k-','LineWidth',2)
    hold on;
    scatter(ExpKprime_bothTD{i_win}(ind_selelecs,1), ExpKprime_bothTD{i_win}(ind_selelecs,2), ...
        200, RatioTD_ExpKprime{i_win}(ind_selelecs,:), '.')
    scatter(nanmean(ExpKprime_bothTD{i_win}(ind_selelecs,1)), ...
        nanmean(ExpKprime_bothTD{i_win}(ind_selelecs,2)), ...
        200, nanmean(RatioTD_ExpKprime{i_win}(ind_selelecs,:)), ...
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
    
    %Save Figure
    if save_poststepFigs == 1
        filename     = ['3DSurf_DiffTD_ExpKprime' num2str(NumRuns) 'runs_' ...
            sub '_' input_DataLabel '_' param.ElecSelect 'elec.png'];
        figfile      = [path_fig filename];
        saveas(gcf, figfile, 'png'); %save png version
        close all;
    end
end

end