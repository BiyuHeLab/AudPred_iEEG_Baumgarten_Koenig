function NASTD_ECoG_HisTrack_PlotExpKvals_AllSubperTW...
    (subs, input_ToneDurLabel, input_DataLabel, ...
    SamplesTW, Label_TW, NumRuns, ...
    param, ...
    plot_SubplotperTW, save_poststepFigs, ...
    paths_NASTD_ECoG)

%Aim: Plot exp Kprime values for all electrodes and all subjects on MNI volume
%per TW and averaged across all TW

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

path_fig_KprimeSurf = [paths_NASTD_ECoG.ECoGdata_HisTrack ...
    'Allsub_n' num2str(length(subs)) '/Exp/' Label_TW '/Figs/Kprime_Surf/'];
if (~exist(path_fig_KprimeSurf, 'dir')); mkdir(path_fig_KprimeSurf); end

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
    
    %Load preprocessed data and Kprime data
    path_inputdata = [paths_NASTD_ECoG.ECoGdata_HisTrack sub '/Exp/' Label_TW '/'];
    load([si.path_preprocdata_sub]);
    fsample = DataClean_AllTrials.fsample; %read this out while preproc data in workspace
    load([path_inputdata sub '_' input_DataLabel '_' tonedur_title ...
        '_HisTrack_ExpKprimeCombFolds_unbalanced_' num2str(NumRuns) 'runs.mat']); %filename: ExpKprime_data, labels_loadedData
    
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
        ind_selelecs_MNI = [ind_selelecs_MNI, ...
            find(strcmp(DataClean_AllTrials.elec.label, ...
            labels_loadedData{ind_selelecs(i_selelec)}))];
    end
    
    data_selElec_allSubs.elecpos = ...
        [data_selElec_allSubs.elecpos; ...
        DataClean_AllTrials.elec.chanpos(ind_selelecs_MNI,:)];    

    %Read out selected elecs per sub,and append selected electrodes
    %to one common summary struct
    label_selelecs_currSub  = {};
    
    for i_selelec = 1:length(ind_selelecs)
        label_selelecs_currSub{i_selelec,1} = ...
            strcat(labels_loadedData{ind_selelecs(i_selelec)},'_', sub); 
        %Add sub-identifier to each elec label
    end
    
    data_selElec_allSubs.label                  = ...
        [data_selElec_allSubs.label; label_selelecs_currSub];
    data_selElec_allSubs.sub                    = ...
        [data_selElec_allSubs.sub; ones(length(label_selelecs_currSub),1)*i_sub];
    data_selElec_allSubs.elecs_persub{i_sub}    = ...
        find(data_selElec_allSubs.sub == i_sub);
           
    for i_win = 1:length(ExpKprime_data.avgRuns.Kprime_Avg)
        if i_sub == 1
            data_selElec_allSubs.kprime{i_win} = [];
        end        
        data_selElec_allSubs.kprime{i_win} = ...
            [data_selElec_allSubs.kprime{i_win}; ...
            ExpKprime_data.avgRuns.Kprime_Avg{i_win}(ind_selelecs,1)];
    end
    
    clear DataClean_AllTrials ExpKprime_data labels_loadedData ind_selelecs_MNI
    disp(['-- ' sub ' finished in ' num2str(round(toc)) ' sec --'])
end

%% 2) Define analysis parameters
%2.1 Define analysis parameters
toneDur_inSecs  = str2num(input_ToneDurLabel);
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

%% 3) Plot exp Kprime values per TW on surface

%3.1 find zlimits for constant plotting across time windows
dv_absmax = 0;
dv_absmin = 16;

for i_win = 1:size(windows,1)
    dv_win_absmax = ceil(max(data_selElec_allSubs.kprime{i_win}));
    dv_win_absmin = floor(min(data_selElec_allSubs.kprime{i_win}));
    
    if dv_win_absmax > dv_absmax
        dv_absmax = dv_win_absmax;
    end
    if dv_win_absmin < dv_absmin
        dv_absmin = dv_win_absmin;
    end
end

if plot_SubplotperTW == 1 %one subplot per TW
    h = figure('visible','off'); %ensures figure doesn't pop up during plotting
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    
    DimSubplot = [4,size(windows,1)/2];
    CounterSubplot = 1;
end

%3.2 Plot data for current time window
for i_win = 1:size(windows,1)
    
    %3.2.1 Set up title information
    w1 = num2str( 1000 * windows(i_win, 1) / fsample , 3 );
    w2 = num2str( 1000 * windows(i_win, 2) / fsample , 3 );
    
    effect_title = ['Exp Kprime values (no threshold - all elecs)'];
    win_title = ['TW = ' w1 ' - ' w2 ' ms'];
   
    if i_win == 1
        stat_title = ['min/max(mean+-SD) Kprime = ' ...
            num2str(min(data_selElec_allSubs.kprime{i_win})) '/' ...
            num2str(max(data_selElec_allSubs.kprime{i_win})) ' (' ...
            num2str(mean(data_selElec_allSubs.kprime{i_win})) ' +- ' ...
            num2str(std(data_selElec_allSubs.kprime{i_win})) ')'];
    else
         stat_title = [num2str(min(data_selElec_allSubs.kprime{i_win})) '/' ...
             num2str(max(data_selElec_allSubs.kprime{i_win})) ' (' ...
            num2str(mean(data_selElec_allSubs.kprime{i_win})) ' +- ' ...
            num2str(std(data_selElec_allSubs.kprime{i_win})) ')'];       
    end
    
    %3.2.2 Set up data struct for plotting
    clear dat
    dat.dimord  = 'chan_time';
    dat.time    = 0;
    dat.label   = data_selElec_allSubs.label;
    
    %3.2.1) Plot electrodes as spheres on MNI brain with each sphere color-coded
    %depending on Kprime-value
    coords      = data_selElec_allSubs.elecpos; %MNI coordinates for selected electrodes
    vals        = data_selElec_allSubs.kprime{i_win}; %parameter of interest for resp. electrodes
    chanSize    = (vals./vals)*2.25; %electrode size (arbitrary, 3 is usually good)
    %Optional: plot higher Kprime values bigger
    for i_elec = 1:length(vals)
        chanSize(i_elec,1) = ...
            (vals(i_elec)./vals(i_elec)) * ...
            (data_selElec_allSubs.kprime{i_win}(i_elec)/2);
    end
             %clims = [0 max([dv_absmax, 1.1])]; %free scaling
    clims       = [0 15]; %fixed scaling
    cmap        = 'jet';
    
    if plot_SubplotperTW == 0 %one separate figure per TW
        
        Figtitle = {['Allsub (n = ' num2str(length(subs)) ') - ' ...
            input_DataLabel ' - ToneDur: ' tonedur_title ' - ' win_title]...
            [effect_title ' - ' stat_title]};        

            NASTD_ECoG_Plot_PlotElecsSurf_SubplotHemis(coords,vals,chanSize,clims,cmap);
            suptitle(Figtitle)
        
        if save_poststepFigs == 1
            filename     = ['3DSurf_ExpKprime_Allsubn' num2str(length(subs)) '_' input_DataLabel ...
                '_' tonedur_title '_' Label_TW num2str(i_win) '_' param.ElecSelect 'elec.png'];
            figfile      = [path_fig_KprimeSurf filename];
            saveas(gcf, [figfile], 'png'); %save png version
            close all;            
        end
        
    elseif plot_SubplotperTW == 1 %one subplot per TW
        
        Figtitle = {['Allsub (n = ' num2str(length(subs)) ') - '  ...
            input_DataLabel ' - ToneDur: ' tonedur_title] [effect_title]};
        
            NASTD_ECoG_Plot_SubplotElecsSurf_SubplotHemis...
                (coords,vals,chanSize,clims,cmap,DimSubplot,CounterSubplot,...
                [win_title ' - ' stat_title]);
            CounterSubplot = CounterSubplot+2;
            
    end
end

if plot_SubplotperTW == 1 %one subplot per TW
    sgtitle(Figtitle)
    if save_poststepFigs == 1
        filename     = ['3DSurf_ExpKprime_Allsubn' num2str(length(subs)) '_' input_DataLabel ...
            '_' tonedur_title '_' Label_TW 'all_' param.ElecSelect 'elec.png'];
        figfile      = [path_fig_KprimeSurf filename];
        saveas(gcf, [figfile], 'png'); %save png version
        close all;
    end
end

%% 4) Plot exp Kprime values averaged across TW on surface
%Average Kprime values across all TW
for i_win = 1:size(windows,1)
    ExpKprime_data_perTW(:,i_win) = ...
        data_selElec_allSubs.kprime{i_win};
end
ExpKprime_data_avgTW = mean(ExpKprime_data_perTW,2);

%Plot avg Kprime
h = figure('visible','off');
set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen

effect_title = ['Exp Kprime values (no threshold - all elecs - ' ...
    num2str(NumRuns) ' runs)'];
win_title = ['Averaged across TW - max/mean(SD) Kprime = ' ...
    num2str(round(max(ExpKprime_data_avgTW),2)) ' / ' ...
    num2str(round(mean(ExpKprime_data_avgTW),2)) ...
    ' (' num2str(round(std(ExpKprime_data_avgTW),2)) ')'];

%Set up data struct for plotting
clear dat
dat.dimord  = 'chan_time';
dat.time    = 0;
dat.label   = data_selElec_allSubs.label;
vals        = ExpKprime_data_avgTW;
%     clims       = [0 max([dv_absmax, 1.1])]; %free scaling
clims       = [0 15]; %fixed scaling
cmap        = 'jet';
chanSize    = (vals./vals)*3; %electrode size (arbitrary)
%Optional: plot higher Kprime values bigger
for i_elec = 1:length(vals)
    chanSize(i_elec,1) = ...
        (vals(i_elec)./vals(i_elec)) * ...
        (ExpKprime_data_avgTW(i_elec)/2);
end

Figtitle = {['Allsub (n = ' num2str(length(subs)) ') - '  ...
    input_DataLabel ' - ToneDur: ' tonedur_title] [win_title] [effect_title]};

NASTD_ECoG_Plot_PlotElecsSurf_SubplotHemis...
    (coords,vals,chanSize,clims,cmap);
sgtitle(Figtitle)

if save_poststepFigs == 1
    filename = ...
        ['3DSurf_ExpKprime_Allsubn' num2str(length(subs)) '_' input_DataLabel ...
        '_' tonedur_title '_' Label_TW 'Avgall_' param.ElecSelect 'elec.png'];
    figfile      = [path_fig_KprimeSurf filename];
    saveas(gcf, figfile, 'png'); %save png version
    close all;    
end

end