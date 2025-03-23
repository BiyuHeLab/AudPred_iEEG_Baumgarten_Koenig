function NASTD_ECoG_Predict_PlotSignElec_AllSub...
    (subs, ...
    InputData_Type, Type_PredEffect, ToneDur_text, ...
    param, ...
    save_poststepFigs, paths_NASTD_ECoG)

%Aim: Plot electrodes showing sign. prediction effects per input data & TD
%condition aggregated over subjects.
%Output: Plot single figure per TW or subplot across TWs containing
%surface plots with marked sign. electrodes. Color-coding indicates t- or p-value.

% set(0, 'DefaultFigureVisible', 'on');

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')

%Add base dir and own script dir
addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

path_fig = ([paths_NASTD_ECoG.ECoGdata_Prediction ...
    '/PredEffects/Allsub_n' num2str(length(subs)) '/Figs/PredEffects_Surf/' Type_PredEffect '/']);
if (~exist(path_fig, 'dir')); mkdir(path_fig); end

%% 1) Load prediction effect data and aggregate relevant info across subjects
coords_allsub = [];
labels_allsub = [];
sub_index = [];

for i_sub = 1:length(subs)
    
    sub = subs{i_sub};
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;

    path_inputdata = [paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
    sub '/Data/' param.Label_TW 'sTW/'];

    load([path_inputdata sub '_PredEffects_' InputData_Type '_' ToneDur_text 'sTD.mat']);

    %Also load ECoG preproc data for channel labels and position
    loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];%path to indiv preprocessed ECoG data
    tic
    disp(['Loading preprocessed data set for sub: ' sub])
    load(loadfile_ECoGpreprocdata);
    disp(['done loading in ' num2str(toc) ' sec'])
    
    %Store coordinates and labels for analyzed electrodes
    for i_elec = 1:length(labels_loadedData)
        used_elecs_chanposIndex(1,i_elec) = ...
            find(strcmp(labels_loadedData{i_elec}, ...
            DataClean_AllTrials.elec.label));
    end
    coords_sub{i_sub} = ...
        DataClean_AllTrials.elec.chanpos(used_elecs_chanposIndex,:);
    coords_allsub = [coords_allsub; coords_sub{i_sub}];
    for i_elec = 1:length(used_elecs_chanposIndex)
        labels_allsub{end+1} = [labels_loadedData{i_elec} '_' sub];
    end
        
    %Select effect and store t and p-values from lin reg prediction analyis
    tval{i_sub} = eval([Type_PredEffect '.tval{1}']);
    sub_index = [sub_index; ones(length(coords_sub{i_sub}),1)*i_sub]; %used to differentiate subjects
    pval{i_sub} = eval([Type_PredEffect '.pval{1}']);
        
    %FDR-correction
    for i_win = 1:length(pval{i_sub})
        [~,~,~,pval_FDR{i_sub}{i_win}] = ...
            fdr_bh(pval{i_sub}{i_win}, param.pval_FDR, 'pdep','no');
    end

    %t- and p-value aggregation over subjects
    for i_win = 1:length(tval{i_sub})
        if i_sub == 1
            tval_allsubs{i_win} = tval{i_sub}{i_win}';
            pval_allsubs{i_win} = pval{i_sub}{i_win}';   
            pval_FDR_allsubs{i_win} = pval_FDR{i_sub}{i_win}';           
        else
            tval_allsubs{i_win} = [tval_allsubs{i_win}; tval{i_sub}{i_win}'];
            pval_allsubs{i_win} = [pval_allsubs{i_win}; pval{i_sub}{i_win}'];
            pval_FDR_allsubs{i_win} = [pval_FDR_allsubs{i_win}; pval_FDR{i_sub}{i_win}'];   
        end
    end
    
    %Cleanup
    used_elecs_chanposIndex = [];
    labels_loadedData = [];  
    clear PredEffect SimplePredErrEffect ComplexPredErrEffect
    
end
labels_allsub = labels_allsub';

%% 2) Define analysis parameters
SampleFreq = DataClean_AllTrials.fsample;
nSamplesPerTone = str2num(ToneDur_text) * SampleFreq;

win_size    = param.SamplesTW;  %25 data points with fs 512 = ~49ms as mentioned in manuscript
s_per_win = (1/SampleFreq)*win_size;
win_overlap = 0; %as mentioned in manuscript

windows = [1 win_size];
while windows(end,end) < nSamplesPerTone
    windows = [windows; windows(end,:) + (win_size - win_overlap)];
end

if windows(end,end) > nSamplesPerTone
    windows(end,:) = [];
end

windows_inms = (windows / SampleFreq) * 1000;

%% 3) find color limits for t-values
dv_absmax = 0;
for i_win = 1:size(windows,1)
    dv_win_absmax = max(tval_allsubs{i_win});
    if dv_win_absmax > dv_absmax
        dv_absmax = dv_win_absmax;
    end
end

%% 4) Plot figure
%4.1 Set up figure for TW-subplotting
set(gcf,'renderer','opengl');

if param.plot_SubplotperTW == 1
    h = figure('visible','off'); %ensures figure doesn't pop up during plotting
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    
    DimSubplot = [4,size(windows,1)/2]; %For both hemispheres
 
    CounterSubplot = 1;
    CounterSurfplot = 1;
    
    SizeFactor = 4;    
end

%4.2 Set up input data for current TW
for i_win = 1:size(windows,1)
    
    if param.plot_SubplotperTW == 0 %1 Figure per TW
        h = figure('visible','off'); %ensures figure doesn't pop up during plotting
        set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
        DimSubplot = [1,2];%For both hemispheres
        CounterSubplot = 1;
        CounterSurfplot = 1;        
        SizeFactor = 4;
    end
    
    %4.3 Determine sign. electrodes    
    if param.FDRcorrect == 1
        FDR_label = 'FDR corr';
        Index_SuprathresElecs{i_win} = pval_FDR_allsubs{i_win} < param.pval_plotting;
        sumSignELec = sum(pval_FDR_allsubs{i_win} < param.pval_plotting);
        pval_inputplot = pval_FDR_allsubs{i_win};
    else
        FDR_label = 'uncorr';
        Index_SuprathresElecs{i_win} = pval_allsubs{i_win} < param.pval_plotting;
        sumSignELec = sum(pval_allsubs{i_win} < param.pval_plotting);
        pval_inputplot = pval_allsubs{i_win};
    end
    
    %4.4 Set up title information
    w1 = num2str( round(1000 * windows(i_win, 1) / SampleFreq) , 3 );
    w2 = num2str( round(1000 * windows(i_win, 2) / SampleFreq) , 3 );
    win_title = ['TW = [' w1 ' - ' w2 ' ms]'];
    
    switch Type_PredEffect
        case 'PredEffect'
            effect_title = ['Prediction - ' InputData_Type ' - ' ...
                ToneDur_text 'sTD - p < ' num2str(param.pval_plotting) ' (' FDR_label ')'];
        case 'SimplePredErrEffect'
            effect_title = ['Simple Prediction error - ' InputData_Type ' - ' ...
                ToneDur_text 'sTD - p < ' num2str(param.pval_plotting) ' (' FDR_label ')'];
        case 'ComplexPredErrEffect'
            effect_title = ['Complex Prediction error - ' InputData_Type ' - ' ...
                ToneDur_text 'sTD - p < ' num2str(param.pval_plotting) ' (' FDR_label ')'];
    end
    
    %4.5 Determine labels of sign. elecs
    index_SignElec = find(Index_SuprathresElecs{i_win});
    labels_SignElec = labels_allsub(index_SignElec);
    if ~isempty(index_SignElec)
        sign_title = [num2str(sumSignELec)...
            ' / ' num2str(length(labels_allsub)) ' supra-thresh. elecs']; %1 line
    else
        sign_title = ['0 / ' num2str(length(labels_allsub)) ' supra-thresh. elecs']; %1 line
    end
    
    %4.6 Set up data struct for plotting
    clear dat
    dat.dimord = 'chan_time';
    dat.time   = 0;
    dat.label  = labels_allsub;
    dat.avg = tval_allsubs{i_win};
    dat.sign_elecs = Index_SuprathresElecs{i_win};   
    
    %4.7 Plot electrodes as spheres on MNI brain with each sphere color-coded
    %depending on t-/p-value
    vals = dat.avg; %t-value lin reg
    clims = [-4 4]; %fixed scaling
    
%     vals = pval_inputplot; %p-value lin reg    
%     clims = [0 0.05]; %fixed scaling

    chanSize = ones(1,length(dat.avg))*SizeFactor; %electrode size (arbitrary)
    %         clims = [max([dv_absmax, 1.1])*-1 max([dv_absmax, 1.1])]; %free scaling
    cmap = 'jet';   
    
    %4.8 Adjust subplot position according to number of subplots
    if param.plot_SubplotperTW == 1
        SubplotPosition = [0 0 0 0];
        ColorbarPosition = [0 0 0 0];
    else
        SubplotPosition = [0 0 0.2 0.2]; %enlarge and place higher
        ColorbarPosition = [0.1 0.3 0 0.3];
    end
    
%     %4.9 Plot surface with sign. electrodes
%     if param.FDRcorrect == 1 %Visible electrode labels for FDR corrected results
%         labels_plotting = labels_allsub;
%     else
        for i_elec = 1:length(labels_allsub)%No electrode labels for uncorrected results
            labels_plotting{i_elec} = '';
        end
%     end
    
    sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_SubplotHemis...
        (coords_allsub, labels_plotting, vals, Index_SuprathresElecs{i_win},...
        chanSize, clims, cmap, ...
        DimSubplot, CounterSubplot, SubplotPosition, ColorbarPosition, []);
    sp_handle_surf{CounterSurfplot} = sp_handle_surf_temp.L;
    sp_handle_surf{CounterSurfplot+1} = sp_handle_surf_temp.R;
    CounterSubplot = CounterSubplot+2;
    CounterSurfplot = CounterSurfplot+2;
      
    %4.10 Adjust header/title
    if param.plot_SubplotperTW == 1 %All TW in 1 plot (2 subplot (surface + 2D plot) per TW)
        title({[win_title] [sign_title]},'FontSize',10,'Interpreter','none')
    else %1 Plot per TW (2 subplot (surface + 2D plot) per TW)
%         title([sub effect_title [win_title  ' - ' sign_title]],'Interpreter','none') %Full title
        sgtitle(['            ' win_title],'FontSize',30,'Interpreter','none'); %TW only
    end
    
    %% 5) Save Figure
    if param.plot_SubplotperTW == 0 %1 plot per TW
        
        %Adjust colorbar
        ColorbarPosition = [0.1 0.35 0 0.5];
        h = colorbar;
        h.Position(1) = ColorbarPosition(1); %sets colorbar to the right
        h.Position(2) = ColorbarPosition(2); %sets colorbar higher
        h.Position(4) = h.Position(4)*ColorbarPosition(4); %makes colorbar shorter
        h.Label.String = ['t-statistic'];
        h.FontSize = 20;
        caxis(clims)
        h.Ticks = [-4, -2, 0, 2, 4];
            
        if save_poststepFigs == 1
            path_fig1 = [path_fig FDR_label '/pval' num2str(param.pval_plotting) '/'];
            if (~exist(path_fig1, 'dir')); mkdir(path_fig1); end
            filename     = ['Allsub_n' num2str(length(subs)) '_' ...
                Type_PredEffect FDR_label '_' InputData_Type ...
                '_' w1 '-' w2 'msTW_' ToneDur_text 'msTD.png'];
            figfile      = [path_fig1 filename];
            saveas(gcf, figfile, 'png'); %save png version
            close all;
        end
    end
    
end %End TW loop

if param.plot_SubplotperTW == 1 %All TW in 1 plot
    %Adjust subplot sizes
    SubplotPosition = [0 -1 0 0.11]; %slightly lower and enlarge
    ColorbarPosition = [0.1 0.25 0 6];
    
    for i_surfplot = 1:length(sp_handle_surf)
        if i_surfplot <= DimSubplot(2)
            pos = get(sp_handle_surf{i_surfplot},'Position');
            newpos = pos;
            newpos(2) = pos(2) -0.1;
            newpos(4) = pos(4) + SubplotPosition(4)/4;
            set(sp_handle_surf{i_surfplot},'Position',newpos)
        else
            pos = get(sp_handle_surf{i_surfplot},'Position');
            pos_above = get(sp_handle_surf{i_surfplot-DimSubplot(2)},'Position');
            newpos = pos;
            newpos(2) = pos_above(2) - 0.22;
            newpos(4) = pos(4)  + SubplotPosition(4)/4;
            set(sp_handle_surf{i_surfplot},'Position',newpos)
        end        
        
        if i_surfplot == 1
            %Adjust title
            %Adjust colorbar
            h = colorbar;
            h.Position(1) = ColorbarPosition(1); %sets colorbar to the right
            h.Position(2) = ColorbarPosition(2); %sets colorbar higher
            h.Position(4) = h.Position(4)*ColorbarPosition(4)/2; %makes colorbar shorter
            %                 h.Label.String = ['t-statistic (lin reg p33 activity on p*34)'];
            h.Label.String = ['t-statistic'];
            h.FontSize = 20;
            caxis(clims)
            h.Ticks = [-4, -2, 0, 2, 4];
        end
    end
    sub_title = ['Allsub_n' num2str(length(subs))];
    Figtitle = [sub_title ' - ' effect_title];
    sgtitle(Figtitle, 'Interpreter','none')
    
    %Save
    if save_poststepFigs == 1
        path_fig1 = [path_fig FDR_label '/pval' num2str(param.pval_plotting) '/'];
        if (~exist(path_fig1, 'dir')); mkdir(path_fig1); end
        filename     = ['Allsub_n' num2str(length(subs)) '_'...
            Type_PredEffect FDR_label '_' InputData_Type ...
            '_allTW_' ToneDur_text 'msTD.png'];
        figfile      = [path_fig1 filename];
        saveas(gcf, figfile, 'png'); %save png version
        close all;
    end
end

end