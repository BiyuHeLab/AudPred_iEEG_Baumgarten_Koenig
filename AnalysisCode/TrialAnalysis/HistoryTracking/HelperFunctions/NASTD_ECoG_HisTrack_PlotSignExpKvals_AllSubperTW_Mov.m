function NASTD_ECoG_HisTrack_PlotSignExpKvals_AllSubperTW_Mov...
    (subs, input_ToneDurLabel, input_DataLabel, ...
    SamplesTW, Label_TW,...
    param, ...
    FDR_correct, pval_plotting,...
    paths_NASTD_ECoG)

%Aim: Plot exp Kprime values for all electrodes and all subjects on MNI volume,per TW

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

path_vid = [paths_NASTD_ECoG.ECoGdata_HisTrack ...
    'Allsub_n' num2str(length(subs)) '/ExpvsShuff/' Label_TW '/Figs/Vid/'];
if (~exist(path_vid, 'dir')); mkdir(path_vid); end
path_frames = ([path_vid 'Frames/']);
if (~exist(path_frames, 'dir')); mkdir(path_frames); end

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
    
    %Determine selected electrode index to read out elec coords
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
    
    %Read out Kprime values for selected electrodes
    for i_win = 1:length(Kprime_data.Exp.avgRuns.Kprime_Avg)
        if i_sub == 1
            data_selElec_allSubs.kprime{i_win} = [];
            data_selElec_allSubs.index_signkprime{i_win} = [];
        end        
        data_selElec_allSubs.kprime{i_win} = ...
            [data_selElec_allSubs.kprime{i_win}; ...
            Kprime_data.Exp.avgRuns.Kprime_Avg{i_win}(ind_selelecs)];
    
    %Thresholding and FDR
         if FDR_correct == 1
            pval_FDRcorrected = mafdr(...
                Kprime_data.Exp.avgRuns.pval_ExpvsShuff{i_win}(ind_selelecs), ...
                'BHFDR', true);
            data_selElec_allSubs.index_signkprime{i_win} = ...
                [data_selElec_allSubs.index_signkprime{i_win}; ...
                (pval_FDRcorrected < pval_plotting)];
            sumSignELec = sum(pval_FDRcorrected < pval_plotting);
            FDR_label = 'FDRcorrected';
        else
            data_selElec_allSubs.index_signkprime{i_win} = ...
                [data_selElec_allSubs.index_signkprime{i_win}; ...
                (Kprime_data.Exp.avgRuns.pval_ExpvsShuff{i_win}(ind_selelecs)...
                < pval_plotting)];
            sumSignELec = sum(Kprime_data.Exp.avgRuns.pval_ExpvsShuff{i_win}(ind_selelecs)...
                < pval_plotting);
            FDR_label = 'uncorrected';
         end         
    end
    
    clear DataClean_AllTrials Kprime_data labels_loadedData 'corr_ttest'
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

%% 3) Plot exp Kprime values on surface
for i_win = 1:size(windows,1)    

    h = figure('visible','off'); %ensures figure doesn't pop up during plotting
    set(gcf,'renderer','opengl');   
    DimSubplot = [1 1];
    CounterSubplot = 1;
    SizeFactor      = 3;
    SubplotPosition = [0.05 -0.02 0 0];

%3.2 Plot data for current time window    
    kprime_signelecs = ...
        data_selElec_allSubs.kprime{i_win}...
        (logical(data_selElec_allSubs.index_signkprime{i_win}));
    index_signelecs = ...
        find(logical(data_selElec_allSubs.index_signkprime{i_win}));
    
    %3.2.1 Set up title information
    w1 = num2str( 1000 * windows(i_win, 1) / fsample , 3 );
    w2 = num2str( 1000 * windows(i_win, 2) / fsample , 3 );
    
    effect_title = ['Thresholded Exp Kprime values for all StimCorr sign. elecs (Exp vs. Shuff; p < ' num2str(pval_plotting) '; ' FDR_label ')'];
    win_title = [w1 ' - ' w2 ' ms'];
          
    if ~isempty(index_signelecs)
        for i_elec = 1:length(data_selElec_allSubs.label)
            label_elecs{i_elec,1} = '';
%             label_elecs{i_elec,1} = data_selElec_allSubs.label{i_elec};
        end
    else
        label_elecs = [];       
    end
       
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
        
    NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH...
        (coords, label_elecs, vals, sign_index,...
        chanSize, clims, cmap, textcolor_rgb, ...
        DimSubplot, CounterSubplot, SubplotPosition, [], []);    
   
    t = title(['Time window: ' win_title]);
        
    %Set Colorbar
    Label_Colorbar      = 'Kprime';
    ColorbarPosition    = [0.1 0.2 0 0.5]; %Left of surface
    h                   = colorbar;
    h.Position(1)       = ColorbarPosition(1); %sets colorbar to the right
    h.Position(2)       = ColorbarPosition(2); %sets colorbar higher
    h.Position(4)       = h.Position(4)*0.75; %makes colorbar shorter
    h.Label.String      = Label_Colorbar;
    h.FontSize          = 12;
    caxis(clims) 
    
    disp(['Processing Frame ' num2str(i_win)]);

 %Save frame    
    filename     = ['Frame' num2str(i_win) '.png'];
    figfile      = [path_frames filename];
    saveas(gcf, figfile, 'png'); %save png version
    close all;    
end

%% Aggregate frames to movie
cd(path_vid);
writerObj = VideoWriter(...
    ['Video_' ...
    input_DataLabel '_' input_ToneDurLabel 'sTD_' ...
    'p' num2str(pval_plotting) FDR_label '_' param.ElecSelect 'elec.avi']);

writerObj.FrameRate = 1;
open(writerObj)

for i_win = 1:size(windows,1)
    frame = sprintf([path_frames 'Frame' num2str(i_win) '.png']);
    input = imread(frame);
    
    writeVideo(writerObj, input);
end
close(writerObj);

end