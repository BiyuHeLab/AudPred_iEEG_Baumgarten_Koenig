function data_ECoGraw = NASTD_ECoG_Preproc_CreateFTstruct(i_sub, save_rawFTstruct, vars, preproc_dir)

%Aim: Load, match (MNI to EDF), and create FT-structure for raw ECoG data
%Based on script from Richard (03/14/19).

%% 0) Specify vars, paths, and setup fieldtrip (done by masterfunction)
paths_NASTD_ECoG = NASTD_ECoG_paths;

%Add base dir and own script dir
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));

%% 1) Determine subject-specific parameters 
sub_list = vars.sub_list(vars.validSubjs); %patients
% i_sub = 1;
sub = sub_list{i_sub};

NASTD_ECoG_subjectinfo %load subject info file (var: si)

%% 2) Match electrode info between MNI(anatomical location) and EDF(electrophysiol. recording) files
%Give output Elec_info, which contains info on electrode label, type, and
%number from both files. This info is then used to find matching
%electrodes. These matching electrodes are selected for further
%preprocessing.

%2.1) Determine paths for input data
    MNIfile = si.path_elecMNIcoord; %load in electrode labels, MNI coords, and type from .txt file
    EDFfile = si.path_rawdata_sub; %determine subject-specific EDF file %can be 512 Hz or 2kHz SF
    T1file  = si.path_elecT1coord; %load in electrode labels, T1 coords, and anatomical localtion label from .txt file
    
%2.2) Read in MNI data
    MNIelec_table = readtable(MNIfile);

%2.3) Read in T1 data incl. anatomical labels    
    [Elec_labels, T1x, T1y, T1z, Anat_labels] = ...
        importElecsfile_AnatLabels(T1file,1, size(MNIelec_table,1));
    T1ElecAnat_table = table(Elec_labels, T1x, T1y, T1z, Anat_labels);

%2.4) Read in EDF data 
    cfg = [];
    cfg.dataset = EDFfile;
    cfg.continuous = 'yes';
    cfg.channel = 'all';
    data_ECoGraw = ft_preprocessing(cfg);
    %Note about Warning OPENEDF: Physical Minimum larger than Maximum.
    %According to the EDF standard, inverted physical min/max values are 
    %used to implement the negative-up convention in clinical
    %neurophysiology.
    %https://mailman.science.ru.nl/pipermail/fieldtrip/2020-November/040486.html
    
    LabelEDF_all = data_ECoGraw.label;

%2.5) Match MNI and EDF channels by label, select strip+grid electrodes, 
    %and present summary electrode info file
    disp(['Found recordings for ' num2str(size(LabelEDF_all,1)) ...
        ' channels and localized ' num2str(size(MNIelec_table,1)) ' electrodes.'])

    %Sort channels for both MNI and EDF 
    for i_chan = 1:length(MNIelec_table{1:end,1})
        MNInames(i_chan,1) = MNIelec_table{i_chan,1};
    end
    for i_chan = 1:length(LabelEDF_all)
        EDFnames(i_chan,1) = LabelEDF_all(i_chan);
        EDFindex{i_chan,1} = num2str(i_chan);
    end
    
    %Present sorted channels   
    maxNumChannels = max(length(MNInames), length(EDFnames));
    if length(MNInames) < maxNumChannels %Fill empty channels with NaNs to get label vectors with equal length
        for i_chan = (length(MNInames)+1):maxNumChannels
            MNInames{i_chan} = 'zzz';%proxy
        end
    end
    if length(EDFnames) < maxNumChannels 
        for i_chan = (length(EDFnames)+1):maxNumChannels
            EDFnames(i_chan) = 'zzz';%proxy
        end
    end
    disp(['TASK1 (manual): Compare Electrode labels and note mismatches'])
    pause
    [sort(MNInames), sort(EDFnames)]
    pause
    %TASK1: Compare alphabetically sorted Channel names.
    %Aim: Note potential missmatch in labeling (e.g., T08 vs TT8)
    %(except for additional zeroes), note the corresponding indices (in
    %unsorted MNInames) and labels
    %in subs_PreProcSettings.(sub).MNI2EDFmismatches.
    %Note: EDFnames might contain channel labels that are not present in
    %MNInames (i.e., channels that could not be located).  
    disp(['TASK2 (manual): Note non-iEEG channels (i.e., ECG, DC; usally at the end of electrode array)'])
    pause
    [EDFnames EDFindex] 
    pause
    %TASK2: Check labels and indices of non-iEEG channels (i.e., ECG/EKG,
    %DC). Those are usually found at the every bottom of the EDFnames array.
    %Note this information in subs_PreProcSettings.(sub).number_nonEEGchan/number_ECoGchan
    %and subs_PreProcSettings.(sub).EKG_chanindex/Trigger_chanindex
     
    [Elec_info, SelectElec_info] = ...
        NASTD_ECoG_Preproc_GetMNIandEDFElecInfo...
        (sub, MNIelec_table, T1ElecAnat_table, LabelEDF_all);
    %Determines and matches indices and labels of non-matching electrodes

%2.5) Import MNI electrode locations/coordinates and store them in FT-compatible struct
    ElecStruct_MNI = ...
        NASTD_ECoG_Preproc_ImportMNILocations(sub, MNIfile, Elec_info.NumMNI_all); 

%2.6) Import T1 electrode locations/coordinates and store them in FT-compatible struct
    ElecStruct_T1 = ...
        NASTD_ECoG_Preproc_ImportT1Locations(T1file, Elec_info.NumMNI_all); 
    
%% 3) Create FT-compatible struct containing 
%Struct contains all electrodes present in EDF set (including non-iEEG channels),
%but spatial location and T1 anatomical label is only specified for 
%EDF-MNI matched electrodes.

    %3.1) Create struct from above input
    elec = NASTD_ECoG_Preproc_CreateElecStruct...
        (ElecStruct_MNI, Elec_info.LabelEDF_all,...
        Elec_info.EDF2MNI_ElecLabelmatch, Elec_info.T1AnatLabel_MNIorder, ...
        Elec_info.index_ECGchan, Elec_info.index_triggerchan);
    
    elec_T1 = NASTD_ECoG_Preproc_CreateElecStruct...
        (ElecStruct_T1, Elec_info.LabelEDF_all,...
        Elec_info.EDF2MNI_ElecLabelmatch, Elec_info.T1AnatLabel_MNIorder, ...
        Elec_info.index_ECGchan, Elec_info.index_triggerchan);
    
    elec.chanpos_T1 = elec_T1.chanpos;
    
    clear elec_T1
    
    %3.2) Copy struct into summary struct
    data_ECoGraw.elec = elec;

    %3.3) Preprocess summary file to select only relevant/chosen
    %channels (otherwise we get problems down the line)
    cfg = [];
    cfg.channel = [data_ECoGraw.elec.label; 'TRIG'];
    data_ECoGraw = ft_preprocessing(cfg, data_ECoGraw);
    
    %3.4) Copy info_elec subfield to final struct
    data_ECoGraw.info_elec.all      = Elec_info;
    data_ECoGraw.info_elec.selected = SelectElec_info;
        
    %3.5) Save summary file to indiv. preproc path
    if save_rawFTstruct == 1
        savefile = [preproc_dir sub '_rawdata_FTstruct.mat'];
        save(savefile, 'data_ECoGraw','-v7.3');    
    end

% %Test: Compare EDF and MNI data (labels and coordinates) in struct - should agree
% a = 100 %channel number
% elec.label(a) %channel label
% elec.chanpos(a,:) %channel position (MNI coordinates) 
% 
% b = find(strcmp(Elec_info.LabelMNI_adjusted_all,elec.label(a))) %find above electrode in MNI file
% ElecStruct_MNI.label(b)%channel label
% ElecStruct_MNI.chanpos(b,:)%channel position (MNI coordinates) 

clear elec Elec_info SelectElec_info EDFfile MNIfile LabelEDF_all MNIelec_table%cleanup

% %% 4) Visually inspect data 
% %Plot ECoG data in databrowser
% cfg = [];
% cfg.channel = data_ECoGraw.info_elec.selected.Label; %selected ECoG channels
% cfg.viewmode = 'vertical';
% cfg = ft_databrowser(cfg,data_ECoGraw);
% 
% %ECG data
% if ~isempty(data_ECoGraw.info_elec.all.index_ECGchan)
%     figure;
%     plot(data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_ECGchan(1),:) ...
%         - data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_ECGchan(2),:))
%     title([sub ' - ECG data'])
% end
% 
% %Trigger data
% figure;
% plot(data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,:))
% title([sub ' - Trigger data (Chan: ' data_ECoGraw.label{data_ECoGraw.info_elec.all.index_triggerchan} ')'])
% 
% % 5) Plot electrodes on MNI volume to inspect visual fit 
% LH = load([paths_NASTD_ECoG.FieldTrip '/template/anatomy/surface_pial_left.mat']);
% RH = load([paths_NASTD_ECoG.FieldTrip '/template/anatomy/surface_pial_right.mat']);
% 
% figure
% % set(gcf,'units','normalized','outerposition',[0 0 1 1])
% title([sub ' Electrode position on MNI standrad brain' ]) 
% ft_plot_mesh(LH.mesh);
% ft_plot_mesh(RH.mesh);
% % ft_plot_sens(ElecStruct_MNI, 'label','label'); %from MNI file
% ft_plot_sens(data_ECoGraw.elec, 'label','label','facecolor', [0 0 1], ...
%     'fontsize',8); %from combined FT-struct
% 
% view([-90 20])
% material dull
% lighting gouraud
% camlight

end
