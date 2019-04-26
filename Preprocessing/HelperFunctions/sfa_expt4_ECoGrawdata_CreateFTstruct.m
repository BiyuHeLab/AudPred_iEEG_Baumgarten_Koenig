function data_ECoGraw = sfa_expt4_ECoGrawdata_CreateFTstruct(i_sub, save_rawFTstruct, vars)

%TJB: Load, match (MNI to EDF), and FT-structure raw ECoG data
%Original script from Richard (03/14/19).

% %% 0) Specify vars, paths, and setup fieldtrip (done by masterfunction)
paths_sfa_expt4 = sfa_expt4_paths(vars.location);

% addpath('/data/gogodisk4/thomas/AuditoryPrediction/iEEG/')
% % addpath('//gogo.sb.nyumc.org/data/gogodisk4/thomas/AuditoryPrediction/iEEG/')
% 
% sfa_expt4_setVars
% paths_sfa_expt4 = sfa_expt4_paths(vars.location);
% 
%Add base dir and own script dir
addpath(genpath(paths_sfa_expt4.BaseDir));
addpath(genpath(paths_sfa_expt4.ScriptsDir));
addpath(genpath(paths_sfa_expt4.SubfunctionsDir));

% addpath(genpath('/data/gogodisk4/thomas/AuditoryPrediction/iEEG/Scripts/FromRichard'));

%% 1) Determine subject-specific parameters 
sub_list = vars.sub_list; %patients
% i_sub = 1;
sub = sub_list{i_sub};

sfa_expt4_subjectinfo %load subject info file (var: si)

% save_rawFTstruct = 0;
if save_rawFTstruct == 1
    %Define directory to save FT-struct raw data
    preproc_dir = [paths_sfa_expt4.Preproc_ECoGdata '/' sub '/'];
    mkdir(preproc_dir);
    cd(preproc_dir);
end

%% 2) Get info on electrodes and mact info from MNI(location) and EDF(recording) files
%Give output Elec_info, which contains info on electrode label, type, and
%number from both files. This info is then used to find matching
%electrodes. These matching electrodes are selected for further
%preprocessing.
%2.1) Determine paths for input data
    MNIfile = si.path_elecMNIcoord; %load in electrode labels & MNI coords from .txt file
    EDFfile = si.path_rawdata_sub; %determine subject-specific EDF file %can be 512 Hz or 2kHz SF

%2.2) Read MNI data from .txt file
    MNIelec_table = readtable(MNIfile);

%2.3) Read EDF data 
    cfg = [];
    cfg.dataset = EDFfile;
    cfg.continuous = 'yes';
    cfg.channel = 'all';
    data_ECoGraw = ft_preprocessing(cfg);
    LabelEDF_all = data_ECoGraw.label;

%2.4) Match channels, select strip+grid electrodes, and give out summary
%electrode info file 
    [Elec_info, SelectElec_info] = sfa_expt4_getMNIandEDFElecInfo(sub,MNIelec_table, LabelEDF_all);

%2.5) Import electrode locations/coordinates and store them in FT-compatible struct
    ElecStruct_MNI = import_MNI_locations(MNIfile,Elec_info.NumMNI_all); %reads in MNI locations for each electrode

%% 3) Create FT-compatible struct containing 
%Struct contains all electrodes present in EDF set (including non-iEEG channels),
%but location is only specified for EDF-MNI matched electrodes.
    %3.1) Create struct from above input
    elec = createElectrodeStruct(ElecStruct_MNI, Elec_info.LabelEDF_all,...
        Elec_info.EDF2MNI_ElecLabelmatch, Elec_info.index_ECGchan,...
        Elec_info.index_triggerchan);
    
    %3.2) Copy struct ato FT-loaded EDF data to get 1 summary
    %file
    data_ECoGraw.elec = elec;

    %3.3) preprocess summary file again to select only relevant/chosen
    %channels (otherwise we get problems down the line)
    cfg = [];
    cfg.channel = [data_ECoGraw.elec.label; 'TRIG'];
    data_ECoGraw = ft_preprocessing(cfg,data_ECoGraw);
    
    %3.4) Copy info_elec subfield to final struct
    data_ECoGraw.info_elec.all = Elec_info;
    data_ECoGraw.info_elec.selected = SelectElec_info;
        
    %3.5) Save summary file to indiv. preproc path
    if save_rawFTstruct == 1
        savefile = [preproc_dir sub '_rawdata_FTstruct.mat'];
        save(savefile, 'data_ECoGraw','-v7.3');    
    end

% %Test
% a = 100
% elec.label(a)
% elec.chanpos(a,:)
% 
% b = find(strcmp(Elec_info.LabelMNI_adjusted_all,elec.label(a)))
% ElecStruct_MNI.label(b)
% ElecStruct_MNI.chanpos(b,:)

clear elec Elec_info SelectElec_info EDFfile MNIfile LabelEDF_all MNIelec_table%cleanup

%% 4) Visually inspect data 
%ECoG data
% cfg = [];
% cfg.channel = data_ECoGraw.info_elec.selected.Label;
% cfg.viewmode = 'vertical';
% cfg = ft_databrowser(cfg,data_ECoGraw);
% 
% %ECG data
% if ~isempty(data_ECoGraw.info_elec.all.index_ECGchan)
%     plot(data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_ECGchan(1),:)-data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_ECGchan(2),:))
% end
% %Trigger data
% plot(data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,:))

%% 5) Plot electrodes on MNI volume to inspect visual fit 
% LH = load([paths_sfa_expt4.FieldTrip '/template/anatomy/surface_pial_left.mat']);
% RH = load([paths_sfa_expt4.FieldTrip '/template/anatomy/surface_pial_right.mat']);
% 
% figure
% title([sub ' Electrode position on MNI standrad brain' ])      
% 
% ft_plot_mesh(LH.mesh);
% ft_plot_mesh(RH.mesh);
% ft_plot_sens(ElecStruct_MNI);
% view([-90 20])
% material dull
% lighting gouraud
% camlight

end
