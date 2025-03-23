%% Set up paths (for data input/output, toolboxes (e.g., fieldtrip) and personal scripts)
function paths_sfa_expt4 = NASTD_paths(location)

restoredefaultpath; %Restores the MATLAB search path to installed products

if strcmp(location,'desktop') %for local matlab (Psychtoolbox installed)
    %% Add main path (gago/gogo)
    addpath('//gogo.sb.nyumc.org/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/') %main project path

    %% Add toolbox paths
    %%% Fieldtrip folders (gago/gogo)
%     addpath('//gogo.sb.nyumc.org/data/gogodisk4/thomas/toolboxes/fieldtrip-20180725/') %older FT - used for C2F and sfa_expt3
    addpath('//gogo.sb.nyumc.org/data/gogodisk4/thomas/toolboxes/fieldtrip-20190314/') %newer FT - used for sfa_expt4
    paths_sfa_expt4.FieldTrip = '//gogo.sb.nyumc.org/data/gogodisk4/thomas/toolboxes/fieldtrip-20190314/';
    paths_sfa_expt4.Freesurfer = '//gogo.sb.nyumc.org/data/gogodisk4/thomas/toolboxes/fieldtrip-20190314/external/freesurfer/';

    global ft_defaults;
    ft_default.trackusage = 'no';
    
    %% Add data paths
    %% Base & Scripts
    paths_sfa_expt4.BaseDir = '//gogo.sb.nyumc.org/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/';
    paths_sfa_expt4.ScriptsDir = '//gogo.sb.nyumc.org/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/NASTD_ECoG_Matlab/';
%     paths_sfa_expt4.SubfunctionsDir = '//gogo.sb.nyumc.org/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Scripts/Subfunctions/';

    %% Raw Data
    %determined via sfa_expt4_subjectinfo.m
    paths_sfa_expt4.RawData = '//gogo.sb.nyumc.org/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Data/';

    %% Analysis
    %Behavioral Hit rates
    paths_sfa_expt4.Analysis_Behavior = '//gogo.sb.nyumc.org/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Behavior/';
    
    %prepocessed ECoG data
    paths_sfa_expt4.Preproc_ECoGdata = '//gogo.sb.nyumc.org/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Preproc_ECoGdata/';
    paths_sfa_expt4.ECoGdata_HisTrack = '//gogo.sb.nyumc.org/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/HisTrack/';
    paths_sfa_expt4.ECoGdata_Prediction = '//gogo.sb.nyumc.org/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Prediction/';
   
    %% Figures
    %Behavior:
    paths_sfa_expt4.Fig_Behavior_SingleSubs = '//gogo.sb.nyumc.org/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Figures/Behavior/SingleSubs/';
    paths_sfa_expt4.Fig_Behavior_GroupAvg = '//gogo.sb.nyumc.org/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Figures/Behavior/GroupAvg/';
    %ECoG recordings
        %raw
        paths_sfa_expt4.Fig_ECoGdataraw_SingleSubs = '//gogo.sb.nyumc.org/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Figures/ECoGdata/raw_data/SingleSubs/';
        paths_sfa_expt4.Fig_ECoGdataraw_GroupAvg = '//gogo.sb.nyumc.org/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Figures/ECoGdata/raw_data/GroupAvg/';
        %History Tracking
        paths_sfa_expt4.Fig_HisTrack_SingleSubs = '//gogo.sb.nyumc.org/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Figures/ECoGdata/HisTrack/SingleSubs/';
        %Prediction
        paths_sfa_expt4.Fig_Prediction_SingleSubs = '//gogo.sb.nyumc.org/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Figures/ECoGdata/Prediction/SingleSubs/';

        
        
elseif strcmp(location,'server') %for gago/gogo matlab (no Psychtoolbox installed)
    %% Add script paths (gago/gogo)
    addpath('/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/'); %Own project folder

    %% Add toolbox paths
    %%% Fieldtrip folders (gago/gogo)
%     addpath('/data/gogodisk4/thomas/toolboxes/fieldtrip-20180725/') %older FT - used for C2F and sfa_expt3
    addpath('/data/gogodisk4/thomas/toolboxes/fieldtrip-20190314/') %newer FT - used for sfa_expt4
    paths_sfa_expt4.FieldTrip = '/data/gogodisk4/thomas/toolboxes/fieldtrip-20190314/';
    paths_sfa_expt4.Freesurfer = '/data/gogodisk4/thomas/toolboxes/fieldtrip-20190314/external/freesurfer/';

    global ft_defaults;
    ft_default.trackusage = 'no';
    %% Add data paths
    %% Base & Scripts
    paths_sfa_expt4.BaseDir = '/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/';
    paths_sfa_expt4.ScriptsDir = '/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/NASTD_ECoG_Matlab/';
%     paths_sfa_expt4.SubfunctionsDir = '/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Scripts/Subfunctions/';

    %% Raw Data
    %determined via sfa_expt4_subjectinfo.m
    paths_sfa_expt4.RawData = '/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Data/';
    
    %% Analysis
    %Behavioral Hit rates
    paths_sfa_expt4.Analysis_Behavior = '/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Behavior/';
    %prepocessed ECoG data
    paths_sfa_expt4.Preproc_ECoGdata = '/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Preproc_ECoGdata/';
    paths_sfa_expt4.ECoGdata_HisTrack = '/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/HisTrack/';
    paths_sfa_expt4.ECoGdata_Prediction = '/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Prediction/';
     
    %% Figures
    %Behavior:
    paths_sfa_expt4.Fig_Behavior_SingleSubs = '/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Figures/Behavior/SingleSubs/';
    paths_sfa_expt4.Fig_Behavior_GroupAvg = '/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Figures/Behavior/GroupAvg/';
     %ECoG recordings
        %raw
        paths_sfa_expt4.Fig_ECoGdataraw_SingleSubs = '/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Figures/ECoGdata/raw_data/SingleSubs/';
        paths_sfa_expt4.Fig_ECoGdataraw_GroupAvg = '/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Figures/ECoGdata/raw_data/GroupAvg/';
        %History Tracking
        paths_sfa_expt4.Fig_HisTrack_SingleSubs = '/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Figures/ECoGdata/HisTrack/SingleSubs/';
        %Prediction
        paths_sfa_expt4.Fig_Prediction_SingleSubs = '/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Figures/ECoGdata/Prediction/SingleSubs/';
        
end