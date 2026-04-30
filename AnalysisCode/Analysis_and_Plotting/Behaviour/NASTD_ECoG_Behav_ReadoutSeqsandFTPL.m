%Aim: Copy stimulus and behavioral data to structure (For Audrey to use
%for behavioral modeling  - April 2020)
%% 0) Specify paths
location = 'server'; %gago/gogo
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')

NASTD_ECoG_setVars
paths_NASTD = NASTD_ECoG_paths(location);

%Add base dir and own script dir
addpath(genpath(paths_NASTD.BaseDir));
addpath(genpath(paths_NASTD.ScriptsDir));

path_data = paths_NASTD.RawData;
path_fig = paths_NASTD.Fig_Behavior_SingleSubs;

i_sub = 1;
subs = vars.sub_list; %patients
sub = subs{i_sub};

tonedur_text = {'0.2' '0.4' 'all'};

 %% 1) load behavioral data, specify paths
for i_sub = 1:length(subs)
    
    sub = subs{i_sub}    
    NASTD_ECoG_subjectinfo %load subject info file (var: si)    
    load([si.path_behavioraldata_sub]);%Load raw behavioral data
    
    BehavData{i_sub}.sub = sub;
    BehavData{i_sub}.stim.betalevel = data.stim.beta;
    BehavData{i_sub}.stim.toneDur = data.stim.toneDur;
    BehavData{i_sub}.stim.logf_presentedfinaltonepitch = data.stim.logf_final;
    BehavData{i_sub}.stim.logf_predictedfinaltonepitch = data.stim.logf_pred;
    BehavData{i_sub}.stim.tonesequenceHz = data.stim.series_f;
    BehavData{i_sub}.stim.nTrials = data.stim.nTrials;
    BehavData{i_sub}.stim.nBlocks = data.stim.nBlocks;
    
    BehavData{i_sub}.resp.FTPLrating = data.resp_prob;

end