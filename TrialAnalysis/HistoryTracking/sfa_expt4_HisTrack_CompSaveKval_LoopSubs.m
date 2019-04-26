%TJB: New version of tone_history_CV_loop.m
%Aim: Loops K-value computation (via sfa_expt3_HisTrack_CompSaveKvalues_SLIDE16.m)
%across subjects, tone durations conditions, and sets (odd vs. even as train vs. test)

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/data/gogodisk4/thomas/AuditoryPrediction/iEEG/')
% addpath('//gogo.sb.nyumc.org/data/gogodisk4/thomas/AuditoryPrediction/iEEG/')

sfa_expt4_setVars
paths_sfa_expt4 = sfa_expt4_paths(vars.location);

%Add base dir and own script dir
addpath(genpath(paths_sfa_expt4.BaseDir));
addpath(genpath(paths_sfa_expt4.ScriptsDir));

% addpath(genpath('/data/gogodisk4/thomas/AuditoryPrediction/iEEG/Scripts/FromRichard'));

%% 0.2) Determine subject-specific parameters (whole-recording)
i_sub = 1; %proxy to load subs_PreProcSettings
sub_list = vars.sub_list; %patients
sub = sub_list{i_sub};

sfa_expt4_subjectinfo %load subject info file (var: si)
subs_PreProcSettings = sfa_expt4_subjectPreProcSettings; %load in file with individual preproc infos

plot_poststepFigs = 0;
save_poststepFigs = 0;

subs = si.sub_list;
tonedur_text = {'0.2' '0.4'};

%% 1) Run K-value computation  for ERF (via sfa_expt4_HisTrack_CompSaveKval.m) across loops
for baseline_correct = 0 %[0 1] %Loop baseline off/on
    for i_sub = 1:length(subs) %Loop subjects
        for i_tonedur = 1:length(tonedur_text) %tone duration condition     
            for i_odd = [0 1] %Loop train set - i_odd = 0 == train on even set
                
                 sfa_expt4_HisTrack_CompSaveKval(subs{i_sub}, baseline_correct, i_odd, tonedur_text{i_tonedur}); %Low Freq
                 
                 sfa_expt4_HisTrack_CompSaveKval_GAMMA(subs{i_sub}, baseline_correct, i_odd, tonedur_text{i_tonedur}); %Gamma
 
            end
        end
    end
end

