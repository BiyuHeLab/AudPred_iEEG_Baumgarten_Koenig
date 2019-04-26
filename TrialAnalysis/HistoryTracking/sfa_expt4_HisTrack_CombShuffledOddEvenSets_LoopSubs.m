%TJB: New version of tone_history_CV_plot_loop.m
%Aim: 1)Loops SHUFFLED k-value combination from even and odd sets and combined-set
%k'-value computation (via sfa_expt3_HisTrack_CombOddEvenSets_SLIDE16_TJB.m)
%across subjects, tone durations conditions, and baseline options
%2) Loops single-subject topo-plotting of k' values

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

subs = si.sub_list;
tonedur_text = {'0.2' '0.4'};

nReps = 100; %determine which files are loaded in and combined
IdentShuffleOrder_AllTineDur = 0; %identical shuffling across tone durs?

% ShuffleSets = 'appended'; %'averaged'
ShuffleSets = 'averaged'; %'averaged'


%% 1) Run SHUFFLED K-value combination 

for baseline_correct = [0] %[1 0 -1] %Loop baseline off/on/sequence beginnig
    for i_sub = 1:2%length(subs) %Loop subjects
        for i_tonedur = 1:length(tonedur_text) %tone duration condition      
            
            disp(['Currently processing: Subject: ' subs{i_sub} ...
                '; Baseline: ' num2str(baseline_correct)...
                '; Tone Duration Condition: ' tonedur_text{i_tonedur}...
                'ms; IdenticalShuffling4ToneOrder: ' num2str(IdentShuffleOrder_AllTineDur)])
            
            tic
   
            sfa_expt4_HisTrack_CombShuffledOddEvenSets...
                (subs{i_sub}, tonedur_text{i_tonedur}, baseline_correct, nReps, IdentShuffleOrder_AllTineDur, ShuffleSets)
%             sfa_expt4_HisTrack_CombShuffledOddEvenSets_GAMMA...
%                 (subs{i_sub}, tonedur_text{i_tonedur}, baseline_correct, nReps, IdentShuffleOrder_AllTineDur, ShuffleSets)
            
            disp(['done combining SHUFFLED k in ' num2str(toc) ' sec'])

        end
    end
end

