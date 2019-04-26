%TJB: New version of tone_history_CV_plot_loop.m
%Aim: 1)Loops k-value combination from even and odd sets and combined-set
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

pval = 0.05;
ShuffleSets = 'averaged'; %'appended'
saveplot = 1;

%% 1) Run K-value combination and k' computation across loops

for baseline_correct = [0] %[1 0 -1] %Loop baseline off/on/sequence beginnig
    for i_sub = 1:length(subs) %Loop subjects
            for i_tonedur = 1:length(tonedur_text) %tone duration condition      
            
            subs{i_sub}
            disp(['Currently processing: Baseline: ' num2str(baseline_correct) '; Subject: ' subs{i_sub} '; Tone Duration Condition: ' tonedur_text{i_tonedur} 'ms'])

%             sfa_expt4_HisTrack_CombOddEvenSets(subs{i_sub}, tonedur_text{i_tonedur}, baseline_correct);
            sfa_expt4_HisTrack_CombOddEvenSets_GAMMA(subs{i_sub}, tonedur_text{i_tonedur}, baseline_correct);

            end
    end
end


%% 2) Run single-subject k' topo-plotting across loops (no statistics)
for baseline_correct = [0] %[1 0 -1] %Loop baseline off/on
    for i_sub = 1:length(subs) %Loop subjects
            for i_tonedur = 1:length(tonedur_text) %tone duration condition      
                       
%             sfa_expt4_HisTrack_PlotCombSets_SingleSubs(subs{i_sub}, tonedur_text{i_tonedur}, baseline_correct, saveplot);
            sfa_expt4_HisTrack_PlotCombSets_SingleSubs_GAMMA(subs{i_sub}, tonedur_text{i_tonedur}, baseline_correct, saveplot);

            end
    end
end

%% 3) Run single-subject k' topo-plotting EXP vs. SHUFFLED across loops
nReps = 100;
for baseline_correct = [0] %[1 0 -1] %Loop baseline off/on
    for i_sub = 1:length(subs) %Loop subjects
            for i_tonedur = 1:length(tonedur_text) %tone duration condition      
                       
%             sfa_expt4_HisTrack_PlotEXPvsSHUFFLEk_SingleSubs(subs{i_sub}, tonedur_text{i_tonedur}, baseline_correct, nReps, pval, ShuffleSets, saveplot)
            sfa_expt4_HisTrack_PlotEXPvsSHUFFLEk_SingleSubs_GAMMA(subs{i_sub}, tonedur_text{i_tonedur}, baseline_correct, nReps, pval, ShuffleSets, saveplot)

            end
    end
end