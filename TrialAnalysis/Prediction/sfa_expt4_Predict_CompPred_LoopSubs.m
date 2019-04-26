%TJB: Based on original script: 'continuous_prediction_loop_jenn.m'

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

ks = [33]; %range of tone pitches based on which prediction should be computed
toneIndex = 33; %tone index indicating for which tone ERF/ERP should be computed, which is then used in regression analysis of prediciton effect

save_poststepFigs = 0;
% set(0, 'DefaultFigureVisible', 'off')

%% 1) Compute prediction effect (neural activity at tone 33 as function of
% p*34) and prediction error for single subjects
for baseline_correct = [0]
    for i_sub = 1:length(subs)
        for i_tonedur = 1:length(tonedur_text) %tone duration condition     
            sfa_expt4_Predict_CompPred_SingleSub(subs{i_sub}, baseline_correct, tonedur_text{i_tonedur}, ks, toneIndex);
        end
    end
end

%% 2) Plot prediction effect
for baseline_correct = [0]
    for i_sub = 1:length(subs)
        for i_tonedur = 1:length(tonedur_text) %tone duration condition     
            sfa_expt4_Predict_PlotPredict_SingleSub(subs{i_sub}, baseline_correct, tonedur_text{i_tonedur}, ks, toneIndex, save_poststepFigs);
        end
    end
end
