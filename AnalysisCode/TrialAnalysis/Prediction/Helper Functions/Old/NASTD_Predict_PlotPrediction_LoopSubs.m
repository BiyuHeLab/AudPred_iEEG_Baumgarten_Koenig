%TJB: Based on original script: 'continuous_prediction_loop_jenn.m'

%% 0.1) Specify vars, paths, and setup fieldtrip
% addpath('/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')

NASTD_setVars
paths_NASTD = NASTD_paths(vars.location);

%Add base dir and own script dir
addpath(genpath(paths_NASTD.BaseDir));
addpath(genpath(paths_NASTD.ScriptsDir));

% addpath(genpath('/data/gogodisk4/thomas/AuditoryPrediction/iEEG/Scripts/FromRichard'));

%% 0.2) Determine subject-specific parameters (whole-recording)
i_sub = 1; %proxy to load subs_PreProcSettings
sub_list = vars.sub_list; %patients
sub = sub_list{i_sub};

NASTD_subjectinfo %load subject info file (var: si)
subs_PreProcSettings = NASTD_SubjectPreProcSettings; %load in file with individual preproc infos


subs = si.sub_list;
tonedur_text = {'0.2' '0.4'};

predictive_sequencerange = [33]; %range of tone pitches based on which prediction should be computed
toneIndex = 33; %tone index indicating for which tone ERF/ERP should be computed, which is then used in regression analysis of prediciton effect

save_poststepFigs = 0;
pval_plotting = 0.05; %p-value for plotting
%  set(0, 'DefaultFigureVisible', 'off') %show surface plots

%% 2) Plot prediction effect
for baseline_correct = [0]
    for i_sub = 1:length(subs)
        for i_tonedur = 1:length(tonedur_text) %tone duration condition 
            
            NASTD_Predict_PlotPrediction_Subs(subs{i_sub}, baseline_correct, ...
                tonedur_text{i_tonedur}, predictive_sequencerange, toneIndex, ...
                pval_plotting, save_poststepFigs); %ERP

            NASTD_Predict_PlotPrediction_Gamma_Subs(subs{i_sub}, baseline_correct, ...
                tonedur_text{i_tonedur}, predictive_sequencerange, toneIndex, ...
                pval_plotting, save_poststepFigs); %Gamma
       end
    end
end