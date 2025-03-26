clear

% - create fGn/fBm from (beta, sigma2) with proper range and final value
% made with series_selection.m
load series_selection_H.01_beta.99_1000.mat

% % - create a mirrored list with smaller sigma
% series_selection_sc 
%%% SCALED with a different sigma (original expt -- we did away with this in sfa_expt3)

% % - create list of final tones based on SD_epsilon differences from expected final tone
% series_selection_tf_SD
%%% spread around 440 Hz

% - create list of final tones based on pitch (log Hz) differences from previous tone
% series_selection_tf_Hz
%%% same fixed spacing of final tone pitch values

% - convert the stimuli to Hz, converging on 440 Hz
series_selection_scale2Hz

% - enumerate trial instantiations