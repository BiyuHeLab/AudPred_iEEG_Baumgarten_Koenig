clear

% - create fGn/fBm from (beta, sigma2) with proper range and final value
% made with series_selection.m
load series_selection_H.01.mat

% - create a mirrored list with smaller sigma
series_selection_sc

% - create list of final tones based on SD_epsilon differences from expected final tone
series_selection_tf_SD

% - create list of final tones based on pitch (log Hz) differences from previous tone
% series_selection_tf_Hz

% - convert the stimuli to Hz, converging on 440 Hz
series_selection_scale2Hz

% - enumerate trial instantiations