function sfa_expt4_HisTrack_CombShuffledOddEvenSets(sub, tonedur_text, baseline_correct, nReps, IdentShuffleOrder_AllTineDur, ShuffleSets)
%TJB: New version of tone_history_CV_plot.m
%Aim: Combine results from both even and odd shuffled traning set regressions 
%1) Compute k-values for oth train and test set (SLIDE16 = last 16 tones, sliding time window)
%2) Save output (matrix k-values)
%TJB: New version of tone_history_CV_shuffle0_plot.m
%Aim: Combine results from both even and odd traning set regressions 
%(i.e., results from K-value computation (via sfa_expt4_HisTrack_CompSaveShuffledKval.m)
%Implementation:
%1) Combine null distribution k-values for both train and test set (SLIDE16 = last 16 tones, sliding time window)
%2) Save output combined output

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/data/gogodisk4/thomas/AuditoryPrediction/iEEG/')
% addpath('//gogo.sb.nyumc.org/data/gogodisk4/thomas/AuditoryPrediction/iEEG/')

sfa_expt4_setVars
paths_sfa_expt4 = sfa_expt4_paths(vars.location);

%Add base dir and own script dir
addpath(genpath(paths_sfa_expt4.BaseDir));
addpath(genpath(paths_sfa_expt4.ScriptsDir));

%% 0.2) Determine subject- and condition- specific parameters (whole-recording)
sfa_expt4_subjectinfo %load subject info file (var: si)
subs_PreProcSettings = sfa_expt4_subjectPreProcSettings; %load in file with individual preproc infos

%Tone duration condition for data load-in
if strcmp(tonedur_text,'0.2')
    tonedur_title = '200';
elseif  strcmp(tonedur_text,'0.4')
    tonedur_title = '400';
end

if baseline_correct == 1
    path_load_EXPk = [paths_sfa_expt4.ECoGdata_HisTrack sub '/baseline_corrected/Exp/'];
    path_load_SHUFFLEDk = [paths_sfa_expt4.ECoGdata_HisTrack sub '/baseline_corrected/Shuffled/'];
    path_save = [paths_sfa_expt4.ECoGdata_HisTrack sub '/baseline_corrected/EXPvsSHUFFLE/'];
    BC_text = 'BC';

elseif baseline_correct == 0
    path_load_EXPk = [paths_sfa_expt4.ECoGdata_HisTrack sub '/not_baseline_corrected/Exp/'];
    path_load_SHUFFLEDk = [paths_sfa_expt4.ECoGdata_HisTrack sub '/not_baseline_corrected/Shuffled/'];
    path_save = [paths_sfa_expt4.ECoGdata_HisTrack sub '/not_baseline_corrected/EXPvsSHUFFLE/'];
    BC_text = 'NBC';
end
    mkdir(path_save); %make new directory for output files

%% 1) Load in shuffled EVEN set
%1.1 load regstats data
if IdentShuffleOrder_AllTineDur == 1
    load([path_load_SHUFFLEDk sub '_HisTrackSLIDE16_regstatsSHUFFLEDidentical' num2str(nReps) 'Reps_' BC_text '_' tonedur_title 'msToneDur_trainEven.mat']);
else
    load([path_load_SHUFFLEDk sub '_HisTrackSLIDE16_regstatsSHUFFLED' num2str(nReps) 'Reps_' BC_text '_' tonedur_title 'msToneDur_trainEven.mat']);
end
% data format: sig2_test_shuffle{i_win, i_k, i_sensor}(1, nReps)
%matrix with sum of squared residuals divided by amount of tested models
%per time window, model order, sensor, repetition


%1.2 Define analysis parameters
fsample = 512;
toneDur_inSecs  = str2num(tonedur_text);
nSamplesPerTone = toneDur_inSecs * fsample;

win_size    = 25;  %25 data points with fs 512 = ~49ms as mentioned in manuscript
s_per_win = (1/fsample)*win_size;
win_overlap = 0; %as mentioned in manuscript

windows = [1 win_size];
while windows(end,end) < nSamplesPerTone
    windows = [windows; windows(end,:) + (win_size - win_overlap)];
end

if windows(end,end) > nSamplesPerTone
    windows(end,:) = [];
end

windows_inms = (windows / fsample) * 1000;

nSensors = size(sig2_test_shuffle,3); %all sensors, not only selected ones
nWindows = length(windows);

n_k      = size(sig2_test_shuffle,2); %number model order;

k_list = 1:n_k;

%1.3 Restructure regstats data and read out relevant stats for model selection from even set
for i_rep = 1:nReps
    for i_win = 1:nWindows
        for i_sensor = 1:nSensors        
            for i_k = 1:n_k
                s(i_k) = sig2_test_shuffle{i_win, i_k, i_sensor}(i_rep); 
                %read out shuffled sum of squared residuals divided by length of model order, per model order
            end            
            
            [min_s, min_s_ind] = min( s ); 
            %read out minimum squared residuals divided by length of model order and the respective position in array
            sig2_test_k_even_shuffle{i_win}(i_sensor, i_rep) = k_list(min_s_ind);
            %k-prime-value (i.e., preferred model order/Number of tones back best explaining MEG data)
             clear s %cleanup
        end
    end
end

clear sig2_test_shuffle
%% 2) Load in ODD set analysis
%2.1 load regstats data
if IdentShuffleOrder_AllTineDur == 1
    load([path_load_SHUFFLEDk sub '_HisTrackSLIDE16_regstatsSHUFFLED_IdentRepsAllToneDur_' nReps 'Reps_' BC_text '_' tonedur_title 'msToneDur_trainOdd.mat']);
else
    load([path_load_SHUFFLEDk sub '_HisTrackSLIDE16_regstatsSHUFFLED' num2str(nReps) 'Reps_' BC_text '_' tonedur_title 'msToneDur_trainOdd.mat']);
end

%1.2 Define analysis parameters
win_size    = 25;  %25 data points with fs 512 = ~49ms as mentioned in manuscript
s_per_win = (1/fsample)*win_size;
win_overlap = 0; %as mentioned in manuscript

fsample = 512;
toneDur_inSecs  = str2num(tonedur_text);
nSamplesPerTone = toneDur_inSecs * fsample;

windows = [1 win_size];
while windows(end,end) < nSamplesPerTone
    windows = [windows; windows(end,:) + (win_size - win_overlap)];
end

if windows(end,end) > nSamplesPerTone
    windows(end,:) = [];
end

windows_inms = (windows / fsample) * 1000;

nSensors = size(sig2_test_shuffle,3); %all sensors, not only selected ones
nWindows = length(windows);
n_k      = size(sig2_test_shuffle,2) %number model order;

k_list = 1:n_k;

%2.3 Restructure regstats data and read out relevant stats for model selection from odd set
for i_rep = 1:nReps
    for i_win = 1:nWindows
        for i_sensor = 1:nSensors        
            for i_k = 1:n_k
                s(i_k) = sig2_test_shuffle{i_win, i_k, i_sensor}(i_rep);
            end
            
            [min_s, min_s_ind] = min( s );
            sig2_test_k_odd_shuffle{i_win}(i_sensor, i_rep) = k_list(min_s_ind);
            clear s %cleanup
            
        end
    end
end

%% 3) Combine EVEN & ODD sets
%i.e., average model selection criteria by averaging across sets
if strcmp(ShuffleSets, 'appended')
    %3.1 Append both shuffled sets
    for i_win = 1:nWindows
        sig2_test_k_all_shuffle{i_win} = [sig2_test_k_even_shuffle{i_win}  sig2_test_k_odd_shuffle{i_win}];
    end %Output: Chan*trial (appended 100+100) matrix for each time window
    
elseif strcmp(ShuffleSets, 'averaged')           
    %3.2 Average across both shuffled sets
    for i_win = 1:nWindows
        for i_rep = 1:size(sig2_test_k_even_shuffle{i_win},2)
            sig2_test_k_all_shuffle{i_win}(:,i_rep) = mean([sig2_test_k_even_shuffle{i_win}(:,i_rep)  sig2_test_k_odd_shuffle{i_win}(:,i_rep)],2); %average across reps
        end
    end %Output: Chan*trial (averaged 100) matrix for each time window
end

%% 4) Load in combined unshuffled k-prime values and compute number of shuffled k-prime higher than original
%4.1 Load in combined unshuffled k-prime values
load([path_load_EXPk sub '_HisTrackSLIDE16_regstatsEXP_' BC_text '_' tonedur_title 'msToneDur_kModelComp_CombSets.mat']);

%4.2 i.e., computes the amount of shuffled k-prime values that are
%equal or bigger to the experimental/original k-prime value and
%divides this by the number of shuffled k-prime values to get the ratio
%for each time window and sensor, in each repetition
for i_win = 1:nWindows    
    for i_sensor = 1:nSensors
        %i.e., computes the amount of shuffled k-prime values that are
        %equal or bigger to the experimental/original k-prime value and
        %divides this by the number of reps/permutations to produce a p-value
        sig2_test_k_pval{i_win}(i_sensor) = sum( sig2_test_k_all_shuffle{i_win}(i_sensor, :) >= sig2_test_k{i_win}(i_sensor) ) ...
                                            / length(sig2_test_k_all_shuffle{i_win}(i_sensor, :));
                                                                              
    end
end 

%% 5) Save output variables
%Combined output file with 1) combined original k-prime values, combined
%2) shuffled k-prime values, 3) amount of shuffled_
if IdentShuffleOrder_AllTineDur == 1
    if strcmp(ShuffleSets, 'appended')
        savefile = [path_save sub '_HisTrackSLIDE16_regstatsSHUFFLEDidentical' num2str(nReps) '_' BC_text '_' tonedur_title 'msToneDur_kModelComp_appSets.mat'];
    elseif strcmp(ShuffleSets, 'averaged')           
        savefile = [path_save sub '_HisTrackSLIDE16_regstatsSHUFFLEDidentical' num2str(nReps) '_' BC_text '_' tonedur_title 'msToneDur_kModelComp_avgSets.mat'];
    end
else
    if strcmp(ShuffleSets, 'appended')
        savefile = [path_save sub '_HisTrackSLIDE16_regstatsSHUFFLED' num2str(nReps) 'Reps_' BC_text '_' tonedur_title 'msToneDur_kModelComp_appSets.mat']; %for appended shufle sets
    elseif strcmp(ShuffleSets, 'averaged')           
        savefile = [path_save sub '_HisTrackSLIDE16_regstatsSHUFFLED' num2str(nReps) 'Reps_' BC_text '_' tonedur_title 'msToneDur_kModelComp_avgSets.mat']; %for averaged shuffle sets
    end
end

save(savefile, 'sig2_test', 'sig2_test_k', 'sig2_test_k_sig', ...
               'sig2_test_k_even_shuffle', 'sig2_test_k_odd_shuffle', 'sig2_test_k_all_shuffle',...
               'sig2_test_k_pval');
        
end