function sfa_expt4_Predict_CompPred_SingleSub(sub, baseline_correct, tonedur_text, ks, toneIndex)
%TJB: Based on original script 'continuous_prediction_jenn.m'
%Aim: How is neural activity at tone 33 modulated by the expected value of
%tone 34, which itself depends on the previous tone sequence input. In
%other words, do tones 1:32 influence the neural activity at tone 33?
%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/data/gogodisk4/thomas/AuditoryPrediction/iEEG/')
% addpath('//gogo.sb.nyumc.org/data/gogodisk4/thomas/AuditoryPrediction/iEEG/')

sfa_expt4_setVars
paths_sfa_expt4 = sfa_expt4_paths(vars.location);

%Add base dir and own script dir
addpath(genpath(paths_sfa_expt4.BaseDir));
addpath(genpath(paths_sfa_expt4.ScriptsDir));

if baseline_correct == 1
    path_save = [paths_sfa_expt4.ECoGdata_Prediction sub '/baseline_corrected/Exp/'];
    BC_text = 'BC';

elseif baseline_correct == 0
    path_save = [paths_sfa_expt4.ECoGdata_Prediction sub '/not_baseline_corrected/Exp/'];
    BC_text = 'NBC';
end

%% 0.2) Determine subject-specific parameters (whole-recording)
sfa_expt4_subjectinfo %load subject info file (var: si)
subs_PreProcSettings = sfa_expt4_subjectPreProcSettings; %load in file with individual preproc infos

%Tone duration condition for data load-in
if strcmp(tonedur_text,'0.2')
    tonedur_title = '200';
elseif  strcmp(tonedur_text,'0.4')
    tonedur_title = '400';
end

%Define directory for Kprime data
Predict_dir = [paths_sfa_expt4.ECoGdata_Prediction sub '/'];
mkdir(Predict_dir);
cd(Predict_dir);

loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];%path to indiv preprocessed ECoG data
loadfile_behav = [si.path_behavioraldata_sub];%path to indiv behavioral/stimulus data

mkdir(path_save); %make new directory for output files


%% 1) load preprocessed ECoG data
tic
disp('loading...')
load(loadfile_ECoGpreprocdata);
load(loadfile_behav);
disp(['done loading in ' num2str(toc) ' sec'])

%% 2) Select input data 
%2.1 Select only trials present in both behavioral and ECoG data and bring
%them into same order
%Determine number of trials present in behav data
numoftrials_behav = num2str(length(data.trialNum)); %Should be around 120
numoftrials_ECoG = num2str(length(data_ECoGfiltref_trials.trial)); %Should be around 120

%Restrict data sets to those valis trials (i.e., without no-response trials)
    %Should be same amount, possibly with different trial indexing
index_goodtrials_behav = subs_PreProcSettings.(sub).goodtrials_afterCAR_indexBehav;
index_goodtrials_ECoG = subs_PreProcSettings.(sub).goodtrials_afterCAR_indexECoG;

index_goodtrials_both = [index_goodtrials_behav; index_goodtrials_ECoG]; %for checking

%Restrict trials to those with good ECoG recording
    %Beware of potential different trial indexing (i.e., make sure to
    %delete same trials for behavior and ECoG)
if isequal(index_goodtrials_behav, index_goodtrials_ECoG) %If trialstructs are identical
    index_goodtrials_ECoG(subs_PreProcSettings.(sub).rejectedTrials) = []; %delete bad trials in ECoG struct filter
    index_goodtrials_behav(subs_PreProcSettings.(sub).rejectedTrials) = []; %delete bad trials in behav struct filter
else %If trialstructs are NOT identitcal
    for i_rejecttrial = 1:length(subs_PreProcSettings.(sub).rejectedTrials) %for each to be removed trial
        index_goodtrials_behav(find(index_goodtrials_ECoG == subs_PreProcSettings.(sub).rejectedTrials(i_rejecttrial))) = []; %delete behav trial
        index_goodtrials_ECoG (find(index_goodtrials_ECoG == subs_PreProcSettings.(sub).rejectedTrials(i_rejecttrial))) = []; %delete ECoG trial
    end
end
%Output: trial-filter with different trial indices for ECoG  behav, but
%same trials rejected in both cases. If filter are used to select trials
%from the respective struts, output matches again (in terms of trialcontent
%to trialorder)

%2.2 Select chosen trials from both behav & ECoG structs

%Behavioral data
dat_fields  = fieldnames(data);
for i = 1:length(dat_fields)
    if eval(['length(data.' dat_fields{i} ') == ' numoftrials_behav])
         eval(['data.' dat_fields{i} ' = data.' dat_fields{i} '(index_goodtrials_behav);']);
    end
end

%Stimulus data
stim_fields = fieldnames(data.stim);
for i = 1:length(stim_fields)
    if eval(['length(data.stim.' stim_fields{i} ') == ' numoftrials_behav])
         eval(['data.stim.' stim_fields{i} ' = data.stim.' stim_fields{i} '(index_goodtrials_behav);']);
    end
end

%ECoG data
ECoG_fields = fieldnames(data_ECoGfiltref_trials);
for i = 1:length(ECoG_fields)
    if eval(['size(data_ECoGfiltref_trials.' ECoG_fields{i} ') == [1 ' numoftrials_ECoG ']'])
        eval(['data_ECoGfiltref_trials.' ECoG_fields{i} ' = data_ECoGfiltref_trials.' ECoG_fields{i} '(index_goodtrials_ECoG);']);
    elseif eval(['size(data_ECoGfiltref_trials.' ECoG_fields{i} ') == [' numoftrials_ECoG ' 2]'])
        eval(['data_ECoGfiltref_trials.' ECoG_fields{i} ' = data_ECoGfiltref_trials.' ECoG_fields{i} '(index_goodtrials_ECoG,:);']);
    end
end

%2.3 Select trials with specific tone length 
    %should be 60 trials per tone dur
    %tone length ID: 0.2 = 1, 0.4 = 2

%Update trial numbers    
% numoftrials = num2str(length(data.trialNum));
numoftrials_behav = num2str(length(data.trialNum)); %Should be around 120
numoftrials_ECoG = num2str(length(data_ECoGfiltref_trials.trial)); %Should be around 120

filt_tonedur = data.stim.toneDur == str2double(tonedur_text); %filter to select trials for certain tone dur 

%Behavioral data
dat_fields  = fieldnames(data);
for i = 1:length(dat_fields)
    if eval(['length(data.' dat_fields{i} ') == ' numoftrials_behav]) 
        eval(['data.' dat_fields{i} ' = data.' dat_fields{i} '(filt_tonedur);']);
    end
end

%Stimulus data
stim_fields = fieldnames(data.stim);
for i = 1:length(stim_fields)
    if eval(['length(data.stim.' stim_fields{i} ') == ' numoftrials_behav]) 
        eval(['data.stim.' stim_fields{i} ' = data.stim.' stim_fields{i} '(filt_tonedur);']);
    end
end

%ECoG data
ECoG_fields = fieldnames(data_ECoGfiltref_trials);
for i = 1:length(ECoG_fields)
    if eval(['size(data_ECoGfiltref_trials.' ECoG_fields{i} ') == [1 ' numoftrials_ECoG ']'])
        eval(['data_ECoGfiltref_trials.' ECoG_fields{i} ' = data_ECoGfiltref_trials.' ECoG_fields{i} '(filt_tonedur);']);
    elseif eval(['size(data_ECoGfiltref_trials.' ECoG_fields{i} ') == [' numoftrials_ECoG ' 2]'])
        eval(['data_ECoGfiltref_trials.' ECoG_fields{i} ' = data_ECoGfiltref_trials.' ECoG_fields{i} '(filt_tonedur,:);']);    
    end
end

%2.4 Create common input struct
data_input = data_ECoGfiltref_trials;
data_input.behav = data;
data_input.stim = data.stim;

%% 3) Extract last two tones per trial, as well as relevant stim/behav info
%3.1) Define time window parameters
nTrials  = length(data_input.trial);
% nSensors = size( data_input.trial{1}, 1 ); %all channels
nSensors = length(data_input.cfg.info_elec.selected.index4EDF); %selected channels

fsample         = data_input.fsample;
toneDur_inSecs  = str2num(tonedur_text);
nSamplesPerTone = toneDur_inSecs * fsample;

nSamplesPerSeries = length(data_input.trial{1}(1,:));
nSecsPerSeries = nSamplesPerSeries/fsample;
nSecsPerSeriesinTheory = toneDur_inSecs*34;

%Define starting point (i.e., first tone)
t_start_ind(1) = 1;
t_end_ind(1)   = nSamplesPerTone;

%Read out start and stop samples for each tone (except the first)
for i_tone = 2:34 %loop across tones
    t_start_ind(i_tone) = t_end_ind(i_tone-1) + 1; %Start index is end index of previous tone +1
    t_end_ind(i_tone)   = t_start_ind(i_tone) + nSamplesPerTone - 1;%End index/final data point of i_tone is start+number of samples per tone -1
end

%Ensure that, if trial has an offset (i.e., trial-definition begins before
%first tone onset), t_start_ind agrees with onset index of first tone
ind_firsttone = find(data_input.time{1} == 0); %find index for start of first tone
    
t_start_ind = ind_firsttone+(t_start_ind);
t_end_ind = ind_firsttone+(t_end_ind);

%round start and stop samples to get nearest integer
t_start_roundind = round(t_start_ind);
t_end_roundind = round(t_end_ind);


%Note: Due to the rounding, some tones are covered by 102/204 samples, whereas other
%tones are covered by 103/205 samples. However, we need identical number of sample per
%tone to compare the tones. Thus, subtract the last additional sample for tones covered 
%by 103/205 samples
for i_tone = 1:length(t_start_ind)
    roundedsamples_per_tone(i_tone,1) = length(t_start_roundind(i_tone):t_end_roundind(i_tone));
    if roundedsamples_per_tone(i_tone) == 103 || roundedsamples_per_tone(i_tone) == 205
       t_end_roundind(i_tone) = t_end_roundind(i_tone)-1;
    end
    roundedsamples_per_tone_corrected(i_tone,1) = t_end_roundind(i_tone) - t_start_roundind(i_tone)';   
end

[t_start_ind',t_end_ind',t_end_ind'-t_start_ind',...
    t_start_roundind',t_end_roundind',t_end_roundind'-t_start_roundind',...
    roundedsamples_per_tone, roundedsamples_per_tone_corrected] 
%summary pre-rounded vs rounded and distance

%3.2) Initialize and fill 3D data arrays
ECoGdata_forCurrTone = zeros(nSensors, length(t_start_roundind(1):t_end_roundind(1)), nTrials); %create empty prox
ECoGdata_forNextTone = ECoGdata_forCurrTone;
%copy trial wise info (i.e., all samples for selected tone for all channels and all trials into data arrays
for i_trial = 1:nTrials
    
    ECoGdata_forCurrTone(:, :, i_trial) = ...
        data_input.trial{i_trial}(data_input.cfg.info_elec.selected.index4EDF, t_start_roundind(toneIndex) : t_end_roundind(toneIndex));
    ECoGdata_forNextTone(:, :, i_trial) = ...
        data_input.trial{i_trial}(data_input.cfg.info_elec.selected.index4EDF, t_start_roundind(toneIndex+1) : t_end_roundind(toneIndex+1));
 
end

%Check for any NaNs in data struct (due to too short trials)
    if any(find(isnan(ECoGdata_forCurrTone))) == 1
                disp(['NaN entries in p33 struct'])
    end
 
    if any(find(isnan(ECoGdata_forNextTone))) == 1
                disp(['NaN entries in p34 struct'])
    end    

%% 4) Optional baseline correct

%4.1) Get the first tone in the sequence (for optional sequence baseline correction)
% define start and end indeces for the final tone
t_first_end_ind   = round(nSamplesPerTone);
t_first_start_ind = 1;

% initialize 3D array
ECoGdata_forFirstTone = zeros(nSensors, round(nSamplesPerTone), nTrials);

for i_trial = 1:nTrials
    
    ECoGdata_forFirstTone(:, :, i_trial) = ...
        data_input.trial{i_trial}...
        (data_input.cfg.info_elec.selected.index4EDF, t_first_start_ind:t_first_end_ind);
    
end

%4.2 Baseline correct
ind_i = 1;
ind_f = round(.02 * fsample); % correct for first 20 ms

if baseline_correct == 1
    for i_sensor = 1:nSensors
        for i_trial = 1:nTrials           
            baseline = mean( ECoGdata_forCurrTone(i_sensor, ind_i:ind_f, i_trial) );
            ECoGdata_forCurrTone(i_sensor, :, i_trial) = ECoGdata_forCurrTone(i_sensor, :, i_trial) - baseline;
            
            baseline = mean( ECoGdata_forNextTone(i_sensor, ind_i:ind_f, i_trial) );
            ECoGdata_forNextTone(i_sensor, :, i_trial) = ECoGdata_forNextTone(i_sensor, :, i_trial) - baseline;            
            
        end
    end    
    
elseif baseline_correct == -1
    for i_sensor = 1:nSensors
        for i_trial = 1:nTrials           
            baseline = mean( ECoGdata_forFirstTone(i_sensor, ind_i:ind_f, i_trial) );
            ECoGdata_forCurrTone(i_sensor, :, i_trial) = ECoGdata_forCurrTone(i_sensor, :, i_trial) - baseline;
            
            ECoGdata_forNextTone(i_sensor, :, i_trial) = ECoGdata_forNextTone(i_sensor, :, i_trial) - baseline;            
            
        end
    end
    
end

%% 5) Get predicted final tones based on each ks ***
%for each trial, compute/read out 
%1) the tone pitch predicted given the data so far
%2) the p*34 (presumably the same as above, only discretized to full tones)
%3) the p34 (i.e., the really presented last tone)
%4) prediction error (i.e., the distance between presented (p34) and predicted tone 
% ks = [33]; %tone number on which corresponding neural data prediction is based?
for i_k = 1:length(ks)
    
    series_start_ind = toneIndex - ks(i_k) + 1; %beginning of tone sequence part used for prediction
    series_end_ind   = toneIndex; %end of tone sequence part used for prediction
    
    for i_trial = 1:nTrials
        series = data_input.stim.series_f{i_trial}(series_start_ind : series_end_ind); %(non-log) selected sequence tones in Hz
        beta = data_input.stim.beta(i_trial); %always 1.5 in ECoG task
        
        if ks(i_k) == 1 %for first tone
            series_pred(i_k, i_trial) = log( data_input.stim.series_f{i_trial}(toneIndex) ); % predict same tone repeating, because no further info for prediciton
        else
             series_pred(i_k, i_trial) = log( find_predicted_tone(series, beta) ); %subfunction: reads out predicted next tone based on beta level and input tone series
        end
              
        series_pitch(i_trial) = log( data_input.stim.series_f{i_trial}(toneIndex) ); %log(p33)
        series_pitch_next(i_trial) = log( data_input.stim.series_f{i_trial}(toneIndex+1) ); %log(p34)

%         series_pred_error(i_k, i_trial) = log( data_input.stim.series_f{i_trial}(toneIndex+1) ) - log( series_pred(i_k, i_trial) );
        series_pred_error(i_k, i_trial) = series_pitch_next(i_trial) - series_pred(i_k, i_trial); 
        %prediction error as difference between actually presented final
        %tone (log(p34)) and computationally determined final tone
        
    end
    
end

%compare predicted tone via subfunction with sequence-defiend p*34 - should
%be very close, p*34 should be same but only scaled to next full tone
[series_pred', data_input.stim.logf_pred', series_pred'-data_input.stim.logf_pred']


%% 6) Compute avg ERF activity per time window for current tone and next tone

%6.1) Define windows
fsample         = data_input.fsample;

win_size    = 25; %25 data points with fs 512 = ~49ms as mentioned in manuscript
s_per_win = (1/fsample)*win_size;
win_overlap = 0; %as mentioned in manuscript

toneDur_inSecs  = str2num(tonedur_text);
nSamplesPerTone = toneDur_inSecs * fsample; %no rounding

%Define number, start and end sample of window per tone 
windows = [1 win_size];
while windows(end,end) < nSamplesPerTone
    windows = [windows; windows(end,:) + (win_size - win_overlap)];
end

if windows(end,end) > nSamplesPerTone %TJB: windows within single tone
    windows(end,:) = [];
end

windows_inms = (windows / fsample) * 1000;

% nSensors = size(data_input.trial{1},1); %all ECoG sensors
nSensors = size(data_input.cfg.info_elec.selected.Label,1); %selected ECoG sensors

%6.2) Compute avg ERF
for i_win = 1:size(windows,1)

    ind_start = windows(i_win, 1);
    ind_end   = windows(i_win, 2);

    erf_win{i_win} = squeeze( mean( ECoGdata_forCurrTone(:, ind_start:ind_end, :), 2) );
    erf_win_next{i_win} = squeeze( mean( ECoGdata_forNextTone(:, ind_start:ind_end, :), 2) );
    %Output: channel*trial matrix per window with 1 ERF values averaged
    %across all samples per window
end


%% 7) Compute measures of association between ERF windows and prediction / prediction error at each sensor

for i_k = 1:length(ks)
    i_k
    for i_win = 1:size(windows,1)

        for i_sensor = 1:nSensors

            if ks(i_k) > 1 %only if the currently focused tone is not the first tone (since for first tone we can't really predict anything)
            
                %%% for prediction at current tone (i.e., with ERF of current tone (33)) %%%
                dv = erf_win{i_win}(i_sensor,:)';  % mean ERF window at this sensor for all trials
                iv_pred = series_pred(i_k,:)'; %predicted tone pitch for all trials
                iv_pitch = series_pitch'; %tone pitch of focused on tone (33)

                iv_pred_quad = iv_pred - log(440); % center prediction on log(440) for quadratic regression

                % linear regression between ERF (per channel, time window)
                % and predicted tone pitch
    %             stats = regstats(dv, [iv_pred, iv_pitch], 'linear');
                stats = regstats(dv, iv_pred, 'linear', 'tstat');

                pred_t1{i_k}{i_win}(i_sensor) = stats.tstat.t(2);
                pred_t1_stats{i_k}{i_win}{i_sensor} = stats;
                pred_t1_p{i_k}{i_win}(i_sensor) = stats.tstat.pval(2);


                % quadratic regression t-stat
    %             stats = regstats(dv, [iv_pred, iv_pitch], 'purequadratic');
%                 stats = regstats(dv, iv_pred_quad, 'purequadratic', 'tstat');
% 
%                 pred_t2{i_k}{i_win}(i_sensor) = stats.tstat.t(3);
%                 pred_t2_stats{i_k}{i_win}{i_sensor} = stats;
%                 pred_t2_p{i_k}{i_win}(i_sensor) = stats.tstat.pval(3);

            end
            
            
            %%% for prediction error at next tone (i.e., with ERF of next tone (34)) %%%
            dv = erf_win_next{i_win}(i_sensor,:)';  % mean ERF window at this sensor for all trials
%             iv_pred_error = series_pred_error(i_k,:)';
            iv_pred_error = abs( series_pred_error(i_k,:)' );         
            %prediction error as difference between actually presented final
            %tone (log(p34)) and computationally determined final tone

            iv_pitch_next = series_pitch_next';
            
            % linear regression t-stat
            stats = regstats(dv, [iv_pred_error, iv_pitch_next], 'linear', 'tstat');
%             stats = regstats(dv, iv_pred_error, 'linear');
            
            pred_err_t1{i_k}{i_win}(i_sensor) = stats.tstat.t(2);
            pred_err_t1_stats{i_k}{i_win}{i_sensor} = stats;
            pred_err_t1_p{i_k}{i_win}(i_sensor) = stats.tstat.pval(2);


%             % quadratic regression t-stat
%             stats = regstats(dv, [iv_pred_error, iv_pitch_next], 'purequadratic', 'tstat');
% %             stats = regstats(dv, iv_pred_error, 'purequadratic');
%             
%             pred_err_t2{i_k}{i_win}(i_sensor) = stats.tstat.t(4);
%             pred_err_t2_stats{i_k}{i_win}{i_sensor} = stats;
%             pred_err_t2_p{i_k}{i_win}(i_sensor) = stats.tstat.pval(4);            
            
            
        end

    end
end

%% 8) Save variables
mkdir([path_save]);
savefile = [path_save sub '_PredictionERF_regstatsEXP_' BC_text '_' tonedur_title 'msToneDur.mat'];

save(savefile, 'pred_t1', 'pred_t1_stats', 'pred_t1_p', ... %'pred_t2', 'pred_t2_stats', 'pred_t2_p', ...
               'pred_err_t1', 'pred_err_t1_stats', 'pred_err_t1_p', ...
               'ks', 'baseline_correct', 'toneIndex', ...
               'erf_win', 'erf_win_next', ...
               'series_pred', 'series_pitch', 'series_pred_error', 'series_pitch_next');
           
end