function sfa_expt4_HisTrack_CompSaveShuffledKval(sub, baseline_correct, train_odd, tonedur_text, nReps)

%TJB: Adapted from sfa_expt3_HisTrack_CompSaveShuffledKval_SLIDE16, which
%is a new version of tone_history_CV_shuffle0.m (by Jenn)
%For each specified subject,train/test set, tone duration:
%1) Compute null-hypothesis k-values (SLIDE16 = last 16 tones, sliding time window)
%by randomly permuting tone pitch sequences, thus destryoing the temporal
%dependency of the tone pitches
%2) Save output (matrix k-values)

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/data/gogodisk4/thomas/AuditoryPrediction/iEEG/')
% addpath('//gogo.sb.nyumc.org/data/gogodisk4/thomas/AuditoryPrediction/iEEG/')

sfa_expt4_setVars
paths_sfa_expt4 = sfa_expt4_paths(vars.location);

%Add base dir and own script dir
addpath(genpath(paths_sfa_expt4.BaseDir));
addpath(genpath(paths_sfa_expt4.ScriptsDir));

% addpath(genpath('/data/gogodisk4/thomas/AuditoryPrediction/iEEG/Scripts/FromRichard'));

%% 0.2) Determine subject- and condition- specific parameters (whole-recording)

sfa_expt4_subjectinfo %load subject info file (var: si)
subs_PreProcSettings = sfa_expt4_subjectPreProcSettings; %load in file with individual preproc infos

%Define directory for Kprime data
HisTrack_dir = [paths_sfa_expt4.ECoGdata_HisTrack sub '/'];
mkdir(HisTrack_dir);
cd(HisTrack_dir);

%Tone duration condition for data load-in
if strcmp(tonedur_text,'0.2')
    tonedur_title = '200';
elseif  strcmp(tonedur_text,'0.4')
    tonedur_title = '400';
end

loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];%path to indiv preprocessed ECoG data
loadfile_behav = [si.path_behavioraldata_sub];%path to indiv behavioral/stimulus data


if baseline_correct == 1
    path_save = [paths_sfa_expt4.ECoGdata_HisTrack sub '/baseline_corrected/Shuffled/'];
    BC_text = 'BC';

elseif baseline_correct == 0
    path_save = [paths_sfa_expt4.ECoGdata_HisTrack sub '/not_baseline_corrected/Shuffled/'];
    BC_text = 'NBC';
end

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

%% 3) Define what is testing and training data set (even vs. odd) and create respective data structures
%train_odd = 1 %manual proxy

if train_odd
    % define training set as odd trials, test set as even trials
    nTrials = length(data_ECoGfiltref_trials.trial);
    f_test  = mod( 1:nTrials, 2 ) == 0;%filter defining test trials
    f_train = ~f_test;%filter defining train trials
    
else
    % define training set as even trials, test set as odd trials
    nTrials = length(data_ECoGfiltref_trials.trial);
    f_train = mod( 1:nTrials, 2 ) == 0;
    f_test  = ~f_train;
end

if sum(f_test) > sum(f_train)
    f_test(max(find(f_test))) = 0;
    disp('UNEVEN TEST TRAIN SETS, DELETING ONE TRIAL');

elseif sum(f_train) > sum(f_test)
    f_train(max(find(f_train))) = 0;
    disp('UNEVEN TEST TRAIN SETS, DELETING ONE TRIAL');
end

%3.1 Define training set
%Copy field struct and input from selected train trials
%3.1.1 For ECoG data
data_ECoGfiltref_trials_train = data_ECoGfiltref_trials;
data_ECoGfiltref_trials_train.time  = data_ECoGfiltref_trials_train.time(f_train);
data_ECoGfiltref_trials_train.trial = data_ECoGfiltref_trials_train.trial(f_train);
data_ECoGfiltref_trials_train.sampleinfo = data_ECoGfiltref_trials_train.sampleinfo(f_train);
data_ECoGfiltref_trials_train.trialinfo  = data_ECoGfiltref_trials_train.trialinfo(f_train);

%3.1.1 For behavioral data
data_train = data;
numoftrials = num2str(length(data.trialNum));

dat_fields  = fieldnames(data_train);
for i = 1:length(dat_fields)
    if eval(['length(data_train.' dat_fields{i} ') == ' numoftrials])
        eval(['data_train.' dat_fields{i} ' = data_train.' dat_fields{i} '(f_train);']);
    end
end

stim_fields = fieldnames(data_train.stim);
for i = 1:length(stim_fields)
    if eval(['length(data_train.stim.' stim_fields{i} ') == ' numoftrials])
        eval(['data_train.stim.' stim_fields{i} ' = data_train.stim.' stim_fields{i} '(f_train);']);
    end
end

%Copy response/behav data and stim data into MEG struct
data_ECoGfiltref_trials_train.behav = data_train;
data_ECoGfiltref_trials_train.stim  = data_train.stim;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3.2 Define test set
%Copy field struct and input from selected train trials
%3.2.1 For ECoG data
data_ECoGfiltref_trials_test = data_ECoGfiltref_trials;
data_ECoGfiltref_trials_test.time  = data_ECoGfiltref_trials_test.time(f_test);
data_ECoGfiltref_trials_test.trial = data_ECoGfiltref_trials_test.trial(f_test);
data_ECoGfiltref_trials_test.sampleinfo = data_ECoGfiltref_trials_test.sampleinfo(f_test);
data_ECoGfiltref_trials_test.trialinfo  = data_ECoGfiltref_trials_test.trialinfo(f_test);

%3.1.1 For behavioral data
data_test = data;
numoftrials = num2str(length(data.trialNum));

dat_fields  = fieldnames(data_test);
for i = 1:length(dat_fields)
    if eval(['length(data_test.' dat_fields{i} ') == ' numoftrials])
        eval(['data_test.' dat_fields{i} ' = data_test.' dat_fields{i} '(f_test);']);
    end
end

stim_fields = fieldnames(data_test.stim);
for i = 1:length(stim_fields)
    if eval(['length(data_test.stim.' stim_fields{i} ') == ' numoftrials])
        eval(['data_test.stim.' stim_fields{i} ' = data_test.stim.' stim_fields{i} '(f_test);']);
    end
end

data_ECoGfiltref_trials_test.behav = data_test;
data_ECoGfiltref_trials_test.stim  = data_test.stim;


%% Note: Check if even/uneven trial split does not result in big differences in response ratings and hit rates
mean(data_ECoGfiltref_trials_train.behav.resp_prob) %Resp rating (1:5)
mean(data_ECoGfiltref_trials_test.behav.resp_prob)

% mean(data_ECoGfiltref_trials_train.behav.correct_uSeq) %Resp correctnes (-1,0,1)
% mean(data_ECoGfiltref_trials_train.behav.resp_prob)

%% 4) Analyse time window defined MEG data
%4.1 Define analysis parameters
nTrials_train = length(data_ECoGfiltref_trials_train.trial);
nTrials_test  = length(data_ECoGfiltref_trials_test.trial);
% nSensors = size( data_ECoGfiltref_trials_train.trial{1}, 1 ); %all sensors in data set
nSensors = subs_PreProcSettings.(sub).number_ECoGchan; %all ECoG channels

fssample              = data_ECoGfiltref_trials_train.fsample;
toneDur_inSecs  = str2num(tonedur_text);
nSamplesPerTone = toneDur_inSecs * fssample; %no rounding
% nSamplesPerTone = floor(toneDur_inSecs * fssample); %rounding down to get integer
nSamplesPerSeries = length(data_ECoGfiltref_trials_train.trial{1}(1,:));

%Define time windows used for analysis
win_size    = 25; %25 data points with fs 512 = ~49ms as mentioned in manuscript
s_per_win = (1/fssample)*win_size;
win_overlap = 0; %as mentioned in manuscript

%Define starting point (i.e., first tone)
t_start_ind(1) = 1;
t_end_ind(1)   = nSamplesPerTone;

%Read out start and stop samples for each tone (except the first)
for i_tone = 2:34 %loop across tones
    t_start_ind(i_tone) = t_end_ind(i_tone-1) + 1; %Start index is end index of previous tone +1
    t_end_ind(i_tone)   = t_start_ind(i_tone) + nSamplesPerTone - 1;%End index/final data point of i_tone is start+number of samples per tone -1
end

%round start and stop samples to get nearest integer
t_start_ind = round(t_start_ind);
t_end_ind = round(t_end_ind);

%Ensure that, if trial has an offset (i.e., trial-definition begins before
%first tone onset), t_start_ind agrees with onset index of first tone
ind_firsttone = find(data_ECoGfiltref_trials_train.time{1} == 0); %find index for start of first tone
    
t_start_ind = ind_firsttone+(t_start_ind);
t_end_ind = ind_firsttone+(t_end_ind);


%Define number, start and end sample of window per tone 
windows = [1 win_size];
while windows(end,end) < nSamplesPerTone
    windows = [windows; windows(end,:) + (win_size - win_overlap)];
end

if windows(end,end) > nSamplesPerTone %TJB: windows within single tone
    windows(end,:) = [];
end

%Compute window start/end time in ms for each time window
windows_inms = (windows / fssample) * 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.2 Compute across-time-point-averaged MEG data for each sensor, time window, trial, and tone for the training set
for toneIndex = 1:33 %loop across tones
    for i_trial = 1:nTrials_train %loop across training trials
        
        % get MEG data (sensors*time points) for current tone
        ECoG_tone_trial = data_ECoGfiltref_trials_train.trial{i_trial}(:, t_start_ind(toneIndex) : t_end_ind(toneIndex));       
        
        % baseline correct (if chosen above)
        ind_i = 1;
        ind_f = .02 * 600; % correct for first 20 ms

        if baseline_correct == 1 %baseline correction for tone-specific baseline
            for i_sensor = 1:nSensors
               baseline = mean( ECoG_tone_trial(i_sensor, ind_i:ind_f) ); %for each sensor, compute mean value for baseline period (i.e., first 20ms/12 entries)
               ECoG_tone_trial(i_sensor, :) = ECoG_tone_trial(i_sensor, :) - baseline; %subtract baseline from current tone
            end
            
        elseif baseline_correct == -1 %baseline correction for first-tone-of-sequence specific baseline
            meg_first_tone = data_ECoGfiltref_trials_train.trial{i_trial}(:, t_start_ind(1) : t_end_ind(1)); %read out all sensor-wise activity for every timepoint during first tone 
            for i_sensor = 1:nSensors
               baseline = mean( meg_first_tone(i_sensor, ind_i:ind_f) ); %Compute sensor-wise average of first 20ms/12 entries of first tone of sequence
               ECoG_tone_trial(i_sensor, :) = ECoG_tone_trial(i_sensor, :) - baseline; %subtract baseline from current tone
            end
        
        end
                
        % get mean ERF activity for all windows for current tone
        for i_win = 1:size(windows,1)
            ind_start = windows(i_win, 1);
            ind_end   = windows(i_win, 2);

            erf_win_train{i_win}{toneIndex}(:, i_trial) = mean(ECoG_tone_trial(:, ind_start:ind_end),2);
            %TJB: output = Time-point-averaged MEG activity for all channels, for current trial, organized per time window, per tone
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.3 Compute across-time-point-averaged MEG data for each sensor, time window, trial, and tone for the test set
for toneIndex = 1:33
    for i_trial = 1:nTrials_test

        % get MEG data for current tone at current trial
        ECoG_tone_trial = data_ECoGfiltref_trials_test.trial{i_trial}(:, t_start_ind(toneIndex) : t_end_ind(toneIndex));
                
        % baseline correct
        ind_i = 1;
        ind_f = .02 * 600; % correct for first 20 ms

        if baseline_correct == 1
            for i_sensor = 1:nSensors
               baseline = mean( ECoG_tone_trial(i_sensor, ind_i:ind_f) );
               ECoG_tone_trial(i_sensor, :) = ECoG_tone_trial(i_sensor, :) - baseline;
            end
            
        elseif baseline_correct == -1
            meg_first_tone = data_ECoGfiltref_trials_test.trial{i_trial}(:, t_start_ind(1) : t_end_ind(1));
            for i_sensor = 1:nSensors
               baseline = mean( meg_first_tone(i_sensor, ind_i:ind_f) );
               ECoG_tone_trial(i_sensor, :) = ECoG_tone_trial(i_sensor, :) - baseline;
            end
        
        end
        
        
        % get mean ERF activity for all windows
        for i_win = 1:size(windows,1)
            ind_start = windows(i_win, 1);
            ind_end   = windows(i_win, 2);

            erf_win_test{i_win}{toneIndex}(:, i_trial) = mean( ECoG_tone_trial(:, ind_start:ind_end), 2 );
            %TJB: output = Time-point-averaged MEG activity for all channels, for current trial, organized per time window, per tone
        end        
    end
end


%% 5) Read out previously-presented tone pitches (i.e., sequence history)

%read out presented tone pitches for training and test set from stim data
for i_trial = 1:nTrials_train
    tone_pitch_train_o(i_trial, :) = log( data_ECoGfiltref_trials_train.stim.series_f{i_trial} );
    %output: trial*tone matrix containing tone pitches (freq) of all tones per sequence
end

for i_trial = 1:nTrials_test
    tone_pitch_test_o(i_trial, :)  = log( data_ECoGfiltref_trials_test.stim.series_f{i_trial} );
end

%SLIDE_16 modifier (take respective tone pitch (no diff)):
tone_pitch_train = tone_pitch_train_o(:, [1:end-1]); %take original matrix, delete last column/tone, because we only use first 33
tone_pitch_test = tone_pitch_test_o(:, [1:end-1]);

%% 6) Shuffle temporal order of tone pitches (tone 1:32) in training set

tone_pitch_train_orig = tone_pitch_train; %store original tone pitch order

for i_rep = 1:nReps %loop for shuffle repetitions
    disp(['--- Current Rep# = ' num2str(i_rep) ' ---'])               
    
    clear tone_pitch_train %delete var for each Rep-run
    for i_trial = 1:size(tone_pitch_train_orig,1)
        t = tone_pitch_train_orig(i_trial, 1:32); %take the first 32 tones and save them separately
        %%Critical Step: Shuffling of original tone order in sequence - dissolves temporal dependencies
        t = Shuffle(t); %Shuffle.m function fro Psychtoolbox3
        
        tone_pitch_train(i_trial,:) = t;
    end

%% 7) Compute k values based on (shuffled) traning set
%compute k as linear regression (across trials) between sensor-wise, window-wise ERF and 
%tone pitch for tones 32:16 (17 tones total)
    ks = 1:16; %TJB: 1:16 regression terms/k-values/models - all potential k vals

    for i_k = 1:length(ks) %index/current k val (i.e., amount of k models)
        i_k
        for i_win = 1:size(windows,1)
            for i_sensor = 1:nSensors            

                clear dv iv_cell
                dv = [];
                for toneIndex = 16:32 %17 entries
                    dv = [dv; erf_win_train{i_win}{toneIndex}(i_sensor, :)']; %Column vector with MEG activity for time window at current tone
                    %Output: Column with sensor-specific, time-window averaged MEG activity/ERF for all trials for specific window and toneIndex
                    %TJB: Column vector (all trials for specific previous tones, then append all trials for next tone): 

                    %Ad hoc change (shuffled) tone_pitch_train so that current tone pitch,
                    %rather than being shuffled, corresponds to true/unshuffled current tone pitch
                    %Ouput: Matrix with trials*tone pitches (1:32) with shuffled tones for all tones except current tone
                    tone_pitch_train_adj = tone_pitch_train; %copy shuffled tone pitch sequence into new var
                    tone_pitch_train_adj(toneIndex) = tone_pitch_train_orig(toneIndex); %replace current tone pitch with true/unshuffled current tone pitch

                    if toneIndex == 16 %for tone index 16 (i.e, first index entry), start new row
                        for kk = 1:ks(i_k)%loop over k values from 1 to current k val
                            iv_cell{kk} = tone_pitch_train_adj(:, toneIndex - kk + 1); %take tone pitch for all trials for curent tone and also for kk tones back 
                            %TJB: Output: column with pitch values from all trials for tone number X (appended by trials), 
                            %cell means amount of tones going back from kk= 1 -> 0 back (only target tone) to kk = 16 -> 15 back (target + 15 previous tones) 
                        end

                    else %for all other tone indices, append above matrix
                        for kk = 1:ks(i_k)
                            iv_cell{kk} = [iv_cell{kk}; tone_pitch_train_adj(:, toneIndex - kk + 1)];
                        end
                    end
                end

                iv = [];
                for kk = 1:ks(i_k)
                    iv(:, kk) = iv_cell{kk}; %conversion from cell to matrix
                    %Output: matrix with tone pitch frequencies;
                    %rows = appended trials per tone (tone 16 to 32),
                    %columns = model order (i.e., tones back (0 to 15 previous
                    %tones)
                end


                % linear regression
                %i.e., for each i_k round/model order, an additional previously presented tone is added 
                %and linear regression is computed
                ss = regstats(dv, iv, 'linear',  {'tstat', 'fstat', 'r', 'rsquare', 'adjrsquare'});
                %beta = regression coefficient = degree of change in outcome var for every 1-unit change in predictor var

                r = ss.r; %residuals
                n = length(r);
                sig2_train = sum(r.^2) / n; %sum of squared residuals divided by length of model order
                K    = length(ss.tstat.beta) + 1; %number of regression coefficients +1

                R2_train = ss.rsquare;

                ss = rmfield(ss, 'r');

                ss.n    = n;
                ss.K    = K;

                ss.R2_train   = R2_train;
                ss.sig2_train = sig2_train;
                ss.AIC_train  = n*log(sig2_train) + 2*K; 
                %AIC = Akaike information criterion - esitmator of relative quality of statistical models for given set of data; means for model selection
                ss.BIC_train  = n*log(sig2_train) + K*log(n);
                %BIC = Bayesian information criterion - criterion for model selection
                %Both model selection criteria introduce a penalty term for the
                %numbers of parameters of a model, thus preventing overfitting

                stats_linear{i_win, i_k, i_sensor} = ss;
                %TJB: Output:cell with linear regression stats
                %dimensions = time window*model order*sensors cell 

                clear ss
            end

        end
    end


    %% 7) Characterize model fit on test set
    %7.1 compute k as linear regression (across trials) between sensor-wise, window-wise ERF and 
    %tone pitch for tones 32:16 (17 tones total), similar to above
    for i_k = 1:length(ks)
        i_k
        for i_win = 1:size(windows,1)
            for i_sensor = 1:nSensors


                clear dv iv_cell
                dv = [];
                for toneIndex = 16:32
                    dv = [dv; erf_win_test{i_win}{toneIndex}(i_sensor, :)'];

                    if toneIndex == 16
                        for kk = 1:ks(i_k)
                            iv_cell{kk} = tone_pitch_test(:, toneIndex - kk +1);
                        end

                    else
                        for kk = 1:ks(i_k)
                            iv_cell{kk} = [iv_cell{kk}; tone_pitch_test(:, toneIndex - kk +1)];
                        end
                    end
                end

                iv = [];
                for kk = 1:ks(i_k)
                    iv(:, kk) = iv_cell{kk};
                end


                %7.2 Find error of training set on test set
                %Take stats from training set (for current time window, model order, and sensor)
                %and apply them to test set
                ss = stats_linear{i_win, i_k, i_sensor};

                beta = ss.tstat.beta';

                try
                for jj = 1:length(dv) %for each trial across all tones from 16:32
                    %compute estimated ERF for each trial across all tones from 16:32
                    %based on the beta weights*tone pitches
                    dv_est(jj) = beta(1) + sum( beta(2:end) .* iv(jj,:) );
                end

                catch
                    keyboard
                end
                r = dv - dv_est';%compute residuals/difference between estimated and real test set ERF activity for each trial across all tones from 16:32

                SS_res = sum(r.^2);%Compute summed (across trials) square of residual
                SS_tot = sum((dv - mean(dv)).^2);%Compute summed square of single trial across-trial average ERF

                R2_test = 1 - SS_res / SS_tot;   %r-square for test set based on estimated data     


                n = length(r);
                sig2_test = sum(r.^2) / n; %%Note: This is the critical k-value by which we later define the k-prime value
                K    = length(ss.tstat.beta) + 1;
                
                sig2_test_shuffle{i_win, i_k, i_sensor}(i_rep) = sig2_test;
                %output: matrix with sum of squared residuals divided by amount of tested models
                %per time window, model order, sensor, repetition                
 
                sig2_test_shuffle{i_win, i_k, i_sensor}(i_rep) = sig2_test;

                clear ss
            end

        end
    end
end %end Reps

%% 8) Save output variables
%8.1 Define Output labels depending on input
%Training/Test set definition
if train_odd
    oddity_text = 'trainOdd';
else
    oddity_text = 'trainEven';
end

%8.2 save output var
savefile = [path_save sub '_HisTrackSLIDE16_regstatsSHUFFLED' num2str(nReps) 'Reps_' BC_text '_' tonedur_title 'msToneDur_' oddity_text '.mat'];
save(savefile, 'sig2_test_shuffle','-v7.3');           

end