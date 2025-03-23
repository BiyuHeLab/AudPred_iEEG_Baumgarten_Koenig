function NASTD_ECoG_HisTrack_CompShuffKvals...
    (sub, input_ToneDurLabel, input_DataLabel, ...
    input_Data, ... 
    SamplesTW, Label_TW, NumFolds, nReps,...
    paths_NASTD_ECoG, ...
    plotFig_CompFolds)

%For each specified subject, train/test set, tone duration, electrode:
%1) Compute shuffled k-values for both train and test set (SLIDE16 = last 16 tones, sliding time window)
%2) Combine shuffled Kprime values across folds for the shuffled condition
%by averaging each rep across folds to receive null distribution of 100
%K-prime values per time window/sensor

%Additions:
%-Read out regression weights (beta weights) for winning k-model/k-prime
%-Plot K-prime and MeanSquareResiduals (parameter on which Kprime selection
%is based) for both sets to allow for comparison across sets

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));

%% 0.2) Determine subject-specific parameters (whole-recording)
%Tone duration condition for data load-in
if strcmp(input_ToneDurLabel,'0.2')
    tonedur_title = '200msTD';
elseif  strcmp(input_ToneDurLabel,'0.4')
    tonedur_title = '400msTD';
end

%Define directory for Kprime data
path_save = [paths_NASTD_ECoG.ECoGdata_HisTrack sub '/Shuff/' Label_TW '/'];
if (~exist(path_save, 'dir')); mkdir(path_save); end

%Define directory for potential figures
path_fig_Kprime = [paths_NASTD_ECoG.ECoGdata_HisTrack sub '/Shuff/' Label_TW '/Figs/Kprime_Folds/'];
if (~exist(path_fig_Kprime, 'dir')); mkdir(path_fig_Kprime); end
path_fig_MSR = [paths_NASTD_ECoG.ECoGdata_HisTrack sub '/Shuff/' Label_TW '/Figs/MinSquarRes_Folds/'];
if (~exist(path_fig_MSR, 'dir')); mkdir(path_fig_MSR); end

%% 1) Select input data 
data_input          = input_Data;
data_input.trial    = input_Data.(input_DataLabel);
data_input.stim     = input_Data.behav.stim;

disp(['Total trial number for TD ' input_ToneDurLabel 's: ' ...
    num2str(length(data_input.trial))])

clear preprocData_perTD

%% 2) Define TP and samples of Tone Presentation
%2.1) Define time window parameters depending on TD
NumSensors          = size(data_input.trial{1}, 1 );
NumTrials           = length(data_input.trial);
fsample             = data_input.fsample;
toneDur_inSecs      = str2num(input_ToneDurLabel);
NumSamplesPerTone   = toneDur_inSecs * fsample;

%2.2 Determine TP/samples for each tone start+end
TP_Tone_StartStop = NaN(35,2);
TP_Tone_StartStop(1,1) = 0;
%Determine Tone Start as TP in recording that is clostest to 'theoretical'
%tone start (doesn't always agree due to sample freq)
for i_tone = 2:35
    Dist = abs(data_input.time{1} - ...
        ((str2num(input_ToneDurLabel)*i_tone)-str2num(input_ToneDurLabel)));
    minDist = min(Dist);
    i_minDist = find(Dist == minDist);
    TP_Tone_StartStop(i_tone,1) = data_input.time{1}(i_minDist);
end
%Determine Tone End by taking start of following tone
for i_tone = 1:34
    i_LastSampleTone = ...
        find(data_input.time{1} == TP_Tone_StartStop(i_tone+1,1));
    TP_Tone_StartStop(i_tone,2) = data_input.time{1}(i_LastSampleTone);
end
TP_Tone_StartStop = TP_Tone_StartStop(1:34,:);

%check if all tones are of equal length, if not then choose min length by deleting the last sample
minSeqLength_sec = min(TP_Tone_StartStop(:,2)-TP_Tone_StartStop(:,1));
for i_tone = 1:34
    if TP_Tone_StartStop(i_tone,2)-TP_Tone_StartStop(i_tone,1) > ...
            minSeqLength_sec
        TP_Tone_StartStop(i_tone,2) = ...
            data_input.time{1}(find(data_input.time{1} == ...
            TP_Tone_StartStop(i_tone,2))-1);
    end
end

%Read out samples corresponding to TPs
Sample_Tone_StartStop = NaN(34,2);
for i_tone = 1:34
    Sample_Tone_StartStop(i_tone,1) = ...
        find(data_input.time{1} == TP_Tone_StartStop(i_tone,1));
    Sample_Tone_StartStop(i_tone,2) = ...
        find(TP_Tone_StartStop(i_tone,2) == data_input.time{1});
end

%% 3) Determine Folds (i.e., trial-to-fold allocation)
%3.1. Determine random trial selection%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3.1.1 Determine number of folds based on trial number
    %(optimally = 20, should be >10, otherwise reduce number of folds)
while NumTrials/NumFolds < 10 
	NumFolds = NumFolds-1;
	disp('-- Number of folds changed due to insufficient number of reps per uSeq --')                       
end
disp([num2str(NumFolds) ' folds used for cross-validation'])

%3.1.2 Determine how many trials are going to be placed in each fold
AvgTrials_perFold = NumTrials/NumFolds;
Remaining_trials = NumTrials;
for i_fold = 1:NumFolds
    if Remaining_trials >= Remaining_trials - (Remaining_trials - ceil(AvgTrials_perFold))
        Trials_perFold(1, i_fold) = Remaining_trials - (Remaining_trials - ceil(AvgTrials_perFold));
        Remaining_trials =  Remaining_trials - ceil(AvgTrials_perFold);
    else
        Trials_perFold(1, i_fold) = Remaining_trials;
    end
end
disp([num2str(Trials_perFold) ' trials per fold used for cross-validation'])

%3.1.2 Shuffle trial order
trial_index_orig = 1:NumTrials;
trial_index_shuffled = Shuffle(trial_index_orig,1);
                 
%3.1.3 Allocate shuffled trials to folds
counter = 1;
for i_fold = 1:NumFolds
    index_Trials_perFold{i_fold} = ...
        trial_index_shuffled(counter:(counter + Trials_perFold(i_fold) -1));
    counter = counter + Trials_perFold(i_fold);
end
    
%3.1.4 Store selected reps in summary file
for i_fold = 1:NumFolds
    Info_TrialSelection.Number_Usedreps_perFold{i_fold} = ...
        length(index_Trials_perFold{i_fold});
    Info_TrialSelection.TrialIndex_ToneDurDataSet_UsedTrialsperFold{i_fold} = ...
        index_Trials_perFold{i_fold};
end
    
%% 4) Loop subsequent computation over folds
for i_fold = 1:NumFolds
    
    disp(['- Computing fold: ' num2str(i_fold) '/' num2str(NumFolds) '-'])
    
    %4.1 for each fold_loop, specify train and test data
    index_train = []; %train data = all other folds
    index_test = [];
    
    for i_fold2 = 1:NumFolds
        index_train = [index_train index_Trials_perFold{i_fold2}];
    end
    for i_trial = 1:length(index_Trials_perFold{i_fold})        
        index_train(find(index_train == index_Trials_perFold{i_fold}(i_trial))) = [];
    end
    index_test = index_Trials_perFold{i_fold}; %test data = current fold
       
    %4.2 construct corresponding filters   
    f_train = [1:length(data_input.trial)]*0; %create empty filter-proxies 
    f_test = [1:length(data_input.trial)]*0;

    f_train(index_train) = 1; %fill set-specific chosen trial indices
    f_test(index_test) = 1;

    f_train = logical(f_train); %convert arrays to logical
    f_test = logical(f_test);


%% 5) Define training set
    %Copy field struct and input from selected train trials
    data_train              = data_input;
    data_train.time         = data_train.time(f_train);
    data_train.sampleinfo   = data_train.sampleinfo(f_train);
    data_train.trialinfo    = data_train.trialinfo(f_train);
    data_train.trial        = data_train.trial(f_train);
    
    data_behav_fields  = fieldnames(data_train.behav);
    for i = 1:length(data_behav_fields)
        if eval(['length(data_train.behav.' data_behav_fields{i} ')']) == NumTrials
            eval(['data_train.behav.' data_behav_fields{i} ...
                ' = data_train.behav.' data_behav_fields{i} '(f_train);']);
        end
    end

    data_stim_fields = fieldnames(data_train.stim);
    for i = 1:length(data_stim_fields)
        if eval(['length(data_train.stim.' data_stim_fields{i} ')']) == NumTrials
            eval(['data_train.stim.' data_stim_fields{i} ...
                ' = data_train.stim.' data_stim_fields{i} '(f_train);']);
        end
    end

%% 6) Define test set
    %Copy field struct and input from selected train trials
    data_test               = data_input;
    data_test.time          = data_test.time(f_test);
    data_test.sampleinfo    = data_test.sampleinfo(f_test);
    data_test.trialinfo     = data_test.trialinfo(f_test);
    data_test.trial         = data_test.trial(f_test);

    data_behav_fields  = fieldnames(data_test.behav);
    for i = 1:length(data_behav_fields)
        if eval(['length(data_test.behav.' data_behav_fields{i} ')']) == NumTrials
            eval(['data_test.behav.' data_behav_fields{i} ...
                ' = data_test.behav.' data_behav_fields{i} '(f_test);']);
        end
    end

    data_stim_fields = fieldnames(data_test.stim);
    for i = 1:length(data_stim_fields)
        if eval(['length(data_test.stim.' data_stim_fields{i} ')']) == NumTrials
            eval(['data_test.stim.' data_stim_fields{i} ...
                ' = data_test.stim.' data_stim_fields{i} '(f_test);']);
        end
    end

%% 7) Analyse time window defined MEG data
    %7.1 Define analysis parameters
    NumTrials_train = length(data_train.trial);
    NumTrials_test  = length(data_test.trial);
%     NumSamplesPerSeries = length(data_train.trial{1}(1,:));

    %Define time windows used for analysis
    win_size    = SamplesTW; %25 data points with fs 512 = 50ms as mentioned in manuscript
    win_overlap = 0; %as mentioned in manuscript

    %Define number, start and end sample of window per tone 
    windows = [1 win_size];
    while windows(end,end) < NumSamplesPerTone
        windows = [windows; windows(end,:) + (win_size - win_overlap)];
    end

    if windows(end,end) > NumSamplesPerTone %TJB: windows within single tone
        windows(end,:) = [];
    end

    %Compute window start/end time in ms for each time window
    windows_inms = (windows / fsample) * 1000;
    
    %Determine number of time windows
    NumWindows = size(windows,1);
    
%% 8) Compute within-window-averaged ECoG data for each electrode, time window, trial, and tone 
%8.1 for the training set
    for toneIndex = 1:33 %loop across tones
        for i_trial = 1:NumTrials_train %loop across training trials

            %Get ECoG data (sensors*time points) for current tone
            ECoGdata_perToneTrial = ...
                data_train.trial{i_trial}...
                (:, Sample_Tone_StartStop(toneIndex,1) : ...
                Sample_Tone_StartStop(toneIndex,2));       

            %Average ECoG data within each window for current tone
            for i_win = 1:size(windows,1)
                ind_start = windows(i_win, 1);
                ind_end   = windows(i_win, 2);

                ECoGdata_Avgwin_train{i_win}{toneIndex}(:, i_trial) = ...
                    mean(ECoGdata_perToneTrial(:, ind_start:ind_end),2);
                %output = Time-point-averaged  activity for 
                %all channels, for current trial, organized per time window, per tone
            end
        end
    end

%8.2 for the test set
    for toneIndex = 1:33
        for i_trial = 1:NumTrials_test

            % get MEG data for current tone at current trial
            ECoGdata_perToneTrial = ...
                data_test.trial{i_trial}...
                (:, Sample_Tone_StartStop(toneIndex,1) : ...
                Sample_Tone_StartStop(toneIndex,2));


            % get mean ERF activity for all windows
            for i_win = 1:size(windows,1)
                ind_start = windows(i_win, 1);
                ind_end   = windows(i_win, 2);

                ECoGdata_Avgwin_test{i_win}{toneIndex}(:, i_trial) = ...
                    mean(ECoGdata_perToneTrial(:, ind_start:ind_end), 2 );
            end        
        end
    end

%% 9) Read out previously-presented tone pitches (i.e., sequence history)
    %read out presented tone pitches for training and test set from stim data
    for i_trial = 1:NumTrials_train
        TonePitch_train_unshuffled(i_trial, :) = ...
            log(data_train.stim.series_f{i_trial});
        %output: trial*tone matrix containing tone pitches (freq) of all tones per sequence
    end

    for i_trial = 1:NumTrials_test
        TonePitch_test_unshuffled(i_trial, :)  = ...
            log(data_test.stim.series_f{i_trial});
    end

    %take original matrix, delete last 2 tones, because we only use first 32 tones
    TonePitch_train_unshuffled = TonePitch_train_unshuffled(:, [1:end-2]); 
    TonePitch_test_unshuffled = TonePitch_test_unshuffled(:, [1:end-2]);

%% 10) Shuffle tone order (1:32) in training set
    %Shuffle details:
	%1) Shuffle tone pitch order anew for each trial
    %2) Shuffle anew for each repetition
	%3) Shuffle tone pitch order anew for every fold      

    for i_rep = 1:nReps %loop for shuffle repetitions
        
        disp(['--- Current Rep# = ' num2str(i_rep) ' ---'])
        
        clear tone_pitch_train %empty var for each Rep-run
        
        for i_trial = 1:size(TonePitch_train_unshuffled,1)
            [ShuffTonePitch, Shuffle_index] = ...
                Shuffle(TonePitch_train_unshuffled(i_trial,:)); %Shuffle.m function fro Psychtoolbox3
            TonePitch_train_shuffled(i_trial,:) = ShuffTonePitch;
            Info_TrialSelection.ShuffleIndex_perTrial{i_fold}{i_rep}(i_trial,:) = ...
                Shuffle_index; %Save ShuffleIndex (i.e., shows original position of tone)
            Shuffle_index = [];
        end
        
        %% 11) Compute k values based on traning set
        %compute k as linear regression (across trials) between ...
        %1)sensor-wise, window-wise ERF and
        %2) tone pitch for tones 32:17 (16 tones total)
        kModels = 1:16; %16 regression terms/k-values/models - all potential k vals
        
        for i_kModel = 1:length(kModels) %index/current k val (i.e., amount of k models)
%             disp(['ToneDur: ' tonedur_title ' - Training Set (Fold: ' ...
%                 num2str(i_fold) ') - Current Rep: ' num2str(i_rep) ...
%                 ' - Current k-model: ' num2str(i_kModel)])
            for i_win = 1:size(windows,1)
                for i_sensor = 1:NumSensors
                    
                    clear dv iv_cell
                    dv = [];
                    for toneIndex = 16:32
                        dv = [dv; ECoGdata_Avgwin_train{i_win}{toneIndex}(i_sensor, :)'];
                        %DV: Neural activity for every trial of tone X in training set
                        %(sensor-specific, TW-specific and averaged across TW), then appended
                        %across tones - i.e., #trials(1:40=40)*#tones(16:32=17)
                        
                        
                        if toneIndex == 16 %for tone index 16 (i.e, first index entry), start new row
                            for kk = 1:kModels(i_kModel)%add tone depending on Kmodel
                                iv_cell{kk} = TonePitch_train_shuffled(:, toneIndex - kk + 1);
                                %IV: Tone pitch of tone X for every trial, plus X-k previous tone pitches for every trial
                                %(depending on Kmodel) - i.e., first only current tone, then current tone + current tone -1, ...
                                %cell means amount of tones going back from kk= 1 -> 0 back (only target tone) to kk = 16 -> 15 back (target + 15 previous tones)
                            end
                            
                        else %for all other tone indices, append above matrix
                            for kk = 1:kModels(i_kModel)
                                iv_cell{kk} = [iv_cell{kk}; TonePitch_train_shuffled(:, toneIndex - kk + 1)];
                            end
                        end
                    end
                    %Note: Same kernel/beta weight is assumed for all tones.
                    %This is why we append trials from all tones, i.e., it
                    %increases our trial number and stat. powerm but relies on
                    %this assumption.
                    
                    
                    iv = [];
                    for kk = 1:kModels(i_kModel)
                        iv(:, kk) = iv_cell{kk}; %conversion from cell to matrix
                        %Output: matrix with tone pitch frequencies;
                        %rows = appended trials per tone (tone 16 to 32),
                        %columns = model order (i.e., tones back (0 to 15 previous
                        %tones)
                    end
                    
                    %Compute  linear regression between neural data and tone pitch
                    %(i.e., for each i_k round/model order, an additional previously presented tone
                    %is added and linear regression is computed
                    stats_LinReg_Train = ...
                        regstats(dv, iv, 'linear',  ...
                        {'tstat', 'fstat', 'r', 'rsquare', 'adjrsquare'});
                    
                    %Read out residuals from linear regression
                    Residuals_Train = stats_LinReg_Train.r; %residuals
                    NumResiduals_Train = length(Residuals_Train);
                    stats_LinReg_Train = rmfield(stats_LinReg_Train, 'r');
                    stats_LinReg_Train.NumResiduals    = NumResiduals_Train;
                    
                    %Compute K-prime defining metric = sum of squared residuals divided by length of model order
                    SumSquaredResiduals_Train = sum(Residuals_Train.^2) / NumResiduals_Train;
                    stats_LinReg_Train.SumSquaredResiduals_TrainOnly = SumSquaredResiduals_Train;
                    
                    %Copy statistical output into format holding loop-vars
                    TempStats_LinReg4Kprime{i_win, i_kModel, i_sensor} = stats_LinReg_Train;
                    
                    %Cleanup (everything except output var TempStats_LinReg4Kprime)
                    clear dv iv_cell iv stats_LinReg_Train Residuals_Train NumResiduals_Train SumSquaredResiduals_Train
                end
                
            end
        end        
        
        %% 11) Characterize model fit on test set
        %compute k as linear regression (across trials) between sensor-wise, window-wise ERF and
        %tone pitch for tones 32:16 (17 tones total), similar to above
        for i_kModel = 1:length(kModels)
%             disp(['ToneDur: ' tonedur_title ' - Test Set (Fold: ' ...
%                 num2str(i_fold) ') - Current Rep: ' num2str(i_rep) ...
%                 ' - Current k-model: ' num2str(i_kModel)])
            for i_win = 1:size(windows,1)
                for i_sensor = 1:NumSensors
                    
                    %11.1 Setup DV and IV
                    clear dv iv_cell
                    dv = [];
                    for toneIndex = 16:32
                        dv = ...
                            [dv; ECoGdata_Avgwin_test{i_win}{toneIndex}(i_sensor, :)'];
                        
                        if toneIndex == 16
                            for kk = 1:kModels(i_kModel)
                                iv_cell{kk} = ...
                                    TonePitch_test_unshuffled(:, toneIndex - kk +1);
                            end
                            
                        else
                            for kk = 1:kModels(i_kModel)
                                iv_cell{kk} = ...
                                    [iv_cell{kk}; ...
                                    TonePitch_test_unshuffled(:, toneIndex - kk +1)];
                            end
                        end
                    end
                    
                    iv = [];
                    for kk = 1:kModels(i_kModel)
                        iv(:, kk) = iv_cell{kk};
                    end
                    
                    %11.2 Find error of training set on test set
                    %Take beta weight from training set (for current time window, model order, and sensor) and apply them to test set
                    stats_LinReg_TestvsTrain = ...
                        TempStats_LinReg4Kprime{i_win, i_kModel, i_sensor};
                    beta = stats_LinReg_TestvsTrain.tstat.beta';
                    %beta = regression coefficient = degree of change in outcome var for every 1-unit change in predictor var
                    
                    for jj = 1:length(dv) %for each trial across all tones from 16:32
                        %compute estimated ERF for each trial across all tones from 16:32
                        %based on the beta weights*tone pitches
                        dv_est(jj) = beta(1) + sum( beta(2:end) .* iv(jj,:) );
                    end
                    
                    Residuals_TestvsTrain = dv - dv_est';
                    %compute residuals/difference between estimated and 
                    %real test set ERF activity for each trial across all tones from 16:32
                    NumResiduals_TestvsTrain = length(Residuals_TestvsTrain);
                    
                    SumSquaredResiduals_TestvsTrain = ...
                        sum(Residuals_TestvsTrain.^2) / ...
                        NumResiduals_TestvsTrain; 
                    %%Note: This is the critical k-value by which we later define the k-prime value
                    stats_LinReg_TestvsTrain.SumSquaredResiduals_TestvsTrain ...
                        = SumSquaredResiduals_TestvsTrain;
                    %%Note: This is the critical k-value by which we later define the k-prime value
                    
                    %11.3 linear regression between neural data and tone pitch  for test set
                    stats_LinReg_Test = ...
                        regstats(dv, iv, 'linear', ...
                        {'tstat', 'fstat', 'r', 'rsquare', 'adjrsquare'});
                    
                    Residuals_Test = stats_LinReg_Test.r;
                    NumResiduals_Test = length(Residuals_Test);
                    SumSquaredResiduals_Test = ...
                        sum(Residuals_Test.^2) / NumResiduals_Test;
                    
                    stats_LinReg_TestvsTrain.SumSquaredResiduals_TestOnly = ...
                        SumSquaredResiduals_Test;
                    
                    %Copy statistical output into format holding loop-vars                    
                    stats_LinReg4Kprime{i_fold}{i_win, i_kModel, i_sensor}.SumSquaredResiduals_TrainOnly(i_rep) = ...
                        stats_LinReg_TestvsTrain.SumSquaredResiduals_TrainOnly; %Add test stats to common stats output var
                    stats_LinReg4Kprime{i_fold}{i_win, i_kModel, i_sensor}.SumSquaredResiduals_TestOnly(i_rep) = ...
                        stats_LinReg_TestvsTrain.SumSquaredResiduals_TestOnly; %Add test stats to common stats output var
                    stats_LinReg4Kprime{i_fold}{i_win, i_kModel, i_sensor}.SumSquaredResiduals_TestvsTrain(i_rep) = ...
                        stats_LinReg_TestvsTrain.SumSquaredResiduals_TestvsTrain; %Add test stats to common stats output var
%                     stats_LinReg4Kprime{i_fold}{i_win, i_kModel, i_sensor}.tstat_train{i_rep} = ...
%                         stats_LinReg_TestvsTrain.tstat;
                    
                    %Cleanup (everything except output var stats_LinReg4Kprime)
                    clear dv iv_cell iv ...
                        stats_LinReg_TestvsTrain beta dv_est ...
                        Residuals_TestvsTrain NumResiduals_TestvsTrain SumSquaredResiduals_TestvsTrain ...
                        Residuals_Test NumResiduals_Test SumSquaredResiduals_Test
                    
                end                
            end
        end
    end
    %Cleanup before each new fold
	clear TonePitch_* ECoGdata_Avgwin_train ECoGdata_Avgwin_test TempStats_LinReg4Kprime
end

% %% 12) Save output variables
% %12.1 save output var
% savefile = [path_save sub '_' input_DataLabel '_' tonedur_title ...
%     '_HisTrack_regstats_ShuffKperfold_unbalanced.mat'];
% save(savefile, 'stats_LinReg4Kprime', 'Info_TrialSelection','-v7.3');  

%% 13) Determine Kprime as 'k-model with minimal sum squared residuals
%13.1 Restructure regstats data and read out relevant stats for model selection from each fold
NumModelOrder = size(stats_LinReg4Kprime{1},2); %number model order;
kModel_list = 1:NumModelOrder;
NumReps = size(stats_LinReg4Kprime{1}{1}.SumSquaredResiduals_TestvsTrain,2);

for i_fold = 1:NumFolds
    for i_win = 1:NumWindows
        for i_sensor = 1:NumSensors
            for i_rep = 1:NumReps
                
                for i_kModel = 1:NumModelOrder
                    SumSquaredResiduals_proxy(i_kModel) = ...
                        stats_LinReg4Kprime{i_fold}...
                        {i_win, i_kModel, i_sensor}.SumSquaredResiduals_TestvsTrain(i_rep);
                    %Criteria for Kprime model selection = sum of squared residuals 
                    %divided by length of model order for each k (per sensor,time win/fold)
                end
                
                %13.2 Read out minimum squared residuals divided by length of model order
                %and the respective position in array per sensor/time win/fold from all folds in common struct
                %(respective position in array = 1:16 = current tone (1) and up to 15 tones back)
                [minSumSquaredResiduals, index_minSumSquaredResiduals] = ...
                    min( SumSquaredResiduals_proxy );
                
                %13.3 Determine Kprime value (i.e., preferred model order/Number of 
                %tones back in sequence history best explaining MEG data)
                %based on index where inimum squared residuals have minimal values
                %and correct value by -1, so that it ranges from 0 (only current tone) to 15 (i.e., current tone + 15 previous tones)
                Kprime_perFoldRep{i_win}{i_sensor}(i_fold, i_rep) = ...
                    kModel_list(index_minSumSquaredResiduals) - 1;
                SumSquaredResiduals_Kprime_perFoldRep{i_win}{i_sensor}(i_fold, i_rep) = ...
                    minSumSquaredResiduals;
                %minimum squared residuals divided by length of model order of K-value (i.e., preferred model order/Number of tones back)
                
                SumSquaredResiduals_perK_perFoldRep{i_win}{i_sensor}(i_fold, i_rep, :) = ...
                    SumSquaredResiduals_proxy;
                %sum of squared residuals divided by length of model order, stored in common matrix across folds
                
                %13.4 Clean-up
                clear SumSquaredResiduals_proxy minSumSquaredResiduals index_minSumSquaredResiduals
                
            end
        end
    end
end

%% 14) Create Kprime null distribution by averaging Shuff Kprime across folds for each rep

%14.1 Compute Average and STD across folds (for each sensor/time win)
for i_win = 1:NumWindows
    for i_sensor = 1:NumSensors
        %Kprime
        Kprime_perRepavgFold{i_win}(i_sensor,:) = ...
            mean(Kprime_perFoldRep{i_win}{i_sensor},1);
        Kprime_perRepSDFold{i_win}(i_sensor,:) = ...
            std(Kprime_perFoldRep{i_win}{i_sensor});
        %SSR Kprime
        SumSquaredResiduals_Kprime_perRepavgFold{i_win}(i_sensor,:) = ...
            mean(SumSquaredResiduals_Kprime_perFoldRep{i_win}{i_sensor},1);
        SumSquaredResiduals_Kprime_perRepSDFold{i_win}(i_sensor,:) = ...
            std(SumSquaredResiduals_Kprime_perFoldRep{i_win}{i_sensor});
        %SSR all K
        SumSquaredResiduals_perK_perRepavgFold{i_win}{i_sensor} = ...
            squeeze(mean(SumSquaredResiduals_perK_perFoldRep{i_win}{i_sensor},1));
        SumSquaredResiduals_perK_perRepSDFold{i_win}{i_sensor} = ...
            squeeze(std(SumSquaredResiduals_perK_perFoldRep{i_win}{i_sensor}));
    end
end

%14.2 Average Kprime and SumSquaredRes across reps within fold
%(only for plotting and comparison between folds, not for analysis)
for i_fold = 1:NumFolds
    for i_win = 1:NumWindows
        for i_sensor = 1:NumSensors
            %Kprime
            Kprime_perFoldavgRep{i_win}(i_sensor,i_fold) = ...
                mean(Kprime_perFoldRep{i_win}{i_sensor}(i_fold,:));
            Kprime_perFoldSDRep{i_win}(i_sensor,i_fold) = ...
                std(Kprime_perFoldRep{i_win}{i_sensor}(i_fold,:));
            %SSR Kprime
            SumSquaredResiduals_Kprime_perFoldavgRep{i_win}(i_sensor,i_fold) = ...
                mean(SumSquaredResiduals_Kprime_perFoldRep{i_win}{i_sensor}(i_fold, :));
            SumSquaredResiduals_Kprime_perFoldSDRep{i_win}(i_sensor,i_fold) = ...
                std(SumSquaredResiduals_Kprime_perFoldRep{i_win}{i_sensor}(i_fold, :));
            %SSR all K
            SumSquaredResiduals_perK_perFoldavgRep{i_win}{i_sensor}(i_fold, :) = ...
                squeeze(mean(SumSquaredResiduals_perK_perFoldRep{i_win}{i_sensor}(i_fold, :, :),2));
            SumSquaredResiduals_perK_perFoldSDRep{i_win}{i_sensor}(i_fold, :) = ...
                squeeze(std(SumSquaredResiduals_perK_perFoldRep{i_win}{i_sensor}(i_fold, :, :)));
        end
    end
end

%14.3 Create summary file containing all relevant vars
ShuffKprime_data = struct;

ShuffKprime_data.perRepavgFold.Kprime_NullDistribution          = Kprime_perRepavgFold;
ShuffKprime_data.perRepavgFold.Kprime_NullDistribution_SD       = Kprime_perRepSDFold;

ShuffKprime_data.perRepavgFold.SumSquaredResiduals_Kprime       = SumSquaredResiduals_Kprime_perRepavgFold;
ShuffKprime_data.perRepavgFold.SumSquaredResiduals_Kprime_SD    = SumSquaredResiduals_Kprime_perRepSDFold;

ShuffKprime_data.perRepavgFold.SumSquaredResiduals_perK         = SumSquaredResiduals_perK_perRepavgFold;
ShuffKprime_data.perRepavgFold.SumSquaredResiduals_perK_SD      = SumSquaredResiduals_perK_perRepSDFold;

ShuffKprime_data.perFoldavgRep.Kprime                           = Kprime_perFoldavgRep;
ShuffKprime_data.perFoldavgRep.Kprime_SD                        = Kprime_perFoldSDRep;
ShuffKprime_data.perFoldavgRep.SumSquaredResiduals_Kprime       = SumSquaredResiduals_Kprime_perFoldavgRep;
ShuffKprime_data.perFoldavgRep.SumSquaredResiduals_Kprime_SD    = SumSquaredResiduals_Kprime_perFoldSDRep;
ShuffKprime_data.perFoldavgRep.SumSquaredResiduals_perK         = SumSquaredResiduals_perK_perFoldavgRep;
ShuffKprime_data.perFoldavgRep.SumSquaredResiduals_perK_SD      = SumSquaredResiduals_perK_perFoldSDRep;

%% 15) Save output variables
labels_loadedData = input_Data.label;

savefile = ...
    [path_save sub '_' input_DataLabel '_' tonedur_title ...
    '_HisTrack_ShuffKprimeCombFolds_unbalanced.mat'];
save(savefile, 'ShuffKprime_data','Info_TrialSelection','labels_loadedData', ...
    '-v7.3');

%% 16) Plot kprime per fold to compare consistency across sets (Optional)
if plotFig_CompFolds == 1
    
    %16.1) Plot shuffled Kprime per fold, averaged across reps within each fold 
    %(per sensor time win)
    figure('visible','off');
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    hold on; 
    sgtitle(['SHUFF Kprime per fold (Avg+STD across reps (' num2str(NumReps) ...
        ')) per sensor for: ' sub '; ' input_DataLabel '; ' ...
        tonedur_title '; ' Label_TW ]);
    
    plot_counter = 0;
    
    for i_win = 1:NumWindows
        for i_Fold = 1:NumFolds
            plot_counter = plot_counter+1;
            subplot(NumWindows, NumFolds+1, plot_counter)
            shadedErrorBar(1:NumSensors, ...
                Kprime_perFoldavgRep{i_win}(:,i_Fold)', ...
                Kprime_perFoldSDRep{i_win}(:,i_Fold)', ...
                'lineprops', 'k-');
            hold on;
            plot(1:NumSensors, ...
                Kprime_perFoldavgRep{i_win}(:,i_Fold)', ...
                'k-','LineWidth',1);
            
            title(['TW =' num2str(i_win) '; Fold: ' num2str(i_Fold)]);
            set(gca,'FontSize',8)
            ylabel('Kprime');
            xlabel('Electrodes');
            ylim([-2 20])
            xlim([0 NumSensors+1])
        end
        
        %Plot average across folds
        plot_counter = plot_counter+1;
        subplot(NumWindows,NumFolds+1,plot_counter)        
        shadedErrorBar(1:NumSensors, ...
            mean(Kprime_perFoldavgRep{i_win},2)', ...
            std(Kprime_perFoldavgRep{i_win}'), ...
            'lineprops', 'b-');
        hold on;
        plot(1:NumSensors, ...
            mean(Kprime_perFoldavgRep{i_win},2)', ...
            'b-','LineWidth',1);
        
        title(['TW =' num2str(i_win) ';Avg (SD) across Folds']);
        set(gca,'FontSize',8)
        ylabel('Kprime');
        xlabel('Electrodes');
        ylim([-2 20])
        xlim([0 NumSensors+1])
    end
    
    %save Fig
    filename = ['SHUFFKprimeAcrossFolds_' sub '_' input_DataLabel '_' ...
        tonedur_title '_' Label_TW '.png'];
    figfile      = [path_fig_Kprime filename];
    saveas(gcf, [figfile], 'png'); %save png version
   
    
    %16.2) Plot difference in minimum squared residuals divided by length of model order
    %Value by which Kprime is determined across odd/even sets (per sensor, time win)
    figure('visible','off');
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    hold on; 
    sgtitle(['MinSumSquaredRes/ModelOrder (defining Kprime - Avg+STD across folds (' ...
        num2str(NumFolds) ')) per sensor for: ' sub '; ' input_DataLabel '; ' ...
        tonedur_title '; ' Label_TW ]);

    plot_counter = 0;    
    for i_win = 1:NumWindows
        for i_Fold = 1:NumFolds
            plot_counter = plot_counter+1;
            subplot(NumWindows, NumFolds+1, plot_counter)            
            shadedErrorBar(1:NumSensors, ...
                SumSquaredResiduals_Kprime_perFoldavgRep{i_win}(:,i_Fold)', ...
                SumSquaredResiduals_Kprime_perFoldSDRep{i_win}(:,i_Fold)', ...
                'lineprops', 'k-');
            hold on;
            plot(1:NumSensors, ...
                SumSquaredResiduals_Kprime_perFoldavgRep{i_win}(:,i_Fold)', ...
                'k-','LineWidth',1);
            
            title(['TW =' num2str(i_win) '; Fold: ' num2str(i_Fold)]);
            set(gca,'FontSize',8)
            ylabel('MinSumSquaredRes/ModelOrder');
            xlabel('Electrodes');
            ylim auto
            xlim([0 NumSensors+1])
        end
        %Plot average across folds
        plot_counter = plot_counter+1;
        subplot(NumWindows, NumFolds+1, plot_counter)        
        shadedErrorBar(1:NumSensors, ...
            mean(SumSquaredResiduals_Kprime_perFoldavgRep{i_win},2)', ...
            std(SumSquaredResiduals_Kprime_perFoldavgRep{i_win}'), ...
            'lineprops', 'b-');
        hold on;
        plot(1:NumSensors, ...
            mean(SumSquaredResiduals_Kprime_perFoldavgRep{i_win},2)', ...
            'b-','LineWidth',1);
        
        title(['TW =' num2str(i_win) ';Avg (SD) across Folds']);
        set(gca,'FontSize',8)
        ylabel('MinSumSquaredRes/ModelOrder');
        xlabel('Electrodes');
        ylim auto
        xlim([0 NumSensors+1])
    end
    
    %save Fig
    filename = ['ShuffminSSRAcrossFolds_' sub '_' input_DataLabel '_' ...
        tonedur_title '_' Label_TW '.png'];
    figfile      = [path_fig_MSR filename];
    saveas(gcf, [figfile], 'png'); %save png version
    
    close all
end

end