function NASTD_ECoG_Predict_CompPredp33_Subs...
    (sub, ToneDur_text, inputData, ...
    predictive_sequencerange, toneIndex, SamplesTW, BLcorrect,...
    paths_NASTD_ECoG, vars)
%Aim: How is neural activity at tone 33 modulated by the expected value of
%tone 34 (Predp34, which itself depends on the previous tone sequence input.
%Linear regression of pitch values p1-32 on neural activity during p33.

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')

%Add base dir and own script dir
addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));

%% 0.2) Determine subject-specific parameters (whole-recording)
NASTD_ECoG_subjectinfo %load subject info file (var: si)
subs_PreProcSettings = NASTD_ECoG_SubPreprocSettings; %load in file with individual preproc infos

%Tone duration condition for data load-in
if strcmp(ToneDur_text,'0.2')
    tonedur_title = '200';
elseif  strcmp(ToneDur_text,'0.4')
    tonedur_title = '400';
elseif  strcmp(ToneDur_text,'all')
    tonedur_title = 'all';   
end

%Load in preprocessed neural and behavioral data
loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];%path to indiv preprocessed ECoG data

tic
disp(['Loading preprocessed data set for sub: ' sub])
load(loadfile_ECoGpreprocdata);
disp(['done loading in ' num2str(toc) ' sec'])

%Define directory for data
Label_TW = num2str(round(SamplesTW/preprocData_AllTrials.fsample,2));

path_save = [paths_NASTD_ECoG.ECoGdata_Prediction sub '/Exp/' Label_TW 'sTW/'];
mkdir(path_save); %make new directory for output files

%% 1) Select input data 
%1.1) Select clean/valid electrodes and trials
Selected_Channels = setdiff(preprocData_AllTrials.cfg.info_elec.selected.index4EDF, subs_PreProcSettings.(sub).rejectedChan_index);
preprocData_CleanTrials = NASTD_ECoG_SelCleanTrials(sub,preprocData_AllTrials);

%1.2) Select only trials with specific Tone Duration (TD)    
preprocData_perTD = NASTD_ECoG_SelTrialsperTD(ToneDur_text, preprocData_CleanTrials);
disp(['Total trial number for TD ' ToneDur_text 's: ' num2str(length(preprocData_perTD.trial))])

% figure;
% plot(preprocData_perTD.trial{1}(3,:))
% hold on;
% plot(preprocData_perTD.trial{2}(3,:))
% hold on;
% plot(preprocData_perTD.trial{3}(3,:))

%1.3 Create common input struct and select input measure
if strcmp(inputData,'LF')
    inputDataFormat = 'trial';
% elseif strcmp(inputData,'GammaAmp')
%     inputDataFormat = 'GammaAmp_Norm';
elseif strcmp(inputData,'LogGammaAmp')
    inputDataFormat = 'LogGammaAmp_Norm';
end

data_input = preprocData_perTD;
data_input.trial = preprocData_perTD.(inputDataFormat);
data_input.stim = preprocData_perTD.behav.stim;

clear preprocData_perTD


%1.3) Baseline-correct (per trial, prestimulus baseline) - optional
if BLcorrect == 1
    BL_label = 'prestimBC';
    
    BLwindow = [-0.5 0]; %window for baseline correction (0 = start of 1st tone)
    data_input = NASTD_ECoG_Predict_BLCdata(BLwindow, data_input);
    disp(['Baseline correcting data (BL window: ' num2str(BLwindow(1)) '-' num2str(BLwindow(2)) 's'])

    % %visually check BC
    % test = data_input;
    % % test.trial = data_input.LogGammaAmp;
    % 
    % cfg             = [];
    % cfg.ylim        = 'maxabs';
    % cfg.viewmode    = 'vertical';
    % cfg.channel     = Selected_Channels;
    % ft_databrowser(cfg,test)
else
     BL_label = 'noBC';   
end

% figure;
% plot(data_input.trial{1}(3,:))
% hold on;
% plot(data_input.trial{2}(3,:))
% hold on;
% plot(data_input.trial{3}(3,:))

%% 2) Extract last two tones per trial, as well as relevant stim/behav info
%2.1) Define time window parameters depending on TD
nTrials  = length(data_input.trial);
nSensors = length(Selected_Channels); %selected channels

fsample         = data_input.fsample;
toneDur_inSecs  = str2num(ToneDur_text);
nSamplesPerTone = toneDur_inSecs * fsample;

%2.2 Determine TP/samples for each tone start+end
TP_Tone_StartStop = NaN(35,2);
TP_Tone_StartStop(1,1) = 0;
for i_tone = 2:35
    Dist = abs(data_input.time{1} - ((str2num(ToneDur_text)*i_tone)-str2num(ToneDur_text)));
    minDist = min(Dist);
    i_minDist = find(Dist == minDist);
    TP_Tone_StartStop(i_tone,1) = data_input.time{1}(i_minDist);
end
for i_tone = 1:34
    i_LastSampleTone = find(data_input.time{1} == TP_Tone_StartStop(i_tone+1,1));
    TP_Tone_StartStop(i_tone,2) = data_input.time{1}(i_LastSampleTone);
end
TP_Tone_StartStop = TP_Tone_StartStop(1:34,:);

%check if all tones are of equal length, if not then choose min length by deleting the last sample
minSeqLength_sec = min(TP_Tone_StartStop(:,2)-TP_Tone_StartStop(:,1));
for i_tone = 1:34
    if TP_Tone_StartStop(i_tone,2)-TP_Tone_StartStop(i_tone,1) > minSeqLength_sec
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

%2.3) Initialize and fill 3D data arrays (nSens, nSamples per tone, nTrials)
ECoGdata_p33 = zeros(nSensors, length(Sample_Tone_StartStop(1,1):Sample_Tone_StartStop(1,2)), nTrials);
ECoGdata_p34 = ECoGdata_p33;
%copy trial wise info (i.e., all samples for selected tone for all channels and all trials into data arrays
for i_trial = 1:nTrials    
    ECoGdata_p33(:, :, i_trial) = ...
        data_input.trial{i_trial}(Selected_Channels, Sample_Tone_StartStop(toneIndex,1) : Sample_Tone_StartStop(toneIndex,2));
    ECoGdata_p34(:, :, i_trial) = ...
        data_input.trial{i_trial}(Selected_Channels, Sample_Tone_StartStop(toneIndex+1,1) : Sample_Tone_StartStop(toneIndex+1,2)); 
end

%% 3) Get predicted final tones based on each predictive_sequencerange ***
%for each trial, compute/read out 
%1) the tone pitch predicted given the data so far (p*34), non-discretized
%2) the p*34 (same as above, only discretized to next full tones)
%3) the p34 (i.e., the really presented last tone)
%4) prediction error (i.e., the distance between presented (p34) and predicted tone 
% predictive_sequencerange = [33]; %tone number on which corresponding neural data prediction is based?

for i_k = 1:length(predictive_sequencerange) %for each tone where predition is analyzed
    
    series_start_ind = toneIndex - predictive_sequencerange(i_k) + 1; %beginning of tone sequence part used for prediction
    series_end_ind   = toneIndex; %end of tone sequence part used for prediction
    
    for i_trial = 1:nTrials
        series = data_input.stim.series_f{i_trial}(series_start_ind : series_end_ind); %(non-log) selected sequence tones in Hz
        beta = data_input.stim.beta(i_trial); %always 1.5 in ECoG task
        
        if predictive_sequencerange(i_k) == 1 %for first tone
            series_predp34_nondiscretized(i_k, i_trial) = ...
                log( data_input.stim.series_f{i_trial}(toneIndex) ); % predict same tone repeating, because no further info for prediciton
        else
            series_predp34_nondiscretized(i_k, i_trial) = ...
                log( find_predicted_tone(series, beta) ); %subfunction: reads out predicted next tone based on beta level and input tone series
        end
              
        series_p33(i_trial) = log( data_input.stim.series_f{i_trial}(toneIndex) ); %log(p33)
        series_p34(i_trial) = log( data_input.stim.series_f{i_trial}(toneIndex+1) ); %log(p34)
        series_predp34_discretized(i_trial) = data_input.stim.logf_pred(i_trial); %log(p*34) discretized             

%         series_pred_error(i_k, i_trial) = log( data_input.stim.series_f{i_trial}(toneIndex+1) ) - log( series_pred(i_k, i_trial) );
%         series_pred_error(i_k, i_trial) = series_p34(i_trial) - series_predp34_nondiscretized(i_k, i_trial); %non-discretized p*34 used
        series_pred_error(i_k, i_trial) = series_p34(i_trial) - series_predp34_discretized (i_k, i_trial);  %discretized p*34 used
        %prediction error as difference between actually presented final tone (log(p34)) and computationally determined final tone
       
        series_simplepred_error(i_k, i_trial) = series_p34(i_trial) - series_p33(i_trial); 
        %simple prediction error as difference between actually presented final tone (log(p34)) and p33
        
    end
    
end

% %compare predicted tone via subfunction with sequence-defiend p*34 - should
% %be very close, p*34 should be same but only scaled to next full tone
% [series_p34_nondiscretized', data_input.stim.logf_pred', series_p34_nondiscretized'-data_input.stim.logf_pred']
% [exp(series_p34_nondiscretized'), exp(data_input.stim.logf_pred'), series_p34_nondiscretized'-data_input.stim.logf_pred']

%% 4) Compute avg ERF activity per time window for current tone and next tone
%4.1) Define windows
fsample         = data_input.fsample;

win_size    = SamplesTW; %25 data points with fs 512 = ~49ms as mentioned in manuscript
disp(['Sample-Size for TimeWindow: ' num2str(SamplesTW) '; = ' num2str(SamplesTW/fsample) 's'])

s_per_win = (1/fsample)*win_size;
win_overlap = 0; %as mentioned in manuscript

toneDur_inSecs  = str2num(ToneDur_text);
nSamplesPerTone = toneDur_inSecs * fsample; %no rounding

%4.2 Define number, start and end sample of window per tone 
windows = [1 win_size];
while windows(end,end) < nSamplesPerTone
    windows = [windows; windows(end,:) + (win_size - win_overlap)];
end

if windows(end,end) > nSamplesPerTone %TJB: windows within single tone
    windows(end,:) = [];
end

windows_inms = (windows / fsample) * 1000;

%4.3) Compute avg activity per window for each sensor and trial
for i_win = 1:size(windows,1)

    ind_start = windows(i_win, 1);
    ind_end   = windows(i_win, 2);

    ERF_p33_win{i_win} = squeeze( mean( ECoGdata_p33(:, ind_start:ind_end, :), 2) );
    ERF_p34_win{i_win} = squeeze( mean( ECoGdata_p34(:, ind_start:ind_end, :), 2) );
    %Output: channel*trial matrix per window with 1 ERF values averaged
    %across all samples per window
end


%% 5) Compute measures of association between ERF windows and prediction / prediction error at each sensor

for i_k = 1:length(predictive_sequencerange)
    for i_win = 1:size(windows,1)
        for i_sensor = 1:nSensors

%%% for prediction at tone 33 (i.e., with ERF of current tone (33)) %%%            
            if predictive_sequencerange(i_k) > 1 %only if the currently focused tone is not the first tone (since for first tone we can't really predict anything)
                dv = ERF_p33_win{i_win}(i_sensor,:)';  % Across window averaged activity at this sensor for all trials
%                 iv_pred = series_predp34_nondiscretized(i_k,:)'; %predicted p34 tone pitch for all trials
                iv_pred = series_predp34_discretized(i_k,:)'; %predicted p34 tone pitch for all trials

                % linear regression between ERF (per channel, time window) and predicted tone pitch
%                 iv_pitch = series_p33'; %tone pitch of focused on tone (33)                
%                 stats = regstats(dv, [iv_pred, iv_pitch], 'linear');
                stats = regstats(dv, iv_pred, 'linear', 'tstat');
                
                
%                 iv_pred_quad = iv_pred - log(440); % center prediction on log(440) for quadratic regression
%                 stats = regstats(dv, iv_pred_quad, 'purequadratic', 'tstat');


                PredEffect.stats{i_k}{i_win}{i_sensor} = stats;
                PredEffect.tval{i_k}{i_win}(i_sensor) = stats.tstat.t(2);
                PredEffect.pval{i_k}{i_win}(i_sensor) = stats.tstat.pval(2);
            end
            
            
%%% for prediction error at tone 34 (i.e., with ERF of next tone (34)) %%%
            dv = ERF_p34_win{i_win}(i_sensor,:)';  % mean ERF window at this sensor for all trials
%             iv_pred_error = series_pred_error(i_k,:)';
            iv_pred_error = abs( series_pred_error(i_k,:)' );         
            %prediction error as difference between actually presented final
            %tone (log(p34)) and predicted final tone (log(predp34))

            iv_p34 = series_p34'; %presented final tone (p34)
            
            % linear regression t-stat
            stats = regstats(dv, [iv_pred_error, iv_p34], 'linear', 'tstat'); 
            %linear regression with 2 predictors (predictio error + p34)
%             stats = regstats(dv, iv_pred_error, 'linear');

            PredErrorEffect.stats{i_k}{i_win}{i_sensor} = stats;            
            PredErrorEffect.tval{i_k}{i_win}(i_sensor) = stats.tstat.t(2);
            PredErrorEffect.pval{i_k}{i_win}(i_sensor) = stats.tstat.pval(2);

%%% for simple prediction error at tone 34 (i.e., with ERF of next tone (34)) %%%
            dv = ERF_p34_win{i_win}(i_sensor,:)';  % mean ERF window at this sensor for all trials
%             iv_pred_error = series_pred_error(i_k,:)';
            iv_simplepred_error = abs( series_simplepred_error(i_k,:)' );         
            %prediction error as difference between actually presented final
            %tone (log(p34)) and computationally determined final tone

            iv_p34 = series_p34';
            
            % linear regression t-stat
            stats = regstats(dv, [iv_simplepred_error, iv_p34], 'linear', 'tstat');
%             stats = regstats(dv, iv_pred_error, 'linear');

            SimplePredErrorEffect.stats{i_k}{i_win}{i_sensor} = stats;            
            SimplePredErrorEffect.tval{i_k}{i_win}(i_sensor) = stats.tstat.t(2);
            SimplePredErrorEffect.pval{i_k}{i_win}(i_sensor) = stats.tstat.pval(2);
           
            
        end

    end
end

%% 8) Save variables
savefile = [path_save sub '_Predp33_regstatsEXP_' inputData '_' tonedur_title 'msTD_' BL_label '.mat'];

save(savefile, 'PredEffect', 'PredErrorEffect', 'SimplePredErrorEffect', ...              
               'predictive_sequencerange', 'toneIndex', ...
               'ERF_p33_win', 'ERF_p34_win', ...
               'series_predp34_nondiscretized', 'series_p33', 'series_pred_error', 'series_p34');
           
end