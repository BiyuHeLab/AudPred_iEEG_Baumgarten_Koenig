%% 1) Specify SPM-input data in FT

%0) Set up paths and input vars%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%0.1) Specify paths and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')

NASTD_setVars
paths_NASTD = NASTD_paths(vars.location);

%Add base dir and own script dir
addpath(genpath(paths_NASTD.BaseDir));
addpath(genpath(paths_NASTD.ScriptsDir));
addpath(paths_NASTD.FieldTrip);


%% 0.2) Determine subject-specific parameters (whole-recording)
i_sub = 1; %proxy to load subs_PreProcSettings
sub_list = vars.sub_list; %patients
sub = sub_list{i_sub}

NASTD_subjectinfo %load subject info file (var: si)
subs_PreProcSettings = NASTD_SubjectPreProcSettings; %load in file with individual preproc infos

subs = si.sub_list;
tonedur_text = '0.4' %{'0.2' '0.4'};

%Define directory for Kprime data
Predict_dir = [paths_NASTD.ECoGdata_Prediction sub '/'];
mkdir(Predict_dir);
cd(Predict_dir);

loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];%path to indiv preprocessed ECoG data
loadfile_behav = [si.path_behavioraldata_sub];%path to indiv behavioral/stimulus data

%% 1) load preprocessed ECoG data
tic
disp('loading...')
load(loadfile_ECoGpreprocdata);
load(loadfile_behav);
disp(['done loading in ' num2str(toc) ' sec'])


%1) Determine input data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.1 Select only trials present in both behavioral and ECoG data and bring
    %them into same order
%Determine number of trials present in behav data
numoftrials_behav = num2str(length(data.trialNum)); %Should be around 120
numoftrials_ECoG = num2str(length(data_ECoGfiltref_trials.trial)); %Should be around 120

%Restrict data sets to those valid trials (i.e., without no-response trials)
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

%1.2 Select chosen trials from both behav & ECoG structs
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

%1.3 Select trials with specific tone length 
    %should be 60 trials per tone dur
    %tone length ID: 0.2 = 1, 0.4 = 2
if strcmp(tonedur_text,'0.2')
    tonedur_title = '200';
elseif  strcmp(tonedur_text,'0.4')
    tonedur_title = '400';
end

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

%% 2) Assure that all trials with same tone dur have equal trial-length
%2.1) Define time window parameters
nTrials  = length(data_ECoGfiltref_trials.trial);
% nSensors = size( data_input.trial{1}, 1 ); %all channels
nSensors = length(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF); %selected channels

fsample         = data_ECoGfiltref_trials.fsample;
toneDur_inSecs  = str2num(tonedur_text);
nSamplesPerTone = toneDur_inSecs * fsample;

nSamplesPerSeries = length(data_ECoGfiltref_trials.trial{1}(1,:));
nSecsPerSeries = nSamplesPerSeries/fsample;
nSecsPerSeriesinTheory = toneDur_inSecs*34;

%2.2 Define starting point (i.e., first tone)
t_start_ind(1) = 1;
t_end_ind(1)   = nSamplesPerTone;

%Read out start and stop samples for each tone (except the first)
for i_tone = 2:34 %loop across tones
    t_start_ind(i_tone) = t_end_ind(i_tone-1) + 1; %Start index is end index of previous tone +1
    t_end_ind(i_tone)   = t_start_ind(i_tone) + nSamplesPerTone - 1;%End index/final data point of i_tone is start+number of samples per tone -1
end

%Ensure that, if trial has an offset (i.e., trial-definition begins before
%first tone onset), t_start_ind agrees with onset index of first tone
ind_firsttone = find(data_ECoGfiltref_trials.time{1} == 0); %find index for start of first tone
    
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

%% 3 Create common input struct (for SPM)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_ECoG = [];
data_ECoG.hdr       = data_ECoGfiltref_trials.hdr;
data_ECoG.time      = data_ECoGfiltref_trials.time;
data_ECoG.fsample   = data_ECoGfiltref_trials.fsample;
data_ECoG.label     = data_ECoGfiltref_trials.label;
data_ECoG.trialinfo = data_ECoGfiltref_trials.trialinfo;
data_ECoG.sampleinfo= data_ECoGfiltref_trials.sampleinfo;
data_ECoG.cfg       = data_ECoGfiltref_trials.cfg;
data_ECoG.behav     = data;
data_ECoG.stim      = data.stim;

%copy trial wise info (i.e., all samples for selected tone for all channels and all trials into data arrays
for i_trial = 1:nTrials
    
    data_ECoG.trial{i_trial} = ...
        data_ECoGfiltref_trials.trial{i_trial}...
        (data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF, t_start_roundind(1) : t_end_roundind(end));
    data_ECoG.time{i_trial} = ...
        data_ECoGfiltref_trials.time{i_trial}(t_start_roundind(1) : t_end_roundind(end));
 
end
    
data_ECoG.label = data_ECoG.label(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF);
   
%3.2) Adjust header%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_hdr = data_ECoG.hdr;
temp_hdr.nChans = length(data_ECoG.label);
temp_hdr.label = data_ECoG.label;
temp_hdr.chantype = temp_hdr.chantype(1:length(data_ECoG.label));
for i_chan = 1:length(data_ECoG.label)
    temp_hdr.chantype{i_chan} = 'eeg';
end
temp_hdr.chanunit = temp_hdr.chanunit(1:length(data_ECoG.label));
temp_hdr.nTrials = length(data_ECoG.trial);

temp_hdr.elec.chanpos = temp_hdr.elec.chanpos(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF,:);
for i_chan = 1:length(data_ECoG.label)
    temp_hdr.elec.chantype{i_chan} = 'other'; %SPM accepts no other specification
end
temp_hdr.elec.chanunit = temp_hdr.elec.chanunit(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF);
temp_hdr.elec.elecpos = temp_hdr.elec.chanpos;
temp_hdr.elec.label = temp_hdr.elec.label(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF);
temp_hdr.elec.tra = diag(ones(1,length(data_ECoG.label)));

data_ECoG.hdr = temp_hdr;

%% 4) Reduce data set to only ROI channels and likely vs unlikely trials
%4.1 Determine trial groups (i.e., large vs. small p*34-p34 distance)
for i_trial = 1:length(data_ECoG.trial)
    logdiff_predvspresent(i_trial) = data_ECoG.stim.logf_pred(i_trial) - data_ECoG.stim.logf_final(i_trial);
end

f_likely = find(abs(logdiff_predvspresent) < median(abs(logdiff_predvspresent))); %lower than median p*34-p34 distance
f_unlikely = find(abs(logdiff_predvspresent) > median(abs(logdiff_predvspresent)));
f_unsure = find(abs(logdiff_predvspresent) == median(abs(logdiff_predvspresent))); %higher than median p*34-p34 distance

f_likely_resp = find(data_ECoG.behav.resp_prob > 3); %rating over 3
f_unsure_resp = find(data_ECoG.behav.resp_prob == 3);
f_unlikely_resp = find(data_ECoG.behav.resp_prob < 3); %rating below 3

data_ECoG.trialinfo(f_likely,3) = 1;
data_ECoG.trialinfo(f_unlikely,3) = -1;
data_ECoG.trialinfo(f_unsure,3) = 0;
data_ECoG.trialinfo(f_likely_resp,4) = 1;
data_ECoG.trialinfo(f_unlikely_resp,4) = -1;
data_ECoG.trialinfo(f_unsure_resp,4) = 0;

chan.A1 = {'G63'};
chan.STG = {'G46'};
chan.IFG = {'IF1'};

%MNI coordinate from: Yale BioImage Suite Package (yale.edu/mni2tal) -
%Lacadie et al. Neuromage (2008)
%A1
    %T1 coordinates:
    % G39 -64.159805 -50.935322 19.088045 100.00% cSTG 
    % G46 -62.884319 -39.909706 27.642714 100.00% cSTG 
    % G47 -65.450333 -40.965710 19.031191 100.00% cSTG 
    % G54 -64.122597 -30.033758 27.670658 100.00% cSTG 
    % G55 -66.271805 -31.305193 18.915991 100.00% cSTG 
    % G62 -63.625118 -20.027098 27.078493 100.00% mSTG 
    % G63 -65.397896 -21.648180 18.049524 100.00% mSTG 
    % T06 -64.231972 -7.734340 11.056163 100.00% mSTG 
    % T07 -65.498558 -15.902032 15.346996 100.00% mSTG 
    % T08 -65.058823 -24.614094 21.837561 100.00% cSTG 

    % A1 for MNI coordinates: -52 -19 7
    %MNI coordinates:
    % G39 -69.500000 -46.500000 15.500000 G 
    % G46 -70.000000 -34.000000 21.000000 G 
    % G47 -72.000000 -37.333333 12.666667 G 
    % G54 -72.000000 -26.000000 17.000000 G 
    % G55 -73.200000 -29.200000 8.400000 G 
    % G62 -72.666667 -17.000000 13.333333 G 
    % G63 -74.000000 -19.333333 5.333333 G 
    % T06 -73.142857 -7.142857 -5.142857 S 
    % T07 -74.000000 -14.666667 0.666667 S 
    % T08 -73.333333 -22.000000 8.666667 S 

%STG
    %T1 coordinates:
    % G39 -64.159805 -50.935322 19.088045 100.00% cSTG 
    % G46 -62.884319 -39.909706 27.642714 100.00% cSTG 
    % G47 -65.450333 -40.965710 19.031191 100.00% cSTG 
    % G54 -64.122597 -30.033758 27.670658 100.00% cSTG 
    % G55 -66.271805 -31.305193 18.915991 100.00% cSTG 
    % G62 -63.625118 -20.027098 27.078493 100.00% mSTG 
    % G63 -65.397896 -21.648180 18.049524 100.00% mSTG 
    % T06 -64.231972 -7.734340 11.056163 100.00% mSTG 
    % T07 -65.498558 -15.902032 15.346996 100.00% mSTG 
    % T08 -65.058823 -24.614094 21.837561 100.00% cSTG 

    %MNI coordinates:
    % G39 -69.500000 -46.500000 15.500000 G 
    % G46 -70.000000 -34.000000 21.000000 G 
    % G47 -72.000000 -37.333333 12.666667 G 
    % G54 -72.000000 -26.000000 17.000000 G 
    % G55 -73.200000 -29.200000 8.400000 G 
    % G62 -72.666667 -17.000000 13.333333 G 
    % G63 -74.000000 -19.333333 5.333333 G 
    % T06 -73.142857 -7.142857 -5.142857 S 
    % T07 -74.000000 -14.666667 0.666667 S 
    % T08 -73.333333 -22.000000 8.666667 S 

%IFG
    %T1
    % IF01 -42.187981 34.914268 34.117516 100.00% parstriangularis 
    % IF02 -46.726345 26.285583 36.303383 100.00% parstriangularis 
    % IF03 -51.203339 16.831656 37.808311 100.00% parsopercularis 
    % IF04 -52.773922 8.375415 41.319572 100.00% parsopercularis 
    % IF05 -54.814228 -1.799271 42.940117 100.00% precentral 
    % IF06 -56.709232 -11.409325 44.624222 100.00% precentral 

    %MNI
    % IF01 -57.500000 41.500000 11.000000 S 
    % IF02 -60.666667 34.000000 15.333333 S 
    % IF03 -64.666667 23.000000 18.666667 S 
    % IF04 -64.500000 15.500000 26.500000 S 
    % IF05 -66.000000 6.666667 31.333333 S 
    % IF06 -66.000000 -2.666667 34.666667 S 

%Change labels of selected channels to ensure right order
data_ECoG.label{find(strcmp(data_ECoG.label, 'G63'))} = '1_A1_G63'; % G63 -74.000000 -19.333333 5.333333 G 
data_ECoG.label{find(strcmp(data_ECoG.label, 'G46'))} = '2_STG_G46'; % G46 -70.000000 -34.000000 21.000000 G 
data_ECoG.label{find(strcmp(data_ECoG.label, 'IF1'))} = '3_IFG_IF1'; % IF01 -57.500000 41.500000 11.000000 S 

cd '/isilon/LFMI/VMdrive/Thomas/toolboxes/fieldtrip-20190314/utilities/' %FT link
% cd '/isilon/LFMI/VMdrive/Thomas/toolboxes/fieldtrip-20190314/fileio/private/' %FT link

cfg = [];
%  cfg.trials = [f_likely f_unlikely];
% cfg.channel = {'G63', 'G46', 'IF1'};
cfg.channel = {'1_A1_G63', '2_STG_G46', '3_IFG_IF1'};
data_ECoG_reduced = ft_selectdata(cfg,data_ECoG);

data_ECoG_reduced.sampleinfo = data_ECoG.sampleinfo;

%Switch order of channels in label and trial subfields
temp = data_ECoG_reduced.label{2};
data_ECoG_reduced.label{2} = data_ECoG_reduced.label{1};
data_ECoG_reduced.label{1} = temp;

for i_trials = 1:length(data_ECoG_reduced.trial)
    temp = data_ECoG_reduced.trial{i_trials}(2,:);
    data_ECoG_reduced.trial{i_trials}(2,:) = data_ECoG_reduced.trial{i_trials}(1,:);
    data_ECoG_reduced.trial{i_trials}(1,:) = temp;
end


%% 5. ERP/Timelock analysis
%determine trials
f_likely = find(data_ECoG_reduced.trialinfo(:,3) == 1);
f_unlikely = find(data_ECoG_reduced.trialinfo(:,3) == -1);
f_likely_resp = find(data_ECoG_reduced.trialinfo(:,4) == 1);
f_unlikely_resp = find(data_ECoG_reduced.trialinfo(:,4) == -1);

%determine samples
samples_p33 = data_ECoG_reduced.time{1}(find(data_ECoG_reduced.time{1} == data_ECoG_reduced.time{1}(end-floor(nSamplesPerTone)*2))...
    :find(data_ECoG_reduced.time{1} == data_ECoG_reduced.time{1}(end-floor(nSamplesPerTone))));
index_p33 = find(data_ECoG_reduced.time{1} == data_ECoG_reduced.time{1}(end-floor(nSamplesPerTone)*2))...
    :find(data_ECoG_reduced.time{1} == data_ECoG_reduced.time{1}(end-floor(nSamplesPerTone)));

samples_p34 = data_ECoG_reduced.time{1}(find(data_ECoG_reduced.time{1} == data_ECoG_reduced.time{1}(end-floor(nSamplesPerTone)))...
    :find(data_ECoG_reduced.time{1} == data_ECoG_reduced.time{1}(end)));
index_p34 = find(data_ECoG_reduced.time{1} == data_ECoG_reduced.time{1}(end-floor(nSamplesPerTone)))...
    :find(data_ECoG_reduced.time{1} == data_ECoG_reduced.time{1}(end));

%5.1 Compute ERP
%Manually
%Likely
    for i_chan = 1:length(data_ECoG_reduced.label)
        for i_trial = f_likely'
            p33_perchan_likely{i_chan}(i_trial,:) = data_ECoG_reduced.trial{i_trial}(i_chan,index_p33);
            p34_perchan_likely{i_chan}(i_trial,:) = data_ECoG_reduced.trial{i_trial}(i_chan,index_p34);
        end
    ERP_p33_likely{i_chan} = mean(p33_perchan_likely{i_chan});
    ERP_p34_likely{i_chan} = mean(p34_perchan_likely{i_chan});

    % STD_p33_likely{i_chan} = std(p33_perchan_likely{i_chan});
    % STD_p34_likely{i_chan} = std(p34_perchan_likely{i_chan});
    SE_p33_likely{i_chan} = std(p33_perchan_likely{i_chan}) / sqrt(length(p33_perchan_likely{i_chan}));
    SE_p34_likely{i_chan} = std(p34_perchan_likely{i_chan}) / sqrt(length(p34_perchan_likely{i_chan}))


    baseline_amp_p33_likely{i_chan} = mean(ERP_p33_likely{i_chan}');
    ERP_p34_BC_likely{i_chan} = ERP_p34_likely{i_chan} - baseline_amp_p33_likely{i_chan};
    end
%Unlikely
    for i_chan = 1:length(data_ECoG_reduced.label)
        for i_trial = f_unlikely'
            p33_perchan_unlikely{i_chan}(i_trial,:) = data_ECoG_reduced.trial{i_trial}(i_chan,index_p33);
            p34_perchan_unlikely{i_chan}(i_trial,:) = data_ECoG_reduced.trial{i_trial}(i_chan,index_p34);
        end
    ERP_p33_unlikely{i_chan} = mean(p33_perchan_unlikely{i_chan});
    ERP_p34_unlikely{i_chan} = mean(p34_perchan_unlikely{i_chan});

    SE_p33_unlikely{i_chan} = std(p33_perchan_unlikely{i_chan}) / sqrt(length(p33_perchan_unlikely{i_chan}));
    SE_p34_unlikely{i_chan} = std(p34_perchan_unlikely{i_chan}) / sqrt(length(p34_perchan_unlikely{i_chan}))

    baseline_amp_p33_unlikely{i_chan} = mean(ERP_p33_unlikely{i_chan}');
    ERP_p34_BC_unlikely{i_chan} = ERP_p34_unlikely{i_chan} - baseline_amp_p33_unlikely{i_chan};
    end
%likely_resp
for i_chan = 1:length(data_ECoG_reduced.label)
    for i_trial = f_likely_resp'
        p33_perchan_likely_resp{i_chan}(i_trial,:) = data_ECoG_reduced.trial{i_trial}(i_chan,index_p33);
        p34_perchan_likely_resp{i_chan}(i_trial,:) = data_ECoG_reduced.trial{i_trial}(i_chan,index_p34);
    end
ERP_p33_likely_resp{i_chan} = mean(p33_perchan_likely_resp{i_chan});
ERP_p34_likely_resp{i_chan} = mean(p34_perchan_likely_resp{i_chan});

% STD_p33_likely_resp{i_chan} = std(p33_perchan_likely_resp{i_chan});
% STD_p34_likely_resp{i_chan} = std(p34_perchan_likely_resp{i_chan});
SE_p33_likely_resp{i_chan} = std(p33_perchan_likely_resp{i_chan}) / sqrt(length(p33_perchan_likely_resp{i_chan}));
SE_p34_likely_resp{i_chan} = std(p34_perchan_likely_resp{i_chan}) / sqrt(length(p34_perchan_likely_resp{i_chan}))


baseline_amp_p33_likely_resp{i_chan} = mean(ERP_p33_likely_resp{i_chan}');
ERP_p34_BC_likely_resp{i_chan} = ERP_p34_likely_resp{i_chan} - baseline_amp_p33_likely_resp{i_chan};
end

for i_chan = 1:length(data_ECoG_reduced.label)
    for i_trial = f_unlikely_resp'
        p33_perchan_unlikely_resp{i_chan}(i_trial,:) = data_ECoG_reduced.trial{i_trial}(i_chan,index_p33);
        p34_perchan_unlikely_resp{i_chan}(i_trial,:) = data_ECoG_reduced.trial{i_trial}(i_chan,index_p34);
    end
ERP_p33_unlikely_resp{i_chan} = mean(p33_perchan_unlikely_resp{i_chan});
ERP_p34_unlikely_resp{i_chan} = mean(p34_perchan_unlikely_resp{i_chan});

SE_p33_unlikely_resp{i_chan} = std(p33_perchan_unlikely_resp{i_chan}) / sqrt(length(p33_perchan_unlikely_resp{i_chan}));
SE_p34_unlikely_resp{i_chan} = std(p34_perchan_unlikely_resp{i_chan}) / sqrt(length(p34_perchan_unlikely_resp{i_chan}))

baseline_amp_p33_unlikely_resp{i_chan} = mean(ERP_p33_unlikely_resp{i_chan}');
ERP_p34_BC_unlikely_resp{i_chan} = ERP_p34_unlikely_resp{i_chan} - baseline_amp_p33_unlikely_resp{i_chan};
end

% figure;
% subplot(1,2,1)
% hold on;
% plot(samples_p34,ERP_p34_BC_likely{1},'b')
% plot(samples_p34,ERP_p34_BC_likely{2},'r')
% plot(samples_p34,ERP_p34_BC_likely{3},'y')
% legend('A1','STG','IFG')
% title('Likely (< median p*34-p34) p34')
% ylim([-15 15])
% subplot(1,2,2)
% hold on;
% plot(samples_p34,ERP_p34_BC_unlikely{1},'b')
% plot(samples_p34,ERP_p34_BC_unlikely{2},'r')
% plot(samples_p34,ERP_p34_BC_unlikely{3},'y')
% legend('A1','STG','IFG')
% title('Unlikely (> median p*34-p34) p34')
% ylim([-15 15])

% BC - likely vs unlikely
    cd '/isilon/LFMI/VMdrive/Thomas'/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/MEG/NASTD_MEG_Matlab/Plotting/
    figure;
    subplot(1,2,1)
    hold on;
    shadedErrorBar(samples_p34,ERP_p34_BC_likely{1},SE_p34_likely{1},'lineProps','b','patchSaturation',0.2)
    plot(samples_p34,ERP_p34_BC_likely{1},'b','lineWidth',2)
    shadedErrorBar(samples_p34,ERP_p34_BC_likely{2},SE_p34_likely{2},'lineProps','r','patchSaturation',0.2)
    plot(samples_p34,ERP_p34_BC_likely{2},'r','lineWidth',2)
    shadedErrorBar(samples_p34,ERP_p34_BC_likely{3},SE_p34_likely{3},'lineProps','y','patchSaturation',0.2)
    plot(samples_p34,ERP_p34_BC_likely{3},'y','lineWidth',2)
    title('Likely (< median p*34-p34) p34')
    legend('A1','','STG','','IFG')
    ylim([-10 15])
    xlabel('Time [ms]')
    subplot(1,2,2)
    hold on;
    shadedErrorBar(samples_p34,ERP_p34_BC_unlikely{1},SE_p34_unlikely{1},'lineProps','b','patchSaturation',0.2)
    plot(samples_p34,ERP_p34_BC_unlikely{1},'b','lineWidth',2)
    shadedErrorBar(samples_p34,ERP_p34_BC_unlikely{2},SE_p34_unlikely{2},'lineProps','r','patchSaturation',0.2)
    plot(samples_p34,ERP_p34_BC_unlikely{2},'r','lineWidth',2)
    shadedErrorBar(samples_p34,ERP_p34_BC_unlikely{3},SE_p34_unlikely{3},'lineProps','y','patchSaturation',0.2)
    plot(samples_p34,ERP_p34_BC_unlikely{3},'y','lineWidth',2)
    title('Unlikely (> median p*34-p34) p34')
    legend('A1','','STG','','IFG')
    ylim([-10 15])
    xlabel('Time [ms]')
    suptitle('Baseline-corrected (by avg amp of p33) ERPs for p34 (shading = SE across trials)')
%Non-BC - likely vs unlikely
    cd '/isilon/LFMI/VMdrive/Thomas'/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/MEG/NASTD_MEG_Matlab/Plotting/
    figure;
    subplot(1,2,1)
    hold on;
    shadedErrorBar(samples_p34,ERP_p34_likely{1},SE_p34_likely{1},'lineProps','b','patchSaturation',0.2)
    plot(samples_p34,ERP_p34_likely{1},'b','lineWidth',2)
    shadedErrorBar(samples_p34,ERP_p34_likely{2},SE_p34_likely{2},'lineProps','r','patchSaturation',0.2)
    plot(samples_p34,ERP_p34_likely{2},'r','lineWidth',2)
    shadedErrorBar(samples_p34,ERP_p34_likely{3},SE_p34_likely{3},'lineProps','y','patchSaturation',0.2)
    plot(samples_p34,ERP_p34_likely{3},'y','lineWidth',2)
    title('Likely (< median p*34-p34) p34')
    legend('A1','','STG','','IFG')
    ylim([-10 15])
    xlabel('Time [ms]')
    subplot(1,2,2)
    hold on;
    shadedErrorBar(samples_p34,ERP_p34_unlikely{1},SE_p34_unlikely{1},'lineProps','b','patchSaturation',0.2)
    plot(samples_p34,ERP_p34_unlikely{1},'b','lineWidth',2)
    shadedErrorBar(samples_p34,ERP_p34_unlikely{2},SE_p34_unlikely{2},'lineProps','r','patchSaturation',0.2)
    plot(samples_p34,ERP_p34_unlikely{2},'r','lineWidth',2)
    shadedErrorBar(samples_p34,ERP_p34_unlikely{3},SE_p34_unlikely{3},'lineProps','y','patchSaturation',0.2)
    plot(samples_p34,ERP_p34_unlikely{3},'y','lineWidth',2)
    title('Unlikely (> median p*34-p34) p34')
    legend('A1','','STG','','IFG')
    ylim([-10 15])
    xlabel('Time [ms]')
    suptitle('Non-Baseline-corrected ERPs for p34 (shading = SE across trials)')

% BC - likely_resp vs unlikely_resp
cd '/isilon/LFMI/VMdrive/Thomas'/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/MEG/NASTD_MEG_Matlab/Plotting/
figure;
subplot(1,2,1)
hold on;
shadedErrorBar(samples_p34,ERP_p34_BC_likely_resp{1},SE_p34_likely_resp{1},'lineProps','b','patchSaturation',0.2)
plot(samples_p34,ERP_p34_BC_likely_resp{1},'b','lineWidth',2)
shadedErrorBar(samples_p34,ERP_p34_BC_likely_resp{2},SE_p34_likely_resp{2},'lineProps','r','patchSaturation',0.2)
plot(samples_p34,ERP_p34_BC_likely_resp{2},'r','lineWidth',2)
shadedErrorBar(samples_p34,ERP_p34_BC_likely_resp{3},SE_p34_likely_resp{3},'lineProps','y','patchSaturation',0.2)
plot(samples_p34,ERP_p34_BC_likely_resp{3},'y','lineWidth',2)
title('Response Likely (> 3) p34')
legend('A1','','STG','','IFG')
ylim([-10 15])
xlabel('Time [ms]')
subplot(1,2,2)
hold on;
shadedErrorBar(samples_p34,ERP_p34_BC_unlikely_resp{1},SE_p34_unlikely_resp{1},'lineProps','b','patchSaturation',0.2)
plot(samples_p34,ERP_p34_BC_unlikely_resp{1},'b','lineWidth',2)
shadedErrorBar(samples_p34,ERP_p34_BC_unlikely_resp{2},SE_p34_unlikely_resp{2},'lineProps','r','patchSaturation',0.2)
plot(samples_p34,ERP_p34_BC_unlikely_resp{2},'r','lineWidth',2)
shadedErrorBar(samples_p34,ERP_p34_BC_unlikely_resp{3},SE_p34_unlikely_resp{3},'lineProps','y','patchSaturation',0.2)
plot(samples_p34,ERP_p34_BC_unlikely_resp{3},'y','lineWidth',2)
title('Unlikely (> median p*34-p34) p34')
legend('A1','','STG','','IFG')
ylim([-10 15])
xlabel('Time [ms]')
suptitle('Baseline-corrected (by avg amp of p33) ERPs for p34 (shading = SE across trials)')

%Non-BC - likely_resp vs unlikely_resp
cd '/isilon/LFMI/VMdrive/Thomas'/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/MEG/NASTD_MEG_Matlab/Plotting/
figure;
subplot(1,2,1)
hold on;
shadedErrorBar(samples_p34,ERP_p34_likely_resp{1},SE_p34_likely_resp{1},'lineProps','b','patchSaturation',0.2)
plot(samples_p34,ERP_p34_likely_resp{1},'b','lineWidth',2)
shadedErrorBar(samples_p34,ERP_p34_likely_resp{2},SE_p34_likely_resp{2},'lineProps','r','patchSaturation',0.2)
plot(samples_p34,ERP_p34_likely_resp{2},'r','lineWidth',2)
shadedErrorBar(samples_p34,ERP_p34_likely_resp{3},SE_p34_likely_resp{3},'lineProps','y','patchSaturation',0.2)
plot(samples_p34,ERP_p34_likely_resp{3},'y','lineWidth',2)
title('Response Unlikely (< 3) p34')
legend('A1','','STG','','IFG')
ylim([-10 15])
xlabel('Time [ms]')
subplot(1,2,2)
hold on;
shadedErrorBar(samples_p34,ERP_p34_unlikely_resp{1},SE_p34_unlikely_resp{1},'lineProps','b','patchSaturation',0.2)
plot(samples_p34,ERP_p34_unlikely_resp{1},'b','lineWidth',2)
shadedErrorBar(samples_p34,ERP_p34_unlikely_resp{2},SE_p34_unlikely_resp{2},'lineProps','r','patchSaturation',0.2)
plot(samples_p34,ERP_p34_unlikely_resp{2},'r','lineWidth',2)
shadedErrorBar(samples_p34,ERP_p34_unlikely_resp{3},SE_p34_unlikely_resp{3},'lineProps','y','patchSaturation',0.2)
plot(samples_p34,ERP_p34_unlikely_resp{3},'y','lineWidth',2)
title('Response Unlikely (< 3) p34')
legend('A1','','STG','','IFG')
ylim([-10 15])
xlabel('Time [ms]')
suptitle('Non-Baseline-corrected ERPs for p34 (shading = SE across trials)')


cd
%FT
cfg = [];
% cfg.latency = [data_ECoG_reduced.time{1}(end-floor(nSamplesPerTone)) data_ECoG_reduced.time{1}(end)]; %only final tone (p34)
cfg.keeptrials = 'no';
cfg.removemean = 'no';


cfg.trials = f_likely;
ERP_likely = ft_timelockanalysis(cfg, data_ECoG_reduced);
cfg.trials = f_unlikely;
ERP_unlikely = ft_timelockanalysis(cfg, data_ECoG_reduced);
cfg.trials = f_likely_resp;
ERP_likely_resp = ft_timelockanalysis(cfg, data_ECoG_reduced);
cfg.trials = f_unlikely_resp;
ERP_unlikely_resp = ft_timelockanalysis(cfg, data_ECoG_reduced);

%5.2 Compute baseline correction of p34 by subtracting average from p33
cfg = [];
cfg.baseline = [data_ECoG_reduced.time{1}(end-floor(nSamplesPerTone)*2) data_ECoG_reduced.time{1}(end-floor(nSamplesPerTone))];
ERP_likely_BC = ft_timelockbaseline(cfg, ERP_likely);
ERP_unlikely_BC = ft_timelockbaseline(cfg, ERP_unlikely);
ERP_likely_resp_BC = ft_timelockbaseline(cfg, ERP_likely_resp);
ERP_unlikely_resp_BC = ft_timelockbaseline(cfg, ERP_unlikely_resp);


%5.3 Plot ERP
% cfg = [];
% cfg.xlim = data_ECoG_reduced.time{1}(find(data_ECoG_reduced.time{1} == data_ECoG_reduced.time{1}(end-floor(nSamplesPerTone)))...
%     :find(data_ECoG_reduced.time{1} == data_ECoG_reduced.time{1}(end)))
% cfg.ylim = 'maxmin';
% cfg.channels = '1_A1_G63';
% ft_singleplotER(cfg,ERP_likely)


data_ECoG_reduced.time{1}(find(data_ECoG_reduced.time{1} == data_ECoG_reduced.time{1}(end-floor(nSamplesPerTone)))...
    :find(data_ECoG_reduced.time{1} == data_ECoG_reduced.time{1}(end)))

figure;
hold on;
plot(data_ECoG_reduced.time{1}(find(data_ECoG_reduced.time{1} == data_ECoG_reduced.time{1}(end-floor(nSamplesPerTone)))...
    :find(data_ECoG_reduced.time{1} == data_ECoG_reduced.time{1}(end))), ...
    ERP_likely_BC.avg(1,find(data_ECoG_reduced.time{1} == data_ECoG_reduced.time{1}(end-floor(nSamplesPerTone)))...
    :find(data_ECoG_reduced.time{1} == data_ECoG_reduced.time{1}(end))))
plot(data_ECoG_reduced.time{1}(find(data_ECoG_reduced.time{1} == data_ECoG_reduced.time{1}(end-floor(nSamplesPerTone)))...
    :find(data_ECoG_reduced.time{1} == data_ECoG_reduced.time{1}(end))), ...
    ERP_likely_BC.avg(2,find(data_ECoG_reduced.time{1} == data_ECoG_reduced.time{1}(end-floor(nSamplesPerTone)))...
    :find(data_ECoG_reduced.time{1} == data_ECoG_reduced.time{1}(end))))
plot(data_ECoG_reduced.time{1}(find(data_ECoG_reduced.time{1} == data_ECoG_reduced.time{1}(end-floor(nSamplesPerTone)))...
    :find(data_ECoG_reduced.time{1} == data_ECoG_reduced.time{1}(end))), ...
    ERP_likely_BC.avg(3,find(data_ECoG_reduced.time{1} == data_ECoG_reduced.time{1}(end-floor(nSamplesPerTone)))...
    :find(data_ECoG_reduced.time{1} == data_ECoG_reduced.time{1}(end))))

legend('A1','STG','IFG')
figure;
plot(1:length(ERP_unlikely.avg), ERP_unlikely.avg)
legend('A1','STG','IFG')

%% 4) Save reduced data
%4.1 FT data 
savefile_matlab = [path_save sub '_reducedData_' tonedur_title 'ms_matlab.mat'];
save(savefile_matlab, 'data_ECoG_reduced','-v7.3');

%4.2 SPM data 
addpath('D:\Programs\MATLAB\spm12\') %path SPM
savefile_spm = [path_save sub '_reducedData_' tonedur_title 'ms.mat'];
data_ECoG_SPM = spm_eeg_ft2spm(data_ECoG_reduced, savefile_spm);

%% 5) Change Specific ins SPM data
cd 'D:\Work\Kongresse&Workshops\2019\OCNC Okinawa Summer School\Project\Data\Analysis\Preproc_ECoGdata\NY688\'
load('NY688_reducedData_200ms.mat')
load('NY688_reducedData_200ms_matlab.mat')

% Apply trial labels to SPM data
for i_trial = 1:length(D.trials)
    
    if data_ECoG_reduced.trialinfo(i_trial,3) == 1
        D.trials(i_trial).label = 'likely';
    elseif data_ECoG_reduced.trialinfo(i_trial,3) == -1
        D.trials(i_trial).label = 'unlikely';
    end
end

% Apply labels to condlist subfield
D.condlist{1} = 'likely';
D.condlist{2} = 'unlikely';

%Save SPM data
path_save = 'D:\Work\Kongresse&Workshops\2019\OCNC Okinawa Summer School\Project\Data\Analysis\Preproc_ECoGdata\NY688\';
savefile_spm = [path_save 'NY688_reducedData_' tonedur_title 'ms.mat'];
save(savefile_spm, 'D','-v7.3');