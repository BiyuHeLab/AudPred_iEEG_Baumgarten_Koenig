function subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings

%Subject-wise settings for data preprocessing

%% NY688
subs_PreProcSettings.NY688.MNI2EDFmismatches.index_mismatchchannels ...
    = [71:78]';
subs_PreProcSettings.NY688.MNI2EDFmismatches.label_mismatchchannels ...
    = [{'T1'}, {'T2'}, {'T3'},{'T4'},{'T5'},{'T6'},{'T7'},{'T8'}];
subs_PreProcSettings.NY688.MNI2EDFmismatches.correctedlabel_mismatchchannels ...
    = [{'TT1'}, {'TT2'}, {'TT3'},{'TT4'},{'TT5'},{'TT6'},{'TT7'},{'TT8'}];
subs_PreProcSettings.NY688.number_ECoGchan= 144; 
%number of channels in EDF set that contain ECoG recordings
subs_PreProcSettings.NY688.number_nonEEGchan= 4; 
%number of empty/useless channels in EDF set (e.g., Pleth, ...) - don't enter DC/ECG/TRIGGER etc here

subs_PreProcSettings.NY688.EKG_chanindex = []; %EKG missing
subs_PreProcSettings.NY688.ApplyPulseArtifactRejection = 0; %Because no EKG

subs_PreProcSettings.NY688.Trigger_chanindex= 145; 
%indexed position related to EDF all chan
subs_PreProcSettings.NY688.Triggerthresh = 0.5; 
%threshhold (%) for determining deflections of trigge channel

subs_PreProcSettings.NY688.index_FirstExpBlock = 2; 
%index of first block of real experiment (i.e., is there a training block in recording?)

subs_PreProcSettings.NY688.bandStopWidth = 2; 
%Width of band-stop filter in Hz

subs_PreProcSettings.NY688.numSTDthresh_forCAR = 2; 
%Threshhold (multiples of STD) to determine electrodes used for for Common Average Referencing 

subs_PreProcSettings.NY688.unsureChan_afterCAR_label ...
    = {'G14', 'G18', 'G21', 'G22', 'G27', 'G31','G44', 'G47', 'G48', ...
    'G50', 'G55','IO1','SO1'}; 
%Identified after CAR
%Note: Whole block 11 seems bad signal quality 

subs_PreProcSettings.NY688.unsureChan_afterCAR_index ...
    = [14, 18, 21, 22, 27, 31, 44, 47, 48, 50, 55, 73, 77];
subs_PreProcSettings.NY688.unsureBlocks = [11]; %whole block bad signal across electrodes

subs_PreProcSettings.NY688.IEDcontamElecs_index = [31]; %Elecs showing >10% of IED-contaminated samples (IED + fixed duration)
subs_PreProcSettings.NY688.IEDcontamElecs_label = {'G31'}; %Elecs showing >10% of IED-contaminated samples (IED + fixed duration)

subs_PreProcSettings.NY688.rejectedChan_label ...
    = {'G14', 'G31', 'G50', 'G55', 'IO1', 'SO1' }; 
subs_PreProcSettings.NY688.rejectedChan_index = [14, 31, 50, 55, 73, 77];
%Decided after ft_databrowser

subs_PreProcSettings.NY688.unsureTrials ...
    = [15, 25, 41, 47, 51, 52, 54, 65, 68, 69, 73, 76, 78, 82, 88, 95, ...
    97, 103, 105, 106, 107, 108, 120]; 
%Trial indices after verification of ECoG-Behav trial indices
subs_PreProcSettings.NY688.rejectedTrials ...
    = [15, 47, 68, 73, 78, 82, 88, 103, 107]; 
%Trial indices after verification of ECoG-Behav trial indices

subs_PreProcSettings.NY688.N1localizer_STDthresh.LF = 4;
subs_PreProcSettings.NY688.N1localizer_STDthresh.GammaAmp = 4;

%% NY704
subs_PreProcSettings.NY704.MNI2EDFmismatches.index_mismatchchannels = []';

subs_PreProcSettings.NY704.number_ECoGchan= 174; 
%number of channels in EDF set that contain ECoG recordings
subs_PreProcSettings.NY704.number_nonEEGchan= 4; 
%number of empty/useless channels in EDF set (e.g., Pleth, ...)  - don't enter DC/ECG/TRIGGER etc here

subs_PreProcSettings.NY704.EKG_chanindex = [176 175]; %Order switched due to polarity
subs_PreProcSettings.NY704.ApplyPulseArtifactRejection = 1; %
subs_PreProcSettings.NY704.EKG_thresh = [96 96 96 96 97 96 96 97 96 96 96 96]; 

subs_PreProcSettings.NY704.Trigger_chanindex= [177]; 
%indexed position related to EDF all chan
subs_PreProcSettings.NY704.Triggerthresh = 0.4; 
%threshhold (%) for determining deflections of trigge channel

subs_PreProcSettings.NY704.index_FirstExpBlock = 2; 
%index of first block of real experiment (i.e., no training)

subs_PreProcSettings.NY704.bandStopWidth = 2; 
%Wifth of band-pass filter in Hz

subs_PreProcSettings.NY704.numSTDthresh_forCAR = 2; 
%Threshhold (multiples of STD) to determine electrodes used for Common Average Referencing 
subs_PreProcSettings.NY704.unsureChan_afterCAR_label = {'G44', 'PF6'}; 
%G44 close to A1, so remains, however flatter than other chan
subs_PreProcSettings.NY704.unsureChan_afterCAR_index = [44, 76];

subs_PreProcSettings.NY704.IEDcontamElecs_index = []; %Elecs showing >10% of IED-contaminated samples (IED + fixed duration)
subs_PreProcSettings.NY704.IEDcontamElecs_label = {''}; %Elecs showing >10% of IED-contaminated samples (IED + fixed duration)

subs_PreProcSettings.NY704.rejectedChan_label = {'G1', 'G39', 'G64', 'PF6'};
%G1, G39, G64 show heartbeat signal, not properly cleaned in G39, G64, PF2,
%PF3, PF4, PF5, AP6, AP7, AP8,
subs_PreProcSettings.NY704.rejectedChan_index = [1, 39, 64, 76];
%G44 shows bursts in some trials, but lies in STG and is fine in the majority of trials

subs_PreProcSettings.NY704.unsureTrials = [8, 34, 91, 104, 119]; 
%Trial indices after verification of ECoG-Behav trial indices
subs_PreProcSettings.NY704.rejectedTrials = [];

subs_PreProcSettings.NY704.N1localizer_STDthresh.LF = 5;
subs_PreProcSettings.NY704.N1localizer_STDthresh.GammaAmp = 3;

%% NY708
subs_PreProcSettings.NY708.MNI2EDFmismatches.index_mismatchchannels = []';

subs_PreProcSettings.NY708.number_ECoGchan= 126; 
%number of channels in EDF set that contain ECoG recordings
subs_PreProcSettings.NY708.number_nonEEGchan= 4; 
%number of empty/useless channels in EDF set (e.g., Pleth, ...)  - don't enter DC/ECG/TRIGGER etc here

subs_PreProcSettings.NY708.EKG_chanindex = [127 128];
subs_PreProcSettings.NY708.ApplyPulseArtifactRejection = 0; 
%Testwise Pulse rejection revealed no pulse artifacts in iEEG data - thus no cleaning
subs_PreProcSettings.NY708.EKG_thresh = [96 97 98 97 97 98 98 97 98 97 96 97]; 

subs_PreProcSettings.NY708.Trigger_chanindex= [129]; 
%indexed position related to EDF all chan

subs_PreProcSettings.NY708.Triggerthresh = 0.25; 
%threshhold (%) for determining deflections of trigge channel
subs_PreProcSettings.NY708.index_FirstExpBlock = 1; 
%index of first block of real experiment (i.e., no training)

subs_PreProcSettings.NY708.bandStopWidth = 2; 
%Wifth of band-pass filter in Hz

subs_PreProcSettings.NY708.numSTDthresh_forCAR = 2; 
%Threshhold (multiples of STD) to determine electrodes used for for Common Average Referencing 
subs_PreProcSettings.NY708.unsureChan_afterCAR_label = {'P1', 'P2'};
subs_PreProcSettings.NY708.unsureChan_afterCAR_index = [87, 88];

subs_PreProcSettings.NY708.IEDcontamElecs_index = [46]; %Elecs showing >10% of IED-contaminated samples (IED + fixed duration)
subs_PreProcSettings.NY708.IEDcontamElecs_label = {'G46'}; %Elecs showing >10% of IED-contaminated samples (IED + fixed duration)

subs_PreProcSettings.NY708.rejectedChan_label = {'G46', 'P1', 'P2'};
subs_PreProcSettings.NY708.rejectedChan_index = [46, 87, 88];
%Notes: O2 very highfreq notes; AT&MT channels show high deflections throughout recording

subs_PreProcSettings.NY708.unsureTrials = []; %Trial indices after verification of ECoG-Behav trial indices
subs_PreProcSettings.NY708.rejectedTrials = [];
%Strong continuous alpha activity in temporal pole channels (AT & MT)


subs_PreProcSettings.NY708.N1localizer_STDthresh.LF = 4;
subs_PreProcSettings.NY708.N1localizer_STDthresh.GammaAmp = 4;


%% NY723
subs_PreProcSettings.NY723.MNI2EDFmismatches.index_mismatchchannels = []';

subs_PreProcSettings.NY723.number_ECoGchan= 162; 
%number of channels in EDF set that contain ECoG recordings
subs_PreProcSettings.NY723.number_nonEEGchan= 4; 
%number of empty/useless channels in EDF set (e.g., Pleth, ...)  - don't enter DC/ECG/TRIGGER etc here

subs_PreProcSettings.NY723.EKG_chanindex = [163 164]; 
subs_PreProcSettings.NY723.ApplyPulseArtifactRejection = 0; %EKG recording very bad quality
subs_PreProcSettings.NY723.EKG_thresh = []; 
%'ECG1', 'ECG2'

subs_PreProcSettings.NY723.Trigger_chanindex= [165]; 
%indexed position related to EDF all chan

subs_PreProcSettings.NY723.Triggerthresh = 0.7; 
%threshhold (%) for determining deflections of trigge channel
subs_PreProcSettings.NY723.index_FirstExpBlock = 1; 
%index of first block of real experiment (i.e., no training)

subs_PreProcSettings.NY723.bandStopFreqs = [60 120 180 240]; 
%band-passs freqs (power line noise)
subs_PreProcSettings.NY723.bandStopWidth = 2; 
%Wifth of band-pass filter in Hz

subs_PreProcSettings.NY723.numSTDthresh_forCAR = 2; 
%Threshhold (multiples of STD) to determine electrodes used for for Common Average Referencing 
subs_PreProcSettings.NY723.unsureChan_afterCAR_label = {'RIF11'}; %Very strong heartbeat artifact
subs_PreProcSettings.NY723.unsureChan_afterCAR_index = [63];
%Note: Freq peaks at ~220Hz in mulitpl channels
%Note: RPT4, RPT5, RPT6 high amplitude low freq oscilations in every trial,
%might be heartbeat

subs_PreProcSettings.NY723.IEDcontamElecs_index = []; %Elecs showing >10% of IED-contaminated samples (IED + fixed duration)
subs_PreProcSettings.NY723.IEDcontamElecs_label = {''}; %Elecs showing >10% of IED-contaminated samples (IED + fixed duration)

subs_PreProcSettings.NY723.rejectedChan_label = {'RIF11'};
subs_PreProcSettings.NY723.rejectedChan_index = [63];

subs_PreProcSettings.NY723.unsureTrials = [29, 69, 97, 103]; 
%Trial indices after verification of ECoG-Behav trial indices
subs_PreProcSettings.NY723.rejectedTrials = [];

subs_PreProcSettings.NY723.N1localizer_STDthresh.LF = 2;
subs_PreProcSettings.NY723.N1localizer_STDthresh.GammaAmp = 3;

%% NY742
%Note: raw data (both 512 and 2kHz) initially appear highly similar across
%channels due to high-level artefacts (presumably movement). Thus, on
%block-level, basically no differences in the recordings are visible across
%channels. This gets much better after CAR and on trial level.

subs_PreProcSettings.NY742.MNI2EDFmismatches.index_mismatchchannels = [137:172]';
subs_PreProcSettings.NY742.MNI2EDFmismatches.label_mismatchchannels = ... 
    [{'DAI1'},{'DAI2'},{'DAI3'},{'DAI4'},{'DAI5'},{'DAI6'},...
        {'DAMT1'},{'DAMT2'},{'DAMT3'},{'DAMT4'},{'DAMT5'},{'DAMT6'},...
        {'DMF1'},{'DMF2'},{'DMF3'},{'DMF4'},{'DMF5'},{'DMF6'},...
        {'DPI1'},{'DPI2'},{'DPI3'},{'DPI4'},{'DPI5'},{'DPI6'},...
        {'DPMT1'},{'DPMT2'},{'DPMT3'},{'DPMT4'},{'DPMT5'},{'DPMT6'},... 
        {'DSMA1'},{'DSMA2'},{'DSMA3'},{'DSMA4'},{'DSMA5'},{'DSMA6'}];
    %old labels in (adjusted) MNIfile
subs_PreProcSettings.NY742.MNI2EDFmismatches.correctedlabel_mismatchchannels = ... 
    [{'dAI1'},{'dAI2'},{'dAI3'},{'dAI4'},{'dAI5'},{'dAI6'},...
        {'dAMT1'},{'dAMT2'},{'dAMT3'},{'dAMT4'},{'dAMT5'},{'dAMT6'},...
        {'dMF1'},{'dMF2'},{'dMF3'},{'dMF4'},{'dMF5'},{'dMF6'},...
        {'dPI1'},{'dPI2'},{'dPI3'},{'dPI4'},{'dPI5'},{'dPI6'},...
        {'dPMT1'},{'dPMT2'},{'dPMT3'},{'dPMT4'},{'dPMT5'},{'dPMT6'},... 
        {'dSMA1'},{'dSMA2'},{'dSMA3'},{'dSMA4'},{'dSMA5'},{'dSMA6'}];
    %corresponding labels in EDF file
    
subs_PreProcSettings.NY742.number_ECoGchan = 172; 
%number of channels in EDF set that contain ECoG recordings
subs_PreProcSettings.NY742.number_nonEEGchan = 0; 
%number of empty/useless channels in EDF set (e.g., Pleth, ...)  - don't enter DC/ECG/TRIGGER etc here

subs_PreProcSettings.NY742.EKG_chanindex = [173 174]; 
subs_PreProcSettings.NY742.ApplyPulseArtifactRejection = 1;
subs_PreProcSettings.NY742.EKG_thresh = [98 98 98 98 98 98 98 98 98 98 98 98];

subs_PreProcSettings.NY742.Trigger_chanindex = [175]; 
%DC1 - indexed position related to EDF all chan

subs_PreProcSettings.NY742.Triggerthresh = 0.325; 
%threshhold (%) for determining deflections of trigger channel
subs_PreProcSettings.NY742.index_FirstExpBlock = 3; 
%index of first block of real experiment (i.e., no training)

subs_PreProcSettings.NY742.bandStopFreqs = [60 120 180 240]; 
%band-passs freqs (power line noise)
subs_PreProcSettings.NY742.bandStopWidth = 2; 
%Wifth of band-pass filter in Hz

subs_PreProcSettings.NY742.numSTDthresh_forCAR = 1.5; 
%Threshhold (multiples of STD) to determine electrodes used for for Common Average Referencing 
subs_PreProcSettings.NY742.unsureChan_afterCAR_label = {'G81', 'G94', 'G96', 'G105', 'G117', 'G126', 'IT1', 'IT4', 'IF4'}; 
subs_PreProcSettings.NY742.unsureChan_afterCAR_index = [81, 94, 96, 105, 117, 126, 129, 132, 136]; 

subs_PreProcSettings.NY742.IEDcontamElecs_index = []; %Elecs showing >10% of IED-contaminated samples (IED + fixed duration)
subs_PreProcSettings.NY742.IEDcontamElecs_label = {''}; %Elecs showing >10% of IED-contaminated samples (IED + fixed duration)

subs_PreProcSettings.NY742.rejectedChan_label = {'G81', 'G96', 'IT1', 'IT4', 'IF4'}; 
subs_PreProcSettings.NY742.rejectedChan_index = [81, 96, 129, 132, 136]; 

subs_PreProcSettings.NY742.unsureTrials = [6 8 10 26 32 46 50 74 84 89 91 95 116 119]; 
%Trial indices after verification of ECoG-Behav trial indices
subs_PreProcSettings.NY742.rejectedTrials = [8 10 26 46 50 89 95 116 119];

subs_PreProcSettings.NY742.N1localizer_STDthresh.LF = 4;
subs_PreProcSettings.NY742.N1localizer_STDthresh.GammaAmp = 4;

%% NY751
subs_PreProcSettings.NY751.MNI2EDFmismatches.index_mismatchchannels = []';
subs_PreProcSettings.NY751.MNI2EDFmismatches.label_mismatchchannels = [];
subs_PreProcSettings.NY751.MNI2EDFmismatches.correctedlabel_mismatchchannels = [];

subs_PreProcSettings.NY751.number_ECoGchan = 148; 
%number of channels in EDF set that contain ECoG recordings
subs_PreProcSettings.NY751.number_nonEEGchan = 4; 
%number of empty/useless channels in EDF set (e.g., Pleth, ...)  - don't enter DC/ECG/TRIGGER etc here
%NOTE: DC5 presumably audio recording - high activity during tone
%presentation 

subs_PreProcSettings.NY751.EKG_chanindex = [147 148]; 
%C148 present but only noise
subs_PreProcSettings.NY751.ApplyPulseArtifactRejection = 0; 
%Very bad EKG signal quality due to C148 being broken, thus no pulse correction
subs_PreProcSettings.NY751.EKG_thresh = [98 98 98 98 98 98 98 98 98 98 98 98]; 

subs_PreProcSettings.NY751.Trigger_chanindex= 149; 
%DC1, indexed position related to EDF all chan

subs_PreProcSettings.NY751.Triggerthresh = 0.5; %threshhold (%) for determining deflections of trigge channel
subs_PreProcSettings.NY751.index_FirstExpBlock = 3; %index of first block of real experiment (i.e., no training)

subs_PreProcSettings.NY751.bandStopFreqs = [60 120 180 240]; %band-passs freqs (power line noise)
subs_PreProcSettings.NY751.bandStopWidth = 2; %Wifth of band-pass filter in Hz

subs_PreProcSettings.NY751.numSTDthresh_forCAR = 2; 
%Threshhold (multiples of STD) to determine electrodes used for for Common Average Referencing 
subs_PreProcSettings.NY751.unsureChan_afterCAR_label = {'AT4', 'TO2', 'AF4'}; 
subs_PreProcSettings.NY751.unsureChan_afterCAR_index = [68 74 96]; 

subs_PreProcSettings.NY751.IEDcontamElecs_index = []; %Elecs showing >10% of IED-contaminated samples (IED + fixed duration)
subs_PreProcSettings.NY751.IEDcontamElecs_label = {''}; %Elecs showing >10% of IED-contaminated samples (IED + fixed duration)

subs_PreProcSettings.NY751.rejectedChan_label = {'AT4', 'TO2', 'AF4'};
subs_PreProcSettings.NY751.rejectedChan_index = [68 74 96];
%Note: AF4 has strong pulse artifact

subs_PreProcSettings.NY751.unsureTrials = []; 
%Trial indices after verification of ECoG-Behav trial indices
subs_PreProcSettings.NY751.rejectedTrials = []; 

subs_PreProcSettings.NY751.N1localizer_STDthresh.LF = 3;
subs_PreProcSettings.NY751.N1localizer_STDthresh.GammaAmp = 4;

%% NY787
subs_PreProcSettings.NY787.MNI2EDFmismatches.index_mismatchchannels = []';
subs_PreProcSettings.NY787.MNI2EDFmismatches.label_mismatchchannels = [];
subs_PreProcSettings.NY787.MNI2EDFmismatches.correctedlabel_mismatchchannels = [];

subs_PreProcSettings.NY787.number_ECoGchan = [172]; 
%number of channels in EDF set that contain ECoG recordings
subs_PreProcSettings.NY787.number_nonEEGchan = 0; 
%number of empty/useless channels in EDF set (e.g., Pleth, ...)  - don't enter DC/ECG/TRIGGER etc here

subs_PreProcSettings.NY787.EKG_chanindex = [174 173]; 
subs_PreProcSettings.NY787.ApplyPulseArtifactRejection = 0; 
subs_PreProcSettings.NY787.EKG_thresh = []; 
%Bad EKG signal quality in channel 174 - thus no pulse correction applied

subs_PreProcSettings.NY787.Trigger_chanindex= [175]; 
%DC1, indexed position related to EDF all chan

subs_PreProcSettings.NY787.Triggerthresh = 0.7; 
%threshhold for determining deflections of trigge channel
subs_PreProcSettings.NY787.index_FirstExpBlock = 1; 
%index of first block of real experiment (i.e., no training)

subs_PreProcSettings.NY787.bandStopFreqs = [60 120 180 240]; 
%band-passs freqs (power line noise)
subs_PreProcSettings.NY787.bandStopWidth = 2; 
%Width of band-pass filter in Hz

subs_PreProcSettings.NY787.numSTDthresh_forCAR = 2; 
%Threshhold (multiples of STD) to determine electrodes used for for Common Average Referencing 
subs_PreProcSettings.NY787.unsureChan_afterCAR_label = {'G9', 'G16', 'G63', 'O8'};
subs_PreProcSettings.NY787.unsureChan_afterCAR_index = [9, 16, 63, 80]; %

subs_PreProcSettings.NY787.IEDcontamElecs_index = []; %Elecs showing >10% of IED-contaminated samples (IED + fixed duration)
subs_PreProcSettings.NY787.IEDcontamElecs_label = {''}; %Elecs showing >10% of IED-contaminated samples (IED + fixed duration)

subs_PreProcSettings.NY787.rejectedChan_label = {'G9', 'G16', 'G18', 'G40', 'G48', 'G56', 'G63', 'G64', 'O8'};
subs_PreProcSettings.NY787.rejectedChan_index = [9 16 18 40 48 56 63 64 80];
%G40, G48, G50, G56 show strong IED contamination, under 10% but still
%rejected. G64 shows strong pulse artifact

subs_PreProcSettings.NY787.unsureTrials = []; 
%Trial indices after verification of ECoG-Behav trial indices
subs_PreProcSettings.NY787.rejectedTrials = []; 

subs_PreProcSettings.NY787.N1localizer_STDthresh.LF = [];
subs_PreProcSettings.NY787.N1localizer_STDthresh.GammaAmp = [];

%Note: Data looks quite irregular/spikey after cleaning (doesn't seem to be
%heartbeat or IEDs)

%% NY794
subs_PreProcSettings.NY794.MNI2EDFmismatches.index_mismatchchannels = []';
subs_PreProcSettings.NY794.MNI2EDFmismatches.label_mismatchchannels = [];
subs_PreProcSettings.NY794.MNI2EDFmismatches.correctedlabel_mismatchchannels = [];

subs_PreProcSettings.NY794.number_ECoGchan = [140]; 
%number of channels in EDF set that contain ECoG recordings
subs_PreProcSettings.NY794.number_nonEEGchan = 0; 
%number of empty/useless channels in EDF set (e.g., Pleth, ...)  - don't enter DC/ECG/TRIGGER etc here

subs_PreProcSettings.NY794.EKG_chanindex = [141 142];
%Channels labeled EKG but no EKG signal or very noisy - thus no pulse
%correction
subs_PreProcSettings.NY794.ApplyPulseArtifactRejection = 0; 
subs_PreProcSettings.NY794.EKG_thresh = []; 

subs_PreProcSettings.NY794.Trigger_chanindex= [143]; 
%DC1, indexed position related to EDF all chan

subs_PreProcSettings.NY794.Triggerthresh = 0.7; 
%threshhold (%) for determining deflections of trigge channel
subs_PreProcSettings.NY794.index_FirstExpBlock = 2; 
%index of first block of real experiment (i.e., no training)

subs_PreProcSettings.NY794.bandStopFreqs = [60 120 180 240]; 
%band-passs freqs (power line noise)
subs_PreProcSettings.NY794.bandStopWidth = 2; 
%Width of band-pass filter in Hz

subs_PreProcSettings.NY794.numSTDthresh_forCAR = 2; 
%Threshhold (multiples of STD) to determine electrodes used for for Common Average Referencing 
subs_PreProcSettings.NY794.unsureChan_afterCAR_label = ...
    {'G5', 'G6', 'G14', 'G24', 'G32', 'G39', 'G44', 'G48', 'LF5', 'LAT3'}; %
subs_PreProcSettings.NY794.unsureChan_afterCAR_index = [5 6 14 24 32 39 44 48 73 109]; %
%Most bad channels due to largepulse artifacts

subs_PreProcSettings.NY794.IEDcontamElecs_index = []; %Elecs showing >10% of IED-contaminated samples (IED + fixed duration)
subs_PreProcSettings.NY794.IEDcontamElecs_label = {''}; %Elecs showing >10% of IED-contaminated samples (IED + fixed duration)

subs_PreProcSettings.NY794.rejectedChan_label = {'G14', 'G32', 'LF5', 'LAT3'};
subs_PreProcSettings.NY794.rejectedChan_index = [14 32 73 109];

subs_PreProcSettings.NY794.unsureTrials = [8 16 24 25 32 35 37 43 46 50 55 63 66 73 76 80 86 90 97 104 107 108 113]; 
%Trial indices after verification of ECoG-Behav trial indices
subs_PreProcSettings.NY794.rejectedTrials = [8 16 24 25 43 46 55 63 66 76 107]; 

subs_PreProcSettings.NY794.N1localizer_STDthresh.LF = [];
subs_PreProcSettings.NY794.N1localizer_STDthresh.GammaAmp = [];

%note: G24, G39, G46, G48 show strong EKG artefacts, but are in/around STG
%LP6, LF4-LF6 show strong EKG artefacts
%% NY798
subs_PreProcSettings.NY798.MNI2EDFmismatches.index_mismatchchannels = []';
subs_PreProcSettings.NY798.MNI2EDFmismatches.label_mismatchchannels = [];
subs_PreProcSettings.NY798.MNI2EDFmismatches.correctedlabel_mismatchchannels = [];

subs_PreProcSettings.NY798.number_ECoGchan = 198; 
%number of channels in EDF set that contain ECoG recordings
subs_PreProcSettings.NY798.number_nonEEGchan = 0; 
%number of empty/useless channels in EDF set (e.g., Pleth, ...)  - don't enter DC/ECG/TRIGGER etc here

subs_PreProcSettings.NY798.EKG_chanindex = [199 200];
%EKG 1 (199) bad - oscillaotry artifact
subs_PreProcSettings.NY798.ApplyPulseArtifactRejection = 0;
%Noisy EKG in 4/12 blocks (2,5,6,10), no correction applied
subs_PreProcSettings.NY798.EKG_thresh = [97 97 97 97 94 98 97 96 97 97 97 98]; 

subs_PreProcSettings.NY798.Trigger_chanindex= [202]; 
%DC2, indexed position related to EDF all chan

subs_PreProcSettings.NY798.Triggerthresh = 0.5; 
%threshhold (%) for determining deflections of trigge channel
subs_PreProcSettings.NY798.index_FirstExpBlock = 2; 
%index of first block of real experiment (i.e., no training)

subs_PreProcSettings.NY798.bandStopFreqs = [60 120 180 240]; 
%band-passs freqs (power line noise)
subs_PreProcSettings.NY798.bandStopWidth = 2; 
%Width of band-pass filter in Hz

subs_PreProcSettings.NY798.numSTDthresh_forCAR = 2; 
%Threshhold (multiples of STD) to determine electrodes used for for Common Average Referencing 
subs_PreProcSettings.NY798.unsureChan_afterCAR_label = {'G64', 'G65', 'G91', 'G94', 'G107', 'PIT4', 'O8'}; %
subs_PreProcSettings.NY798.unsureChan_afterCAR_index = [64 65 91 94 107 136 144]; 

subs_PreProcSettings.NY798.IEDcontamElecs_index = []; %Elecs showing >10% of IED-contaminated samples (IED + fixed duration)
subs_PreProcSettings.NY798.IEDcontamElecs_label = {''}; %Elecs showing >10% of IED-contaminated samples (IED + fixed duration)
%Many electrodes close to 10% (>5%). Automatic detection suboptimal,
%however, as many IEDs are not found and added duration often too short. 

subs_PreProcSettings.NY798.unsureTrials = []; %IEDs in almost all trials

%Keep heavily IED-contaminated Channels (mostly temporal) and use auto+manual cleaning 
subs_PreProcSettings.NY798.rejectedChan_label = {'G65', 'G91', 'G94', 'PIT4', 'O8', 'IF2'};
subs_PreProcSettings.NY798.rejectedChan_index = [65 91 94 136 144 166];
   subs_PreProcSettings.NY798.rejectedTrials = ...
       [3 20 31 32 35 37 39 42 48 53 58 72 73 75 76 83 84 88 90 94 97 ...
       98 102 104 105 109 113 115 117]; 
   
%    %Remove heavily IED-contaminated Channels (mostly temporal)
%    subs_PreProcSettings.NY798.IEDcontaminatedChan_index = ...
%        [33 34 36 41:46 49 50 51 52 53 54 55 57 58 62 63 ...
%        112 116:119 124:126 128 130:132 134:136 142];
%    subs_PreProcSettings.NY798.rejectedChan_index = ...
%        [33 34 36 41:46 49 50 51 52 53 54 55 57 58 62 63 65 91 94 ...
%        112 116:119 124:126 128 130:132 134:136 142 144 166];
%    subs_PreProcSettings.NY798.rejectedChan_label = ...
%        {'G33','G34','G36','G41','G42','G43','G44','G45','G46','G49',...
%        'G50','G51','G52','G53','G54','G55','G57','G58','G62','G63',...
%        'G65', 'G91', 'G94', 'G112','G116','G117','G118','G119','G124','G125',...
%        'G126','G128','AIT2','AIT3','AIT4','PIT2','PIT3','PIT4','O6','O8','IF2'};
%    subs_PreProcSettings.NY798.rejectedTrials = [15 34 90];


subs_PreProcSettings.NY798.N1localizer_STDthresh.LF = [];
subs_PreProcSettings.NY798.N1localizer_STDthresh.GammaAmp = [];

%Note: Subject very problematic: Many IEDs throughout the trials in ~35
%electrodes, all lying in inferior frontal and superior temporal regions.
%~43 trials per TD after preprocessing and ~32 trials per TD after
%NaN-rejection and filtering (in these trials ~5-10 with averaged values in
%relevant TOI). Unclear whast is best: Keep IED-elecs and deal with very
%low data or reject IED-elecs and loose most of ROI...

