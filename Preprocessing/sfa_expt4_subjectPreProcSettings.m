function subs_PreProcSettings = sfa_expt4_subjectPreProcSettings

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NY688
subs_PreProcSettings.NY688.MNI2EDFmismatches.index_mismatchchannels = [71:78]';
subs_PreProcSettings.NY688.MNI2EDFmismatches.label_mismatchchannels = [{'T1'}, {'T2'}, {'T3'},{'T4'},{'T5'},{'T6'},{'T7'},{'T8'}];
subs_PreProcSettings.NY688.MNI2EDFmismatches.correctedlabel_mismatchchannels = [{'TT1'}, {'TT2'}, {'TT3'},{'TT4'},{'TT5'},{'TT6'},{'TT7'},{'TT8'}];
subs_PreProcSettings.NY688.number_ECoGchan= 144; %number of channels in EDF set that contain ECoG recordings (no DC, Trigger,...)
subs_PreProcSettings.NY688.number_nonEEGchan= 4; %number of empty/useless channels in EDF set (e.g., DC, Pleth, ...)
subs_PreProcSettings.NY688.EKG_chanindex = []; %EKG missing
subs_PreProcSettings.NY688.Trigger_chanindex= 145; %indexed position related to EDF all chan
subs_PreProcSettings.NY688.Triggerthresh = 0.5; %threshhold (%) for determining deflections of trigge channel

subs_PreProcSettings.NY688.index_FirstExpBlock = 2; %index of first block of real experiment (i.e., is there a training block in recording?)

subs_PreProcSettings.NY688.bandStopFreqs = [60 120 180 240]; %band-passs freqs (power line noise+harmonics)
subs_PreProcSettings.NY688.bandStopWidth = 2; %Wifth of band-pass filter in Hz
subs_PreProcSettings.NY688.HighPassFreq = 0.01;

subs_PreProcSettings.NY688.numSTDthresh_forCAR = 2; %Threshhold (multiples of STD) to determine electrodes used for for Common Average Referencing 
subs_PreProcSettings.NY688.unsureChan_afterCAR_label = {'G14', 'G18', 'G21', 'G27', 'G31','G44', 'G47', 'G48', 'G50', 'G55','IO1','SO1'};
subs_PreProcSettings.NY688.unsureChan_afterCAR_index = [14, 18, 21, 27, 31, 44, 47, 48, 50, 55, 73, 77];
subs_PreProcSettings.NY688.unsureBlocks = [11];

subs_PreProcSettings.NY688.rejectedChan_label = {'G14','G31','G50', 'G55', 'IO1', 'SO1' };
subs_PreProcSettings.NY688.rejectedChan_index = [14, 31, 50, 55, 73, 77];
subs_PreProcSettings.NY688.rejectedTrials = [15, 47, 68, 78, 82, 88, 97, 103];
subs_PreProcSettings.NY688.unsureTrials = [54, 69, 70, 73, 76, 95, 107, 108, 109, 120];

subs_PreProcSettings.NY688.goodtrials_afterCAR_indexECoG = 1:120; 
subs_PreProcSettings.NY688.goodtrials_afterCAR_indexBehav = 1:120;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NY704
subs_PreProcSettings.NY704.MNI2EDFmismatches.index_mismatchchannels = []';

subs_PreProcSettings.NY704.number_ECoGchan= 174; %number of channels in EDF set that contain ECoG recordings (no DC, Trigger,...)
subs_PreProcSettings.NY704.number_nonEEGchan= 4; %number of empty/useless channels in EDF set (e.g., DC, Pleth, ...) 

subs_PreProcSettings.NY704.EKG_chanindex = [175 176];
subs_PreProcSettings.NY704.Trigger_chanindex= [177]; %indexed position related to EDF all chan
subs_PreProcSettings.NY704.Triggerthresh = 0.4; %threshhold (%) for determining deflections of trigge channel

subs_PreProcSettings.NY704.index_FirstExpBlock = 2; %index of first block of real experiment (i.e., no training)

subs_PreProcSettings.NY704.bandStopFreqs = [60 120 180 240]; %band-passs freqs (power line noise)
subs_PreProcSettings.NY704.bandStopWidth = 2; %Wifth of band-pass filter in Hz
subs_PreProcSettings.NY704.HighPassFreq = 0.01;

subs_PreProcSettings.NY704.numSTDthresh_forCAR = 2; %Threshhold (multiples of STD) to determine electrodes used for for Common Average Referencing 
subs_PreProcSettings.NY704.unsureChan_afterCAR_label = {'G44', 'PF6'};
subs_PreProcSettings.NY704.unsureChan_afterCAR_index = [44, 76];


subs_PreProcSettings.NY704.rejectedChan_label = {'G1', 'G39', 'G44', 'PF6', 'G64'};%G1, G39, G64 show slow oscillatroy signal
subs_PreProcSettings.NY704.rejectedChan_index = [1, 39, 44, 64, 76];
subs_PreProcSettings.NY704.rejectedTrials = [];

subs_PreProcSettings.NY704.goodtrials_afterCAR_indexECoG = [1:48,50,54:123]; %trial 49 in block 5 missed and repeated 3 times - taken missed trial and repetitions out
subs_PreProcSettings.NY704.goodtrials_afterCAR_indexBehav = [1:48,50:120]; %trial 49 in block 5 missed and repeated 3 times - taken missed trial out

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NY708
subs_PreProcSettings.NY708.MNI2EDFmismatches.index_mismatchchannels = []';

subs_PreProcSettings.NY708.number_ECoGchan= 126; %number of channels in EDF set that contain ECoG recordings (no DC, Trigger,...)
subs_PreProcSettings.NY708.number_nonEEGchan= 4; %number of empty/useless channels in EDF set (e.g., DC, Pleth, ...) 
subs_PreProcSettings.NY708.EKG_chanindex = [127 128];
subs_PreProcSettings.NY708.Trigger_chanindex= [129]; %indexed position related to EDF all chan
subs_PreProcSettings.NY708.Triggerthresh = 0.25; %threshhold (%) for determining deflections of trigge channel
subs_PreProcSettings.NY708.index_FirstExpBlock = 1; %index of first block of real experiment (i.e., no training)

subs_PreProcSettings.NY708.bandStopFreqs = [60 120 180 240]; %band-passs freqs (power line noise)
subs_PreProcSettings.NY708.bandStopWidth = 2; %Wifth of band-pass filter in Hz
subs_PreProcSettings.NY708.HighPassFreq = 0.01;

subs_PreProcSettings.NY708.numSTDthresh_forCAR = 3; %Threshhold (multiples of STD) to determine electrodes used for for Common Average Referencing 
subs_PreProcSettings.NY708.unsureChan_afterCAR_label = {'G44', 'P1', 'P2', };
subs_PreProcSettings.NY708.unsureChan_afterCAR_index = [44, 87, 88];

subs_PreProcSettings.NY708.goodtrials_afterCAR_indexECoG = [1:11, 14:20, 23:122]; %trial 12+13 in block 2 missed and repeated 1 time each - taken missed trial and repetitions out
subs_PreProcSettings.NY708.goodtrials_afterCAR_indexBehav = [1:11, 14:120]; %trial 12+13 in block 2 missed and repeated 1 time each - taken missed trial and repetitions out

subs_PreProcSettings.NY708.rejectedChan_label = {'G44', 'P1', 'P2', };
subs_PreProcSettings.NY708.rejectedChan_index = [44, 87, 88];
%Notes: O2 very highfreq notes; AT&MT channels show igh deflection
%trhoughout the recording; 
subs_PreProcSettings.NY708.rejectedTrials = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NY723
subs_PreProcSettings.NY723.MNI2EDFmismatches.index_mismatchchannels = []';

subs_PreProcSettings.NY723.number_ECoGchan= 162; %number of channels in EDF set that contain ECoG recordings (no DC, Trigger,...)
subs_PreProcSettings.NY723.number_nonEEGchan= 8; %number of empty/useless channels in EDF set (e.g., DC, Pleth, ...) 

subs_PreProcSettings.NY723.EKG_chanindex = [163 164];
subs_PreProcSettings.NY723.Trigger_chanindex= [165]; %indexed position related to EDF all chan
subs_PreProcSettings.NY723.Triggerthresh = 0.5; %threshhold (%) for determining deflections of trigge channel

subs_PreProcSettings.NY723.index_FirstExpBlock = 1; %index of first block of real experiment (i.e., no training)

subs_PreProcSettings.NY723.bandStopFreqs = [60 120 180 240]; %band-passs freqs (power line noise)
subs_PreProcSettings.NY723.bandStopWidth = 2; %Wifth of band-pass filter in Hz
subs_PreProcSettings.NY723.HighPassFreq = 0.01;

subs_PreProcSettings.NY723.numSTDthresh_forCAR = 2; %Threshhold (multiples of STD) to determine electrodes used for for Common Average Referencing 
subs_PreProcSettings.NY723.unsureChan_afterCAR_label = {};
subs_PreProcSettings.NY723.unsureChan_afterCAR_index = [];
%Note: Freq peaks at ~220Hz in mulitpl channels
%Note: RPT4, RPT5, RPT6 high amplitude low freq oscilations in every trial,
%however too irregular for heartbeat
subs_PreProcSettings.NY723.goodtrials_afterCAR_indexECoG = [1:120]; %trial 12+13 in block 2 missed and repeated 1 time each - taken missed trial and repetitions out
subs_PreProcSettings.NY723.goodtrials_afterCAR_indexBehav = [1:120]; %trial 12+13 in block 2 missed and repeated 1 time each - taken missed trial and repetitions out

subs_PreProcSettings.NY723.rejectedChan_label = {};
subs_PreProcSettings.NY723.rejectedChan_index = [];
subs_PreProcSettings.NY723.rejectedTrials = []; %75,97,98,109
