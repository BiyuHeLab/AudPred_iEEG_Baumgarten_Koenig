function param = sfa_expt4_getParams(window, freq, isBB, isMEG)

if ismac
addpath(genpath('/Users/ecog/Documents/MATLAB/sfa_expt4/')) %iEEG Mac path
else
% addpath(genpath('C:/Users/localbiyu/Desktop/Thomas/Experiments/Auditory History Tracking/sfa_expt4_Adeen_July18/')); %Path experimental laptop
% addpath(genpath('\\gogo.sb.nyumc.org\data\gogodisk4\thomas\AuditoryPrediction\iEEG\Experiment_Versions\sfa_expt4_Adeen_Current\')); %Path desktop/gogo
addpath('\\gogo.sb.nyumc.org\data\gogodisk4\thomas\NaturalisticAuditorySequences_ToneDuration(NAS_TD)\ECoG\Experiment_Versions\sfa_expt4_Adeen_Current\stim_creation\');
end


%% stimulus parameters
load series_selection_scale2Hz_jenn %load in stimuli (i.e., sound series, beta levels, indices)
%loads series_selection_scale2Hz_jenn.mat from subfolder stim_creation
param.f_sample_inHz    = freq;%audio sample rate from PsychPortAudio

param.nTonesInSeries   = length(series_Hz{1}(1,:)) + 1; %34 tones in total
param.toneDur_inSecs   = [0.2 0.4]; %Tone Duration
param.seriesDur_inSecs = param.toneDur_inSecs * param.nTonesInSeries;

param.paBuffer_inSecs       = .3; % time between loading stim in PA buffer and playing it
param.fixationBuffer_inSecs = .5; %.3; % time between offset of last trial's Feedback (FB) and onset of this trial's fixation dot

param.FBoffset2stimOnset_inSecs = 1.2; %time between fixation buffer and auditory stim onset
param.trialStartBuffer_inSecs = param.FBoffset2stimOnset_inSecs - param.paBuffer_inSecs; %time between start of a trial and auditory stim of that trial

%% display parameters

param.BGcolor   = 0; %255; %127;
param.font      = 'Helvetica';
param.fontColor = 255; 0; %0;
param.textSizeTrial = 30; %50
param.textSizeBlock = 30; %50
param.textSize  = param.textSizeTrial;
param.sx        = 200;
param.vSpacing  = 1.4;

% fixation point
param.fixationRadius_inDegrees = .5;
if isMEG
    param.distFromScreen_inCm = 60;
else
    param.distFromScreen_inCm = 60;
end

[midW, midH] = getScreenMidpoint(window);
fixationRadius_inPixels = degrees2pixels(param.fixationRadius_inDegrees, param.distFromScreen_inCm);
param.fixationRect = [midW - fixationRadius_inPixels, ...
                      midH - fixationRadius_inPixels, ...
                      midW + fixationRadius_inPixels, ...
                      midH + fixationRadius_inPixels];

%% response

% input device
if IsOSX %Is Mac OS X
    param.keyboardNumber = getKeyboardNumber;
else
    param.keyboardNumber = [];
end

param.isBB = isBB;

% response time allowed for each response before triggering next trial automatically
param.respDur_inSecs = 5;

% input keys
param.exitKey      = 'ESCAPE';

if isBB
    param.readyKey = 'Return';
    param.r1_validKeys = {'1!' '2@' '3#' '4$' '5%' param.exitKey}; %for response 1
    param.r2_validKeys = {'2@' '3#' '4$' param.exitKey}; %for response 2
    
else
%     param.readyKey = [];
    param.readyKey = 'Return';
    param.r1_validKeys = {'1!' '2@' '3#' '4$' '5%' param.exitKey}; %for response 1
    param.r2_validKeys = {'8*' '9(' '0)' param.exitKey};%for response 2
    
end

param.r3_validKeys = param.r1_validKeys;  %for response 3 (irrelevant,no resp3)





%% block structure

param.nTrialsPerBlock     = 10; % should be divisible by 120
%TJB: strange, in outline 30 trials per block
param.breakDur_inSecs     = Inf; % if Inf, break is ended by subject
param.breakWarning_inSecs = 10;  % point at which time left starts to count down
param.breakBuffer_inSecs  = 1;   % buffer b/t end of break and start of next trial


%% triggers
%This is all for MEG only, ECoG triggers are auditory
% potential port address #3
param.io_address = 16376;
param.io_address = hex2dec('3000');

% potential port address #2
% param.io_address = 16376;
% param.io_address = hex2dec('3FF8');

param.triggerDur_inSecs = .02;

param.trig_start = 255;
param.trig_stim  = [11 12 13 14 15];
param.trig_r1    = [21 22 23 24 25 29];
param.trig_r2    = [31 32 33 34 35 39];
param.trig_r3    = [41 42 43 44 45 49];

param.isMEG = isMEG;