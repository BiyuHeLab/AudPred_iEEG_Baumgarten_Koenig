clear


%% startup
%Determine if Lab-StimLaptop (Win) or ECoG-StimLaptop (Mac)
valid = 0;
while ~valid
    isMac = input('Is the experiment run on the iEEG laptop (Mac)? (y/n) ...   ','s');
    
    valid = any( strcmp(isMac, {'y' 'n'}) );
    if ~valid
        fprintf('invalid input\n\n');
        
    else
        switch isMac
            case 'y', isMac = 1;
            case 'n', isMac = 0;
        end
        
        if isMac
            fprintf('Experiment is run on the iEEG laptop (Mac).\n\n');            
        else
            fprintf('Experiment is run on labs own laptop (Win).\n\n');            
        end
    end
end

% set up directory structure
if isMac
addpath('/Users/ecog/Documents/MATLAB/sfa_expt4/') %Path Laptop Adeen
cd('/Users/ecog/Documents/MATLAB/sfa_expt4/') %Path Laptop Adeen
addpath([pwd '/functions/']);

else   
% addpath('\\gogo.sb.nyumc.org\data\gogodisk4\thomas\AuditoryPrediction\iEEG\Experiment_Versions\sfa_expt4_Adeen_Current\');
% cd('\\gogo.sb.nyumc.org\data\gogodisk4\thomas\AuditoryPrediction\iEEG\Experiment_Versions\sfa_expt4_Adeen_Current\');
addpath('\\gogo.sb.nyumc.org\data\gogodisk4\thomas\NaturalisticAuditorySequences_ToneDuration(NAS_TD)\ECoG\Experiment_Versions\sfa_expt4_Adeen_Current\');
cd '\\gogo.sb.nyumc.org\data\gogodisk4\thomas\NaturalisticAuditorySequences_ToneDuration(NAS_TD)\ECoG\Experiment_Versions\sfa_expt4_Adeen_Current\';
addpath([pwd '\functions\']);
end

expName = 'sfa_expt4';

isPractice = 1; %determine training session

startup = startExperiment(expName, isPractice); %Give out startup file with dir, filename, start time

%% check for triggers and button box
isMEG = 0; %Specify that this is no MEG session
isBB = 0; %Specify that input will be given via keyboard 

WaitSecs(2); %psychtoolbox command
% pause(2) %matlab equivalent

%% launch PTB
% initialize PsychPortAudio
fprintf('Initializing PsychPortAudio...\n');
fprintf('------------------------------\n\n');

AssertOpenGL;

lowLatency = 1;
InitializePsychSound( lowLatency );

deviceID       = [];
mode           = 1;  % sound playback only
timingPriority = 2;  % take aggressive control of audio device
freq           = []; % allow PsychPortAudio to set the determined optimal value, depending on other settings 
nChannels      = 2;
bufferSize     = 128*2^3;
pahandle       = PsychPortAudio('Open', deviceID, mode, timingPriority, freq, nChannels, bufferSize);

s    = PsychPortAudio('GetStatus', pahandle)
freq = s.SampleRate; %44100

%% Start training?
fprintf('\n\nReady to go?\n');
proceed = 0;
while ~proceed
    control = input('y --> proceed\nd --> debug mode\nq --> quit\n\n...   ','s');
    
    valid = any( strcmp( control, {'y' 'd' 'q'} ) );
    if ~valid
        fprintf('invalid input\n\n');
        
    else
        switch control
            case 'y', proceed = 1;
            case 'q', PsychPortAudio('Close', pahandle); return;
            case 'd', keyboard;
        end
    end
end

% launch PTB window
[window keyboardNumber] = openScreen;
%Proxy: window = 10, keyboardNumber = []

%% get experiment parameters
%I.e., get stimulus parameters, display parameters, response
%parameters,block structure, trigger parameters
param = sfa_expt4_getParamsTraining(window, freq, isBB, isMEG);
param.restarting = 0; %new start, no restart

Screen('FillRect',  window, param.BGcolor);
Screen('TextSize',  window, param.textSize);
Screen('TextColor', window, param.fontColor);
Screen('TextFont',  window, param.font);

DrawFormattedText(window,'getting ready...','center','center');
Screen('Flip',window);

%% create stimuli for the training 
    stim  = sfa_expt4_makeTrainingStim(param, pahandle);

%% run the experiment
sx='center';
sy='center';
    
    %On-screen instructions (screen 1)
        InstructionText = ['Thank you for participating in our experiment.\n\n' ...
                     'You will listen to brief sequences of multiple auditory tones.\n' ...
                     'These tones will differ in pitch and vary throughout the tone sequence.\n\n' ...
                     'Each tone sequence has an overall trend in how the tone pitches change over time. \n' ...
                     'This trend can be either gradually falling, remain stable, or gradually rising, \n' ...
                     'and it determines which tone pitch is likely to come next. \n' ...
                     'For example, if you hear a falling trend, \n' ...
                     'it is most probable that the next tones will have a low pitch.' ...
                     '\n\n' ...
                     '\n\n' ...
                     'Press any key to continue.'];    

    wrect = Screen('Rect',window);
    sx = round(wrect(3)*.1);
    Screen('TextSize',  window, param.textSizeInstructions);
    DrawFormattedText(window, InstructionText, sx, sy, param.fontColor,[],[],[],param.vSpacingInstructions);
    Screen('Flip',window);
    Screen('TextSize',  window, param.textSizeInstructions);

    WaitSecs(.5);
    KbWait;
    
    %On-screen instructions (screen 2)
        InstructionText = ['Your task is to tell us how well the final tone pitch fits in with the preceding tone sequence. \n\n' ...
                     'The final tone pitch can either fit in well with the preceding tone sequence and continue the trend, \n' ...
                     'or it can not fit in well with the preceding tone sequence and break the trend.' ...
                     '\n\n' ...
                     '\n\n' ...
                     'Press any key to continue.'];    

    wrect = Screen('Rect',window);
    sx = round(wrect(3)*.1);
    Screen('TextSize',  window, param.textSizeInstructions);
    DrawFormattedText(window, InstructionText, sx, sy, param.fontColor,[],[],[],param.vSpacingInstructions);
    Screen('Flip',window);
    Screen('TextSize',  window, param.textSizeInstructions);

    WaitSecs(.5);
    KbWait;

    %On-screen instructions (screen 3)
        InstructionText = ['To tell us how well the final tone pitch fits in with the preceding tone sequence, \n' ...
                     'you will rate how likely it is that the present tone sequence culminates with the presented final tone. \n\n' ...
                     'The rating is done on a scale from 1 to 5, where: \n' ...
                     '1 means that the final tone was very unlikely. It did not fit in at all with the preceding tone sequence. \n' ... ...
                     '5 means that the final tone was very likely. It fit in well with the preceding tone sequence. \n\n' ...
                     'Please try to use all 5 options of the rating scale. Do not limit yourself to using only 2 or 3 responses.' ...
                     '\n\n' ...
                     '\n\n' ...
                     'Press any key to continue.'];    

    wrect = Screen('Rect',window);
    sx = round(wrect(3)*.1);
    Screen('TextSize',  window, param.textSizeInstructions);
    DrawFormattedText(window, InstructionText, sx, sy, param.fontColor,[],[],[],param.vSpacingInstructions);
    Screen('Flip',window);
    Screen('TextSize',  window, param.textSizeInstructions);

    WaitSecs(.5);
    KbWait;

    %On-screen instructions (screen 4)
        InstructionText = ['For your rating, you have up to 5 seconds to respond. \n' ...
                     'If you are completely unsure of what to respond, it is OK to provide your best guess. \n\n' ...                 
                     'During the task, please look directly on the dot when it is on the screen. \n' ...
                     'Please move as little as possible and please do not close your eyes. \n\n' ...
                     'At certain times during the task, you will have the option to take break periods where you can rest. \n' ...
                     'During the break periods, it is OK to close your eyes and move. \n\n' ...
                     'We will start the with a short practice session, so that you get familiar with the task. \n' ...
                     'If anything is unclear, please do not hesitate to ask the experimenter.' ...
                     '\n\n' ...
                     '\n\n' ...
                     'Press any key to continue.'];    

    wrect = Screen('Rect',window);
    sx = round(wrect(3)*.1);
    Screen('TextSize',  window, param.textSizeInstructions);
    DrawFormattedText(window, InstructionText, sx, sy, param.fontColor,[],[],[],param.vSpacingInstructions);
    Screen('Flip',window);
    Screen('TextSize',  window, param.textSizeInstructions);

    WaitSecs(.5);
    KbWait;    
    
    %On-screen instructions (screen )
    DrawFormattedText(window,'If you have no more questions, we will now begin with the practice session.\n\nPress any key to continue.','center','center');
    Screen('Flip',window);

    WaitSecs(.5);
    KbWait;
    
param.isPractice = isPractice;
    param.nTrialsPerBlock = 9;
    stim.nTrials = 18;
    
[data exitNow] = sfa_expt4_runBlockTraining(window, param, stim);

Screen('CloseAll');
PsychPortAudio('Close', pahandle);

startup.endTime    = clock;
startup.endTimeStr = datestr(startup.endTime);
stim.soundwave  = [];
save(startup.dataFile, 'startup', 'param', 'stim', 'data');

