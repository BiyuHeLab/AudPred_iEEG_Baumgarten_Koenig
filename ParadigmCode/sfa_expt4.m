clear

%% restart?
% if restarting, 
% save data in online backup to new folder first
restarting = 0;


%% startup
nPracticeTrials = 0;

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
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Experiment_Versions/sfa_expt4_Adeen_Current/');
cd '/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Experiment_Versions/sfa_expt4_Adeen_Current/';
addpath([pwd '/functions/']);

% addpath('\\gogo.sb.nyumc.org\data\gogodisk4\thomas\NaturalisticAuditorySequences_ToneDuration(NAS_TD)\ECoG\Experiment_Versions\sfa_expt4_Adeen_Current\');
% cd '\\gogo.sb.nyumc.org\data\gogodisk4\thomas\NaturalisticAuditorySequences_ToneDuration(NAS_TD)\ECoG\Experiment_Versions\sfa_expt4_Adeen_Current\';
% addpath('\\gogo.sb.nyumc.org\data\gogodisk4\thomas\AuditoryPrediction\iEEG\Experiment_Versions\sfa_expt4_Adeen_Current\');
% cd('\\gogo.sb.nyumc.org\data\gogodisk4\thomas\AuditoryPrediction\iEEG\Experiment_Versions\sfa_expt4_Adeen_Current\');
% addpath([pwd '\functions\']);
end

expName = 'sfa_expt4';

% %Determine if practice session or experiment session
% valid = 0;
% while ~valid
%     isPractice = input('Is this a practice session? (y/n) ...   ','s');
%     
%     valid = any( strcmp(isPractice, {'y' 'n'}) );
%     if ~valid
%         fprintf('invalid input\n\n');
%         
%     else
%         switch isPractice
%             case 'y', isPractice = 1;
%             case 'n', isPractice = 0;
%         end
%         
%         if isPractice
%             fprintf('This is a practice session of %d trials.\n\n', nPracticeTrials);            
%         else
%             fprintf('This is a full experimental session.\n\n');            
%         end
%     end
% end
isPractice = 0;
startup = startExperiment(expName, isPractice); %Give out startup file with dir, filename, start time

%% check for triggers and button box
%Determine if MEG or ECog session, which in turn changes trigger
%functionality
% valid = 0;
% while ~valid
%     isMEG = input('Is this an MEG session? (y/n) ...   ','s');
%     
%     valid = any( strcmp(isMEG, {'y' 'n'}) );
%     if ~valid
%         fprintf('invalid input\n\n');
%         
%     else
%         switch isMEG
%             case 'y', isMEG = 1;
%             case 'n', isMEG = 0;
%         end
%         
%         if isMEG
%             fprintf('This is an MEG session. Trigger functionality is turned on.\n\n');            
%         else
%             fprintf('This is not an MEG session. Trigger functionality is shut off.\n\n');            
%         end
%     end
% end
% 
% %Determine input method (keyboard vs. buttonbox)
% valid = 0;
% while ~valid
%     isBB = input('Are we using button boxes for input? (y/n) ...   ','s');
%     
%     valid = any( strcmp(isBB, {'y' 'n'}) );
%     if ~valid
%         fprintf('invalid input\n\n');
%         
%     else
%         switch isBB
%             case 'y', isBB = 1;
%             case 'n', isBB = 0;
%         end        
%         
%         if isBB
%             fprintf('Input will be collected through the button boxes.\n\n');
%         else
%             fprintf('Input will be collected through the keyboard.\n\n');
%         end
%     end
% end

isMEG = 0;
isBB = 0;
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
freq = s.SampleRate; 

if isMEG
    fprintf('Initializing port I/O...\n');
    fprintf('-------------------------\n\n');

    config_io; %installs a kernel-level driver to access low-level hardware 
end

%% Start experiment?
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
%Proxy: Window = 10, keyboardNumber = []

%% get experiment parameters
%I.e., get stimulus parameters, display parameters, response
%parameters,block structure, trigger parameters
param = sfa_expt4_getParams(window, freq, isBB, isMEG);
param.restarting = restarting;

Screen('FillRect',  window, param.BGcolor);
Screen('TextSize',  window, param.textSize);
Screen('TextColor', window, param.fontColor);
Screen('TextFont',  window, param.font);

DrawFormattedText(window,'getting ready...','center','center');
Screen('Flip',window);



%% create stimuli for the main experiment
if ~restarting %if not restarting, create new stimuli
    stim  = sfa_expt4_makeStim(param, pahandle);
    stim_ind_order = stim.ind_order; %Individual stim order
    save([pwd '\data\online_backup\stim_ind_order.mat'], 'stim_ind_order');

else
    % load previously created stim if we're restarting
    load([pwd '\data\online_backup\stim_ind_order.mat']);
    param.stim_ind_order = stim_ind_order;
    stim  = sfa_expt4_makeStim(param, pahandle);
    
end

%% run the experiment

if param.isMEG
    DrawFormattedText(window,'Please wait a moment while\n\nthe experimenters get everything set up.','center','center');
    Screen('Flip',window);
    
    % wait indefinitely, to allow experimenter to get MEG computer
    % ready for next recording
    bkey = recordValidKeys(GetSecs, Inf, param.keyboardNumber, {param.readyKey param.exitKey});
    if strcmp(bkey, param.exitKey)
        exitNow = 1;
        return
    end            

    sx = 'center';
    sy='center';
    blockText = 'Press any key to continue.'; 
    DrawFormattedText(window, blockText, sx, sy, param.fontColor);
    Screen('Flip',window);

    WaitSecs(1);
    KbWait;

else            
    DrawFormattedText(window,'Ready to start.\n\nPress any key to continue.','center','center');
    Screen('Flip',window);

    WaitSecs(.5);
    KbWait;

end

param.isPractice = isPractice;
if isPractice
    param.nTrialsPerBlock = nPracticeTrials;
    stim.nTrials = nPracticeTrials;
end
    
[data exitNow] = sfa_expt4_runBlock(window, param, stim);

Screen('CloseAll');
PsychPortAudio('Close', pahandle);

startup.endTime    = clock;
startup.endTimeStr = datestr(startup.endTime);
stim.soundwave  = [];
save(startup.dataFile, 'startup', 'param', 'stim', 'data');

