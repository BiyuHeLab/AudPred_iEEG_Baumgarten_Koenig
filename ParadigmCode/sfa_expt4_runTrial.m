function [data exitNow] = sfa_expt4_runTrial(window, param, stim, trialNum, timingInd, data)


data.timing.trialStart(timingInd) = GetSecs;

exitNow = 0;

data.trialNum(trialNum) = trialNum;

sx = 'center';
sy = 'center';
vSpacing = 1.4;


%% show fixation dot

% give some buffer time between the offset of last trial's FB and the onset
% of this trial's fixation
WaitSecs(param.fixationBuffer_inSecs);

Screen('FillOval', window, param.fontColor, param.fixationRect);
data.timing.fixationOnset(timingInd) = Screen('Flip', window);

% % % image = Screen ('GetImage', window);
% % % file = 'sfa2_fixation.png';
% % % imwrite(image,file,'png');


%% buffer time so that auditory stim onset is 1 sec after FB offset of last trial

if trialNum == 1
    WaitSecs(param.trialStartBuffer_inSecs - param.fixationBuffer_inSecs);
else
    waitUntil(data.timing.FBoffset(timingInd-1), param.trialStartBuffer_inSecs);   
end


%% play the stimulus

t0 = GetSecs + param.paBuffer_inSecs; %determine timepoint when sequence is to be played

% trial_text = ['trial ' num2str(trialNum) ' of 360'];
% DrawFormattedText(window, trial_text, sx, sy, param.fontColor, [], [], [], vSpacing);
% Screen('Flip', window);

% load buffer
%create trigger, set it on top of stimulus sequence, then play both
%Trigger: long square wave (2*4410 data points = 0.2 sec because audio sampling freq = 44.100) at the beginning of sequence (i.e., tone 1), 
%short square wave (1*4410 data points = 0.1 sec) at end of sequence (i.e., final tone)
w  = stim.soundwave{trialNum};
%w2 = [w; w];
trigger = zeros(size(w));
trigger(1:2205) = 1;% add square wave beginning of trial
trigger(2206:4410) = -1;% add square wave beginning of trial
trigger(4411:6615) = 1;% add square wave beginning of trial
trigger(6616:8820) = -1;% add square wave beginning of trial

trigger(end-4410:end-2205) = 1;% add square wave end of trial
trigger(end-2205:end) = -1;% add square wave end of trial

w2 = [trigger; w];
PsychPortAudio('FillBuffer', stim.pahandle, w2);

% initialize playback
data.timing.stimOnset(timingInd)     = PsychPortAudio('Start', stim.pahandle, [], t0);
data.timing.stimOnset_req(timingInd) = t0;

% send the trigger
% % % if param.isMEG
% % %     trigger = param.trig_stim( stim.betaID(trialNum) );
% % %     data.triggers.stim_trigger(trialNum) = trigger;
% % %     data.timing.stim_trigger(trialNum)   = send_trigger(trigger, param);
% % % end

% check pahandle status until it's done
WaitSecs(param.paBuffer_inSecs + param.seriesDur_inSecs(stim.toneDurID(trialNum)) - 2);

s = PsychPortAudio('GetStatus', stim.pahandle);
while s.Active == 1
    s = PsychPortAudio('GetStatus', stim.pahandle);
end

% % remove fixation
% Screen('Flip', window);

% estimate duration of series playback with two methods
data.timing.dur_loop(timingInd) = GetSecs - t0;%current time - time play sequence 
data.timing.dur_est(timingInd)  = s.EstimatedStopTime - t0;
data.pahandle_status(timingInd) = s;


%% get behavioral responses

%%%%%%%%%%% start probability rating %%%%%%%%%%%

% prompt for final tone probability rating
% % % r1_text = 'Final tone likelihood?\n\n1 = very unlikely,  2 = somewhat unlikely,  3 = somewhat likely,  4 = very likely\n(use LEFT HAND)';

if param.isBB
    r1_text = 'Final tone likelihood?\n\n1 = very unlikely   5 = very likely\nfinger mapping: pinky for 1 up through thumb for 5';
else
    r1_text = 'Final tone likelihood?\n\n1 = very unlikely   5 = very likely\nuse keys 1-5';
end

DrawFormattedText(window, r1_text, sx, sy, param.fontColor, [], [], [], vSpacing);

WaitSecs(.4);
data.timing.r1_prompt(timingInd) = Screen('Flip',window);

% % % image = Screen ('GetImage', window);
% % % file = 'sfa2_q1.png';
% % % imwrite(image,file,'png');
 
% collect tone probability rating
[r1_key data.r1_RT(trialNum) data.timing.r1_time(timingInd)] = ... 
    recordValidKeys(data.timing.r1_prompt(timingInd), param.respDur_inSecs, param.keyboardNumber, param.r1_validKeys);

Screen('Flip',window);

% assess response
id = 6;
rk = param.r1_validKeys;
switch r1_key
    case rk{1},         data.resp_prob(trialNum) =  1;  id = 1;
    case rk{2},         data.resp_prob(trialNum) =  2;  id = 2;
    case rk{3},         data.resp_prob(trialNum) =  3;  id = 3;
    case rk{4},         data.resp_prob(trialNum) =  4;  id = 4;
    case rk{5},         data.resp_prob(trialNum) =  5;  id = 5;
    case 'noanswer',    data.resp_prob(trialNum) = -1;
    case 'invalid',     data.resp_prob(trialNum) = -2;
    case 'cell',        data.resp_prob(trialNum) = -3;
    case param.exitKey, data.resp_prob(trialNum) = -4; exitNow = 1;
    otherwise,          data.resp_prob(trialNum) = -5;
end

% send trigger
if param.isMEG
    trigger = param.trig_r1(id);
    data.timing.r1_trigger(timingInd)   = send_trigger(trigger, param);
    data.triggers.r1_trigger(trialNum) = trigger;
end

% give 'too slow' warning if necessary
if strcmp(r1_key,'noanswer')
    data.timing.r1_tooSlow(timingInd) = 1;
    DrawFormattedText(window, 'Too slow!', sx, sy, param.fontColor)
    Screen('Flip', window);
    WaitSecs(2);
end

Screen('Flip',window);

%Give response feedback
% sort out response 
if data.resp_prob(trialNum) == 1
   feedback_response = 'very unlikely';
elseif data.resp_prob(trialNum) == 2
   feedback_response = 'relatively unlikely';
elseif data.resp_prob(trialNum) == 3
   feedback_response = 'unsure';
elseif data.resp_prob(trialNum) == 4
   feedback_response = 'relatively likely';
elseif data.resp_prob(trialNum) == 5
   feedback_response = 'very likely';
end

switch data.resp_prob(trialNum)
    case 1,  FB_text = ['Response = 1\n\n' ...
             'Likelihood rating = ' feedback_response '.'] ;
    case 2,  FB_text = ['Response = 2\n\n' ...
             'Likelihood rating = ' feedback_response '.'] ;
    case 3,  FB_text = ['Response = 3\n\n' ...
             'Likelihood rating = ' feedback_response '.'] ;
    case 4,  FB_text = ['Response = 4\n\n' ...
             'Likelihood rating = ' feedback_response '.'] ;
    case 5,  FB_text = ['Response = 5\n\n' ...
             'Likelihood rating = ' feedback_response '.'] ;

    case -1, FB_text = ['No response recorded.'];
    case -2, FB_text = ['Invalid response.'];
    case -3, FB_text = ['Invalid response.'];
    case -4, FB_text = ['Exit key.'];
    case -5, FB_text = ['Invalid response.'];
        
end
DrawFormattedText(window, FB_text, sx, sy, param.fontColor, [], [], [], vSpacing);
Screen('Flip', window);
WaitSecs(2);

%%%%%%%%%%% start useq estimation %%%%%%%%%%%
% 
% % prompt for useq estimation
% if param.isBB
%     r2_text = '\n\nMatch the sequence to its picture.\n\nring finger --> select A \nmiddle finger --> select B \nindex finger --> select C';
% else
%     r2_text = '\n\nMatch the sequence to its picture.\n\n[A] --> select 8 key          [B] --> 9 key          [C] -->0 key';
% end
% 
% DrawFormattedText(window, r2_text, sx, 5, param.fontColor, [], [], [], vSpacing);
% 
% %     % Option1 (original version): draw 2 random images, and the 1 correct one
% %     uSeq_corrID = stim.uSeqID(trialNum);
% % 
% %     uSeq_all = [1:15];
% %     uSeq_all = uSeq_all(~(uSeq_corrID == uSeq_all)); %Take out correct
% %     uSeq_all = Shuffle(uSeq_all); %Shuffle
% % 
%     % Option2 (easy version): draw 2 random images with other trend/p*34 (low, medium, high), and the correct one
%     uSeq_corrID = stim.uSeqID(trialNum);
%     uSeq_predID = stim.predID(trialNum); %Determine trend/p*34
% 
%     diff_predID = stim.predID ~= uSeq_predID;
%     uSeq_all = unique(stim.uSeqID(diff_predID));
% 
%     uSeq_all = Shuffle(uSeq_all); %Shuffle
% 
% uSeq_decoy1 = uSeq_all(1);
% uSeq_decoy2 = uSeq_all(2);
% 
% uSeq_disp3 = [uSeq_corrID uSeq_decoy1 uSeq_decoy2];
% uSeq_disp3 = Shuffle(uSeq_disp3);
% 
% 
% data.uSeq.disp3(trialNum,:) = uSeq_disp3;
% data.uSeq.corrID(trialNum) = uSeq_corrID;
% data.uSeq.corrABC(trialNum) = find(uSeq_disp3 == uSeq_corrID);
% 
% %TJB: add mac/gogo paths
% % addpath('C:\Users\localbiyu\Desktop\Brian\sfa\sfa_expt4\useq_pics')
% if ismac
% addpath('/Users/ecog/Documents/MATLAB/sfa_expt4/useq_pics') %Path Laptop Adeen
% else
% addpath('\\gogo.sb.nyumc.org\data\gogodisk4\thomas\AuditoryPrediction\iEEG\Experiment_Versions\sfa_expt4_Adeen_Current\useq_pics'); %Path desktop/gogo
% end
% 
% imaA = imread(['uSeq_' num2str(uSeq_disp3(1))], 'png');
% imaB = imread(['uSeq_' num2str(uSeq_disp3(2))], 'png');
% imaC = imread(['uSeq_' num2str(uSeq_disp3(3))], 'png');
% 
% %%% scale down images
% scalefactor = 0.3;
% imaA = imresize(imaA, scalefactor);
% imaB = imresize(imaB, scalefactor);
% imaC = imresize(imaC, scalefactor);
% 
% imgarray = [imaA zeros(size(imaA,1), 150, 3) imaB zeros(size(imaA,1), 150, 3) imaC]; 
% Screen('PutImage', window, imgarray)
% 
% WaitSecs(.4);
% data.timing.r2_prompt(timingInd) = Screen('Flip',window);
% 
% % % % image = Screen ('GetImage', window);
% % % % file = 'sfa2_q2.png';
% % % % imwrite(image,file,'png');
% 
% % % collect useq estimation
% % [r2_key data.r2_RT(trialNum) data.timing.r2_time(timingInd)] = ... 
% %     recordValidKeys(data.timing.r2_prompt(timingInd), param.respDur_inSecs, param.keyboardNumber, param.r2_validKeys); %original command
% 
% [r2_key data.r2_RT(trialNum) data.timing.r2_time(timingInd)] = ... 
%     recordValidKeys(data.timing.r2_prompt(timingInd), param.respDur_inSecs*2, param.keyboardNumber, param.r2_validKeys); %double resp time for task 2
% 
% 
% Screen('Flip',window);
% 
% % assess response
% id = 6;
% rk = param.r2_validKeys;
% switch r2_key
% %     case rk{1},         data.uSeq.respABC(trialNum) =  0;     id = 1;   uSeq_resp = 0;
%     case rk{1},         data.uSeq.respABC(trialNum) =  1;    id = 1;   data.uSeq.respID(trialNum) = uSeq_disp3(1);
%     case rk{2},         data.uSeq.respABC(trialNum) =  2;    id = 2;   data.uSeq.respID(trialNum) = uSeq_disp3(2);
%     case rk{3},         data.uSeq.respABC(trialNum) =  3;    id = 3;   data.uSeq.respID(trialNum) = uSeq_disp3(3);
% %     case rk{5},         data.uSeq.respABC(trialNum) =  2;     id = 5;   uSeq_resp = 4;
%     case 'noanswer',    data.uSeq.respABC(trialNum) = -1;
%     case 'invalid',     data.uSeq.respABC(trialNum) = -2;
%     case 'cell',        data.uSeq.respABC(trialNum) = -3;
%     case param.exitKey, data.uSeq.respABC(trialNum) = -4; exitNow = 1;
%     otherwise,          data.uSeq.respABC(trialNum) = -5;
% end
% 
% % send trigger
% if param.isMEG
%     trigger = param.trig_r2(id);
%     data.timing.r2_trigger(timingInd)   = send_trigger(trigger, param);
%     data.triggers.r2_trigger(trialNum) = trigger;
% end
% 
% % give 'too slow' warning if necessary
% if strcmp(r2_key,'noanswer')
%     data.timing.r2_tooSlow(timingInd) = 1;    
%     DrawFormattedText(window, 'Too slow!', sx, sy, param.fontColor)
%     Screen('Flip', window);
%     WaitSecs(2);
% end
% 
% Screen('Flip',window);
% 
% 
% %% sort out response 
% 
% % sort out correctness
% if data.uSeq.respABC(trialNum) == data.uSeq.corrABC(trialNum)
%     data.correct_uSeq(trialNum) = 1;
% 
% % N/A if response < 1
% elseif data.uSeq.respABC(trialNum) < 0
%     data.correct_uSeq(trialNum) = -1;
% 
% % otherwise, just incorrect
% else
%     data.correct_uSeq(trialNum) = 0;
% end
% 
% if exitNow, return; end
% 
% 
% %% show feedback
% 
% ABC = {'A' 'B' 'C'};
% corrABC = ABC{data.uSeq.corrABC(trialNum)};
% if data.correct_uSeq(trialNum) == 1 || data.correct_uSeq(trialNum) == 0
%     respABC = ABC{data.uSeq.respABC(trialNum)};
% end
% 
% switch data.correct_uSeq(trialNum)
%     case 1,  FB_text = 'Correct!';
%     case 0,  FB_text = ['Incorrect! You answered ' respABC '. The correct answer was ' corrABC '.'];
%     case -1, FB_text = ['No response recorded. The correct answer was ' corrABC];
% end
% DrawFormattedText(window, FB_text, sx, sy, param.fontColor, [], [], [], vSpacing);
% Screen('Flip', window);

data.timing.FBoffset(timingInd) = Screen('Flip', window);

data.timing.trialEnd(timingInd) = GetSecs;
data.timing.trialDur(timingInd) = data.timing.trialEnd(timingInd) - data.timing.trialStart(timingInd);

end