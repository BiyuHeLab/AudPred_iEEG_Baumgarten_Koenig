function [block] = NASTD_ECoG_Preproc_DefineBlocks...
    (data_ECoGraw, offset_inSamples)
%Aim: Determine how to define blocks

% load([preproc_dir sub '_triggerinfo.mat'])

blockCounter = 0;

trigger_info = data_ECoGraw.info_trigger;
fsample = data_ECoGraw.fsample;

if fsample == 512
    additionalsamples4trialend = 52; %52 samples with 512 Hz = 101.6ms
elseif fsample == 2048
    additionalsamples4trialend = 205; %52 samples with 2048 Hz = 100.1ms
end

%Block determination
for counter_blocktrigger = 1:(length(trigger_info.Exp.BlockStart_index)) 
    
    %Count blocks by means of block counter
       blockCounter = blockCounter+1;     
       blockbegin_index = trigger_info.Exp.BlockStart_index(counter_blocktrigger) ...
           - offset_inSamples;  
       %block begin is defined by blockstart trigger - X data points
       
       %Block end is defined by the first sample of the end trigger 
       %of the final trial of that block. However, the last sample would 
       %agree with the real tone sequence end. Since the whole trial-end 
       %trigger (1 pos deflection, 1 negativ deflection)is 100ms long, 
       %we now add these 100+X ms here to cover the whole tone sequence.       
       blockend_index = ...
           trigger_info.Exp.TrialEnd_index{counter_blocktrigger}...
           (end)+ additionalsamples4trialend + offset_inSamples;  
       %block end is defined by end of last trial + X data points/ 1sec 
       
       blockoffset_index = ...
           blockbegin_index - ...
           trigger_info.Exp.BlockStart_index(counter_blocktrigger); 
       %offset determines when relative to the start of the defined block
       %the 'block-start'trigger came (e.g., after 1 sec/ 512 samples). Sample
       %X+1 is then defined as t0.
       block(blockCounter,:) = ...
           [blockbegin_index blockend_index blockoffset_index blockCounter];
    
end

  block_avgdur_insamples = mean(block(:,2)-block(:,1));
  block_avgdur_insec = block_avgdur_insamples/data_ECoGraw.fsample;
  fprintf(['Found ', num2str(length(block)), ' of 12 blocks', ...
      '\nwith ' num2str(block_avgdur_insec/60) ' min average block duration.' ...
      '\n']); 
    
%Estimate of whole block Dur (stim + resp) = 152 s = 2.5 min
%All Tone seq duration in s =  5*(34*0.2) + 5*(34*0.4) = 102 s
%Response time (max 5s per trial) = 10*5 = 50 s    
    
    
       