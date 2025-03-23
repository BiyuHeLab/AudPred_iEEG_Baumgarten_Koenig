function [trial] = NASTD_ECoG_Preproc_DefineTrials(data_ECoGfiltref_blocks,offset_inSamples)
% load([preproc_dir sub '_triggerinfo.mat'])

BlockCounter = 0;
TrialCounter = 0;

trigger_info = data_ECoGfiltref_blocks.cfg.info_trigger;
fsample = data_ECoGfiltref_blocks.fsample;
if fsample == 512
    additionalsamples4trialend = 52; %52 samples with 512 Hz = 101.6ms
elseif fsample == 2048
    additionalsamples4trialend = 205; %52 samples with 2048 Hz = 100.1ms
end

%Block determination
for i_blocktrigger = 1:length(trigger_info.Exp.BlockStart_index)
    BlockCounter = BlockCounter+1;  
    for i_trialtrigger = 1:length(trigger_info.Exp.TrialStart_index{i_blocktrigger})   
    TrialCounter = TrialCounter+1;
       trialbegin_index = ...
           trigger_info.Exp.TrialStart_index{i_blocktrigger}(i_trialtrigger) ...
           - offset_inSamples;  %trial begin is defined by blockstart trigger - offset 
       
       %!! Note: During trigger read-out, we defined the first sample in
       %the trial-end trigger deflection as trial end signal. However, the
       %last sample would agree with the real tone sequence end. Since the
       %whole trial-end trigger (1 pos deflection, 1 negativ deflection)is
       %100ms long, we now add these 100+X ms here to cover the whole tone
       %sequence.
       trialend_index = ...
           trigger_info.Exp.TrialEnd_index{i_blocktrigger}(i_trialtrigger) ...
           + additionalsamples4trialend + offset_inSamples;  
       %trial end is defined by beginning of the trial end trigger + additionalsamples4trialend data
       %points + offset
       trialoffset_index =  trialbegin_index - trigger_info.Exp.TrialStart_index{i_blocktrigger}(i_trialtrigger); 
       %offset determines when relative to the start of the defined trial
       %the 'trial-start'trigger came
       trial(TrialCounter,:) =  [trialbegin_index trialend_index trialoffset_index TrialCounter BlockCounter];
    end
end
    
  trial_avgdur_insec = mean((trial(:,2)-trial(:,1))./fsample);
  fprintf(['Found ', num2str(length(trial)), ' of 120 trials', ...
      '\nwith ' num2str(round(trial_avgdur_insec,1)) ...
      ' sec average trial duration (theory ~ 12 sec).' '\n']); 
    
  %Theoretical trial length:
  %short TD*number tones + max response time = 0.2*34 + 2 = 8.8 sec
  %long TD*number tones + max response time = 0.4*34 +2 = 15.6 sec 
  %Average across TD = 12.2s
  
end        