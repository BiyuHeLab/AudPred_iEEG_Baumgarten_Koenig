function DataClean_AllTrials = NASTD_ECoG_Preproc_SelCleanTrials...
    (sub,DataClean_AllTrials)
 
%Select clean trials from preprocessed data

%load in file with individual preproc infos
subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings; 

Num_AllTrials = length(DataClean_AllTrials.trial);

if ~isempty(subs_PreProcSettings.(sub).rejectedTrials)
    
    disp([num2str(length(subs_PreProcSettings.(sub).rejectedTrials)) ' rejected trials found - selecting clean trials only'])   
    Index_CleanTrials = setdiff(1:Num_AllTrials, subs_PreProcSettings.(sub).rejectedTrials);
  
%Select clean/valid trials in ECoG data
    dat_fields  = fieldnames(DataClean_AllTrials);
    for i_field = 1:length(dat_fields)
        if eval(['length(DataClean_AllTrials.' dat_fields{i_field} ') == ' num2str(Num_AllTrials)]) %if missed trials are encoded in current field
%             disp(['Adjusting subfield: ' dat_fields{i_field}])
            if size(eval(['DataClean_AllTrials.' dat_fields{i_field}]),1) == 1
            eval(['DataClean_AllTrials.' dat_fields{i_field} ' = DataClean_AllTrials.' dat_fields{i_field} '(Index_CleanTrials);']);
            else 
            eval(['DataClean_AllTrials.' dat_fields{i_field} ' = DataClean_AllTrials.' dat_fields{i_field} '(Index_CleanTrials,:);']);
            end
        end

    end
%Select clean/valid trials in behav data    
    dat_fields  = fieldnames(DataClean_AllTrials.behav);
    for i_field = 1:length(dat_fields)
        if eval(['length(DataClean_AllTrials.behav.' dat_fields{i_field} ') == ' num2str(Num_AllTrials)]) %if missed trials are encoded in current field
%             disp(['Adjusting subfield: ' dat_fields{i_field}])
            eval(['DataClean_AllTrials.behav.' dat_fields{i_field} ' = DataClean_AllTrials.behav.' dat_fields{i_field} '(Index_CleanTrials);']);
        end
    end
%Select clean/valid trials in stim data    
    dat_fields  = fieldnames(DataClean_AllTrials.behav.stim);
    for i_field = 1:length(dat_fields)
        if eval(['length(DataClean_AllTrials.behav.stim.' dat_fields{i_field} ') == ' num2str(Num_AllTrials)]) %if missed trials are encoded in current field
%             disp(['Adjusting subfield: ' dat_fields{i_field}])
            eval(['DataClean_AllTrials.behav.stim.' dat_fields{i_field} ' = DataClean_AllTrials.behav.stim.' dat_fields{i_field} '(Index_CleanTrials);']);
        end
    end
        
%Select clean/valid trials in behav.timing data    
    dat_fields  = fieldnames(DataClean_AllTrials.behav.timing);
    for i_field = 1:length(dat_fields)
        if eval(['length(DataClean_AllTrials.behav.timing.' dat_fields{i_field} ') == ' num2str(Num_AllTrials)]) %if missed trials are encoded in current field
%             disp(['Adjusting subfield: ' dat_fields{i_field}])
            eval(['DataClean_AllTrials.behav.timing.' dat_fields{i_field} ' = DataClean_AllTrials.behav.timing.' dat_fields{i_field} '(Index_CleanTrials);']);
        end
    end
    
end

%Final check if ECoG data length, stimulus info, behav data length agree for each trial
Check_SeqLengthperTrial = NaN(length(DataClean_AllTrials.trial),3);
for i_trial = 1:length(DataClean_AllTrials.trial)
    
    if length(DataClean_AllTrials.trial{i_trial}) < 5000 %short seqs ~ 4500 samples, long seqs ~ 8000 samples
        Check_SeqLengthperTrial(i_trial,1) = 1;
    elseif length(DataClean_AllTrials.trial{i_trial}) > 5000
        Check_SeqLengthperTrial(i_trial,1) = 2;
    end
    
    if DataClean_AllTrials.behav.stim.toneDur(i_trial) == 0.2
        Check_SeqLengthperTrial(i_trial,2) = 1;
    elseif DataClean_AllTrials.behav.stim.toneDur(i_trial) == 0.4
        Check_SeqLengthperTrial(i_trial,2) = 2;
    end
    
%     if DataClean_AllTrials.behav.timing.trialDur(i_trial) < 15 %short seqs ~ 11-14 sec, long seqs ~ 17-20 sec
%         Check_SeqLengthperTrial(i_trial,3) = 1;
%     elseif DataClean_AllTrials.behav.timing.trialDur(i_trial) > mean(DataClean_AllTrials.behav.timing.trialDur)
%         Check_SeqLengthperTrial(i_trial,3) = 2;
%     end
    
    if isequal(Check_SeqLengthperTrial(i_trial,1),Check_SeqLengthperTrial(i_trial,2)) ~= 1        
        disp(['Behav-ECoG data mismatch in Trial: ' num2str(i_trial) ' !'])
%         pause
    end

end

clear Check_SeqLengthperTrial