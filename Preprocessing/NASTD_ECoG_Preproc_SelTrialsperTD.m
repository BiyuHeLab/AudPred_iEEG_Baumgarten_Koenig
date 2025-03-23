function DataClean_CleanTrialsElecs = NASTD_ECoG_Preproc_SelTrialsperTD...
    (ToneDur_text, DataClean_CleanTrialsElecs)

%Select only trials where tone sequence was presented with specific tone
%duration (TD)

%1. Read out TD-specific trials and adjust data set
% toneDur = {'0.2' '0.4', 'all'};
if ~strcmp(ToneDur_text,'all')
    TrialFilt_ToneDur = ...
        DataClean_CleanTrialsElecs.behav.stim.toneDur == ...
        str2double(ToneDur_text); %filter to select trials for certain tone dur
    
    %ECoG data
    Num_AllTrials = length(DataClean_CleanTrialsElecs.trial);
    
    dat_fields  = fieldnames(DataClean_CleanTrialsElecs);
    for i_field = 1:length(dat_fields)
        if eval(['length(DataClean_CleanTrialsElecs.' ...
                dat_fields{i_field} ') == ' num2str(Num_AllTrials)]) 
            %if missed trials are encoded in current field
            if size(eval(['DataClean_CleanTrialsElecs.' ...
                    dat_fields{i_field}]),1) == 1
                eval(['DataClean_CleanTrialsElecs.' ...
                    dat_fields{i_field} ' = DataClean_CleanTrialsElecs.' ...
                    dat_fields{i_field} '(TrialFilt_ToneDur);']);
            else
                eval(['DataClean_CleanTrialsElecs.' ...
                    dat_fields{i_field} ' = DataClean_CleanTrialsElecs.' ...
                    dat_fields{i_field} '(TrialFilt_ToneDur,:);']);
            end
        end
    end
    
    %Behav data
    dat_fields  = fieldnames(DataClean_CleanTrialsElecs.behav);
    for i_field = 1:length(dat_fields)
        if eval(['length(DataClean_CleanTrialsElecs.behav.' ...
                dat_fields{i_field} ') == ' num2str(Num_AllTrials)]) 
            %if missed trials are encoded in current field
            eval(['DataClean_CleanTrialsElecs.behav.' ...
                dat_fields{i_field} ' = DataClean_CleanTrialsElecs.behav.' ...
                dat_fields{i_field} '(TrialFilt_ToneDur);']);
        end
    end
    
    %Behav.stim data
    dat_fields  = fieldnames(DataClean_CleanTrialsElecs.behav.stim);
    for i_field = 1:length(dat_fields)
        if eval(['length(DataClean_CleanTrialsElecs.behav.stim.' ...
                dat_fields{i_field} ') == ' num2str(Num_AllTrials)]) 
            %if missed trials are encoded in current field
            eval(['DataClean_CleanTrialsElecs.behav.stim.' ...
                dat_fields{i_field} ' = DataClean_CleanTrialsElecs.behav.stim.' ...
                dat_fields{i_field} '(TrialFilt_ToneDur);']);
        end
    end
    
    %Behav.timing data
    dat_fields  = fieldnames(DataClean_CleanTrialsElecs.behav.timing);
    for i_field = 1:length(dat_fields)
        if eval(['length(DataClean_CleanTrialsElecs.behav.timing.' ...
                dat_fields{i_field} ') == ' num2str(Num_AllTrials)]) 
            %if missed trials are encoded in current field
            eval(['DataClean_CleanTrialsElecs.behav.timing.' ...
                dat_fields{i_field} ' = DataClean_CleanTrialsElecs.behav.timing.' ...
                dat_fields{i_field} '(TrialFilt_ToneDur);']);
        end
    end
    
end

end
