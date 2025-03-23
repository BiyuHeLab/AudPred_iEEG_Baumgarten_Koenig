function preprocData_perTD = NASTD_ECoG_SelTrialsperPredp34(label_predP34, preprocData_perTD)

% label_predP34 = [-1 0 1];% low, medium, high;

TrialFilt_Predp34 = preprocData_perTD.behav.stim.predID == label_predP34; %filter to select trials for certain tone dur

%ECoG data
Num_AllTrials = length(preprocData_perTD.trial);

dat_fields  = fieldnames(preprocData_perTD);
for i_field = 1:length(dat_fields)
    if eval(['length(preprocData_perTD.' dat_fields{i_field} ') == ' num2str(Num_AllTrials)]) %if missed trials are encoded in current field
        %             disp(['Adjusting subfield: ' dat_fields{i_field}])
        if size(eval(['preprocData_perTD.' dat_fields{i_field}]),1) == 1
            eval(['preprocData_perTD.' dat_fields{i_field} ' = preprocData_perTD.' dat_fields{i_field} '(TrialFilt_Predp34);']);
        else
            eval(['preprocData_perTD.' dat_fields{i_field} ' = preprocData_perTD.' dat_fields{i_field} '(TrialFilt_Predp34,:);']);
        end
    end
end

%Behav data
dat_fields  = fieldnames(preprocData_perTD.behav);
for i_field = 1:length(dat_fields)
    if eval(['length(preprocData_perTD.behav.' dat_fields{i_field} ') == ' num2str(Num_AllTrials)]) %if missed trials are encoded in current field
        %             disp(['Adjusting subfield: ' dat_fields{i_field}])
        eval(['preprocData_perTD.behav.' dat_fields{i_field} ' = preprocData_perTD.behav.' dat_fields{i_field} '(TrialFilt_Predp34);']);
    end
end

%Behav.stim data
dat_fields  = fieldnames(preprocData_perTD.behav.stim);
for i_field = 1:length(dat_fields)
    if eval(['length(preprocData_perTD.behav.stim.' dat_fields{i_field} ') == ' num2str(Num_AllTrials)]) %if missed trials are encoded in current field
        %             disp(['Adjusting subfield: ' dat_fields{i_field}])
        eval(['preprocData_perTD.behav.stim.' dat_fields{i_field} ' = preprocData_perTD.behav.stim.' dat_fields{i_field} '(TrialFilt_Predp34);']);
    end
end

%Behav.timing data
dat_fields  = fieldnames(preprocData_perTD.behav.timing);
for i_field = 1:length(dat_fields)
    if eval(['length(preprocData_perTD.behav.timing.' dat_fields{i_field} ') == ' num2str(Num_AllTrials)]) %if missed trials are encoded in current field
        %             disp(['Adjusting subfield: ' dat_fields{i_field}])
        eval(['preprocData_perTD.behav.timing.' dat_fields{i_field} ' = preprocData_perTD.behav.timing.' dat_fields{i_field} '(TrialFilt_Predp34);']);
    end
end

end
