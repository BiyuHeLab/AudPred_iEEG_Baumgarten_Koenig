function input_data = ...
    NASTD_ECoG_FiltNaNinterp_RemoveNaNtrials...
    (sub, TD_label, ...
    input_data)

%Aim: Reject trials where some electrodes have only NaN-entries (in
%particular NY798)

index_excludedtrials = [];

for i_trial = 1:length(input_data.trial)
    for i_elec = 1:size(input_data.trial{i_trial},1)
        ind_noNaNsamples{i_trial}{i_elec} = ...
            find(~isnan(input_data.trial{i_trial}(i_elec,:)));
        
        if isempty(ind_noNaNsamples{i_trial}{i_elec})
            index_excludedtrials = [index_excludedtrials, i_trial];
        end
    end
end
index_excludedtrials                            = unique(index_excludedtrials);
index_remainingtrials                           = 1:length(input_data.trial);
index_remainingtrials(index_excludedtrials)     = [];
filt_remainingtrial                             = true(1,length(input_data.trial));
filt_remainingtrial(index_excludedtrials)       = false;
Num_AllTrials                                   = length(input_data.trial);

if ~isempty(index_excludedtrials)
    %MEG data
    cfg = [];
    cfg.trials = index_remainingtrials;
    input_data = ft_selectdata(cfg,input_data);
    
    %Behav data
    dat_fields  = fieldnames(input_data.behav);
    for i_field = 1:length(dat_fields)
        if eval(['length(input_data.behav.' dat_fields{i_field} ...
                ') == ' num2str(Num_AllTrials)]) %if missed trials are encoded in current field
            eval(['input_data.behav.' dat_fields{i_field} ...
                ' = input_data.behav.' dat_fields{i_field} '(filt_remainingtrial);']);
        end
    end
    
    %Behav.stim data
    dat_fields  = fieldnames(input_data.behav.stim);
    for i_field = 1:length(dat_fields)
        if eval(['length(input_data.behav.stim.' dat_fields{i_field} ...
                ') == ' num2str(Num_AllTrials)]) %if missed trials are encoded in current field
            eval(['input_data.behav.stim.' dat_fields{i_field} ...
                ' = input_data.behav.stim.' dat_fields{i_field} '(filt_remainingtrial);']);
        end
    end
    
    %Behav.timing data
    dat_fields  = fieldnames(input_data.behav.timing);
    for i_field = 1:length(dat_fields)
        if eval(['length(input_data.behav.timing.' dat_fields{i_field} ...
                ') == ' num2str(Num_AllTrials)]) %if missed trials are encoded in current field
            eval(['input_data.behav.timing.' dat_fields{i_field} ...
                ' = input_data.behav.timing.' dat_fields{i_field} '(filt_remainingtrial);']);
        end
    end
    input_data.info.FiltInterp.rejectedtrialind ...
        = index_excludedtrials;
    disp([num2str(length(index_excludedtrials)) ...
        ' trials removed due to all-NaN entries for '...
        sub '; TD:' TD_label])
else
    input_data.info.FiltInterp.rejectedtrialind ...
        = [];
    disp(['No trials removed due to all-NaN entries for '...
        sub '; TD:' TD_label])
end