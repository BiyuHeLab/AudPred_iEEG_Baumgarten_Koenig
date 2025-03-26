%% Aim:
%Save tone sequences used in NASTD_ECoG study as 1) Matlab file and 2) .wav
%file and provide a plot that shows each sequence

%% 0) Specify paths and analysis options
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
NASTD_ECoG_setVars
paths_NASTD_ECoG = NASTD_ECoG_paths;

% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));

path_outputdata = '/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/NASTD_ECoG_Matlab/CreateWav/Sequences/';
if (~exist(path_outputdata, 'dir')); mkdir(path_outputdata); end
path_fig = '/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/NASTD_ECoG_Matlab/CreateWav/Figs/';
if (~exist(path_fig, 'dir')); mkdir(path_fig); end


%Determine subjects
sub_list = vars.sub_list;

%Load in file with individual preproc infos
subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;

ToneDur_text = {'0.2' '0.4'};


for i_sub = 2%:length(subs) %Loop subjects
    for i_tonedur = 1:length(ToneDur_text) %tone duration condition
        
        %Load in behavioral and tone sequence data
        sub = sub_list{i_sub};
        NASTD_ECoG_subjectinfo %load subject info file (var: si)

        loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];    
        load(loadfile_ECoGpreprocdata);
        
        
        %1.3) Filter trials with selected ToneDur
        %tone length ID: 0.15 = 1, 0.3 = 2, 0.6 = 3
        numoftrials = num2str(length(DataClean_AllTrials.trial));
        
        f = DataClean_AllTrials.behav.stim.toneDur == str2double(ToneDur_text{i_tonedur}); %filter to select trials for certain tone dur
        
        dat_fields  = fieldnames(DataClean_AllTrials);
        for i = 1:length(dat_fields)
            if eval(['length(DataClean_AllTrials.' dat_fields{i} ') ==' numoftrials]) %360
                eval(['DataClean_AllTrials.' dat_fields{i} ' = DataClean_AllTrials.' dat_fields{i} '(f);']);
            end
        end
        
        stim_fields = fieldnames(DataClean_AllTrials.behav.stim);
        for i = 1:length(stim_fields)
            if eval(['length(DataClean_AllTrials.behav.stim.' stim_fields{i} ') ==' numoftrials]) %360
                eval(['DataClean_AllTrials.behav.stim.' stim_fields{i} ' = DataClean_AllTrials.behav.stim.' stim_fields{i} '(f);']);
            end
        end
        
        %2. Place tone sequence information in summary file
        
        %2.1 Determine unique tone sequences based on tone 1-33
        %2.1.1 Create seqID proxy field
        DataClean_AllTrials.behav.stim.seqID = (1:length(DataClean_AllTrials.behav.stim.series_f));
        DataClean_AllTrials.behav.stim.seqID(1:end) = NaN;
        seqID_perUniqueSum = [];
        
        %2.1.2 Read out all trials and determine individual sequence ID
        for i_trial = 1:length(DataClean_AllTrials.behav.stim.series_f)
            All_trials_block1(i_trial,:) = DataClean_AllTrials.behav.stim.series_f{i_trial}(1:33); %only first 33 tone, since 34th varies
        end
        seqID_perUniqueSum = unique(sum(All_trials_block1')); %use sum across tone frequency values as metric determinig SeqID
        %should result in 15 unique sequences (independent of final tone pitch/p34)
        %Identical across tone duration conditions
        
        %2.1.3 Individually label all trials based on their across-tone-freq-sum
        for i_trial = 1:length(DataClean_AllTrials.behav.stim.series_f) %loop across all trials
            for i_seqID = 1:length(seqID_perUniqueSum) %loop across uSeqs
                if sum(DataClean_AllTrials.behav.stim.series_f{i_trial}(1:33)) == seqID_perUniqueSum(i_seqID)
                    DataClean_AllTrials.behav.stim.seqID(i_trial) = i_seqID;
                end
            end
        end
        
        %2.2 place all tone sequences in summary file
        for i_seq = 1:length(DataClean_AllTrials.behav.stim.series_f)
            ToneSeq.SequenceID(i_seq,1) = DataClean_AllTrials.behav.stim.seqID(i_seq);
            ToneSeq.Beta(i_seq,1) = DataClean_AllTrials.behav.stim.beta(i_seq);
            ToneSeq.ToneDur_Sec(i_seq,1) = DataClean_AllTrials.behav.stim.toneDur(i_seq);
            ToneSeq.PresentedFinalTonePitch_logF(i_seq,1) = DataClean_AllTrials.behav.stim.logf_final(i_seq);
            ToneSeq.PredictedFinalTonePitch_logF(i_seq,1) = DataClean_AllTrials.behav.stim.logf_pred(i_seq);
            ToneSeq.TonePitch_Hz(i_seq,:) = DataClean_AllTrials.behav.stim.series_f{i_seq};
        end
        
        %2.3 place unique tone sequences in summary file        
        for i_uSeq = unique(ToneSeq.SequenceID)'
            selected_uSeq = find(ToneSeq.SequenceID == i_uSeq);
            
            uniqueToneSeq{i_tonedur}.Beta(i_uSeq,1) = ToneSeq.Beta(selected_uSeq(1),1);
            uniqueToneSeq{i_tonedur}.ToneDur_Sec(i_uSeq,1) = ToneSeq.ToneDur_Sec(selected_uSeq(1),1);
            uniqueToneSeq{i_tonedur}.PresentedFinalTonePitch_logF(i_uSeq,:) = unique(ToneSeq.PresentedFinalTonePitch_logF(selected_uSeq))';
            uniqueToneSeq{i_tonedur}.PredictedFinalTonePitch_logF(i_uSeq,1) = ToneSeq.PredictedFinalTonePitch_logF(selected_uSeq(1),1);
            
            if uniqueToneSeq{i_tonedur}.PredictedFinalTonePitch_logF(i_uSeq,1) < 6
                uniqueToneSeq{i_tonedur}.PredictedFinalTonePitch{i_uSeq,1} = 'low';
            elseif uniqueToneSeq{i_tonedur}.PredictedFinalTonePitch_logF(i_uSeq,1) > 6.0868
                uniqueToneSeq{i_tonedur}.PredictedFinalTonePitch{i_uSeq,1} = 'high';
            else uniqueToneSeq{i_tonedur}.PredictedFinalTonePitch_logF(i_uSeq,1)
                uniqueToneSeq{i_tonedur}.PredictedFinalTonePitch{i_uSeq,1} = 'med';
            end
            
            uniqueToneSeq{i_tonedur}.TonePitch_Hz(i_uSeq,:) = ToneSeq.TonePitch_Hz(selected_uSeq(1),1:33);
            uniqueToneSeq{i_tonedur}.ToneSeries_Soundwave(i_uSeq,:) = ...
                series2soundwave(ToneSeq.TonePitch_Hz(selected_uSeq(1),1:33), str2double(ToneDur_text{i_tonedur}), 44100);
        end
        
        %Sort struct based on beta level
        T = struct2table(uniqueToneSeq{i_tonedur});
        sortedT = sortrows(T,'Beta');
        uniqueToneSeq{i_tonedur} = table2struct(sortedT,'ToScalar',true);
        
        %Plot uSeq and convert to .wav file'
        for i_uSeq = unique(ToneSeq.SequenceID)'
            
            if i_tonedur == 1 %do this only once
                figure;
                scatter(1:length(uniqueToneSeq{i_tonedur}.TonePitch_Hz(i_uSeq,:)),log(uniqueToneSeq{i_tonedur}.TonePitch_Hz(i_uSeq,:)),'o', 'filled', 'k')  %plot all tones including p34
                hold on;
                plot(1:length(uniqueToneSeq{i_tonedur}.TonePitch_Hz(i_uSeq,:)),log(uniqueToneSeq{i_tonedur}.TonePitch_Hz(i_uSeq,:)),'k') %plot line connecting tones
                hold on;
                scatter(34,uniqueToneSeq{i_tonedur}.PredictedFinalTonePitch_logF(i_uSeq,:),'o', 'filled', 'r')  %plot p*34
                for i_predp34 = 1:size(uniqueToneSeq{i_tonedur}.PresentedFinalTonePitch_logF,2)
                    hold on
                    scatter(34,uniqueToneSeq{i_tonedur}.PresentedFinalTonePitch_logF(i_uSeq,i_predp34)','o', 'filled', 'b')  %plot p34
                end
                hold on
                plot(1:length(DataClean_AllTrials.behav.stim.series_f{i_seq}),ones(1,34)*log(440),'k--') %plot middle line
                xlim([1 34])
                ylim([log(220) log(880)])
                set(gca,'YTick',[log(220) log(440) log(880)])
                yticklabels({'220 Hz', '440 Hz', '880 Hz'})
                title(['Beta: ' num2str(uniqueToneSeq{i_tonedur}.Beta(i_uSeq)) ...
                    '; predicted final tone pitch = ' uniqueToneSeq{i_tonedur}.PredictedFinalTonePitch{i_uSeq}])
                
                filename     = [path_fig ...
                    'ToneSeq' num2str(i_uSeq) '_ToneDur' ToneDur_text{i_tonedur} ...
                    's_Beta' num2str(uniqueToneSeq{i_tonedur}.Beta(i_uSeq)) ...
                    '_predp34' uniqueToneSeq{i_tonedur}.PredictedFinalTonePitch{i_uSeq} '.png'];
                saveas(gcf, [filename], 'png'); %save png version
                close;
            end
                        
            filename = [path_outputdata ...
                'ToneSeq' num2str(i_uSeq) '_ToneDur' ToneDur_text{i_tonedur} ...
                's_Beta' num2str(uniqueToneSeq{i_tonedur}.Beta(i_uSeq)) ...
                '_predp34' uniqueToneSeq{i_tonedur}.PredictedFinalTonePitch{i_uSeq}];
            audiowrite([filename '.wav'],uniqueToneSeq{i_tonedur}.ToneSeries_Soundwave(i_uSeq,:),44100)
        end
    end
end

save([path_outputdata, 'AllToneSeq.mat'], 'uniqueToneSeq', '-v7.3');

