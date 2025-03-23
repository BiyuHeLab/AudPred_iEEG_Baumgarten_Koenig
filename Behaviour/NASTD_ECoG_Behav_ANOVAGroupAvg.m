function [stats_anova2] = NASTD_ECoG_Behav_ANOVAGroupAvg(subs, tonedur_text, vars)

%% Group-level ANOVAs for FTPLestimation
%Script computes 1) across subject average ANOVAs for FTPLestimates
    %1way ANOVA (dep var: FTPLrating, factors/IV: presented FTP (p34))
    %2way ANOVA (dep var: FTPLrating, factors/IV: presented FTP (p34), predicted FTP (p*34))

%% 0) Specify paths
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')

NASTD_ECoG_setVars
paths_NASTD_ECoG = NASTD_ECoG_paths;

%Add base dir and own script dir
addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));


%% 1) load behavioral data, specify paths
for i_sub = 1:length(subs)
    
    sub = subs{i_sub};
    
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    
    load([si.path_behavioraldata_sub]);%Load raw behavioral data
    
    for i_tonedur = 1:length(tonedur_text)
        
        tonedur = tonedur_text{i_tonedur};
        
        %% 2) check counterbalancing (i.e., if all stim combinations present per block)
        
        trials_by_block = reshape(1:120, [10 12]);
        for i_block = 1:12
            
            trials = trials_by_block(:, i_block);
            
            count = zeros(2,3);
            
            for i_trial = trials
                count(data.stim.toneDurID(i_trial),... %tone duration conditions
                    data.stim.predID(i_trial)+2) = ... %label for math. expected FTP (-1 = low, 0 = medium, 1 = high)
                    count(data.stim.toneDurID(i_trial),...
                    data.stim.predID(i_trial)+2) + 1;
            end
            
            if all(count(:) == 1)
                block_check(i_block) = 1;
            else
                block_check(i_block) = 0;
                disp('bad counterbalancing')
            end
        end
        
        
        %% 3.1) define filters for trial selection
        
        % 3.1.1 general filter - ensure all behavioral responses were entered properly
        filter_resp = data.resp_prob >= 0 & data.trialNum > 0;
        %This takes care of the -1 responses in the final tone likelihood estimates
        filter_rt   = data.r1_RT > 0;
        
        % 3.1.2 Tone duration filter
        switch tonedur
            case '0.2'
                filter_toneDur = data.stim.toneDur == 0.2;
                
            case '0.4'
                filter_toneDur = data.stim.toneDur == 0.4;
                
            otherwise
                filter_toneDur = ones(1,120);
        end
        
        % 3.1.3 manually defined filters
        filter_manual = ones(1,120);
        
        % 3.1.4 combined filter
        filter = filter_resp & filter_rt & filter_toneDur & filter_manual; %filter var for trial selection
        
        %% 3.2) apply filters
        % filtered data will be applied to newly defined structs "dataf" and
        % "stimf", since these are the only relevant subfields in behavioral data
        
        %define new struct 'dataf' containing all data subfields with behavioral responses but only for above-filtered trials
        %from data subfields
        dfn = fieldnames(data);
        for j = 1:length(dfn)
            if eval(['length(data.' dfn{j} ') == 120'])
                eval(['dataf.' dfn{j} ' = data.' dfn{j} '(filter);']);
            end
        end
        
%         if i_sub < 4
%             %from data.uSeq (i.e., task2) subfields
%             dfn = fieldnames(data.uSeq);
%             for j = 1:length(dfn)
%                 if eval(['length(data.uSeq.' dfn{j} ') == 120'])
%                     if strcmp(dfn{j},'disp3') %special because 3 columns instead of one
%                         eval(['dataf.' dfn{j} ' = data.uSeq.' dfn{j} '(filter,:);']);
%                     else
%                         eval(['dataf.' dfn{j} ' = data.uSeq.' dfn{j} '(filter);']);
%                     end
%                 end
%             end
%             dataf = rmfield(dataf,{'resp_beta','conf_beta', 'correct_beta', 'diff_beta', 'r3_RT'});
%         end
        %define new struct 'stimf' containing all data subfields with all stimulus-wise parameters but only for above-filtered trials
        sfn = fieldnames(data.stim);
        for j = 1:length(sfn)
            if eval(['length(data.stim.' sfn{j} ') == 120'])
                eval(['stimf.' sfn{j} ' = data.stim.' sfn{j} '(filter);']);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% 4. Create data struct for ANOVA
        %% 4.1 data struct for 1way ANOVA (dep var: FTPLrating, factors/IV: presented FTP (p34))
        %empty placeholder vectors, later filled single subject data
        vec_FTPLrating = []; %subject response/final tone likelihood rating (how well does the
        %final tone pitch fit in given the previously presented sequence
        %(i.e., distance between presented and math. expected final tone;
        %1 - very unlikely, 5 - very likely))
        vec_p34 = []; %actually presented final tone frequency
        %([-3:3],negative values = < 440 Hz, positive values = > 440 Hz)
        
        if strcmp(subs(i_sub),'NY688')
            fprintf('Stop - SubNY688 has different experimental parameters - not comparable with other subjects')
            pause
        else
            
            p34 = {'-2' '-1' '1' '2'}; %actually presented final tone pitches
            if strcmp(tonedur,'all')
                %rowcount determined by [1 - (number of trials per condition))]
                rowCount = -29; % all tone dur = 30 trials/condition (4 conditions)
            else
                rowCount = -14; %for 0.2/0.4s tone dur = 15 trials/condition
            end
            
        end
        
        for i_p34 = 1:length(p34) %p34 loop
            actual = p34{i_p34};
            
            switch tonedur %filter tone duration
                case '0.2'
                    filter_toneDur = data.stim.toneDur == 0.2;
                case '0.4'
                    filter_toneDur = data.stim.toneDur == 0.4;
                otherwise
                    filter_toneDur = ones(1,120);
            end
            
            switch actual %filter actually presented final tone
                %Note: Only NY688 has 6 options, all other have ony 4
                case '-3'
                    filter_actual = data.stim.finalID == -3;
                case '-2'
                    filter_actual = data.stim.finalID == -2;
                case '-1'
                    filter_actual = data.stim.finalID == -1;
                case '1'
                    filter_actual = data.stim.finalID == 1;
                case '2'
                    filter_actual = data.stim.finalID == 2;
                case '3'
                    filter_actual = data.stim.finalID == 3;
            end
            
            filter = filter_toneDur & filter_actual;
            
            %Place trial selection in new structure used for MATLAB ANOVA
            vec_FTPLrating = [vec_FTPLrating; data.resp_prob(filter ==1)'];
            %vector reading out the subject's final tone likelihood rating for each trial in
            %order determined by above loops (i.e., hierarchical order)
            
            nTrialsPerCond = length(data.resp_prob(filter ==1));
            
            %Update row
            rowCount = rowCount + nTrialsPerCond;
            lastRow = rowCount + nTrialsPerCond - 1;
            
            %Copy order of loop vars into vectors (i.e., hierarchical sorting for
            %1) tone duration, 2) math. expected final tone, 3) actually presented
            %final tone
            [vec_ToneDur(rowCount: lastRow)] = deal(i_tonedur);
            [vec_p34(rowCount: lastRow)] = deal(i_p34);
        end
        
        %Output: single-subject cell containing vectors with all trials,
        %hiearchically sorted by 1) tone dur (1:3), 2) Expected tone (1:3), 3) Actual tone (1:6)
        FTPLrating_anova1{i_sub,i_tonedur}.FTPLrating = vec_FTPLrating;
        FTPLrating_anova1{i_sub,i_tonedur}.ToneDur = vec_ToneDur';
        FTPLrating_anova1{i_sub,i_tonedur}.p34 = vec_p34';
        
        %Question: vecFTPLrating should only contain [1:5], but in some enties
        %there are -1. This is mentioned in th Group-level ANOVA script,
        %where -1 entries are then replaced with 1. Presumably, -1 entries
        %resemble missed responses. Therefore, I now replace them with NaNs.
        ii = find(FTPLrating_anova1{i_sub,i_tonedur}.FTPLrating == -1);
        FTPLrating_anova1{i_sub,i_tonedur}.FTPLrating(ii) = NaN;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% 4.2 Data Struct for 2way ANOVA (dep var: FTPLrating, factors/IV: presented FTP (p34),
        %math. expected FTP (p*34))
        vec_FTPLrating = []; %subject response/final tone likelihood rating (how well does the
        %final tone pitch fit in given the previously presented sequence
        %(i.e., distance between presented and math. expected final tone;
        %1 - very unlikely, 5 - very likely))
        vec_Predp34 = []; %math. expected final tone (-1 = low, 0 = medium, 1 = high)
        vec_p34 = []; %actually presented final tone frequency
        %([-3:3],negative values = < 440 Hz, positive values = > 440 Hz)
        
        %empty placeholder vectors, later filled single subject data
        Predp34 = {'low' 'med' 'high'}; %math. expected final tone pitches
        if strcmp(subs(i_sub),'NY688')
            fprintf('Stop - SubNY688 has different experimental parameters - not comparable with other subjects')
            pause
            
        else
            p34 = {'-2' '-1' '1' '2'}; %actually presented final tone pitches
            
            if strcmp(tonedur,'all')
                rowCount = -9; %for all tone dur = 10 trials/condition
            else
                rowCount = -4; %for 0.2/0.4s tone dur = 5 trials/condition
            end
            
        end
        
        for k = 1:length(Predp34) %p*34 loop
            expected = Predp34{k};
            
            for i_p34 = 1:length(p34) %p34 loop
                actual = p34{i_p34};
                
                switch tonedur %filter tone duration
                    case '0.2'
                        filter_toneDur = data.stim.toneDur == 0.2;
                    case '0.4'
                        filter_toneDur = data.stim.toneDur == 0.4;
                    otherwise
                        filter_toneDur = ones(1,120);
                end
                
                switch expected %filter math. expected final tone
                    case 'low'
                        filter_expected = data.stim.predID == -1;
                    case 'med'
                        filter_expected = data.stim.predID == 0;
                    case 'high'
                        filter_expected = data.stim.predID == 1;
                end
                
                switch actual %filter actually presented final tone
                    %Note: Only NY688 has 6 options, all other have ony 4
                    case '-3'
                        filter_actual = data.stim.finalID == -3;
                    case '-2'
                        filter_actual = data.stim.finalID == -2;
                    case '-1'
                        filter_actual = data.stim.finalID == -1;
                    case '1'
                        filter_actual = data.stim.finalID == 1;
                    case '2'
                        filter_actual = data.stim.finalID == 2;
                    case '3'
                        filter_actual = data.stim.finalID == 3;
                end
                
                filter = filter_toneDur & filter_expected & filter_actual;
                
                %Place trial selection in new structure used for MATLAB ANOVA
                vec_FTPLrating = [vec_FTPLrating; data.resp_prob(filter ==1)'];
                %vector reading out the subject's final tone likelihood rating for each trial in
                %order determined by above loops (i.e., hierarchical order)
                
                nTrialsPerCond = length(data.resp_prob(filter ==1));
                
                %Update row
                rowCount = rowCount + nTrialsPerCond;
                lastRow = rowCount + nTrialsPerCond - 1;
                
                %Copy order of loop vars into vectors (i.e., hierarchical sorting for
                %1) tone duration, 2) math. expected final tone, 3) actually presented
                %final tone
                [vec_ToneDur(rowCount: lastRow)] = deal(i_tonedur);
                [vec_Predp34(rowCount: lastRow)] = deal(k);
                [vec_p34(rowCount: lastRow)] = deal(i_p34);
            end
            
        end
        
        %Output: single-subject cell containing vectors with all trials,
        %hiearchically sorted by 1) tone dur (1:3), 2) Expected tone (1:3), 3) Actual tone (1:6)
        FTPLrating_anova2{i_sub,i_tonedur}.FTPLrating = vec_FTPLrating;
        FTPLrating_anova2{i_sub,i_tonedur}.ToneDur = vec_ToneDur';
        FTPLrating_anova2{i_sub,i_tonedur}.Predp34 = vec_Predp34';
        FTPLrating_anova2{i_sub,i_tonedur}.p34 = vec_p34';
        
        %Question: vecFTPLrating should only contain [1:5], but in some enties
        %there are -1. This is mentioned in th Group-level ANOVA script,
        %where -1 entries are then replaced with 1. Presumably, -1 entries
        %resemble missed responses. Therefore, I now replace them with NaNs.
        ii = find(FTPLrating_anova2{i_sub,i_tonedur}.FTPLrating == -1);
        FTPLrating_anova2{i_sub,i_tonedur}.FTPLrating(ii) = NaN;
        
        clear vec* subDurTbl rowCount2 lastRow2
        
    end
end

%% 5 Concatenate trials across subjects per tone duration condition and compute respective ANOVA
% 5.1 1way ANOVA (dep var: FTPLrating, factors/IV: presented FTP (p34))
for i_tone = 1:length(tonedur_text) %Separate RM ANOVA for each tone duration
    
    %Create empty proxies
    GroupVec_FTPLrating{i_tone} = [];
    GroupVec_p34{i_tone} = [];
    GroupVec_sub{i_tone} = [];
    
    for i_sub = 1:length(subs) %Concatenate all trials for all subs per tone condition and fill proxies
        GroupVec_FTPLrating{i_tone}  = [GroupVec_FTPLrating{i_tone};  FTPLrating_anova1{i_sub,i_tone}.FTPLrating];
        GroupVec_p34{i_tone}      = [GroupVec_p34{i_tone};      FTPLrating_anova1{i_sub,i_tone}.p34];
        ntrials                 = length(FTPLrating_anova2{i_sub,i_tone}.FTPLrating);
        GroupVec_sub{i_tone}          = [GroupVec_sub{i_tone};          repelem(i_sub,ntrials)'];
    end
    
    %compute average FTPL rating per p43 for plotting
    for i_sub = 1:length(subs) %Concatenate all trials for all subs per tone condition and fill proxies
        i_p34_1 = find(FTPLrating_anova1{i_sub,i_tone}.p34 == 1);
        AvgFTPLrating_perp34{i_tone}(i_sub,1) = mean(FTPLrating_anova1{i_sub,i_tone}.FTPLrating(i_p34_1));
        i_p34_2 = find(FTPLrating_anova1{i_sub,i_tone}.p34 == 2);
        AvgFTPLrating_perp34{i_tone}(i_sub,2) = mean(FTPLrating_anova1{i_sub,i_tone}.FTPLrating(i_p34_2));  
        i_p34_3 = find(FTPLrating_anova1{i_sub,i_tone}.p34 == 3);
        AvgFTPLrating_perp34{i_tone}(i_sub,3) = mean(FTPLrating_anova1{i_sub,i_tone}.FTPLrating(i_p34_3)); 
        i_p34_4 = find(FTPLrating_anova1{i_sub,i_tone}.p34 == 4);
        AvgFTPLrating_perp34{i_tone}(i_sub,4) = mean(FTPLrating_anova1{i_sub,i_tone}.FTPLrating(i_p34_4));         
    end
    
    figure; plot(1:4,AvgFTPLrating_perp34{i_tone},'.','MarkerSize',15)
    hold on;
    plot(1:4,AvgFTPLrating_perp34{i_tone})
    xlim([0 5])
    xlabel('p34')
    ylabel('Avg FTLrating')
    title(['Average FTPLrating per sub across p34 for TD: ' tonedur_text{i_tone} 's'])
    %Repeated Measures Two-way Analysis of Variance Test
    %dependent var: final tone likelihood rating;
    %factors: math. expected FTP, presented FTP
    %CAVE: results don't agree with SPSS rmANOVA
%     RMAOV1([GroupVec_FTPLrating{i_tone}, GroupVec_p34{i_tone}, GroupVec_sub{i_tone}-1]);
        
end

% 5.2 2way ANOVA (dep var: FTPLrating, factors/IV: expected FTP (p*34), presented FTP (p34),
for i_tone = 1:length(tonedur_text) %Separate RM ANOVA for each tone duration
    
    %Create empty proxies
    GroupVec_FTPLrating{i_tone} = [];
    GroupVec_Predp34{i_tone} = [];
    GroupVec_p34{i_tone} = [];
    GroupVec_sub{i_tone} = [];
    
    for i_sub = 1:length(subs) %Concatenate all trials for all subs per tone condition and fill proxies
        GroupVec_FTPLrating{i_tone}  = [GroupVec_FTPLrating{i_tone};  FTPLrating_anova2{i_sub,i_tone}.FTPLrating];
        GroupVec_Predp34{i_tone} 	= [GroupVec_Predp34{i_tone};    FTPLrating_anova2{i_sub,i_tone}.Predp34];
        GroupVec_p34{i_tone}      = [GroupVec_p34{i_tone};      FTPLrating_anova2{i_sub,i_tone}.p34];
        ntrials                 = length(FTPLrating_anova2{i_sub,i_tone}.FTPLrating);
        GroupVec_sub{i_tone}          = [GroupVec_sub{i_tone};          repelem(i_sub,ntrials)'];
    end
    
    %compute average FTPL rating per p43 and p*34 for plotting
    for i_sub = 1:length(subs) %Concatenate all trials for all subs per tone condition and fill proxies
        for i_Predp34 = 1:3
            for i_p34 = 1:4
                filt_Predp = FTPLrating_anova2{i_sub,i_tone}.Predp34 == i_Predp34;
                filt_p = FTPLrating_anova2{i_sub,i_tone}.p34 == i_p34;
                filt = filt_Predp & filt_p;
                    AvgFTPLrating_perPredp34_p34{i_tone}{i_Predp34}(i_sub,i_p34) = mean(FTPLrating_anova2{i_sub,i_tone}.FTPLrating(filt));
            end
        end      
    end
    
    test_AvgFTPLrating_perPredp34_p34{i_tone} = [];
    filt_sub = [];
    filt_p34 = [];
    filt_predp34 = [];
    for i_Predp34 = 1:3
        for i_p34 = 1:4
            test_AvgFTPLrating_perPredp34_p34{i_tone} = ...
                [test_AvgFTPLrating_perPredp34_p34{i_tone}; ...
                AvgFTPLrating_perPredp34_p34{i_tone}{i_Predp34}(:,i_p34)];
            filt_sub = [filt_sub; [1:8]'];
            filt_p34 = [filt_p34; ones(8,1)*i_p34];
            filt_predp34 = [filt_predp34; ones(8,1)*i_Predp34];
        end
    end
        
    figure;
    for i_Predp34 = 1:3
        hold on;
        plot(1:4,mean(AvgFTPLrating_perPredp34_p34{i_tone}{i_Predp34}),'.','MarkerSize',15)
        hold on;
        plot(1:4,mean(AvgFTPLrating_perPredp34_p34{i_tone}{i_Predp34}))
    end
    
    xlim([0 5])
    xlabel('p34')
    ylabel('Avg FTLrating')
    title(['Average FTPLrating per sub across Predp34 and p34 for TD: ' tonedur_text{i_tone} 's'])

    
    %Repeated Measures Two-way Analysis of Variance Test
    %dependent var: final tone likelihood rating;
    %factors: math. expected FTP, presented FTP
%     RMAOV2([GroupVec_FTPLrating{i_tone}, GroupVec_Predp34{i_tone}, GroupVec_p34{i_tone}, GroupVec_sub{i_tone}-1]);
%     stats_anova2{i_tone} = rm_anova2(GroupVec_FTPLrating{i_tone}, GroupVec_sub{i_tone}, GroupVec_Predp34{i_tone}, GroupVec_p34{i_tone}, {'math.expected tone','presented tone'});
    stats_anova2{i_tone} = rm_anova2(test_AvgFTPLrating_perPredp34_p34{i_tone}, ...
        filt_sub, filt_predp34, filt_p34, ...
        {'p*34','p34'})
    %IV1: math. expected FTP
    %IV2: actual FTP
    
end

end
