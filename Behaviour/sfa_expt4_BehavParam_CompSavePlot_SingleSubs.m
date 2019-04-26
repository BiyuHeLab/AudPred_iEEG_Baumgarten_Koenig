%1) Script computes single-subject behavioral data 
%(final tone likelihood estimate, uSeqID) 
%for all subjects ( only the finally selected ones)
%2) Script saves behavioral output vars
%3) Script plots singe-subject behavioral parameters

%%Note: Think about inclding 3way rmANOVA (tone dur, p*34, p34)
clear 
%% 0) Specify paths
% addpath('//gogo.sb.nyumc.org/data/gogodisk4/thomas/AuditoryPrediction/iEEG/')
location = 'server'; %gago/gogo
addpath('/data/gogodisk4/thomas/AuditoryPrediction/iEEG/')

% location = 'desktop'; %local on office desktop
% addpath('//gogo.sb.nyumc.org/data/gogodisk4/thomas/AuditoryPrediction/iEEG/')
sfa_expt4_setVars
paths_sfa_expt4 = sfa_expt4_paths(location);

%Add base dir and own script dir
addpath(genpath(paths_sfa_expt4.BaseDir));
addpath(genpath(paths_sfa_expt4.ScriptsDir));

    global ft_defaults;
     ft_default.trackusage = 'no';

path_data = paths_sfa_expt4.RawData;
path_fig = paths_sfa_expt4.Fig_Behavior_SingleSubs;

i_sub = 1;
subs = vars.sub_list; %patients
sub = subs{i_sub};

tonedur_text = {'0.2' '0.4' 'all'};

saveplot = 0; %Save Figures?

 %% 1) load behavioral data, specify paths
for i_sub = 1:length(subs)
     
    sub = subs{i_sub}

    sfa_expt4_subjectinfo %load subject info file (var: si)

load([si.path_behavioraldata_sub]);%Load raw behavioral data 
% load([si.path_stimorder_sub]);%Load stimulus presentation order 

BehavData_SingleSub = []; %Create empty proxy file
BehavData_SingleSub.tonedur = tonedur_text; %specify current tone dur condition


for i_tonedur = 1:length(tonedur_text)

tonedur = tonedur_text{i_tonedur};
path_fig_sub = [path_fig sub '/'];

% switch tonedur
%     case '0.2'
%         path_fig_sub = [path_fig_sub 'toneDur=0.2/'];
% 
%     case '0.4'
%         path_fig_sub = [path_fig_sub 'toneDur=0.4/'];
% 
%     case 'all'
%         path_fig_sub = [path_fig_sub 'toneDur=all/'];
% end

mkdir(path_fig_sub); %make new tone-duration-specific-directory where figures are saved

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
    if i_sub < 4 %only first 3 patients had task2
    filter_resp = data.uSeq.respABC >= 0 & data.resp_prob >= 0 & data.trialNum > 0;
    %This takes care of the -1 responses in the final tone likelihood estimates
    filter_rt   = data.r1_RT > 0 & data.r2_RT > 0;
    else
    filter_resp = data.resp_prob >= 0 & data.trialNum > 0;
    %This takes care of the -1 responses in the final tone likelihood estimates
    filter_rt   = data.r1_RT > 0;
    end        

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
    
    %from data.uSeq (i.e., task2) subfields
    if i_sub < 4 %only first 3 patients had task2    
        dfn = fieldnames(data.uSeq);
        for j = 1:length(dfn)       
            if eval(['length(data.uSeq.' dfn{j} ') == 120'])
                if strcmp(dfn{j},'disp3') %special because 3 columns instead of one
                eval(['dataf.' dfn{j} ' = data.uSeq.' dfn{j} '(filter,:);']);            
                else   
                eval(['dataf.' dfn{j} ' = data.uSeq.' dfn{j} '(filter);']);
                end
            end  
        end
        dataf = rmfield(dataf,{'resp_beta','conf_beta', 'correct_beta', 'diff_beta', 'r3_RT'});
    end
    %define new struct 'stimf' containing all data subfields with all stimulus-wise parameters but only for above-filtered trials
    sfn = fieldnames(data.stim);
    for j = 1:length(sfn)
        if eval(['length(data.stim.' sfn{j} ') == 120'])
            eval(['stimf.' sfn{j} ' = data.stim.' sfn{j} '(filter);']);
        end
    end
    
    BehavData_SingleSub.TrialNum(1,i_tonedur) = length(dataf.trialNum);

    %% 4) Final Tone Probability Rating - Task1
    %I.e., rate probability of presented final tone pitch (p34) based on previous sequence
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %4.1 Hit rate computation
    possible_ftp = unique(stimf.logf_final); %possible log freq of final tone (p34)
    predicted_ftp = unique(stimf.logf_pred); %possible log freq of predicted final tone (p*34)
    %3 different options per tone dur
    %6 different options across tone dur
    BehavData_SingleSub.presentedFTP{i_tonedur} = possible_ftp;
    BehavData_SingleSub.expectedFTP{i_tonedur} = predicted_ftp;

    predIDs = [-1 0 1]; %TJB: numerical label for predicted final tone (p*34; low/middle/high)

    for i_pred = 1:3 %index for predicted final tone

        f_exp = stimf.predID == predIDs(i_pred); %filter to select those trials with specific math. expected final tone (low, medium, high)

        for i_tf = 1:length(possible_ftp)

            f_tone  = stimf.logf_final == possible_ftp(i_tf); %filter to select only trials with specific final tone pitch
            ff      = f_exp & f_tone; %filter to select only trials with specific final tone pitch and specific p*34 trend

            %TJB: resp_prob = final tone likelihood rating (1 = very unlikely, 5 = very likely)
            rp_tot(i_tf) = mean(dataf.resp_prob(f_tone)); %mean likelihood rating for trials with specific final tone pitches
            rp_tot_sem(i_tf) = std(dataf.resp_prob(f_tone)) / sqrt(sum(f_tone));%response probability SEM

            % get response probabilities by expected final tone
            %TJB: i.e., response probabilities for specific tone depending on predicted final tone (low/middle/high)
                    %This is whats depicted in Fig. 3C: expectedness rating as a function of final tone pitch (tf) and predicted final tone (predIDs)
            rp_exp{i_pred}(i_tf) = mean(dataf.resp_prob(f_exp & f_tone));%%mean likelihood rating for trials with specific final tone pitches and specific p*34 trend (cell)
            rp_exp_sem{i_pred}(i_tf) = std(dataf.resp_prob(f_exp & f_tone)) / sqrt(sum(f_exp & f_tone));%response probability SEM
        
        end        
  
    end
    
    if strcmp(sub,'NY688') %for 1st subject NY688 (who had ony 3 actual tones per condition)       
        %Because in subject N688, p34 and p*34 are slightly different between tone duration
        %conditions, they are computed separately above, resulting in 6 values. 
        %Subsequently, we average within the low/medium/high conditions to
        %receive 3 values again.
        if strcmp(tonedur,'all')
            possible_ftp2(1,1) = mean([possible_ftp(1,1),possible_ftp(1,2)]);
            possible_ftp2(1,2) = mean([possible_ftp(1,3),possible_ftp(1,4)]);
            possible_ftp2(1,3) = mean([possible_ftp(1,5),possible_ftp(1,6)]);

            rp_tot2(1,1) = mean([rp_tot(1,1),rp_tot(1,2)]);
            rp_tot2(1,2) = mean([rp_tot(1,3),rp_tot(1,4)]);
            rp_tot2(1,3) = mean([rp_tot(1,5),rp_tot(1,6)]);

            rp_tot_sem2(1,1) = mean([rp_tot_sem(1,1),rp_tot_sem(1,2)]);
            rp_tot_sem2(1,2) = mean([rp_tot_sem(1,3),rp_tot_sem(1,4)]);
            rp_tot_sem2(1,3) = mean([rp_tot_sem(1,5),rp_tot_sem(1,6)]);

            for i_pred = 1:3 
            rp_exp2{i_pred}(1,1) = mean([rp_exp{i_pred}(1,1),rp_exp{i_pred}(1,2)]);
            rp_exp2{i_pred}(1,2) = mean([rp_exp{i_pred}(1,3),rp_exp{i_pred}(1,4)]);
            rp_exp2{i_pred}(1,3) = mean([rp_exp{i_pred}(1,5),rp_exp{i_pred}(1,6)]);

            rp_exp_sem2{i_pred}(1,1) = mean([rp_exp_sem{i_pred}(1,1),rp_exp_sem{i_pred}(1,2)]);
            rp_exp_sem2{i_pred}(1,2) = mean([rp_exp_sem{i_pred}(1,3),rp_exp_sem{i_pred}(1,4)]);
            rp_exp_sem2{i_pred}(1,3) = mean([rp_exp_sem{i_pred}(1,5),rp_exp_sem{i_pred}(1,6)]);
            end

            possible_ftp = possible_ftp2;
            rp_tot = rp_tot2;
            rp_tot_sem = rp_tot_sem2;
            rp_exp = rp_exp2;
            rp_exp_sem = rp_exp_sem2;

            clear possible_ftp2 rp_tot2 rp_tot_sem2 rp_exp2 rp_exp_sem2
        end
    end

    
BehavData_SingleSub.FTPLrating_total{i_tonedur} = rp_tot;
BehavData_SingleSub.FTPLrating_total_SEM{i_tonedur} = rp_tot_sem;
BehavData_SingleSub.FTPLrating_exp{i_tonedur} = rp_exp;
BehavData_SingleSub.FTPLrating_total_SEM{i_tonedur} = rp_exp_sem;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %4.2 1way ANOVA (dep var: FTPLrating, factors: presented FTP (p34))
    %empty placeholder vectors, later filled single subject data
    vecPredRating = []; %subject response/final tone likelihood rating (how well does the
    %final tone pitch fit in given the previously presented sequence 
    %(i.e., distance between presented and math. expected final tone;
    %1 - very unlikely, 5 - very likely))
    vecActual = []; %actually presented final tone frequency 
    %([-3:3],negative values = < 440 Hz, positive values = > 440 Hz)
    
    if strcmp(subs(i_sub),'NY688')
        actual_FTP = {'-3' '-2' '-1' '1' '2' '3'}; %actually presented final tone pitches
            if strcmp(tonedur,'all')
            %rowcount determined by [1 - (number of trials per condition))] 
            rowCount = -19; % all tone dur = 20 trials/condition (6 conditions)
            else 
            rowCount = -19; %for 0.2/0.4s tone dur = 20 trials/condition (6 conditions)
            end
    else
        
        actual_FTP = {'-2' '-1' '1' '2'}; %actually presented final tone pitches
            if strcmp(tonedur,'all')
            %rowcount determined by [1 - (number of trials per condition))] 
            rowCount = -29; % all tone dur = 30 trials/condition (4 conditions)
            else 
            rowCount = -14; %for 0.2/0.4s tone dur = 15 trials/condition
            end
        
    end
               
            for L = 1:length(actual_FTP) %p34 loop
                actual = actual_FTP{L}   
        
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
            vecPredRating = [vecPredRating; data.resp_prob(filter ==1)']; 
            %vector reading out the subject's final tone likelihood rating for each trial in
            %order determined by above loops (i.e., hierarchical order)

            nTrialsPerCond = length(data.resp_prob(filter ==1));

            %Update row
            rowCount = rowCount + nTrialsPerCond;
            lastRow = rowCount + nTrialsPerCond - 1;

            %Copy order of loop vars into vectors (i.e., hierarchical sorting for
            %1) tone duration, 2) math. expected final tone, 3) actually presented
            %final tone
            [vecToneDur{rowCount: lastRow}] = deal(tonedur);
            [vecActual{rowCount: lastRow}] = deal(actual);
            end       
        
            %Output: single-subject cell containing vectors with all trials,
            %hiearchically sorted by 1) tone dur (1:3), 2) Expected tone (1:3), 3) Actual tone (1:6)
            FTPLrating_anova1{i_sub}{i_tonedur}.PredRating = vecPredRating;
            FTPLrating_anova1{i_sub}{i_tonedur}.ToneDur = vecToneDur';
            FTPLrating_anova1{i_sub}{i_tonedur}.Actual = vecActual';

            %Question: vecPredRating should only contain [1:5], but in some enties
            %there are -1. This is mentioned in th Group-level ANOVA script,
            %where -1 entries are then replaced with 1. Presumably, -1 entries
            %resemble missed responses. Therefore, I now replace them with NaNs.
            ii = find(FTPLrating_anova1{i_sub}{i_tonedur}.PredRating == -1);
            FTPLrating_anova1{i_sub}{i_tonedur}.PredRating(ii) = NaN;
            
            varnames = {'ActualFTP'}; %variable labeling for ANOVA2
            
            [~, subDurTbl{i_sub}] = anovan(FTPLrating_anova1{i_sub}{i_tonedur}.PredRating, {FTPLrating_anova1{i_sub}{i_tonedur}.Actual},...
            'model', 'interaction', 'varnames', varnames, 'display', 'off')
    
%             FsubDur(i_sub) = subDurTbl{i_sub}(4,6); 
%             %Subject- & ToneDuration-wise F value for math.expected*presented final tone pitch interaction effect

            BehavData_SingleSub.FTPLrating_anova1{i_tonedur} = subDurTbl{i_sub};
            clear vec* subDurTbl

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %4.3 2way ANOVA (dep var: FTPLrating, factors: presented FTP (p34),
    %math. expected FTP (p*34))
    vecPredRating = []; %subject response/final tone likelihood rating (how well does the
        %final tone pitch fit in given the previously presented sequence 
        %(i.e., distance between presented and math. expected final tone;
        %1 - very unlikely, 5 - very likely))
    vecExpected = []; %math. expected final tone (-1 = low, 0 = medium, 1 = high)
    vecActual = []; %actually presented final tone frequency 
    %([-3:3],negative values = < 440 Hz, positive values = > 440 Hz)

    %empty placeholder vectors, later filled single subject data
    expected_FTP = {'low' 'med' 'high'}; %math. expected final tone pitches
    if strcmp(subs(i_sub),'NY688')
        actual_FTP = {'-3' '-2' '-1' '1' '2' '3'}; %actually presented final tone pitches 
    
        rowCount = -7; %for 0.2/0.4s tone dur = 8 trials/condition, but 
        
    else
        actual_FTP = {'-2' '-1' '1' '2'}; %actually presented final tone pitches
        
        if strcmp(tonedur,'all')
        rowCount = -9; %for all tone dur = 10 trials/condition
        else 
        rowCount = -4; %for 0.2/0.4s tone dur = 5 trials/condition
        end    
    
    end


       
        for k = 1:length(expected_FTP) %p*34 loop
            expected = expected_FTP{k}
            
            for L = 1:length(actual_FTP) %p34 loop
                actual = actual_FTP{L}   
        
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
            vecPredRating = [vecPredRating; data.resp_prob(filter ==1)']; 
            %vector reading out the subject's final tone likelihood rating for each trial in
            %order determined by above loops (i.e., hierarchical order)

            nTrialsPerCond = length(data.resp_prob(filter ==1));

            %Update row
            if rowCount < 0 %for first iteration
            rowCount = rowCount + nTrialsPerCond;
            else 
            rowCount = lastRow +1; %first entry of new iteration based on last entry of previous iteration
            %- necessary when conditions have different amount of trials,
            %as with sub NY688
            end
             rowCount2(k,L) = rowCount
            lastRow = rowCount + nTrialsPerCond - 1;
             lastRow2(k,L) = lastRow

            %Copy order of loop vars into vectors (i.e., hierarchical sorting for
            %1) tone duration, 2) math. expected final tone, 3) actually presented
            %final tone
            [vecToneDur{rowCount: lastRow}] = deal(tonedur);
            [vecExpected{rowCount: lastRow}] = deal(expected);
            [vecActual{rowCount: lastRow}] = deal(actual);
            end

        end
        
            %Output: single-subject cell containing vectors with all trials,
            %hiearchically sorted by 1) tone dur (1:3), 2) Expected tone (1:3), 3) Actual tone (1:6)
            FTPLrating_anova2{i_sub}{i_tonedur}.PredRating = vecPredRating;
            FTPLrating_anova2{i_sub}{i_tonedur}.ToneDur = vecToneDur';
            FTPLrating_anova2{i_sub}{i_tonedur}.Expected = vecExpected';
            FTPLrating_anova2{i_sub}{i_tonedur}.Actual = vecActual';

            %Question: vecPredRating should only contain [1:5], but in some enties
            %there are -1. This is mentioned in th Group-level ANOVA script,
            %where -1 entries are then replaced with 1. Presumably, -1 entries
            %resemble missed responses. Therefore, I now replace them with NaNs.
            ii = find(FTPLrating_anova2{i_sub}{i_tonedur}.PredRating == -1);
            FTPLrating_anova2{i_sub}{i_tonedur}.PredRating(ii) = NaN;
            
            varnames = {'ExpectedFTP'; 'ActualFTP'}; %variable labeling for ANOVA2
            
            [~, subDurTbl{i_sub}] = anovan(FTPLrating_anova2{i_sub}{i_tonedur}.PredRating, {FTPLrating_anova2{i_sub}{i_tonedur}.Expected FTPLrating_anova2{i_sub}{i_tonedur}.Actual},...
            'model', 'interaction', 'varnames', varnames, 'display', 'off')
    
%             FsubDur(i_sub) = subDurTbl{i_sub}(4,6); 
%             %Subject- & ToneDuration-wise F value for math.expected*presented final tone pitch interaction effect

            BehavData_SingleSub.FTPLrating_anova2{i_tonedur} = subDurTbl{i_sub};
            clear vec* subDurTbl rowCount2 lastRow2
            
%% Sequence Identification - Task 2
%i.e., select 1 of 3 presented images that visually represents the
%previously presented tone sequence
if i_sub < 4 %only first 3 patients had task2
    addpath(paths_sfa_expt4.ScriptsDir) %necessary for dprime function

    filter_resp_uSeq = data.uSeq.respABC >= 0;
    %analysis aims:
    % % - hit rate in total (120 trials)
    % Hitrate_SeqID_total = sum(data.correct_uSeq(filter_resp_uSeq))/length(data.correct_uSeq(filter_resp_uSeq));
    % [dval_SeqID_total, cval_SeqID_total] = dprime(Hitrate_SeqID_total,1-Hitrate_SeqID_perTonedur,1,2)

    % - hit rate per tone dur (60 trials)
    Hitrate_SeqID_perTonedur = sum(dataf.correct_uSeq)/length(dataf.correct_uSeq);
    % [dval_SeqID_perTonedur, cval_SeqID_perTonedur] = dprime(Hitrate_SeqID_perTonedur,1-Hitrate_SeqID_perTonedur,1,2);           
                dval_SeqID_perTonedur = norminv(Hitrate_SeqID_perTonedur) - norminv(1-Hitrate_SeqID_perTonedur); 
                %norminv = gives the area under the curve of the normal distribution, equivalent of z-transformation

    % - hit rate per individual sequence (8 trials/4 per tone dur)
        possible_uSeqID = unique(stimf.uSeqID); %possible individual sequences
        %15 different options with 8 runs in total/4 runs per tone dur

        for i_uSeqID = 1:length(possible_uSeqID) %index for unique SedIDs
            f_uSeqID = stimf.uSeqID == possible_uSeqID(i_uSeqID); %filter to select those trials with specific SeqID
            Hitrate_SeqID_peruSeqID(1,i_uSeqID) = mean(dataf.correct_uSeq(f_uSeqID));%mean hit rate
    %         [dval_SeqID_peruSeqID(1,i_uSeqID), cval_SeqID_peruSeqID(1,i_uSeqID)] = dprime(Hitrate_SeqID_peruSeqID(1,i_uSeqID),1-Hitrate_SeqID_peruSeqID(1,i_uSeqID),1,2);
        end

    BehavData_SingleSub.SeqID(1,i_tonedur) = Hitrate_SeqID_perTonedur;
    BehavData_SingleSub.SeqID_dprime(1,i_tonedur) = dval_SeqID_perTonedur;
    % BehavData.SeqID_crit(1,i_tonedur) = cval_SeqID_perTonedur;
    BehavData_SingleSub.SeqID_peruSeq{i_tonedur} = Hitrate_SeqID_peruSeqID;
end   
%% 6) Plotting single-subject behavioral results
addpath(paths_sfa_expt4.ScriptsDir) %necessary for figd function
%6.1 plot responses over time
h = figd(15,1); hold on;
set(gcf,'units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1);
plot(dataf.trialNum, dataf.resp_prob, 'b.');
xlabel('trial')
ylabel('Final Tone Pitch Likelihood Rating')
xlim([1 120])
set(gca,'YTick',[1:5])
ylim([.5 5.5])
title(['Responses over time - ' sub '- ' num2str(length(dataf.trialNum)) ' of 120 trials included - ToneDur: ' num2str(tonedur)])
if i_sub < 4 %only first 3 patients had task2
    subplot(2,1,2);
    plot(dataf.trialNum, dataf.respABC, 'b.');
    xlabel('trial')
    ylabel('Sequence Identification Response')
    xlim([1 120])
    set(gca,'YTick',[1:3])
    ylim([0.9 3.1])
end

if saveplot
    filename = [path_fig_sub sub '_' tonedur 's_' 'TimecourseResponses.png'];
    saveas(gcf, filename, 'png'); %save png version
    delete(h);
end


%6.2 plot RTs over time
h = figd(15,1); hold on;
set(gcf,'units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1);
plot(dataf.trialNum, dataf.r1_RT, 'b-');
xlabel('trial')
ylabel('RT for FTPL rating [s]')
xlim([1 120])
title(['Response Time - ' sub '- ' num2str(length(dataf.trialNum)) ' of 120 trials included - ToneDur: ' num2str(tonedur)])
    ylim([0 5])
    
if i_sub < 4 %only first 3 patients had task2
    subplot(2,1,2);
    plot(dataf.trialNum, dataf.r2_RT, 'b-');
    xlabel('trial')
    ylabel('RT for Sequence Identification [s]')
    xlim([1 120])
    if i_sub < 3 %from patient 3 onward we increase response time for task2 to 10sec
        ylim([0 5])
    else
        ylim([0 10])
    end
end

if saveplot
    filename = [path_fig_sub sub '_' tonedur 's_ResponseTimes.png'];
    saveas(gcf, filename, 'png'); %save png version
    delete(h);
end


%6.3 plot final tone probability rating independent of and dependent on math. expected tone
% xlab = {'539' '586' '655'};
if strcmp(sub,'NY688') %for 1st subject NY688 (who had ony 3 actual tones per condition)
    
 a = {num2str(exp(possible_ftp(1))), num2str(exp(possible_ftp(2))), num2str(exp(possible_ftp(3)))};
    xlab = {[a{1} ' Hz/low'],[a{2} ' Hz/medium'],[a{3} ' Hz/high']};
    pred_colors = {'rs', 'ks', 'gs'};

    h = figd(15,2,5);
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    subplot(2,1,1);
    errorbar(possible_ftp, rp_tot, rp_tot_sem, 'bo-','Linewidth',2);
    xlabel('Presented Final Tone Pitch')
    ylabel('Final Tone Pitch Likelihood Rating')
    ylim([1 5])
    % xlim([log(200) log(900)])
    xlim([possible_ftp(1)-0.2 possible_ftp(3)+0.2])
    set(gca,'XTick',possible_ftp)
    set(gca,'XTickLabel', xlab)
    title(['FTPL rating - ' sub '- ' num2str(length(dataf.trialNum)) ' of 120 trials included - ToneDur: ' tonedur]);

    %include inset with ANOVA results:
    inset_axis = axes('position', [0.15    0.625    0.5    0.05]);
    axis off;
    box on
    results_anova1 = ['Main effect p34: F = ' num2str(BehavData_SingleSub.FTPLrating_anova1{i_tonedur}{2,6}) '; p = ' num2str(BehavData_SingleSub.FTPLrating_anova1{i_tonedur}{2,7})];
    text(0.15,0.625,results_anova1)

    set(0,'CurrentFigure',h);
    subplot(2,1,2);
    etitle = {'math.expected FTP = low', 'm.e. FTP = med', 'm.e. FTP  = high'};
    colors = {'ro-' 'ko-' 'go-'};

    for i_pred = 1:3
        hold on;
        errorbar(possible_ftp, rp_exp{i_pred}, rp_exp_sem{i_pred}, colors{i_pred},'Linewidth',2);
        xlabel('Presented Final Tone Pitch')
        ylabel('Final Tone Pitch Likelihood Rating')
        ylim([1 5])
    %     xlim([log(200) log(900)])
        xlim([possible_ftp(1)-0.2 possible_ftp(3)+0.2])
       set(gca,'XTick',possible_ftp)
        set(gca,'XTickLabel', xlab)
        title(['FTPL rating as function of math. expected FTP - ' sub '- ' num2str(length(dataf.trialNum)) ' of 120 trials included - ToneDur: ' tonedur]);

        hl = legend(etitle); 
        set(hl, 'box', 'off','Location','SouthEast')   

    end
        hold on
        for i_pred = 1:3
            plot(predicted_ftp(i_pred), 2, pred_colors{i_pred},'LineWidth', 3,'MarkerSize', 10)
            hold on
        end

        %include inset with ANOVA results:
        inset_axis = axes('position', [0.15    0.125    0.5    0.05]);
        axis off;
        box on
        MainEffect = ['Main effect p*34: F = ' num2str(BehavData_SingleSub.FTPLrating_anova2{i_tonedur}{2,6}) ...
            '; p = ' num2str(BehavData_SingleSub.FTPLrating_anova2{i_tonedur}{2,7})]

        InteractionEffect = ['Interaction effect p34*p*34: F = ' num2str(BehavData_SingleSub.FTPLrating_anova2{i_tonedur}{4,6}) ...
            '; p = ' num2str(BehavData_SingleSub.FTPLrating_anova2{i_tonedur}{4,7})];
        text(0.15,0.625,{MainEffect,InteractionEffect})
        
else %for every other subject (who had 4 p34 tones)    
    a = {num2str(exp(possible_ftp(1))), num2str(exp(possible_ftp(2))), num2str(exp(possible_ftp(3))),num2str(exp(possible_ftp(4)))};
    xlab = {[a{1} ' Hz/low'],[a{2} ' Hz/medium1'],[a{3} ' Hz/medium2'], [a{4} ' Hz/high']};
    pred_colors = {'rs', 'ks', 'gs'};

    h = figd(15,2,5);
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    subplot(2,1,1);
    errorbar(possible_ftp, rp_tot, rp_tot_sem, 'bo-','Linewidth',2);
    xlabel('Presented Final Tone Pitch')
    ylabel('Final Tone Pitch Likelihood Rating')
    ylim([1 5])
    % xlim([log(200) log(900)])
    xlim([possible_ftp(1)-0.2 possible_ftp(4)+0.2])
    set(gca,'XTick',possible_ftp)
    set(gca,'XTickLabel', xlab)
    title(['FTPL rating - ' sub '- ' num2str(length(dataf.trialNum)) ' of 120 trials included - ToneDur: ' tonedur]);

    %include inset with ANOVA results:
    inset_axis = axes('position', [0.15    0.625    0.5    0.05]);
    axis off;
    box on
    results_anova1 = ['Main effect p34: F = ' num2str(BehavData_SingleSub.FTPLrating_anova1{i_tonedur}{2,6}) '; p = ' num2str(BehavData_SingleSub.FTPLrating_anova1{i_tonedur}{2,7})];
    text(0.15,0.625,results_anova1)

    set(0,'CurrentFigure',h);
    subplot(2,1,2);
    etitle = {'math.expected FTP = low', 'm.e. FTP = med', 'm.e. FTP  = high'};
    colors = {'ro-' 'ko-' 'go-'};

    for i_pred = 1:3
        hold on;
        errorbar(possible_ftp, rp_exp{i_pred}, rp_exp_sem{i_pred}, colors{i_pred},'Linewidth',2);
        xlabel('Presented Final Tone Pitch')
        ylabel('Final Tone Pitch Likelihood Rating')
        ylim([1 5])
    %     xlim([log(200) log(900)])
        xlim([possible_ftp(1)-0.2 possible_ftp(4)+0.2])
       set(gca,'XTick',possible_ftp)
        set(gca,'XTickLabel', xlab)
        title(['FTPL rating as function of math. expected FTP - ' sub '- ' num2str(length(dataf.trialNum)) ' of 120 trials included - ToneDur: ' tonedur]);

        hl = legend(etitle); 
        set(hl, 'box', 'off','Location','SouthEast')   

    end
        hold on
        for i_pred = 1:3
            plot(predicted_ftp(i_pred), 2, pred_colors{i_pred},'LineWidth', 3,'MarkerSize', 10)
            hold on
        end

        %include inset with ANOVA results:
        inset_axis = axes('position', [0.15    0.125    0.5    0.05]);
        axis off;
        box on
        MainEffect = ['Main effect p*34: F = ' num2str(BehavData_SingleSub.FTPLrating_anova2{i_tonedur}{2,6}) ...
            '; p = ' num2str(BehavData_SingleSub.FTPLrating_anova2{i_tonedur}{2,7})]

        InteractionEffect = ['Interaction effect p34*p*34: F = ' num2str(BehavData_SingleSub.FTPLrating_anova2{i_tonedur}{4,6}) ...
            '; p = ' num2str(BehavData_SingleSub.FTPLrating_anova2{i_tonedur}{4,7})];
        text(0.15,0.625,{MainEffect,InteractionEffect})
end    
       
if saveplot
    filename = [path_fig_sub sub '_' tonedur '_FTPLrating_Overview.png'];
    saveas(gcf, filename, 'png'); %save png version
    delete(h);
end

end
% %6.4 plot final tone probability rating as function of math. expected final
% %tone
% h = figd(25,2,5);
% etitle = {'math.expected FTP = low', 'm.e. FTP = med', 'm.e. FTP  = high'};
% colors = {'ro-' 'ko-' 'go-'};
% 
% 
% for i_pred = 1:3
%     hold on;
%     errorbar(possible_ftp, rp_exp{i_pred}, rp_exp_sem{i_pred}, colors{i_pred},'Linewidth',2);
%     xlabel('Final Tone Pitch (Hz)')
%     ylabel('Final Tone Pitch Likelihood Rating')
%     ylim([1 5])
%     xlim([log(200) log(900)])
%     set(gca,'XTick',possible_ftp)
%     set(gca,'XTickLabel', xlab)
%     title(['FTPL rating as function of math. expected FTP - ' sub '- ' num2str(length(dataf.trialNum)) ' of 120 trials included - ToneDur: ' tonedur]);
% 
%     hl = legend(etitle); 
%     set(hl, 'box', 'off','Location','SouthEast')
% end
% 
% if saveplot
%     filename = [path_fig_sub sub '_' tonedur 'FTPLErating_permathexpecFTP.png'];
%     saveas(gcf, filename, 'png'); %save png version
%     delete(h);
% end
%6.5 plot Sequence Identification (task2) for all tone dur
if i_sub < 4 %only first 3 patients had task2
    xlab = {'0.2 ms' '0.4 ms' 'all'};

    h = figd(15,2,5);
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    subplot(4,1,1);
    plot(BehavData_SingleSub.SeqID, '.b','MarkerSize',40,'Linewidth',2);   
        xlabel('Tone Duration')
        ylabel('HitRate')
        ylim([-0.1 1.1])
        xlim([0.85 3.15])
        set(gca,'XTick',1:3)
        set(gca,'XTickLabel', xlab)
        title(['Hit Rate (ratio correct responses/all responses) for Sequence Identification  - ' sub]);
        hold on; line([0.85:3.85],[0.33,0.33,0.33,0.33],'LineStyle','--','Color','k','Linewidth',1);

        inset_axis = axes('position', [0.115    0.86    0.5    0.05]);
        axis off;
        box on
        dprime1 = ['d-prime: ' num2str(BehavData_SingleSub.SeqID_dprime(1))];
        text(0.068,0.9,dprime1)
        dprime2 = ['d-prime: ' num2str(BehavData_SingleSub.SeqID_dprime(2))];
        text(0.745,0.9,dprime2)    
        dprime3 = ['d-prime: ' num2str(BehavData_SingleSub.SeqID_dprime(3))];
        text(1.41, 0.9,dprime3)     

    xlab = {'Seq1','Seq2','Seq3','Seq4','Seq5','Seq6','Seq7','Seq8','Seq9','Seq10','Seq11','Seq12','Seq13','Seq14','Seq15'};   
    subplot(4,1,2);
    plot(BehavData_SingleSub.SeqID_peruSeq{1},'.b','MarkerSize',20,'Linewidth',2);   
    %     xlabel('Individual Sequence - 0.2s ToneDur')
        ylabel('HitRate')
        ylim([-0.1 1.1])
        xlim([0.85 15.15])
        set(gca,'XTick',1:15)
        set(gca,'XTickLabel', xlab)
        title(['Hit Rate for unique sequence SeqID - 0.2s ToneDur - AvgTrialsperSeq: ' num2str(ceil(BehavData_SingleSub.TrialNum(1)/15))]);
        hold on; line([0.85:15.85],[0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33],'LineStyle','--','Color','k','Linewidth',1);   

    xlab = {'','','','','','','','','','','','','','',''};      
    subplot(4,1,3);
    plot(BehavData_SingleSub.SeqID_peruSeq{2},'.b','MarkerSize',20,'Linewidth',2);   
    %     xlabel('Individual Sequence - 0.4s ToneDur')
        ylabel('HitRate')
        ylim([-0.1 1.1])
        xlim([0.85 15.15])
    %     set(gca,'XTick',1:15)
         set(gca,'XTickLabel', xlab)
        title(['Hit Rate for unique sequence SeqID - 0.4s ToneDur - AvgTrialsperSeq: ' num2str(ceil(BehavData_SingleSub.TrialNum(2)/15))]);
        hold on; line([0.85:15.85],[0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33],'LineStyle','--','Color','k','Linewidth',1);   

    subplot(4,1,4);
    plot(BehavData_SingleSub.SeqID_peruSeq{3},'.b','MarkerSize',20,'Linewidth',2);   
    %     xlabel('Individual Sequence - All ToneDur')
        ylabel('HitRate')
        ylim([-0.1 1.1])
        xlim([0.85 15.15])
    %     set(gca,'XTick',1:15)
         set(gca,'XTickLabel', xlab)
        title(['Hit Rate for unique sequence SeqID - All ToneDur - AvgTrialsperSeq: ' num2str(ceil(BehavData_SingleSub.TrialNum(3)/15))]);
        hold on; line([0.85:15.85],[0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33],'LineStyle','--','Color','k','Linewidth',1);   


    if saveplot
        filename = [path_fig_sub sub '_SeqIDHitRates_Overview.png'];
        saveas(gcf, filename, 'png'); %save png version
        delete(h);
    end    
end
% 7) Save single-subject hit rates (both tasks) for all tone duration conditions 
mkdir([paths_sfa_expt4.Analysis_Behavior sub '/'])
save([paths_sfa_expt4.Analysis_Behavior sub '/' 'BehavData.mat'], 'BehavData_SingleSub')
end