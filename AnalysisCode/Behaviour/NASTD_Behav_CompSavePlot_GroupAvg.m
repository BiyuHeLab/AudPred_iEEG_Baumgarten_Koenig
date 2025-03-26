%1) Script computes single-subject behavioral data 
%(final tone likelihood estimate,uSeq identification hit rates) 
%for all subjects ( only the finally selected ones)
%2) Script saves behavioral output vars (parameters for correlation
%estimated to presented beta) in 
%3) Script plots across-subject averaged behavioral parameters

%% 0) Setup analysis
% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')

NASTD_ECoG_setVars
paths_NASTD_ECoG = NASTD_ECoG_paths;

addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));

%Determine subjects
sub_list = vars.sub_list;
%Load in file with individual preproc infos
subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;

path_data = paths_NASTD_ECoG.RawData;
path_fig = paths_NASTD_ECoG.Analysis_Behavior;

i_sub = 1;
subs = vars.sub_list; %patients
sub = subs{i_sub};

tonedur_text = {'0.2' '0.4' 'all'};

saveplot = 0; %Save Figures?

BehavData_Group = []; %Create empty proxy file

%% 1) load single sub behavioral data, specify paths
subs = vars.sub_list(vars.validSubjs(2:end)); %don't use sub1 because it has different p34 allocation
% subs = vars.sub_list;

for i_sub = 1:length(subs)
    
    sub = subs{i_sub}
    
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    
    load([si.path_behavioraldata_sub]);%Load raw behavioral data
    % load([si.path_stimorder_sub]);%Load stimulus presentation order
    
    BehavData_Group{i_sub}.tonedur = tonedur_text; %specify current tone dur condition
    
    
    for i_tonedur = 1:length(tonedur_text)
        
        tonedur = tonedur_text{i_tonedur};
        
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
        %     filter_resp = data.uSeq.respABC >= 0 & data.resp_prob >= 0 & data.trialNum > 0;
        filter_resp = data.resp_prob >= 0 & data.trialNum > 0;
        %This takes care of the -1 responses in the final tone likelihood estimates
        %     filter_rt   = data.r1_RT > 0 & data.r2_RT > 0;
        filter_rt   = data.r1_RT > 0 ;
        
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
        
        %     %from data.uSeq (i.e., task2) subfields
        %     if i_sub < 4
        %         dfn = fieldnames(data.uSeq);
        %         for j = 1:length(dfn)
        %             if eval(['length(data.uSeq.' dfn{j} ') == 120'])
        %                 if strcmp(dfn{j},'disp3') %special because 3 columns instead of one
        %                 eval(['dataf.' dfn{j} ' = data.uSeq.' dfn{j} '(filter,:);']);
        %                 else
        %                 eval(['dataf.' dfn{j} ' = data.uSeq.' dfn{j} '(filter);']);
        %                 end
        %             end
        %         end
        %         dataf = rmfield(dataf,{'resp_beta','conf_beta', 'correct_beta', 'diff_beta', 'r3_RT'});
        %     end
        %define new struct 'stimf' containing all data subfields with all stimulus-wise parameters but only for above-filtered trials
        sfn = fieldnames(data.stim);
        for j = 1:length(sfn)
            if eval(['length(data.stim.' sfn{j} ') == 120'])
                eval(['stimf.' sfn{j} ' = data.stim.' sfn{j} '(filter);']);
            end
        end
        
        BehavData_Group{i_sub}.TrialNum(1,i_tonedur) = length(dataf.trialNum);
        
        %% 4) Final Tone Probability Rating - Task1
        %I.e., rate probability of presented final tone pitch (p34) based on previous sequence
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %4.1 Hit rate computation
        possible_ftp = unique(stimf.logf_final); %possible log freq of final tone (p34)
        predicted_ftp = unique(stimf.logf_pred); %possible log freq of predicted final tone (p*34)
        %3 different options per tone dur
        %6 different options across tone dur
        BehavData_Group{i_sub}.presentedFTP{i_tonedur} = possible_ftp;
        BehavData_Group{i_sub}.expectedFTP{i_tonedur} = predicted_ftp;
        
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
        
        BehavData_Group{i_sub}.FTPLrating_total{i_tonedur} = rp_tot;
        BehavData_Group{i_sub}.FTPLrating_total_SEM{i_tonedur} = rp_tot_sem;
        BehavData_Group{i_sub}.FTPLrating_exp{i_tonedur} = rp_exp;
        BehavData_Group{i_sub}.FTPLrating_exp_SEM{i_tonedur} = rp_exp_sem;
        
        %% 5) Sequence Identification - Task 2
        %i.e., select 1 of 3 presented images that visually represents the
        %previously presented tone sequence
        % if i_sub < 4
        %     addpath(paths_NASTD_ECoG.ScriptsDir) %necessary for dprime function
        %
        %     filter_resp_uSeq = data.uSeq.respABC >= 0;
        %     %analysis aims:
        %     % % - hit rate in total (120 trials)
        %     % Hitrate_SeqID_total = sum(data.correct_uSeq(filter_resp_uSeq))/length(data.correct_uSeq(filter_resp_uSeq));
        %     % [dval_SeqID_total, cval_SeqID_total] = dprime(Hitrate_SeqID_total,1-Hitrate_SeqID_perTonedur,1,2)
        %
        %     % - hit rate per tone dur (60 trials)
        %     Hitrate_SeqID_perTonedur = sum(dataf.correct_uSeq)/length(dataf.correct_uSeq);
        %     % [dval_SeqID_perTonedur, cval_SeqID_perTonedur] = dprime(Hitrate_SeqID_perTonedur,1-Hitrate_SeqID_perTonedur,1,2);
        %                 dval_SeqID_perTonedur = norminv(Hitrate_SeqID_perTonedur) - norminv(1-Hitrate_SeqID_perTonedur);
        %                 %norminv = gives the area under the curve of the normal distribution, equivalent of z-transformation
        %
        %     % - hit rate per individual sequence (8 trials/4 per tone dur)
        %         possible_uSeqID = unique(stimf.uSeqID); %possible individual sequences
        %         %15 different options with 8 runs in total/4 runs per tone dur
        %
        %         for i_uSeqID = 1:length(possible_uSeqID) %index for unique SedIDs
        %             f_uSeqID = stimf.uSeqID == possible_uSeqID(i_uSeqID); %filter to select those trials with specific SeqID
        %             Hitrate_SeqID_peruSeqID(1,i_uSeqID) = mean(dataf.correct_uSeq(f_uSeqID));%mean hit rate
        %     %         [dval_SeqID_peruSeqID(1,i_uSeqID), cval_SeqID_peruSeqID(1,i_uSeqID)] = dprime(Hitrate_SeqID_peruSeqID(1,i_uSeqID),1-Hitrate_SeqID_peruSeqID(1,i_uSeqID),1,2);
        %         end
        %
        %     BehavData_Group{i_sub}.SeqID_HitRate(1,i_tonedur) = Hitrate_SeqID_perTonedur;
        %     BehavData_Group{i_sub}.SeqID_dprime(1,i_tonedur) = dval_SeqID_perTonedur;
        %     % BehavData.SeqID_crit(1,i_tonedur) = cval_SeqID_perTonedur;
        %     BehavData_Group{i_sub}.SeqID_HitRateperuSeq{i_tonedur} = Hitrate_SeqID_peruSeqID;
        % end
        
        %% 6) Copy data from cell to matrix to enable later averaging
        proxy.TrialNum(i_sub, i_tonedur) = BehavData_Group{i_sub}.TrialNum(i_tonedur);
        proxy.r1_RT(i_sub,i_tonedur) = mean(dataf.r1_RT);
        % proxy.r2_RT(i_sub,i_tonedur) = mean(dataf.r2_RT);
        
        for i = 1:length(possible_ftp)
            proxy.FTPLrating_total(i_sub,i,i_tonedur) = BehavData_Group{i_sub}.FTPLrating_total{i_tonedur}(i);
            proxy.FTPLrating_total_SEM(i_sub,i,i_tonedur) = BehavData_Group{i_sub}.FTPLrating_total_SEM{i_tonedur}(i);
            for j = 1:length(predicted_ftp)
                proxy.FTPLrating_exp(i_sub,i,j,i_tonedur)  = BehavData_Group{i_sub}.FTPLrating_exp{i_tonedur}{j}(i);
                %organized by: line: subjects, colums = p34, 3D = p*34, 4D = tonedur
                proxy.FTPLrating_exp_SEM(i_sub,i,j,i_tonedur)  = BehavData_Group{i_sub}.FTPLrating_exp_SEM{i_tonedur}{j}(i);
            end
        end
        
        % if i_sub < 4
        %     proxy.SeqID_HitRate(i_sub,i_tonedur) = BehavData_Group{i_sub}.SeqID_HitRate(i_tonedur);
        %     proxy.SeqID_dprime(i_sub,i_tonedur) = BehavData_Group{i_sub}.SeqID_dprime(i_tonedur);
        %
        %     proxy.SeqID_HitRateperuSeq(i_sub,:,i_tonedur) = BehavData_Group{i_sub}.SeqID_HitRateperuSeq{i_tonedur};
        % end
        
    end
end

%% 7) Average across groups
BehavData_GroupAvg = [];
BehavData_GroupAvg.tonedur = BehavData_Group{2}.tonedur;
BehavData_GroupAvg.TrialNum = squeeze(mean(proxy.TrialNum));
BehavData_GroupAvg.r1_RT = (mean(proxy.r1_RT));
% BehavData_GroupAvg.r2_RT = (mean(proxy.r2_RT));

for i_tonedur = 1:length(tonedur)   
    BehavData_GroupAvg.FTPLrating_total(i_tonedur,:) = mean(proxy.FTPLrating_total(:,:,i_tonedur)); 
    %organized by line = tonedur, column = p34
    BehavData_GroupAvg.FTPLrating_total_SEM(i_tonedur,:) = mean(proxy.FTPLrating_total_SEM(:,:,i_tonedur));
    %organized by line = tonedur, column = p34       
        for i_predFTP = 1:length(predicted_ftp)
        BehavData_GroupAvg.FTPLrating_exp(i_tonedur,:,i_predFTP) = mean(proxy.FTPLrating_exp(:,:,i_predFTP,i_tonedur));
        %organized by line = tonedur, column = p34, 3D = p*34     
        BehavData_GroupAvg.FTPLrating_exp_SEM(i_tonedur,:,i_predFTP) = mean(proxy.FTPLrating_exp_SEM(:,:,i_predFTP,i_tonedur));    
        end
%     BehavData_GroupAvg.SeqID_HitRate(i_tonedur) = mean(proxy.SeqID_HitRate(:,i_tonedur));
%     BehavData_GroupAvg.SeqID_HitRate_SEM(i_tonedur) = std(proxy.SeqID_HitRate(:,i_tonedur)) / sqrt(length(proxy.SeqID_HitRate));
%     BehavData_GroupAvg.SeqID_HitRateperuSeq(i_tonedur,:) = mean(proxy.SeqID_HitRateperuSeq(:,:,i_tonedur));
%     BehavData_GroupAvg.SeqID_HitRateperuSeq_SEM(i_tonedur,:) = std(proxy.SeqID_HitRateperuSeq(:,:,i_tonedur)) / sqrt(size(proxy.SeqID_HitRateperuSeq(:,:,i_tonedur),1));
%     BehavData_GroupAvg.SeqID_dprime(i_tonedur) = mean(proxy.SeqID_dprime(:,i_tonedur));
end

%% 8) Group level statistics

%8.1 Group-Level ANOVA
%8.1.1 1way ANOVA & 8.1.2 2way ANOVA

BehavData_GroupAvg.FTPLrating_stats = NASTD_ECoG_Behav_ANOVAGroupAvg(subs, tonedur_text, vars);
%output: stats_anova2 table per tone duration{0.2, 0.4.all}
%Note: Without NY688 because different paradigm

% %8.2. Test if SeqID Hit Rates are sign. different from chance
% for i_tonedur = 1:length(tonedur)   
%     [h,p,ci, stats] = ttest(proxy.SeqID_HitRate(:,i_tonedur),0.33);
%     BehavData_GroupAvg.SeqID_stats.p_tTestvsChance(i_tonedur) = p;    
%     BehavData_GroupAvg.SeqID_stats.tstat_tTestvsChance(i_tonedur) = stats.tstat;    
% end
% clear h p ci stats
%% 9) Plotting group-average behavioral results
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
for i_tonedur = 1:length(tonedur_text)

%9.1 plot final tone probability rating independent of and dependent on math. expected tone
% xlab = {'539' '586' '655'};
a = {num2str(exp(possible_ftp(1))), num2str(exp(possible_ftp(2))), num2str(exp(possible_ftp(3))),num2str(exp(possible_ftp(4)))};
xlab = {[a{1} ' Hz/low'],[a{2} ' Hz/medium'],[a{3} ' Hz/high'], [a{4} ' Hz/high']};
pred_colors = {'rs', 'ks', 'gs'};

h = figd(15,2,5);
set(gcf,'units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1);
errorbar(possible_ftp, BehavData_GroupAvg.FTPLrating_total(i_tonedur,:), BehavData_GroupAvg.FTPLrating_total_SEM(i_tonedur,:), 'bo-','Linewidth',2);
xlabel('Presented Final Tone Pitch')
ylabel('Final Tone Pitch Likelihood Rating')
ylim([1 5])
% xlim([log(200) log(900)])
xlim([possible_ftp(1)-0.2 possible_ftp(4)+0.2])
set(gca,'XTick',possible_ftp)
set(gca,'XTickLabel', xlab)
title(['FTPL rating - GroupAvg with n = ' num2str(length(subs)) ' - ' num2str(BehavData_GroupAvg.TrialNum(i_tonedur)) ' of 120 trials included - ToneDur: ' tonedur_text{i_tonedur}]);

%include inset with ANOVA results:
inset_axis = axes('position', [0.15    0.625    0.5    0.05]);
axis off;
box on
results_anova1 = ['Main effect presented FTP (p34): F = ' num2str(BehavData_GroupAvg.FTPLrating_stats{i_tonedur}{3,5}) '; p = ' num2str(BehavData_GroupAvg.FTPLrating_stats{i_tonedur}{3,6})];
text(0.15,0.625,results_anova1)

set(0,'CurrentFigure',h);
subplot(2,1,2);
etitle = {'math.expected FTP = low', 'm.e. FTP = med', 'm.e. FTP  = high'};
colors = {'ro-' 'ko-' 'go-'};

for i_pred = 1:3
    hold on;
    errorbar(possible_ftp, BehavData_GroupAvg.FTPLrating_exp(i_tonedur,:,i_pred), BehavData_GroupAvg.FTPLrating_exp_SEM(i_tonedur,:,i_pred), colors{i_pred},'Linewidth',2);
    xlabel('Presented Final Tone Pitch')
    ylabel('Final Tone Pitch Likelihood Rating')
    ylim([1 5])
%     xlim([log(200) log(900)])
    xlim([possible_ftp(1)-0.2 possible_ftp(4)+0.2])
   set(gca,'XTick',possible_ftp)
    set(gca,'XTickLabel', xlab)
    title(['FTPL rating as function of math. expected FTP - GroupAvg with n = ' ...
        num2str(length(subs)) ' - ' num2str(BehavData_GroupAvg.TrialNum(i_tonedur)) ' of 120 trials included - ToneDur: ' tonedur_text{i_tonedur}]);

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
    MainEffect = ['Main effect math. expected FTP (p*34): F = ' num2str(BehavData_GroupAvg.FTPLrating_stats{i_tonedur}{2,5}) ...
        '; p = ' num2str(BehavData_GroupAvg.FTPLrating_stats{i_tonedur}{2,6})]
    
    InteractionEffect = ['Interaction effect p34*p*34: F = ' num2str(BehavData_GroupAvg.FTPLrating_stats{i_tonedur}{4,5}) ...
        '; p = ' num2str(BehavData_GroupAvg.FTPLrating_stats{i_tonedur}{4,6})];
    text(0.15,0.625,{MainEffect,InteractionEffect})

if saveplot
    filename = [path_fig 'GroupAvg_' tonedur_text{i_tonedur} '_FTPLrating_Overview.png'];
    saveas(gcf, filename, 'png'); %save png version
    delete(h);
end

end

%%9.2 Add behavioral results from MEG study on healthy subjects
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
load([paths_NASTD_ECoG.Analysis_Behavior 'GroupAvg/BehavFTPLData_Group_MEGstudy'])
common_TF = unique(sort([possible_ftp,tf]));
tonedur_text_subjects = {'0.15', '0.3', '0.6'};
%taken from 2way anova per TD (factors: p34, p*34)
FTPLstat_Subjects = struct;
FTPLstat_Subjects.mainp34.F = [9.01 11.0 13.33];
FTPLstat_Subjects.mainp34.p = [4.9134e-07 2.2848e-08 7.6279e-10];
FTPLstat_Subjects.mainpredp34.F = [3.93 1.49 1.99];
FTPLstat_Subjects.mainpredp34.p = [0.03 0.24 0.15];
FTPLstat_Subjects.interp34predp34.F = [11.77 20.58 25.25];
FTPLstat_Subjects.interp34predp34.p = [1.1102e-15 0 0];

for i_tonedur = 1:2%length(tonedur_text)

%9.2.1 plot final tone probability rating independent of math. expected tone
% xlab = {'539' '586' '655'};
a = {num2str(exp(possible_ftp(1))), num2str(exp(possible_ftp(2))), num2str(exp(possible_ftp(3))),num2str(exp(possible_ftp(4)))};
xlab = {[a{1}], '277', num2str(round(str2num(a{2}))),'349','554',num2str(round(str2num(a{3}))),'698', [a{4}]};

h = figd(15,2,5);
set(gcf,'units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1);
%Patients
errorbar(possible_ftp, BehavData_GroupAvg.FTPLrating_total(i_tonedur,:), BehavData_GroupAvg.FTPLrating_total_SEM(i_tonedur,:),...
    'ko-','Linewidth',2,'MarkerSize', 10, 'MarkerFaceColor','k');
%Healthy subjects
hold on;
e = errorbar(tf, rp_tot{i_tonedur}, rp_tot_sem{i_tonedur},...
    'ko--','Linewidth',2,'MarkerSize', 10, 'MarkerFaceColor','white');

xlabel('P34 [Hz]')
ylabel('FTPL Rating')
ylim([1 5])
% xlim([log(200) log(900)])
xlim([log(200) log(900)])
set(gca,'XTick',common_TF)
set(gca,'XTickLabel', xlab)
legend('Patients','Subjects')
title(['FTPL rating for Patients/Healthy Subjects - GroupAvg n = ' num2str(length(subs)) '/20 - ToneDur: ' tonedur_text{i_tonedur} '/' tonedur_text_subjects{i_tonedur} 'sec']);

%include inset with ANOVA results:
inset_axis = axes('position', [0.15    0.625    0.5    0.05]);
axis off;
box on
results_anova1 = {['Main effect p34 Patients: F = ' num2str(BehavData_GroupAvg.FTPLrating_stats{i_tonedur}{3,5}) '; p = ' num2str(BehavData_GroupAvg.FTPLrating_stats{i_tonedur}{3,6})],...
    ['Main effect p34 Subjects: F = ' num2str(FTPLstat_Subjects.mainp34.F(i_tonedur)) '; p = ' num2str(FTPLstat_Subjects.mainp34.p(i_tonedur))]};
text(0.15,0.625,results_anova1) 

%9.2.2 Plot  final tone probability rating dependent on math. expected tone
set(0,'CurrentFigure',h);
subplot(2,1,2);
% etitle = {'math.expected FTP = low', 'm.e. FTP = med', 'm.e. FTP  = high'};
pred_colors = {'bs', 'cs', 'ys'};

%Patients
colors = {'bo-' 'co-' 'yo-'};
for i_pred = 1:3
    hold on;
    errorbar(possible_ftp, BehavData_GroupAvg.FTPLrating_exp(i_tonedur,:,i_pred), BehavData_GroupAvg.FTPLrating_exp_SEM(i_tonedur,:,i_pred), colors{i_pred},...
        'Linewidth',2, 'MarkerSize', 10, 'MarkerFaceColor','k');
    xlabel('P34 [Hz]')
    ylabel('FTPL Rating')
    ylim([1 5])
    xlim([log(200) log(900)])
    set(gca,'XTick',common_TF)
    set(gca,'XTickLabel', xlab)
end
    hold on
    for i_pred = 1:3
        plot(predicted_ftp(i_pred), 2, pred_colors{i_pred},'LineWidth', 3,'MarkerSize', 10, 'MarkerFaceColor','k')
        hold on
    end
    
%Subjects
colors = {'bo--' 'co--' 'yo--'};
for i_pred = 1:3
    hold on;
    e = errorbar(tf, rp_exp{i_tonedur}(i_pred,:), rp_exp_sem{i_tonedur}(i_pred,:), colors{i_pred},...
        'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor','white');
    e.LineWidth = 3;
    xlabel('P34 [Hz]')
    ylabel('FTPL Rating')
    ylim([1 5])
    xlim([log(200) log(900)])
    set(gca,'XTick',common_TF)
    set(gca,'XTickLabel', xlab)
    title(['GAvg - FTPLrating per p*34 - ToneDur: ' tonedur_text{i_tonedur} ' s']);
end

hold on
for i_pred = 1:3
    plot(meanExpPitch(i_pred), 2, pred_colors{i_pred},'LineWidth', 3,'MarkerSize', 10, 'MarkerFaceColor','white')
    hold on
end

    title(['FTPL rating as function of p*34 for Patients/Healthy Subjects - GroupAvg n = ' num2str(length(subs)) '/20 - ToneDur: ' tonedur_text{i_tonedur} '/' tonedur_text_subjects{i_tonedur} ' sec']);
    etitle = {'p*_{34} = low', 'p*_{34} = med', 'p*_{34} = high'};
%     etitle = {'p*_{34} = low', 'p*_{34} = med', 'p*_{34} = high','E[pitch] = low', 'E[pitch] = med', 'E[pitch] = high'};
    hl = legend(etitle, 'Location', 'Best'); 
    hl.FontSize = 10;
    
    %include inset with ANOVA results:
    inset_axis = axes('position', [0.15    0.125    0.5    0.05]);
    axis off;
    box on
    MainEffect1 = ['Main effect p*34 Patients: F = ' num2str(BehavData_GroupAvg.FTPLrating_stats{i_tonedur}{2,5}) ...
        '; p = ' num2str(BehavData_GroupAvg.FTPLrating_stats{i_tonedur}{2,6})];    
    InteractionEffect1 = ['Interaction p34*p*34 Patients: F = ' num2str(BehavData_GroupAvg.FTPLrating_stats{i_tonedur}{4,5}) ...
        '; p = ' num2str(BehavData_GroupAvg.FTPLrating_stats{i_tonedur}{4,6})];
    MainEffect2 = ['Main effect p*34 Subjects: F = ' num2str(FTPLstat_Subjects.mainpredp34.F(i_tonedur)) ...
        '; p = ' num2str(FTPLstat_Subjects.mainpredp34.p(i_tonedur))];    
    InteractionEffect2 = ['Interaction p34*p*34 Subjects: F = ' num2str(FTPLstat_Subjects.interp34predp34.F(i_tonedur)) ...
        '; p = ' num2str(FTPLstat_Subjects.interp34predp34.p(i_tonedur))];    
    text(0.15,0.625,{MainEffect1,MainEffect2,InteractionEffect1,InteractionEffect2})

if saveplot
    filename = [path_fig 'PatientsSubjects_GroupAvg_' tonedur_text{i_tonedur} '_FTPLrating_Overview.png'];
    saveas(gcf, filename, 'png'); %save png version
    delete(h);
end

end


% %8.2 plot Sequence Identification (task2) for all tone dur
% xlab = {'0.2 ms' '0.4 ms' 'all'};
% 
% h = figd(15,2,5);
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
% subplot(4,1,1);
% errorbar(BehavData_GroupAvg.SeqID_HitRate,BehavData_GroupAvg.SeqID_HitRate_SEM, '.b','MarkerSize',40,'Linewidth',2);   
%     xlabel('Tone Duration')
%     ylabel('HitRate')
%     ylim([-0.1 1.1])
%     xlim([0.85 3.15])
%     set(gca,'XTick',1:3)
%     set(gca,'XTickLabel', xlab)
%     title(['Hit Rate (ratio correct responses/all responses) for Sequence Identification  - GroupAvg with n = ' num2str(length(subs))]);
%     hold on; line([0.85:3.85],[0.33,0.33,0.33,0.33],'LineStyle','--','Color','k','Linewidth',1);
% 
%     inset_axis = axes('position', [0.115    0.86    0.5    0.05]);
%     axis off;
%     box on
%     dprime1 = ['d-prime: ' num2str(BehavData_GroupAvg.SeqID_dprime(1))];
%     tstat1 = ['t = ' num2str(BehavData_GroupAvg.SeqID_stats.tstat_tTestvsChance(1)) '; p = ' num2str(BehavData_GroupAvg.SeqID_stats.p_tTestvsChance(1))];
%     text(0.068,0.7,{dprime1 tstat1})
%     dprime2 = ['d-prime: ' num2str(BehavData_GroupAvg.SeqID_dprime(2))];
%     tstat2 = ['t = ' num2str(BehavData_GroupAvg.SeqID_stats.tstat_tTestvsChance(2)) '; p = ' num2str(BehavData_GroupAvg.SeqID_stats.p_tTestvsChance(2))];
%     text(0.745,0.7,{dprime2 tstat2})    
%     dprime3 = ['d-prime: ' num2str(BehavData_GroupAvg.SeqID_dprime(3))];
%     tstat3 = ['t = ' num2str(BehavData_GroupAvg.SeqID_stats.tstat_tTestvsChance(3)) '; p = ' num2str(BehavData_GroupAvg.SeqID_stats.p_tTestvsChance(3))];
%    text(1.41, 0.7,{dprime3 tstat3})     
% 
% xlab = {'Seq1','Seq2','Seq3','Seq4','Seq5','Seq6','Seq7','Seq8','Seq9','Seq10','Seq11','Seq12','Seq13','Seq14','Seq15'};   
% subplot(4,1,2);
% errorbar(BehavData_GroupAvg.SeqID_HitRateperuSeq(1,:),BehavData_GroupAvg.SeqID_HitRateperuSeq_SEM(1,:), '.b','MarkerSize',20,'Linewidth',2);   
% %     xlabel('Individual Sequence - 0.2s ToneDur')
%     ylabel('HitRate')
%     ylim([-0.1 1.1])
%     xlim([0.85 15.15])
%     set(gca,'XTick',1:15)
%     set(gca,'XTickLabel', xlab)
%     title(['Hit Rate for unique sequence SeqID - 0.2s ToneDur']);
%     hold on; line([0.85:15.85],[0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33],'LineStyle','--','Color','k','Linewidth',1);   
% 
% xlab = {'','','','','','','','','','','','','','',''};      
% subplot(4,1,3);
% errorbar(BehavData_GroupAvg.SeqID_HitRateperuSeq(2,:),BehavData_GroupAvg.SeqID_HitRateperuSeq_SEM(2,:), '.b','MarkerSize',20,'Linewidth',2);   
% %     xlabel('Individual Sequence - 0.4s ToneDur')
%     ylabel('HitRate')
%     ylim([-0.1 1.1])
%     xlim([0.85 15.15])
% %     set(gca,'XTick',1:15)
%      set(gca,'XTickLabel', xlab)
%     title(['Hit Rate for unique sequence SeqID - 0.4s ToneDur']);
%     hold on; line([0.85:15.85],[0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33],'LineStyle','--','Color','k','Linewidth',1);   
% 
% subplot(4,1,4);
% errorbar(BehavData_GroupAvg.SeqID_HitRateperuSeq(3,:),BehavData_GroupAvg.SeqID_HitRateperuSeq_SEM(3,:),'.b','MarkerSize',20,'Linewidth',2);   
% %     xlabel('Individual Sequence - All ToneDur')
%     ylabel('HitRate')
%     ylim([-0.1 1.1])
%     xlim([0.85 15.15])
% %     set(gca,'XTick',1:15)
%      set(gca,'XTickLabel', xlab)
%     title(['Hit Rate for unique sequence SeqID - All ToneDur']);
%     hold on; line([0.85:15.85],[0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33],'LineStyle','--','Color','k','Linewidth',1);   
%  
%     
% if saveplot
%     filename = [path_fig 'GroupAvg_SeqIDHitRates_Overview.png'];
%     saveas(gcf, filename, 'png'); %save png version
%     delete(h);
% end 

   
% 7) Save single-subject hit rates (both tasks) for all tone duration conditions 
mkdir([paths_NASTD_ECoG.Analysis_Behavior 'GroupAvg/'])
save([paths_NASTD_ECoG.Analysis_Behavior 'GroupAvg/' 'BehavData_Group.mat'], 'BehavData_GroupAvg')
