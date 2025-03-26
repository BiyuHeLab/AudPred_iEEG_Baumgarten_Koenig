%% Project: NASTD_ECoG
%Compute behavioral prediction effects

%1. Script computes single-subject behavioral data
%(final tone pitch likelihood (FTPL)nrating
%2. Average across (selected) subjects and plot Grup-Level results

%% 0) Setup analysis
% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')

NASTD_ECoG_setVars
paths_NASTD_ECoG = NASTD_ECoG_paths;

% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));

%Determine subjects
sub_list = vars.sub_list(vars.validSubjs(2:end));
%Load in file with individual preproc infos
subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;

path_outputdata = ...
    [paths_NASTD_ECoG.Analysis_Behavior 'FTPLratings/Allsub_n' num2str(length(sub_list)) '/Data/'];
if (~exist(path_outputdata, 'dir')); mkdir(path_outputdata); end
path_fig = ...
    [paths_NASTD_ECoG.Analysis_Behavior 'FTPLratings/Allsub_n' num2str(length(sub_list)) '/Figs/'];
if (~exist(path_fig, 'dir')); mkdir(path_fig); end

tonedur_text = {'0.2' '0.4' 'all'};

saveplot = 0; %Save Figures?

BehavData.Sub = []; %Create empty proxy file

%% 1) load single sub behavioral data, specify paths
for i_sub = 1:length(sub_list)
    %Don't use NY688 because of different p*34s
    
    sub = sub_list{i_sub}
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    load([si.path_behavioraldata_sub]);%Load raw behavioral data
    
    BehavData.Sub{i_sub}.tonedur = tonedur_text; %specify current tone dur condition
    
    %% 2) check counterbalancing (i.e., if all stim combinations present per block)
    for i_tonedur = 1:length(tonedur_text)
        CurrTonedur = tonedur_text{i_tonedur};
        
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
        
        %% 3) Differentiate and read out different trial types
        %3.1) define filters for trial selection
        %General filter - ensure all behavioral responses were entered properly
        filter_resp = data.resp_prob >= 0 & data.trialNum > 0;
        filter_rt   = data.r1_RT > 0 ;
        
        %Tone duration filter
        switch CurrTonedur
            case '0.2'
                filter_toneDur = data.stim.toneDur == 0.2;
            case '0.4'
                filter_toneDur = data.stim.toneDur == 0.4;
            otherwise
                filter_toneDur = ones(1,120);
        end
        
        %Combined filter
        filter = filter_resp & filter_rt & filter_toneDur; %filter var for trial selection
        
        %3.2) Select trials based on filters
        dfn = fieldnames(data);
        for j = 1:length(dfn)
            if eval(['length(data.' dfn{j} ') == 120'])
                eval(['data_filt.' dfn{j} ' = data.' dfn{j} '(filter);']);
            end
        end
        
        sfn = fieldnames(data.stim);
        for j = 1:length(sfn)
            if eval(['length(data.stim.' sfn{j} ') == 120'])
                eval(['stim_filt.' sfn{j} ' = data.stim.' sfn{j} '(filter);']);
            end
        end
        
        BehavData.Sub{i_sub}.TrialNum(1,i_tonedur) = length(data_filt.trialNum);
        
        %% 4) Final Tone Probability Rating
        %Compute FTPL rate for each p*34 and p34 instantiation
        p34 = unique(stim_filt.logf_final);
        pred_p34 = unique(stim_filt.logf_pred);
        BehavData.Sub{i_sub}.p34{i_tonedur} = p34;
        BehavData.Sub{i_sub}.pred_p34{i_tonedur} = pred_p34;
        
        predp34_IDs = [-1 0 1];
        
        for i_predp34 = 1:3 %index for predicted final tone
            filter_predp34 = stim_filt.predID == predp34_IDs(i_predp34);
            
            for i_p34 = 1:length(p34)%index for presented final tone
                filter_p34  = stim_filt.logf_final == p34(i_p34);
                
                avg_FTPLrating_perp34(i_p34) = ...
                    mean(data_filt.resp_prob(filter_p34));
                SEM_FTPLrating_perp34(i_p34) = ...%SEM across trials
                    std(data_filt.resp_prob(filter_p34)) / sqrt(sum(filter_p34));
                
                avg_FTPLrating_perpredp34{i_predp34}(i_p34) = ...
                    mean(data_filt.resp_prob(filter_predp34 & filter_p34));
                SEM_FTPLrating_perpredp34{i_predp34}(i_p34) = ...%SEM across trials
                    std(data_filt.resp_prob(filter_predp34 & filter_p34)) / sqrt(sum(filter_predp34 & filter_p34));
                
            end
        end
        
        BehavData.Sub{i_sub}.avgFTPLrating_perp34{i_tonedur} = avg_FTPLrating_perp34;
        BehavData.Sub{i_sub}.semFTPLrating_perp34{i_tonedur} = SEM_FTPLrating_perp34;%SEM across trials
        BehavData.Sub{i_sub}.avgFTPLrating_perp34predp34{i_tonedur} = avg_FTPLrating_perpredp34;
        BehavData.Sub{i_sub}.semFTPLrating_perp34predp34{i_tonedur} = SEM_FTPLrating_perpredp34;%SEM across trials
        
        %% 5) Copy data from cell to matrix to enable later averaging
        proxy.TrialNum(i_sub, i_tonedur) = BehavData.Sub{i_sub}.TrialNum(i_tonedur);
        proxy.r1_RT(i_sub,i_tonedur) = mean(data_filt.r1_RT);
        
        for i = 1:length(p34)
            proxy.avgFTPLrating_perp34(i_sub,i,i_tonedur) = ...
                BehavData.Sub{i_sub}.avgFTPLrating_perp34{i_tonedur}(i);
            proxy.semFTPLrating_perp34(i_sub,i,i_tonedur) = ...
                BehavData.Sub{i_sub}.semFTPLrating_perp34{i_tonedur}(i);
            for j = 1:length(pred_p34)
                
                %organized by: line: subjects, colums = p34, 3D = p*34, 4D = tonedur
                proxy.avgFTPLrating_perp34predp34(i_sub,i,j,i_tonedur) = ...
                    BehavData.Sub{i_sub}.avgFTPLrating_perp34predp34{i_tonedur}{j}(i);
                proxy.semFTPLrating_perp34predp34(i_sub,i,j,i_tonedur) = ...
                    BehavData.Sub{i_sub}.semFTPLrating_perp34predp34{i_tonedur}{j}(i);
            end
        end
        
    end
end

%% 6) Average FTPL ratings across subjects
BehavData.Group = [];
BehavData.Group.tonedur = tonedur_text;
BehavData.Group.avgTrialNum = squeeze(mean(proxy.TrialNum));
BehavData.Group.avgRT = (mean(proxy.r1_RT));

for i_tonedur = 1:length(CurrTonedur)
    BehavData.Group.avgFTPLrating_perp34(i_tonedur,:) = ...
        mean(proxy.avgFTPLrating_perp34(:,:,i_tonedur));
    BehavData.Group.std(i_tonedur,:) = ...%STD across subjects
        std(proxy.avgFTPLrating_perp34(:,:,i_tonedur));
    BehavData.Group.semFTPLrating_perp34(i_tonedur,:) = ...%SEM across subjects
        std(proxy.avgFTPLrating_perp34(:,:,i_tonedur))/sqrt(length(sub_list));
    
    for i_predFTP = 1:length(pred_p34)
        BehavData.Group.avgFTPLrating_perp34predp34(i_tonedur,:,i_predFTP) = ...
            mean(proxy.avgFTPLrating_perp34predp34(:,:,i_predFTP,i_tonedur));
        BehavData.Group.stdFTPLrating_perp34predp34(i_tonedur,:,i_predFTP) = ...%STD across subjects
            std(proxy.avgFTPLrating_perp34predp34(:,:,i_predFTP,i_tonedur));
        BehavData.Group.semFTPLrating_perp34predp34(i_tonedur,:,i_predFTP) = ...%SEM across subjects
            std(proxy.avgFTPLrating_perp34predp34(:,:,i_predFTP,i_tonedur))/sqrt(length(sub_list));
    end
end

%% 7) Group level statistics
%7.1 rm ANOVA
BehavData.Group.FTPLrating_stats = ...
    NASTD_ECoG_Behav_GroupStat(sub_list, tonedur_text);

%% 8) Plotting group-average behavioral results
for i_tonedur = 1:length(tonedur_text)
    
    p34_Hz = {num2str(round(exp(p34(1)))), num2str(round(exp(p34(2)))), ...
        num2str(round(exp(p34(3)))),num2str(round(exp(p34(4))))};
    xlabel1 = {[p34_Hz{1} ' Hz'],[p34_Hz{2} ' Hz'],[p34_Hz{3} ' Hz'], [p34_Hz{4} ' Hz']};
    
    h = figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    sgtitle(['FTPL rating (mean/SEM) for Group Level (n = ' num2str(length(BehavData.Sub)) ...
        '); ' num2str(BehavData.Group.avgTrialNum(i_tonedur)) ...
        '/120 trials included; ToneDur: ' tonedur_text{i_tonedur}]);
    
    %8.1 plot final tone probability rating as function of p34
    subplot(2,1,1);
    errorbar(p34, BehavData.Group.avgFTPLrating_perp34(i_tonedur,:), ...
        BehavData.Group.semFTPLrating_perp34(i_tonedur,:), 'ko-','Linewidth',2);
    xlabel('Presented Final Tone Pitch (p34)')
    xlim([p34(1)-0.2 p34(4)+0.2])
    ylabel('Final Tone Pitch Likelihood Rating')
    ylim([1 5])
    yticks([1 2 3 4 5])
    yticklabels({'1 (unlikely)', '2', '3 (unsure)', '4', '5 (likely)'})
    set(gca,'XTick',p34)
    set(gca,'XTickLabel', xlabel1)
    title(['FTPL rating (mean/SEM) per p34']);
    
    %add single subject results
    hold on;
    for i_sub = 1:length(sub_list)
        s = scatter(p34 + 0.025, BehavData.Sub{i_sub}.avgFTPLrating_perp34{i_tonedur},'o','filled');
        s.MarkerFaceColor = [0.7, 0.7, 0.7];
%         plot(p34 + 0.025, BehavData.Sub{i_sub}.avgFTPLrating_perp34{i_tonedur}, 'Color', [0.7, 0.7, 0.7]);
    end
    
%     %add inset with ANOVA results:
%     inset_axis = axes('position', [0.15    0.625    0.5    0.05]);
%     axis off;
%     box on
%     results_anova1 = ['Main effect presented FTP (p34): F = ' ...
%         num2str(BehavData.Group.FTPLrating_stats{i_tonedur}{3,5}) ...
%         '; p = ' num2str(BehavData.Group.FTPLrating_stats{i_tonedur}{3,6})];
%     text(0.15,0.625,results_anova1)
    
    %8.2 plot final tone probability rating as function of p34 and p*34
    subplot(2,1,2);
    legend_entry = {'p*34 = low', 'p*34 = med', 'p*34  = high'};
%     predp34_colors = {'b-' 'c-' 'y-'};
    predp34_colors = {[0,0.45,0.7],[0, 0.6,0.5],[0.8,0.4,0]};
    
    for i_predp34 = 1:3
        hold on;
        errorbar(p34 - (0.02*(i_predp34 -2)), ...
            BehavData.Group.avgFTPLrating_perp34predp34(i_tonedur,:,i_predp34), ...
            BehavData.Group.semFTPLrating_perp34predp34(i_tonedur,:,i_predp34), ...
            'Color', predp34_colors{i_predp34},'Linewidth',2);
    end
    xlabel('Presented Final Tone Pitch (p34)')
    xlim([p34(1)-0.2 p34(4)+0.2])
    ylabel('Final Tone Pitch Likelihood Rating')
    ylim([1 5])
    yticks([1 2 3 4 5])
    yticklabels({'1 (unlikely)', '2', '3 (unsure)', '4', '5 (likely)'})
    set(gca,'XTick',p34)
    set(gca,'XTickLabel', xlabel1)
    title(['FTPL rating (mean/SEM) per p34 as a function of p*34']);
    
    %Add marks for p*34 freqs
    hold on
%     predp34_colors = {'bs', 'cs', 'ys'};
    for i_predp34 = 1:3
        plot(pred_p34(i_predp34), 2, '-s', ...
            'MarkerEdgeColor', predp34_colors{i_predp34}, ...
            'MarkerFaceColor', predp34_colors{i_predp34}, ...
            'LineWidth', 3, ...
            'MarkerSize', 10)
        hold on
    end
    
    %Add single subject results
    hold on;
%     predp34_colors = {'b', 'c', 'y'};

    for i_predp34 = 1:3
        for i_sub = 1:length(sub_list)
            s = scatter(p34 - (0.02*(i_predp34) - 0.05), ...
                BehavData.Sub{i_sub}.avgFTPLrating_perp34predp34{i_tonedur}{i_predp34}, ...
                'o','filled');
            s.MarkerFaceColor = predp34_colors{i_predp34};
%             plot(p34 - ((0.035*i_predp34) - 0.06), ...
%                 BehavData.Sub{i_sub}.avgFTPLrating_perp34predp34{i_tonedur}{i_predp34}, ...
%                 'Color', [0.7, 0.7, 0.7]);
        end
    end
    
    hl = legend(legend_entry);
    set(hl, 'box', 'off','Location','NorthEast')
    
    
%     %include inset with ANOVA results:
%     inset_axis = axes('position', [0.15    0.125    0.5    0.05]);
%     axis off;
%     box on
%     MainEffect = ['Main effect math. expected FTP (p*34): F = ' ...
%         num2str(BehavData.Group.FTPLrating_stats{i_tonedur}{2,5}) ...
%         '; p = ' num2str(BehavData.Group.FTPLrating_stats{i_tonedur}{2,6})];
%     
%     InteractionEffect = ['Interaction effect p34*p*34: F = ' ...
%         num2str(BehavData.Group.FTPLrating_stats{i_tonedur}{4,5}) ...
%         '; p = ' num2str(BehavData.Group.FTPLrating_stats{i_tonedur}{4,6})];
%     text(0.15,0.625,{MainEffect,InteractionEffect})
    
    if saveplot
        filename = [path_fig 'FTPLrating_GroupLeveln' ...
            num2str(length(sub_list)) '_' tonedur_text{i_tonedur} 's.png'];
        saveas(gcf, filename, 'png'); %save png version
        delete(h);
    end
    
end

%8.3 Add behavioral results from MEG study on healthy subjects
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
load([paths_NASTD_ECoG.Analysis_Behavior 'GroupAvg/BehavFTPLData_Group_MEGstudy'])
common_TF = unique(sort([p34,tf]));
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
    p34_Hz = {num2str(exp(p34(1))), num2str(exp(p34(2))), num2str(exp(p34(3))),num2str(exp(p34(4)))};
    xlabel1 = {[p34_Hz{1}], '277', num2str(round(str2num(p34_Hz{2}))),'349','554',num2str(round(str2num(p34_Hz{3}))),'698', [p34_Hz{4}]};
    
    h = figd(15,2,5);
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    subplot(2,1,1);
    %Patients
    errorbar(p34, BehavData.Group.avgFTPLrating_perp34(i_tonedur,:), BehavData.Group.semFTPLrating_perp34(i_tonedur,:),...
        'ko-','Linewidth',2,'MarkerSize', 10, 'MarkerFaceColor','k');
    %Healthy subjects
    hold on;
    e = errorbar(tf, avg_FTPLrating_perp34{i_tonedur}, SEM_FTPLrating_perp34{i_tonedur},...
        'ko--','Linewidth',2,'MarkerSize', 10, 'MarkerFaceColor','white');
    
    xlabel('P34 [Hz]')
    ylabel('FTPL Rating')
    ylim([1 5])
    % xlim([log(200) log(900)])
    xlim([log(200) log(900)])
    set(gca,'XTick',common_TF)
    set(gca,'XTickLabel', xlabel1)
    legend('Patients','Subjects')
    title(['FTPL rating for Patients/Healthy Subjects - GroupAvg n = ' num2str(length(subs)) '/20 - ToneDur: ' tonedur_text{i_tonedur} '/' tonedur_text_subjects{i_tonedur} 'sec']);
    
    %include inset with ANOVA results:
    inset_axis = axes('position', [0.15    0.625    0.5    0.05]);
    axis off;
    box on
    results_anova1 = {['Main effect p34 Patients: F = ' num2str(BehavData.Group.FTPLrating_stats{i_tonedur}{3,5}) '; p = ' num2str(BehavData.Group.FTPLrating_stats{i_tonedur}{3,6})],...
        ['Main effect p34 Subjects: F = ' num2str(FTPLstat_Subjects.mainp34.F(i_tonedur)) '; p = ' num2str(FTPLstat_Subjects.mainp34.p(i_tonedur))]};
    text(0.15,0.625,results_anova1)
    
    %9.2.2 Plot  final tone probability rating dependent on math. expected tone
    set(0,'CurrentFigure',h);
    subplot(2,1,2);
    % etitle = {'math.expected FTP = low', 'm.e. FTP = med', 'm.e. FTP  = high'};
    predp34_colors = {'bs', 'cs', 'ys'};
    
    %Patients
    colors = {'bo-' 'co-' 'yo-'};
    for i_predp34 = 1:3
        hold on;
        errorbar(p34, BehavData.Group.avgFTPLrating_perp34predp34(i_tonedur,:,i_predp34), BehavData.Group.semFTPLrating_perp34predp34(i_tonedur,:,i_predp34), colors{i_predp34},...
            'Linewidth',2, 'MarkerSize', 10, 'MarkerFaceColor','k');
        xlabel('P34 [Hz]')
        ylabel('FTPL Rating')
        ylim([1 5])
        xlim([log(200) log(900)])
        set(gca,'XTick',common_TF)
        set(gca,'XTickLabel', xlabel1)
    end
    hold on
    for i_predp34 = 1:3
        plot(pred_p34(i_predp34), 2, predp34_colors{i_predp34},'LineWidth', 3,'MarkerSize', 10, 'MarkerFaceColor','k')
        hold on
    end
    
    %Subjects
    colors = {'bo--' 'co--' 'yo--'};
    for i_predp34 = 1:3
        hold on;
        e = errorbar(tf, avg_FTPLrating_perpredp34{i_tonedur}(i_predp34,:), SEM_FTPLrating_perpredp34{i_tonedur}(i_predp34,:), colors{i_predp34},...
            'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor','white');
        e.LineWidth = 3;
        xlabel('P34 [Hz]')
        ylabel('FTPL Rating')
        ylim([1 5])
        xlim([log(200) log(900)])
        set(gca,'XTick',common_TF)
        set(gca,'XTickLabel', xlabel1)
        title(['GAvg - FTPLrating per p*34 - ToneDur: ' tonedur_text{i_tonedur} ' s']);
    end
    
    hold on
    for i_predp34 = 1:3
        plot(meanExpPitch(i_predp34), 2, predp34_colors{i_predp34},'LineWidth', 3,'MarkerSize', 10, 'MarkerFaceColor','white')
        hold on
    end
    
    title(['FTPL rating as function of p*34 for Patients/Healthy Subjects - GroupAvg n = ' num2str(length(subs)) '/20 - ToneDur: ' tonedur_text{i_tonedur} '/' tonedur_text_subjects{i_tonedur} ' sec']);
    legend_entry = {'p*_{34} = low', 'p*_{34} = med', 'p*_{34} = high'};
    %     etitle = {'p*_{34} = low', 'p*_{34} = med', 'p*_{34} = high','E[pitch] = low', 'E[pitch] = med', 'E[pitch] = high'};
    hl = legend(legend_entry, 'Location', 'Best');
    hl.FontSize = 10;
    
    %include inset with ANOVA results:
    inset_axis = axes('position', [0.15    0.125    0.5    0.05]);
    axis off;
    box on
    MainEffect1 = ['Main effect p*34 Patients: F = ' num2str(BehavData.Group.FTPLrating_stats{i_tonedur}{2,5}) ...
        '; p = ' num2str(BehavData.Group.FTPLrating_stats{i_tonedur}{2,6})];
    InteractionEffect1 = ['Interaction p34*p*34 Patients: F = ' num2str(BehavData.Group.FTPLrating_stats{i_tonedur}{4,5}) ...
        '; p = ' num2str(BehavData.Group.FTPLrating_stats{i_tonedur}{4,6})];
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

%% 9) Save single-subject hit rates (both tasks) for all tone duration conditions
save([path_data 'FTPLratings_Groupn' num2str(length(sub_list)) '.mat'], ...
    'BehavData')
