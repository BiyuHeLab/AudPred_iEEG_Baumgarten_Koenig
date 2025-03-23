%NAST_Behav_ECoG_CompareMEG2FTPLeffects_Histogram
%Aim: To compare behavioral FTPL rating between helahty subjects (MEG study) 
%and patients (ECoG study).
%Plot pval (log10) for single-subject main effect p34 and interaction effect 
%p34-p*34 from healthy subs as histogram, add patient pvalues

%% 1: Load in behavioral statistics
%1.1 ECoG patients
location = 'server'; %gago/gogo
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
NASTD_ECoG_setVars
paths_NASTD = NASTD_ECoG_paths(location);

%Add base dir and own script dir
addpath(genpath(paths_NASTD.BaseDir));
addpath(genpath(paths_NASTD.ScriptsDir));
path_data = paths_NASTD.RawData;

i_sub = 1;
subs = vars.sub_list; %patients
tonedur_text = {'0.2' '0.4'};

for i_sub = 1:length(subs)
    for i_tonedur = 1:2%length(tonedur_text)
        sub = subs{i_sub};
        load([paths_NASTD.Analysis_Behavior sub '/' 'BehavData.mat']) %BehavData_SingleSub

            %From ANOVA1_p34
%             FTPL_Stat.ECoG.Main_p34{i_tonedur}(i_sub,:) = cell2mat(BehavData_SingleSub.FTPLrating_anova1{i_tonedur}(2,6:7)); %Factors: p34
            %From ANOVA2_p34predp34
            FTPL_Stat.ECoG.Main_p34_ANOVA2{i_tonedur}(i_sub,:) = cell2mat(BehavData_SingleSub.FTPLrating_anova2{i_tonedur}(3,6:7)); %Factors: p34, predp34
            FTPL_Stat.ECoG.Main_predp34_ANOVA2{i_tonedur}(i_sub,:) = cell2mat(BehavData_SingleSub.FTPLrating_anova2{i_tonedur}(2,6:7)); %Factors: p34, predp34
            FTPL_Stat.ECoG.Interaction_p34predp34_ANOVA2{i_tonedur}(i_sub,:) = cell2mat(BehavData_SingleSub.FTPLrating_anova2{i_tonedur}(4,6:7)); %Factors: p34, predp34
    end
    
    %From ANOVA3_TDp34predp34
    FTPL_Stat.ECoG.Main_TD_ANOVA3{3}(i_sub,:) = cell2mat(BehavData_SingleSub.FTPLrating_anova3(2,6:7));
    FTPL_Stat.ECoG.Main_p34_ANOVA3{3}(i_sub,:) = cell2mat(BehavData_SingleSub.FTPLrating_anova3(4,6:7));
    FTPL_Stat.ECoG.Main_predp34_ANOVA3{3}(i_sub,:) = cell2mat(BehavData_SingleSub.FTPLrating_anova3(3,6:7));    
    FTPL_Stat.ECoG.Interaction_p34predp34_ANOVA3{3}(i_sub,:) = cell2mat(BehavData_SingleSub.FTPLrating_anova3(7,6:7));
end

%1.2 healthy subjects
location = 'local'; %'cluster'
server = 'gogo'; 
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/MEG/') %main project path
paths_NASTD_MEG = NASTD_MEG_paths(location, server);
addpath(genpath(paths_NASTD_MEG.ScriptsDir));

path_analysis   = paths_NASTD_MEG.Analysis.Behavior;
path_outputdata = ([path_analysis 'FTPLrating/']); %Output: Fvals

load([path_outputdata 'Fpstat_FTPLrating_AllEffects.mat']); %FTPLstat

for i_sub = 1:length(FTPLstat.SingleSub.ANOVA2.FP_Main_p34)
    for i_tonedur = 1:size(FTPLstat.SingleSub.ANOVA2.FP_Main_p34,2)
        
        %From ANOVA2_p34predp34
        FTPL_Stat.MEG.Main_p34_ANOVA2{i_tonedur}(i_sub,:) = cell2mat(FTPLstat.SingleSub.ANOVA2.FP_Main_p34{i_sub,i_tonedur});
        FTPL_Stat.MEG.Main_predp34_ANOVA2{i_tonedur}(i_sub,:) = cell2mat(FTPLstat.SingleSub.ANOVA2.FP_Main_predp34{i_sub,i_tonedur});
        FTPL_Stat.MEG.Interaction_p34predp34_ANOVA2{i_tonedur}(i_sub,:) = cell2mat(FTPLstat.SingleSub.ANOVA2.FP_Interaction_p34predp34{i_sub,i_tonedur});
    end
    %From ANOVA3_TDp34predp34
    FTPL_Stat.MEG.Main_TD_ANOVA3{3}(i_sub,:) = cell2mat(FTPLstat.SingleSub.ANOVA3.FP_Main_TD(i_sub,:));
    FTPL_Stat.MEG.Main_p34_ANOVA3{3}(i_sub,:) = cell2mat(FTPLstat.SingleSub.ANOVA3.FP_Main_p34(i_sub,:));
    FTPL_Stat.MEG.Main_predp34_ANOVA3{3}(i_sub,:) = cell2mat(FTPLstat.SingleSub.ANOVA3.FP_Main_predp34(i_sub,:));
    FTPL_Stat.MEG.Interaction_p34predp34_ANOVA3{3}(i_sub,:) = cell2mat(FTPLstat.SingleSub.ANOVA3.FP_Interaction_p34predp34(i_sub,:));
end

%% Plot p-values in histogramm
common_ToneDur = {'short (150/200ms)','long (300/400ms)','all'};
effect_perTD = {'Main_p34_ANOVA2', 'Main_predp34_ANOVA2', 'Interaction_p34predp34_ANOVA2'};
effect_acrossTD = {'Main_p34_ANOVA3', 'Main_predp34_ANOVA3', 'Interaction_p34predp34_ANOVA3'};
cmap_patients = jet(length(vars.sub_list));

h = figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
CounterSubplot = 0;

for i_tonedur = 1:length(common_ToneDur);
    if i_tonedur < 3
        effect = effect_perTD;
    else
        effect = effect_acrossTD;
    end
    
    DimSubplot = [length(effect),length(common_ToneDur)]; 
    
    for i_effect = 1:length(effect)  
            
        CounterSubplot = CounterSubplot+1;
        subplot(DimSubplot(1),DimSubplot(2),CounterSubplot)
        histogram(log10(FTPL_Stat.MEG.(effect{i_effect}){i_tonedur}(:,2)),10)
        hold on; 
        for i_patients = 1:length(vars.sub_list)
            l = line([log10(FTPL_Stat.ECoG.(effect{i_effect}){i_tonedur}(i_patients,2)) log10(FTPL_Stat.ECoG.(effect{i_effect}){i_tonedur}(i_patients,2))],[0 10]);
            l.Color = cmap_patients(i_patients,:);
            l.LineWidth = 2;
        end
        l1 = legend('Healthy Subjects',...
            ['NY688, p = ' num2str(round(FTPL_Stat.ECoG.(effect{i_effect}){i_tonedur}(1,2),2))],...
            ['NY704, p = ' num2str(round(FTPL_Stat.ECoG.(effect{i_effect}){i_tonedur}(2,2),2))],...
            ['NY708, p = ' num2str(round(FTPL_Stat.ECoG.(effect{i_effect}){i_tonedur}(3,2),2))],...
            ['NY723, p = ' num2str(round(FTPL_Stat.ECoG.(effect{i_effect}){i_tonedur}(4,2),2))],...
            ['NY742, p = ' num2str(round(FTPL_Stat.ECoG.(effect{i_effect}){i_tonedur}(5,2),2))],...
            ['NY751, p = ' num2str(round(FTPL_Stat.ECoG.(effect{i_effect}){i_tonedur}(6,2),2))],...
            'Location','NorthWest');
        xlabel('log10(p-value)')
        ylabel('Frequency')
        title([effect{i_effect} ' ; TD = ' common_ToneDur{i_tonedur}],'Interpreter', 'none')
    end
end
suptitle('FTPLrating p-value comparison for healthy subjects (MEG) vs. patients (ECoG)')

%% Plot p-values in scatterplot
h = figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
CounterSubplot = 0;

for i_tonedur = 1:length(common_ToneDur)
    if i_tonedur < 3
        effect = effect_perTD;
    else
        effect = effect_acrossTD;
    end
    
    DimSubplot = [length(effect),length(common_ToneDur)]; 
    
    for i_effect = 1:length(effect)  
            
        CounterSubplot = CounterSubplot+1;
        subplot(DimSubplot(1),DimSubplot(2),CounterSubplot)
        
        s1 = scatter(log10(FTPL_Stat.MEG.(effect{i_effect}){i_tonedur}(:,2))',ones(1,20),...
            ones(1,20)*50,...
            'b','d','filled');

        hold on; 
        s2 = scatter(log10(FTPL_Stat.ECoG.(effect{i_effect}){i_tonedur}(:,2)),ones(1,length(vars.sub_list))*2,...
            ones(1,length(vars.sub_list))*50,...
                        cmap,'s','filled');    
  
%         legend('Healthy Subjects','ECoG patients','Location','NorthWest')
        
        hold on;
        line([-3 -3],[0 3],'LineStyle','--','Color','k')
        hold on;
        line([-2 -2],[0 3],'LineStyle','--','Color','k')
        hold on;
        line([-1.301 -1.301],[0 3],'LineStyle','--','Color','k')
        
        xlabel('log10(p-value)')
        ylabel('Sample (1 = sub, 2 = pat)')
        ylim([0.5 2.5])
        title([effect{i_effect} ' ; TD = ' common_ToneDur{i_tonedur}],'Interpreter', 'none')
    end
end
suptitle('FTPLrating p-value comparison for healthy subjects (MEG) vs. patients (ECoG)')

    