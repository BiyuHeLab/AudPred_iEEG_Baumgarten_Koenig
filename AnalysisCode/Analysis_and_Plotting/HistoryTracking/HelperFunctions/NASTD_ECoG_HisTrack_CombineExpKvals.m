function NASTD_ECoG_HisTrack_CombineExpKvals...
    (sub, ToneDur_text, inputData, ...
    SamplesTW, Label_TW,...
    plotFig_CompFolds, ...
    paths_NASTD_ECoG, vars)
%Aim: Combine output from all k-fold crossvalidation sets

%Additions:
%-Read out regression weights (beta weights) for winning k-model/k-prime
%-Plot K-prime and MeanSquareResiduals (parameter on which Kprime selection
%is based) for both sets to allow for comparison across sets

%% 0.1) Specify vars & paths
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')

%Add base dir and own script dir
addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));

%% 0.2) Determine subject-specific parameters (whole-recording)
NASTD_ECoG_subjectinfo %load subject info file (var: si)
subs_PreProcSettings = NASTD_ECoG_SubPreprocSettings; %load in file with individual preproc infos

%Tone duration condition for data load-in
if strcmp(ToneDur_text,'0.2')
    tonedur_title = '200msTD';
elseif  strcmp(ToneDur_text,'0.4')
    tonedur_title = '400msTD';
end

%Define directory for combined Kprime data
path_save = [paths_NASTD_ECoG.ECoGdata_HisTrack sub '/Exp/' Label_TW '/'];
mkdir(path_save); %make new directory for output files

%Define directory for fold-wise uncombined Kprime data
path_load = [paths_NASTD_ECoG.ECoGdata_HisTrack sub '/Exp/' Label_TW '/'];

%Define directory for potential figures
path_fig_Kprime = [paths_NASTD_ECoG.Fig_HisTrack_SingleSubs sub '/Exp/KprimeAcrossFolds/'];
path_fig_MSR = [paths_NASTD_ECoG.Fig_HisTrack_SingleSubs sub '/Exp/MinSquaredResidualsAcrossFolds/'];

%% 1) Load preprocessed data and fold-wise uncombined Kprime data
tic
disp([' -- Loading preprocessed data set for sub: ' sub])
load([si.path_preprocdata_sub]);
disp([' -- Done loading in ' num2str(toc) ' sec'])

tic
disp([' -- Loading fold-wise uncombined Kprime data for sub: ' sub])
load([path_load sub '_' inputData '_' tonedur_title '_HisTrack_regstats_ExpKperfold_unbalanced.mat']); %var name: stats_LinReg4Kprime, Info_TrialSelection
disp([' -- Done loading in ' num2str(toc) ' sec'])

%% 2) Define analysis parameters
fsample         = preprocData_AllTrials.fsample;
toneDur_inSecs  = str2num(ToneDur_text);

NumSamplesPerTone = toneDur_inSecs * fsample;
NumFolds = length(stats_LinReg4Kprime); %Dimensions stats_LinReg4Kprime = {i_fold}{i_win, i_k, i_sensor}
NumSensors = size(stats_LinReg4Kprime{1},3);
IndexSelElecs = setdiff(preprocData_AllTrials.cfg.info_elec.selected.index4EDF, subs_PreProcSettings.(sub).rejectedChan_index);

%Define time windows used for analysis
win_size    = SamplesTW;
win_overlap = 0;

%Define number, start and end sample of window per tone
windows = [1 win_size];
while windows(end,end) < NumSamplesPerTone
    windows = [windows; windows(end,:) + (win_size - win_overlap)];
end

if windows(end,end) > NumSamplesPerTone %TJB: windows within single tone
    windows(end,:) = [];
end
NumWindows = size(windows,1);

NumModelOrder = size(stats_LinReg4Kprime{1},2); %number model order;
k_list = 1:NumModelOrder;

%% 3) Determine Kprime as 'k-model with minimal sum squared residuals
%3.1 Restructure regstats data and read out relevant stats for model selection from each fold
for i_win = 1:NumWindows
    for i_sensor = 1:NumSensors
        for i_fold = 1:NumFolds
            
            for i_kModel = 1:NumModelOrder
                SumSquaredResiduals_proxy(i_kModel) = stats_LinReg4Kprime{i_fold}{i_win, i_kModel, i_sensor}.SumSquaredResiduals_TestvsTrain;
                %Criteria for Kprime model selection = sum of squared residuals divided by length of model order for each k (per sensor,time win/fold)
            end
            
            %3.2 Read out minimum squared residuals divided by length of model order
            %and the respective position in array per sensor/time win/fold from all folds in common struct
            %(respective position in array = 1:16 = current tone (1) and up to 15 tones back)
            [minSumSquaredResiduals, index_minSumSquaredResiduals] = min( SumSquaredResiduals_proxy );
            
            %3.3 Determine Kprime value (i.e., preferred model order/Number of tones back in sequence history best explaining MEG data)
            %based on index where minimum squared residuals have minimal values
            %and correct value by -1, so that it ranges from 0 (only current tone) to 15 (i.e., current tone + 15 previous tones)
            kprime_allFolds{i_win}(i_fold,i_sensor) = k_list(index_minSumSquaredResiduals) - 1;
            kprime_SumSquaredResiduals_allFolds{i_win}(i_fold,i_sensor) = minSumSquaredResiduals;
            %minimum squared residuals divided by length of model order of K-value (i.e., preferred model order/Number of tones back)
            
            SumSquaredResiduals_allK_allFolds{i_win}{i_sensor}(i_fold, :) = SumSquaredResiduals_proxy;
            %sum of squared residuals divided by length of model order, stored in common matrix across folds
            
            %3.4 Average beta regression weights for each k-model across sets
            for i_kModel = 1:NumModelOrder
                BetaRegWeights_allFolds{i_win, i_kModel, i_sensor}(i_fold,:) =  stats_LinReg4Kprime{i_fold}{i_win, i_kModel, i_sensor}.tstat.beta';
            end
            
            %3.5 Clean-up
            clear SumSquaredResiduals_proxy minSumSquaredResiduals index_minSumSquaredResiduals
            
        end
    end
end

%% 4) Combine Kprime across folds by averaging values across cross-validation folds
%i.e., average model selection criteria by averaging across folds
for i_win = 1:NumWindows
    for i_sensor = 1:NumSensors
        
        %4.1 Compute Average and STD across folds (for each sensor/time win)
        %Kprime
        kprime_AVGacrossFolds{i_win}(i_sensor) = mean(kprime_allFolds{i_win}(:,i_sensor));
        kprime_STDacrossFolds{i_win}(i_sensor) = std(kprime_allFolds{i_win}(:,i_sensor));
        %SSR Kprime
        kprime_sumsquaredRes_AVGacrossFolds{i_win}(i_sensor) = mean(kprime_SumSquaredResiduals_allFolds{i_win}(:,i_sensor));
        kprime_sumsquaredRes_STDacrossFolds{i_win}(i_sensor) = std(kprime_SumSquaredResiduals_allFolds{i_win}(:,i_sensor));
        %SSR all K
        sumsquaredRes_allK_AVGacrossFolds{i_win}(i_sensor, :) = mean(SumSquaredResiduals_allK_allFolds{i_win}{i_sensor});
        sumsquaredRes_allK_STDacrossFolds{i_win}(i_sensor, :) = std(SumSquaredResiduals_allK_allFolds{i_win}{i_sensor});
        
        %Beta Regression Weights
        %all K
        for i_kModel = 1:NumModelOrder
            BetaRegWeights_AVGacrossFolds{i_win, i_kModel, i_sensor} = mean(BetaRegWeights_allFolds{i_win, i_kModel, i_sensor});
            BetaRegWeights_STDacrossFolds{i_win, i_kModel, i_sensor} = std(BetaRegWeights_allFolds{i_win, i_kModel, i_sensor});
        end
        %K-prime
        kprime_BetaRegWeight_AVGacrossFolds{i_win}{i_sensor} = ...
            mean(BetaRegWeights_AVGacrossFolds{i_win, round(kprime_AVGacrossFolds{i_win}(i_sensor) + 1), i_sensor});
        kprime_BetaRegWeight_STDacrossFolds{i_win}{i_sensor} = ...
            std(BetaRegWeights_AVGacrossFolds{i_win, round(kprime_AVGacrossFolds{i_win}(i_sensor) + 1), i_sensor});
        %Note: Beta Regression Weights are selected based on 1)
        %across-fold-averaged beta weights and 2) across-fold-averaged Kprime
        %values (here -1 lin reg offset term correction is undone, since we
        %need position in 1:16 aray, not Kprime value)
    end
end

%4.2 Create summary file containing all relevant vars
ExpKprime_data = struct;
ExpKprime_data.perFold = struct;
ExpKprime_data.perFold.Kprime = kprime_allFolds;
ExpKprime_data.perFold.SSR_perK = SumSquaredResiduals_allK_allFolds;
ExpKprime_data.perFold.SSR_Kprime = kprime_SumSquaredResiduals_allFolds;
ExpKprime_data.perFold.BetaRegWeights_perK = BetaRegWeights_allFolds;
ExpKprime_data.avgFolds = struct;
ExpKprime_data.avgFolds.Kprime_Avg = kprime_AVGacrossFolds;
ExpKprime_data.avgFolds.Kprime_STD = kprime_STDacrossFolds;
ExpKprime_data.avgFolds.SSR_Kprime_Avg = kprime_sumsquaredRes_AVGacrossFolds;
ExpKprime_data.avgFolds.SSR_Kprime_STD = kprime_sumsquaredRes_STDacrossFolds;
ExpKprime_data.avgFolds.SSR_perK_Avg = sumsquaredRes_allK_AVGacrossFolds;
ExpKprime_data.avgFolds.SSR_perK_STD = sumsquaredRes_allK_STDacrossFolds;
ExpKprime_data.avgFolds.BetaRegWeights_perK_Avg = BetaRegWeights_AVGacrossFolds;
ExpKprime_data.avgFolds.BetaRegWeights_perK_STD = BetaRegWeights_STDacrossFolds;
ExpKprime_data.avgFolds.BetaRegWeights_Kprime_Avg = kprime_BetaRegWeight_AVGacrossFolds;
ExpKprime_data.avgFolds.BetaRegWeights_Kprime_STD = kprime_BetaRegWeight_STDacrossFolds;


%% 5) Plot kprime per fold to compare consistency across sets (Optional)
if plotFig_CompFolds == 1
    mkdir(path_fig_Kprime)
    mkdir(path_fig_MSR)
    
    %OPTIONAL: plot difference in K' across folds (per sensor (selected sensors only), time win)
    figure;hold on; set(gcf,'units','normalized','outerposition',[0 0 1 1])
    suptitle(['EXP Kprime (Avg+STD across folds (' num2str(NumFolds) ')) per sensor for: ' sub '; ' inputData '; ' tonedur_title '; ' Label_TW ]);
    plot_counter = 0;
    
    for i_win = 1:NumWindows
        plot_counter = plot_counter+1;
        subplot(NumWindows,1,plot_counter)
        
        shadedErrorBar(1:length(IndexSelElecs), kprime_AVGacrossFolds{i_win}(IndexSelElecs),kprime_STDacrossFolds{i_win}(IndexSelElecs), 'lineprops', 'k-');
        hold on;
        plot(1:length(IndexSelElecs),kprime_AVGacrossFolds{i_win}(IndexSelElecs),'k','LineWidth',2)       
        
        title(['TW =' num2str(i_win)]);
        ElecLabels = [];
        for i_chan = 1:length(IndexSelElecs)
            ElecLabels = [ElecLabels, preprocData_AllTrials.label(IndexSelElecs(i_chan))];
        end
        set(gca,'xtick',1:length(IndexSelElecs))
        set(gca,'FontSize',8,'XTickLabelRotation',45)
        xticklabels(ElecLabels);
        ylabel('Kprime');
        
        ylim([-2 20])
        xlim([0 length(IndexSelElecs)+1])
        hold on;
    end
    
    %save Fig
    filename = ['EXPKprimeAcrossFolds_' sub '_' inputData '_' tonedur_title '_' Label_TW '.png'];
    figfile      = [path_fig_Kprime filename];
    saveas(gcf, [figfile], 'png'); %save png version
    
    
    
    %OPTIONAL: plot difference in minimum squared residuals divided by length of model order
    %Value by which Kprime is determined across odd/even sets (per sensor, time win)
    figure;hold on; set(gcf,'units','normalized','outerposition',[0 0 1 1])
    suptitle(['MinSumSquaredRes/ModelOrder (defining Kprime - Avg+STD across folds (' num2str(NumFolds) ')) per sensor for: ' sub '; ' inputData '; ' tonedur_title '; ' Label_TW ]);
    
    plot_counter = 0;
    
    for i_win = 1:NumWindows
        plot_counter = plot_counter+1;
        subplot(NumWindows,1,plot_counter)
        
        shadedErrorBar(1:length(IndexSelElecs), kprime_sumsquaredRes_AVGacrossFolds{i_win}(IndexSelElecs),kprime_sumsquaredRes_STDacrossFolds{i_win}(IndexSelElecs), 'lineprops', 'k-');
        hold on;
        plot(1:length(IndexSelElecs),kprime_sumsquaredRes_AVGacrossFolds{i_win}(IndexSelElecs),'k','LineWidth',2)
                   
        title(['TW =' num2str(i_win)]);
        ElecLabels = [];
        for i_chan = 1:length(IndexSelElecs)
            ElecLabels = [ElecLabels, preprocData_AllTrials.label(IndexSelElecs(i_chan))];
        end
        set(gca,'xtick',1:length(IndexSelElecs))
        set(gca,'FontSize',8,'XTickLabelRotation',45)
        xticklabels(ElecLabels);
        ylabel('MinSumSquaredRes/ModelOrder');
        
        ylim auto
        xlim([0 length(IndexSelElecs)+1])
        hold on;
        
    end
    
    %save Fig
    filename = ['EXPminSSRAcrossFolds_' sub '_' inputData '_' tonedur_title '_' Label_TW '.png'];
    figfile      = [path_fig_MSR filename];
    saveas(gcf, [figfile], 'png'); %save png version
    
    close all
end

%% 6) Save output variables
savefile = [path_save sub '_' inputData '_' tonedur_title '_HisTrack_ExpKprimeCombFolds_unbalanced.mat'];
save(savefile, 'ExpKprime_data','-v7.3');

end