function NASTD_ECoG_HisTrack_CombineShuffKvals...
    (sub, ToneDur_text, inputData, ...
    SamplesTW, Label_TW,...
    plotFig_CompFolds, ...
    paths_NASTD_ECoG, vars)

%Aim: Combine shuffled Kprime values across folds for the shuffled condition
%by averaging each rep across folds to receive null distribution of 100
%K-prime values per time window/sensor

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
path_save = [paths_NASTD_ECoG.ECoGdata_HisTrack sub '/Shuff/' Label_TW '/'];
mkdir(path_save); %make new directory for output files

%Define directory for fold-wise uncombined Kprime data
path_load = [paths_NASTD_ECoG.ECoGdata_HisTrack sub '/Shuff/' Label_TW '/'];
%Define directory for potential figures
path_fig_Kprime = [paths_NASTD_ECoG.Fig_HisTrack_SingleSubs sub '/Shuff/KprimeAcrossFolds/'];
path_fig_MSR = [paths_NASTD_ECoG.Fig_HisTrack_SingleSubs sub '/Shuff/MinSquaredResidualsAcrossFolds/'];

%% 1) Load preprocessed data and fold-wise uncombined Kprime data
tic
disp([' -- Loading preprocessed data set for sub: ' sub])
load([si.path_preprocdata_sub]);
disp([' -- Done loading in ' num2str(toc) ' sec'])

tic
disp([' -- Loading fold-wise uncombined shuff Kprime data for sub: ' sub])
load([path_load sub '_' inputData '_' tonedur_title '_HisTrack_regstats_ShuffKperfold_unbalanced.mat']); %var name: stats_LinReg4Kprime, Info_TrialSelection
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
kModel_list = 1:NumModelOrder;
NumReps = size(stats_LinReg4Kprime{1}{1}.SumSquaredResiduals_TestvsTrain,2);

%% 3) Determine Kprime as K-model with minimal sum squared residuals
%3.1 Restructure regstats data and read out relevant stats for model selection per fold and rep
for i_fold = 1:NumFolds
    for i_win = 1:NumWindows
        for i_sensor = 1:NumSensors
            for i_rep = 1:NumReps
                
                for i_kModel = 1:NumModelOrder
                    SumSquaredResiduals_proxy(i_kModel) = stats_LinReg4Kprime{i_fold}{i_win, i_kModel, i_sensor}.SumSquaredResiduals_TestvsTrain(i_rep);
                    %Criteria for Kprime model selection = sum of squared residuals divided by length of model order for each k (per sensor,time win/fold)
                end
                
                %3.2 Read out minimum squared residuals divided by length of model order
                %and the respective position in array per sensor/time win/fold from all folds in common struct
                %(respective position in array = 1:16 = current tone (1) and up to 15 tones back)
                [minSumSquaredResiduals, index_minSumSquaredResiduals] = min( SumSquaredResiduals_proxy );
                
                %3.3 Determine Kprime value (i.e., preferred model order/Number of tones back in sequence history best explaining MEG data)
                %based on index where inimum squared residuals have minimal values
                %and correct value by -1, so that it ranges from 0 (only current tone) to 15 (i.e., current tone + 15 previous tones)
                Kprime_perFoldRep{i_win}{i_sensor}(i_fold, i_rep) = kModel_list(index_minSumSquaredResiduals) - 1;
                SumSquaredResiduals_Kprime_perFoldRep{i_win}{i_sensor}(i_fold, i_rep) = minSumSquaredResiduals;
                %minimum squared residuals divided by length of model order of K-value (i.e., preferred model order/Number of tones back)
                
                SumSquaredResiduals_perK_perFoldRep{i_win}{i_sensor}(i_fold, i_rep, :) = SumSquaredResiduals_proxy;
                %sum of squared residuals divided by length of model order, stored in common matrix across folds
                
                %3.4 Clean-up
                clear SumSquaredResiduals_proxy minSumSquaredResiduals index_minSumSquaredResiduals
                
            end
        end
    end
end

%% 4) Create Kprime null distribution by averaging Shuff Kprime across folds for each rep

%4.1 Compute Average and STD across folds (for each sensor/time win)
for i_win = 1:NumWindows
    for i_sensor = 1:NumSensors
        %Kprime
        Kprime_perRepavgFold{i_win}(i_sensor,:) = mean(Kprime_perFoldRep{i_win}{i_sensor},1);
        Kprime_perRepSDFold{i_win}(i_sensor,:) = std(Kprime_perFoldRep{i_win}{i_sensor});
        %SSR Kprime
        SumSquaredResiduals_Kprime_perRepavgFold{i_win}(i_sensor,:) = mean(SumSquaredResiduals_Kprime_perFoldRep{i_win}{i_sensor},1);
        SumSquaredResiduals_Kprime_perRepSDFold{i_win}(i_sensor,:) = std(SumSquaredResiduals_Kprime_perFoldRep{i_win}{i_sensor});
        %SSR all K
        SumSquaredResiduals_perK_perRepavgFold{i_win}{i_sensor} = squeeze(mean(SumSquaredResiduals_perK_perFoldRep{i_win}{i_sensor},1));
        SumSquaredResiduals_perK_perRepSDFold{i_win}{i_sensor} = squeeze(std(SumSquaredResiduals_perK_perFoldRep{i_win}{i_sensor}));
    end
end

%4.2 Average Kprime and SumSquaredRes across reps within fold
%(only for plotting and comparison between folds)
for i_fold = 1:NumFolds
    for i_win = 1:NumWindows
        for i_sensor = 1:NumSensors
            %Kprime
            Kprime_perFoldavgRep{i_win}(i_sensor,i_fold) = mean(Kprime_perFoldRep{i_win}{i_sensor}(i_fold,:));
            Kprime_perFoldSDRep{i_win}(i_sensor,i_fold) = std(Kprime_perFoldRep{i_win}{i_sensor}(i_fold,:));
            %SSR Kprime
            SumSquaredResiduals_Kprime_perFoldavgRep{i_win}(i_sensor,i_fold) = mean(SumSquaredResiduals_Kprime_perFoldRep{i_win}{i_sensor}(i_fold, :));
            SumSquaredResiduals_Kprime_perFoldSDRep{i_win}(i_sensor,i_fold) = std(SumSquaredResiduals_Kprime_perFoldRep{i_win}{i_sensor}(i_fold, :));
            %SSR all K
            SumSquaredResiduals_perK_perFoldavgRep{i_win}{i_sensor}(i_fold, :) = squeeze(mean(SumSquaredResiduals_perK_perFoldRep{i_win}{i_sensor}(i_fold, :, :),2));
            SumSquaredResiduals_perK_perFoldSDRep{i_win}{i_sensor}(i_fold, :) = squeeze(std(SumSquaredResiduals_perK_perFoldRep{i_win}{i_sensor}(i_fold, :, :)));
        end
    end
end

%4.2 Create summary file containing all relevant vars
ShuffKprime_data = struct;

ShuffKprime_data.perRepavgFold.Kprime_NullDistribution = Kprime_perRepavgFold;
ShuffKprime_data.perRepavgFold.Kprime_NullDistribution_SD = Kprime_perRepSDFold;
ShuffKprime_data.perRepavgFold.SumSquaredResiduals_Kprime = SumSquaredResiduals_Kprime_perRepavgFold;
ShuffKprime_data.perRepavgFold.SumSquaredResiduals_Kprime_SD = SumSquaredResiduals_Kprime_perRepSDFold;
ShuffKprime_data.perRepavgFold.SumSquaredResiduals_perK = SumSquaredResiduals_perK_perRepavgFold;
ShuffKprime_data.perRepavgFold.SumSquaredResiduals_perK_SD = SumSquaredResiduals_perK_perRepSDFold;

ShuffKprime_data.perFoldavgRep.Kprime = Kprime_perFoldavgRep;
ShuffKprime_data.perFoldavgRep.Kprime_SD = Kprime_perFoldSDRep;
ShuffKprime_data.perFoldavgRep.SumSquaredResiduals_Kprime = SumSquaredResiduals_Kprime_perFoldavgRep;
ShuffKprime_data.perFoldavgRep.SumSquaredResiduals_Kprime_SD = SumSquaredResiduals_Kprime_perFoldSDRep;
ShuffKprime_data.perFoldavgRep.SumSquaredResiduals_perK = SumSquaredResiduals_perK_perFoldavgRep;
ShuffKprime_data.perFoldavgRep.SumSquaredResiduals_perK_SD = SumSquaredResiduals_perK_perFoldSDRep;

%% 5) Plot kprime per set to compare consistency across sets (Optional)
if plotFig_CompFolds == 1
    mkdir(path_fig_Kprime)
    mkdir(path_fig_MSR)
    
    %OPTIONAL: plot shuffled Kprime per fold, averaged across reps within each fold (per sensor (selected sensors only), time win)
    figure;hold on; set(gcf,'units','normalized','outerposition',[0 0 1 1])
    suptitle(['Shuff Kprime per Fold (Avg+STD across reps (' num2str(NumReps) ')) per sensor for: ' sub '; ' inputData '; ' tonedur_title '; ' Label_TW ]);
    plot_counter = 0;
    
    for i_win = 1:NumWindows
        for i_Fold = 1:NumFolds
            plot_counter = plot_counter+1;
            subplot(NumWindows,NumFolds+1,plot_counter)
            
            shadedErrorBar(1:length(IndexSelElecs), Kprime_perFoldavgRep{i_win}(IndexSelElecs,i_Fold)' ,Kprime_perFoldSDRep{i_win}(IndexSelElecs,i_Fold)', 'lineprops', 'k-');
            hold on;
            plot(1:length(IndexSelElecs), Kprime_perFoldavgRep{i_win}(IndexSelElecs,i_Fold)', 'k-','LineWidth',1);
            
            title(['TW =' num2str(i_win) '; Fold: ' num2str(i_Fold)]);
            set(gca,'FontSize',8)
            ylabel('Kprime');
            xlabel('Electrodes');
            ylim([-2 20])
            xlim([0 length(IndexSelElecs)+1])
            hold on;
        end
        
        %Plot average across folds
        plot_counter = plot_counter+1;
        subplot(NumWindows,NumFolds+1,plot_counter)
        
        shadedErrorBar(1:length(IndexSelElecs), mean(Kprime_perFoldavgRep{i_win}(IndexSelElecs,:),2)' ,std(Kprime_perFoldavgRep{i_win}(IndexSelElecs,:)'), 'lineprops', 'b-');
        hold on;
        plot(1:length(IndexSelElecs), mean(Kprime_perFoldavgRep{i_win}(IndexSelElecs,:),2)','b-','LineWidth',1);
        
        title(['TW =' num2str(i_win) ';Avg (SD) across Folds']);
        set(gca,'FontSize',8)
        ylabel('Kprime');
        xlabel('Electrodes');
        ylim([-2 20])
        xlim([0 length(IndexSelElecs)+1])
        hold on;
    end
    
    %save Fig
    filename = ['ShuffKprimeAcrossFolds_' sub '_' inputData '_' tonedur_title '_' Label_TW '.png'];
    figfile      = [path_fig_Kprime filename];
    saveas(gcf, [figfile], 'png'); %save png version
   
    
    %OPTIONAL: plot difference in minimum squared residuals divided by length of model order
    %Value by which Kprime is determined across odd/even sets (per sensor, time win)
    figure;hold on; set(gcf,'units','normalized','outerposition',[0 0 1 1])
    suptitle(['MinSumSquaredRes/ModelOrder (defining Kprime - Avg+STD across folds (' num2str(NumFolds) ')) per sensor for: ' sub '; ' inputData '; ' tonedur_title '; ' Label_TW ]);
    plot_counter = 0;
    
    for i_win = 1:NumWindows
        for i_Fold = 1:NumFolds
            plot_counter = plot_counter+1;
            subplot(NumWindows,NumFolds+1,plot_counter)
            
            shadedErrorBar(1:length(IndexSelElecs), SumSquaredResiduals_Kprime_perFoldavgRep{i_win}(IndexSelElecs,i_Fold)', SumSquaredResiduals_Kprime_perFoldSDRep{i_win}(IndexSelElecs,i_Fold)', 'lineprops', 'k-');
            hold on;
            plot(1:length(IndexSelElecs), SumSquaredResiduals_Kprime_perFoldavgRep{i_win}(IndexSelElecs,i_Fold)', 'k-','LineWidth',1);
            
            title(['TW =' num2str(i_win) '; Fold: ' num2str(i_Fold)]);
            set(gca,'FontSize',8)
            ylabel('MinSumSquaredRes/ModelOrder');
            xlabel('Electrodes');
            ylim auto
            xlim([0 length(IndexSelElecs)+1])
            hold on;
        end
        %Plot average across folds
        plot_counter = plot_counter+1;
        subplot(NumWindows,NumFolds+1,plot_counter)
        
        shadedErrorBar(1:length(IndexSelElecs), mean(SumSquaredResiduals_Kprime_perFoldavgRep{i_win}(IndexSelElecs,:),2)' ,std(SumSquaredResiduals_Kprime_perFoldavgRep{i_win}(IndexSelElecs,:)'), 'lineprops', 'b-');
        hold on;
        plot(1:length(IndexSelElecs), mean(SumSquaredResiduals_Kprime_perFoldavgRep{i_win}(IndexSelElecs,:),2)','b-','LineWidth',1);
        
        title(['TW =' num2str(i_win) ';Avg (SD) across Folds']);
        set(gca,'FontSize',8)
        ylabel('MinSumSquaredRes/ModelOrder');
        xlabel('Electrodes');
        ylim auto
        xlim([0 length(IndexSelElecs)+1])
        hold on;
    end
    
    %save Fig
    filename = ['ShuffminSSRAcrossFolds_' sub '_' inputData '_' tonedur_title '_' Label_TW '.png'];
    figfile      = [path_fig_MSR filename];
    saveas(gcf, [figfile], 'png'); %save png version
    
    close all
end


%% 6) Save output variables
%Combined output file with 1) combined original k-prime values, combined
%2) shuffled k-prime values, 3) amount of shuffled_
savefile = [path_save sub '_' inputData '_' tonedur_title '_HisTrack_ShuffKprimeCombFolds_unbalanced.mat'];
save(savefile, 'ShuffKprime_data','-v7.3');
end

