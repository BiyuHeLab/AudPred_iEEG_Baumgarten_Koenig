function NASTD_ECoG_HisTrack_CombineExpShuffK...
    (sub, input_ToneDurLabel, input_DataLabel, ...
    SamplesTW, Label_TW, NumRuns, ...
    plotFig_CompFolds, ...
    paths_NASTD_ECoG)

%Aim: Combine Exp and Shuffled Kprime distribution to determine significance for Exp Kprime

%% 0.1) Specify vars & paths
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
%Add base dir and own script dir
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));

%% 0.2) Determine subject-specific parameters (whole-recording)
NASTD_ECoG_subjectinfo %load subject info file (var: si)
subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;

%Tone duration condition for data load-in
if strcmp(input_ToneDurLabel,'0.2')
    tonedur_title = '200msTD';
elseif  strcmp(input_ToneDurLabel,'0.4')
    tonedur_title = '400msTD';
end

%Define directory for combined Kprime data
path_save = [paths_NASTD_ECoG.ECoGdata_HisTrack sub '/ExpvsShuff/' Label_TW '/'];
if (~exist(path_save, 'dir')); mkdir(path_save); end

%Define directory for fold-wise uncombined Kprime data
path_load_exp = [paths_NASTD_ECoG.ECoGdata_HisTrack sub '/Exp/' Label_TW '/'];
path_load_shuff = [paths_NASTD_ECoG.ECoGdata_HisTrack sub '/Shuff/' Label_TW '/'];

%Define directory for potential figures
path_fig_avgFolds = [paths_NASTD_ECoG.ECoGdata_HisTrack sub '/ExpvsShuff/' Label_TW '/Figs/Kprime_Folds/'];
if (~exist(path_fig_avgFolds, 'dir')); mkdir(path_fig_avgFolds); end
path_fig_avgRuns = [paths_NASTD_ECoG.ECoGdata_HisTrack sub '/ExpvsShuff/' Label_TW '/Figs/Kprime_Runs/'];
if (~exist(path_fig_avgRuns, 'dir')); mkdir(path_fig_avgRuns); end

%% 1) Load preprocessed data and fold-wise uncombined Kprime data
tic
disp([' -- Loading preprocessed data set for sub: ' sub])
load([si.path_preprocdata_sub]);
disp([' -- Done loading in ' num2str(toc) ' sec'])

tic
disp([' -- Loading combined exp Kprime data for sub: ' sub])
load([path_load_exp sub '_' input_DataLabel '_' tonedur_title ...
    '_HisTrack_ExpKprimeCombFolds_unbalanced_' num2str(NumRuns) 'runs.mat'], ...
    'ExpKprime_data', 'labels_loadedData'); 
disp([' -- Done loading in ' num2str(toc) ' sec'])

tic
disp([' -- Loading combined shuff Kprime data for sub: ' sub])
load([path_load_shuff sub '_' input_DataLabel '_' tonedur_title ...
    '_HisTrack_ShuffKprimeCombFolds_unbalanced.mat'], ...
    'ShuffKprime_data'); 
disp([' -- Done loading in ' num2str(toc) ' sec'])


%% 2) Define analysis parameters
fsample         = DataClean_AllTrials.fsample;
toneDur_inSecs  = str2num(input_ToneDurLabel);

NumSamplesPerTone   = toneDur_inSecs * fsample;
NumSensors          = size(ExpKprime_data.avgRuns.Kprime_Avg{1},1);

clear DataClean_AllTrials

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
NumWindows  = size(windows,1);
NumReps     = size(ShuffKprime_data.perRepavgFold.Kprime_NullDistribution{1},2);
NumFolds    = size(ExpKprime_data.perFold{1}.Kprime{1},1);
NumRuns     = length(ExpKprime_data.avgFolds);


%% 3) Average repetitions in shuffled data across number of runs in exp data
%Exp Kprime esteimates were averaged across 5 runs to increase reliability.
%Perform the same averaging for shuffled Kprime estimates to reach a 
%comparable level of noise.

disp([' -- Averaging shuffled Kprime values across ' num2str(NumRuns) ' randomly selected repetitions --'])

for i_win = 1:NumWindows
    
    %Create empty struct
    ShuffKprime_data.perRepavgRuns.Kprime_NullDistribution{i_win} = ...
        zeros(NumSensors, NumReps);
    
    for i_sensor = 1:NumSensors        
        for i_reps = 1:NumReps
            
            %Randomly select NumRuns repetitions
            selected_reps = randi(NumReps, [1,NumRuns]);
            ShuffKprime_data.perRepavgRuns.selected_reps{i_win}{i_sensor, i_reps} ...
                = selected_reps;
            %Average shuffled Kprime estimates across randomly selected repetitions
            ShuffKprime_data.perRepavgRuns.Kprime_NullDistribution{i_win}(i_sensor, i_reps) = ...
                mean(ShuffKprime_data.perRepavgFold.Kprime_NullDistribution{i_win}(i_sensor, selected_reps));
            
        end
    end
end

% %Optional: Plot averaged and non averaged shuffled k' estimates
% i_sensor = randi(NumSensors, [1,1]);
% 
% figure;
% set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
% 
% subplot(1,2,1)
% plot(ShuffKprime_data.perRepavgFold.Kprime_NullDistribution{1}(i_sensor,:),'b-','LineWidth',1)
% hold on
% plot(ShuffKprime_data.perRepavgRuns.Kprime_NullDistribution{1}(i_sensor,:),'r-','LineWidth',2)
% ylim([0 15])
% legend('non-averaged shuff Kprime', 'averaged shuff Kprime')
% ylabel('Kprime')
% xlabel('Repetition')
% title('Kprime estimate per repetition')
% 
% subplot(1,2,2)
% bin_edges = [0:1:15];
% histogram(ShuffKprime_data.perRepavgFold.Kprime_NullDistribution{1}(i_sensor,:),bin_edges,'FaceColor','b','FaceAlpha',0.5)
% hold on;
% histogram(ShuffKprime_data.perRepavgRuns.Kprime_NullDistribution{1}(i_sensor,:),bin_edges,'FaceColor','r','FaceAlpha',0.5)
% xlim([0 15])
% xlabel('Kprime')
% ylabel('Frequency')
% title('Histogram: Kprime estimate per repetition')
% 
% sgtitle(['Shuff Kprime comparison between non-averaged and across-rep-averaged computation for channel ' num2str(i_sensor)])


%% 4) Compute number of reps in shuffled Kprime distribution higher than in experimental correlate
%Computes the amount of shuffled k-prime values that are
%equal or bigger to the experimental/original k-prime value and
%divides this by the number of reps
%for each time window and sensor, in each repetition
Kprime_data.Exp                 = ExpKprime_data;
% Kprime_data.Shuff.perRepavgRuns = ShuffKprime_data.perRepavgFold; %not averaged across reps
Kprime_data.Shuff.perRepavgRuns = ShuffKprime_data.perRepavgRuns; %averaged across reps

for i_win = 1:NumWindows    
    for i_sensor = 1:NumSensors

        Kprime_data.Exp.avgRuns.pval_ExpvsShuff{i_win}(i_sensor,1) = ...
            sum(ShuffKprime_data.perRepavgRuns.Kprime_NullDistribution{i_win}(i_sensor, :) >= ...
            ExpKprime_data.avgRuns.Kprime_Avg{i_win}(i_sensor)) ...
            / length(ShuffKprime_data.perRepavgRuns.Kprime_NullDistribution{i_win}(i_sensor, :));
        
        %Add binary significance coding for plotting
        if Kprime_data.Exp.avgRuns.pval_ExpvsShuff{i_win}(i_sensor,1) < 0.05                              
            Kprime_data.Exp.avgRuns.mask_ExpvsShuff{i_win}(i_sensor,1) = 1;
        else
            Kprime_data.Exp.avgRuns.mask_ExpvsShuff{i_win}(i_sensor,1) = 0;  
        end
        
%         pval_ExpvsShuff_testnonavg{i_win}(i_sensor,1) = ...
%             sum(ShuffKprime_data.perRepavgFold.Kprime_NullDistribution{i_win}(i_sensor, :) >= ...
%             ExpKprime_data.avgRuns.Kprime_Avg{i_win}(i_sensor)) ...
%             / length(ShuffKprime_data.perRepavgFold.Kprime_NullDistribution{i_win}(i_sensor, :));       
%         if pval_ExpvsShuff_testnonavg{i_win}(i_sensor,1) < 0.05                              
%             mask_ExpvsShuff_testnonavg{i_win}(i_sensor,1) = 1;
%         else
%             mask_ExpvsShuff_testnonavg{i_win}(i_sensor,1) = 0;  
%         end
                                                                              
    end
end                                                                                                                                                                                       

% %Optional: Plot p-values for averaged and non-averaged shuffled Kprime estimates
% signelecs_nonavg = sum(pval_ExpvsShuff_testnonavg{1} < 0.05);
% signelecs_avg = sum(Kprime_data.Exp.avgRuns.pval_ExpvsShuff{1} < 0.05);
% 
% figure;
% set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
% 
% subplot(2,2,1)
% plot(pval_ExpvsShuff_testnonavg{1},'b-','LineWidth',2)
% hold on
% plot(Kprime_data.Exp.avgRuns.pval_ExpvsShuff{1},'r-','LineWidth',1)
% hold on
% plot(1:NumSensors,ones(1,NumSensors)*0.05,'k-','LineWidth',3)
% ylim([0 1])
% legend('p-value for non-averaged shuff Kprime', 'p-value for averaged shuff Kprime')
% ylabel('p-value')
% xlabel('Sensor')
% title('p-values (uncorrected) per sensor')
% title({['p-values (uncorrected) per sensor'], ...
%     ['Sign. elecs non-avg = ' num2str(signelecs_nonavg) '; Sign. elecs avg = ' num2str(signelecs_avg)]})
% 
% subplot(2,2,2)
% bin_edges = [0:0.05:1];
% histogram(pval_ExpvsShuff_testnonavg{1},bin_edges,'FaceColor','b','FaceAlpha',0.5)
% hold on;
% histogram(Kprime_data.Exp.avgRuns.pval_ExpvsShuff{1},bin_edges,'FaceColor','r','FaceAlpha',0.5)
% hold on
% xline(0.05,'k-','LineWidth',3)
% xlim([0 1])
% xlabel('p-value')
% ylabel('Frequency')
% title({['Histogram: p-values (uncorrected) across sensors'], ...
%     ['Sign. elecs non-avg = ' num2str(signelecs_nonavg) '; Sign. elecs avg = ' num2str(signelecs_avg)]})
% 
% 
% pval_FDRcorrected_noavg = mafdr(pval_ExpvsShuff_testnonavg{1},'BHFDR', true);
% pval_FDRcorrected_avg = mafdr(Kprime_data.Exp.avgRuns.pval_ExpvsShuff{1},'BHFDR', true);
% 
% signelecs_nonavg = sum(pval_FDRcorrected_noavg < 0.05);
% signelecs_avg = sum(pval_FDRcorrected_avg < 0.05);
% 
% subplot(2,2,3)
% plot(pval_FDRcorrected_noavg,'b-','LineWidth',2)
% hold on
% plot(pval_FDRcorrected_avg,'r-','LineWidth',1)
% hold on
% plot(1:NumSensors,ones(1,NumSensors)*0.05,'k-','LineWidth',3)
% ylim([0 1])
% ylabel('p-value')
% xlabel('Sensor')
% title({['p-values (FDR-corrected) per sensor'], ...
%     ['Sign. elecs non-avg = ' num2str(signelecs_nonavg) '; Sign. elecs avg = ' num2str(signelecs_avg)]})
%      
% subplot(2,2,4)
% bin_edges = [0:0.05:1];
% histogram(pval_FDRcorrected_noavg,bin_edges,'FaceColor','b','FaceAlpha',0.5)
% hold on;
% histogram(pval_FDRcorrected_avg,bin_edges,'FaceColor','r','FaceAlpha',0.5)
% hold on
% xline(0.05,'k-','LineWidth',3)
% xlim([0 1])
% xlabel('p-value')
% ylabel('Frequency')
% title({['Histogram: p-values (FDR-corrected) across sensors'], ...
%     ['Sign. elecs non-avg = ' num2str(signelecs_nonavg) '; Sign. elecs avg = ' num2str(signelecs_avg)]})
% 
% sgtitle(['p-value comparison between non-averaged and across-rep-averaged computation for ' sub])
        

%% 5) Save output variables
%Combined output file
savefile = [path_save sub '_' input_DataLabel '_' tonedur_title ...
    '_HisTrack_ExpvsShuffKprimeCombFolds_unbalanced.mat'];
save(savefile, 'Kprime_data', 'labels_loadedData','-v7.3');

%% 6) Optional: Plot and compare EXP vs SHUFF
%6.1 plot Exp vs. Shuff Kprime per fold (per sensor (selected sensors only), time win)
if plotFig_CompFolds == 1
    
    %Not feasible, since now 5 folds per run for exp data
          
%     figure;hold on; set(gcf,'units','normalized','outerposition',[0 0 1 1])
%     sgtitle(['Exp vs. Shuff Kprime per Fold (Avg+SD across reps (' num2str(NumReps) ')) per sensor for: ' ...
%         sub '; ' input_DataLabel '; ' tonedur_title '; ' Label_TW]);
%     
%     plot_counter = 0;
%     
%     for i_win = 1:NumWindows
%         for i_Fold = 1:NumFolds
%             
%             plot_counter = plot_counter+1;
%             subplot(NumWindows,NumFolds+1,plot_counter)
%             
%             plot(1:NumSensors, ExpKprime_data.perFold.Kprime{i_win}(i_Fold,:), 'r-','LineWidth',2);
%             hold on;
%             shadedErrorBar(1:NumSensors, ...
%                 ShuffKprime_data.perFoldavgRep.Kprime{i_win}(:,i_Fold)',...
%                 ShuffKprime_data.perFoldavgRep.Kprime_SD{i_win}(:,i_Fold)', 'lineprops', 'k-');
%             hold on;
%             plot(1:NumSensors, ShuffKprime_data.perFoldavgRep.Kprime{i_win}(:,i_Fold)', 'k-','LineWidth',2);
%             
%             title(['TW =' num2str(i_win) '; Fold: ' num2str(i_Fold)]);
%             set(gca,'FontSize',8)
%             ylabel('Kprime');
%             xlabel('Electrodes');
%             ylim([-2 20])
%             xlim([0 NumSensors+1])
%             
%             if plot_counter == 1
%                 if NumWindows < 5
%                     legend('Exp Kprime','Shuff Kprime (Avg+SD across reps)')
%                 end
%             end
%         end
%         
%         %Plot average across folds
%         plot_counter = plot_counter+1;
%         subplot(NumWindows,NumFolds+1,plot_counter)
%         
%         shadedErrorBar(1:NumSensors, ...
%             ExpKprime_data.avgFolds.Kprime_Avg{i_win}(:),...
%             ExpKprime_data.avgFolds.Kprime_STD{i_win}(:), 'lineprops', 'r-');        
%         hold on; 
%         plot(1:NumSensors, ExpKprime_data.avgFolds.Kprime_Avg{i_win}(:), 'r-','LineWidth',2);
%         hold on;
%         
%         shadedErrorBar(1:NumSensors, ...
%             mean(ShuffKprime_data.perRepavgFold.Kprime_NullDistribution{i_win}(:,:),2)',...
%             std(ShuffKprime_data.perRepavgFold.Kprime_NullDistribution{i_win}(:,:)'), 'lineprops', 'k-');
%         hold on;
%         plot(1:NumSensors, mean(ShuffKprime_data.perRepavgFold.Kprime_NullDistribution{i_win}(:,:),2)', 'k-','LineWidth',2);
%         hold on;
%         
%         area(1:NumSensors, Kprime_data.Exp.avgFolds.mask_ExpvsShuff{i_win}(:)*20,...
%         'basevalue',0,'FaceColor',[0.1, 0.1, 0.8],'FaceAlpha', 0.5,'LineStyle','none');
%         
%         title(['TW =' num2str(i_win) '; Avg (SD) across Folds']);
%         set(gca,'FontSize',8)
%         ylabel('Kprime');
%         xlabel('Electrodes');
%         ylim([-2 20])
%         xlim([0 NumSensors+1])
%         hold on;
%         
%     end
%     %save Fig
%     filename = ['ExpvsShuffKprime_perFold_' sub '_' input_DataLabel '_' tonedur_title '_' Label_TW '.png'];
%     figfile      = [path_fig filename];
%     saveas(gcf, [figfile], 'png'); %save png version
%     close;    
       
    %Not feasible, since now shuff Kprime averaged across runs

% %6.2 plot Exp vs. Shuff Kprime avg across folds / runs (for Exp) for each TW
%    figure;hold on; set(gcf,'units','normalized','outerposition',[0 0 1 1])
%     sgtitle(['Exp vs. Shuff Kprime across Folds/Runs (Shuff: Avg+SD across folds & Reps);' num2str(NumReps) ')) per sensor for: ' ...
%         sub '; ' input_DataLabel '; ' tonedur_title '; ' Label_TW]);
%     
%     plot_counter = 0; 
%     
%     if NumWindows < 5
%         Plot_rows = 1;
%     else 
%         Plot_rows = 2;
%     end
%     
%     for i_win = 1:NumWindows
%         
%         plot_counter = plot_counter+1;
%         
%         if NumWindows < 5   
%             subplot(NumWindows,Plot_rows,plot_counter)
%         else
%             subplot(NumWindows/2,Plot_rows,plot_counter)
%         end
%             
%         shadedErrorBar(1:NumSensors, ...
%             ExpKprime_data.avgRuns.Kprime_Avg{i_win}(:),...
%             ExpKprime_data.avgRuns.Kprime_STD{i_win}(:), 'lineprops', 'r-');        
%         hold on; 
%         plot(1:NumSensors, ExpKprime_data.avgRuns.Kprime_Avg{i_win}(:), 'r-','LineWidth',2);
%         hold on;
%         
%         shadedErrorBar(1:NumSensors, ...
%             mean(ShuffKprime_data.perRepavgFold.Kprime_NullDistribution{i_win}(:,:),2)',...
%             std(ShuffKprime_data.perRepavgFold.Kprime_NullDistribution{i_win}(:,:)'), 'lineprops', 'k-');
%         hold on;
%         plot(1:NumSensors, mean(ShuffKprime_data.perRepavgFold.Kprime_NullDistribution{i_win}(:,:),2)', 'k-','LineWidth',2);
%         hold on;
%         
%         area(1:NumSensors, Kprime_data.Exp.avgRuns.mask_ExpvsShuff{i_win}(:)*20,...
%         'basevalue',0,'FaceColor',[0.1, 0.1, 0.8],'FaceAlpha', 0.5,'LineStyle','none');
%         
%         title(['TW =' num2str(i_win) '; ' num2str(sum(Kprime_data.Exp.avgRuns.mask_ExpvsShuff{i_win}(:))) ' sign. electrodes']);
%         set(gca,'FontSize',8)
%         ylabel('Kprime');
%         xlabel('Electrodes');
%         ylim([-2 20])
%         xlim([0 NumSensors+1])
%         hold on;
%         
%         if plot_counter == 1
%             legend('Exp Kprime','','Shuff Kprime','','p < 0.05','Location','northwest')            
%         end
%         
%     end
%     %save Fig
%     filename = ['ExpvsShuffKprime_AcrossFolds_' sub '_' input_DataLabel '_' tonedur_title '_' Label_TW '.png'];
%     figfile      = [path_fig_avgFolds filename];
%     saveas(gcf, [figfile], 'png'); %save png version
%     close
%    
% end

%6.3 plot Exp vs. Shuff Kprime both averaged across runs for each TW
   figure;hold on; set(gcf,'units','normalized','outerposition',[0 0 1 1])
    sgtitle(['Exp vs. Shuff (' num2str(NumReps) ' reps) Kprime across Runs (Avg+SD across Runs ' num2str(NumRuns) ') per sensor for: ' ...
        sub '; ' input_DataLabel '; ' tonedur_title '; ' Label_TW]);
    
    plot_counter = 0; 
    
    if NumWindows < 5
        Plot_rows = 1;
    else 
        Plot_rows = 2;
    end
    
    for i_win = 1:NumWindows
        
        plot_counter = plot_counter+1;
        
        if NumWindows < 5   
            subplot(NumWindows,Plot_rows,plot_counter)
        else
            subplot(NumWindows/2,Plot_rows,plot_counter)
        end
            
        shadedErrorBar(1:NumSensors, ...
            ExpKprime_data.avgRuns.Kprime_Avg{i_win}(:),...
            ExpKprime_data.avgRuns.Kprime_STD{i_win}(:), 'lineprops', 'r-');        
        hold on; 
        plot(1:NumSensors, ExpKprime_data.avgRuns.Kprime_Avg{i_win}(:), 'r-','LineWidth',2);
        hold on;
        
        shadedErrorBar(1:NumSensors, ...
            mean(ShuffKprime_data.perRepavgRuns.Kprime_NullDistribution{i_win}(:,:),2)',...
            std(ShuffKprime_data.perRepavgRuns.Kprime_NullDistribution{i_win}(:,:)'), 'lineprops', 'k-');
        hold on;
        plot(1:NumSensors, mean(ShuffKprime_data.perRepavgRuns.Kprime_NullDistribution{i_win}(:,:),2)', 'k-','LineWidth',2);
        hold on;
        
        area(1:NumSensors, Kprime_data.Exp.avgRuns.mask_ExpvsShuff{i_win}(:)*20,...
        'basevalue',0,'FaceColor',[0.1, 0.1, 0.8],'FaceAlpha', 0.5,'LineStyle','none');
        
        title(['TW =' num2str(i_win) '; ' num2str(sum(Kprime_data.Exp.avgRuns.mask_ExpvsShuff{i_win}(:))) ' sign. electrodes']);
        set(gca,'FontSize',8)
        ylabel('Kprime');
        xlabel('Electrodes');
        ylim([-2 20])
        xlim([0 NumSensors+1])
        hold on;
        
        if plot_counter == 1
            legend('Exp Kprime','','Shuff Kprime','','p < 0.05','Location','northwest')            
        end
        
    end
    %save Fig
    filename = ['ExpvsShuffKprime_AcrossRuns_' sub '_' input_DataLabel '_' tonedur_title '_' Label_TW '.png'];
    figfile      = [path_fig_avgRuns filename];
    saveas(gcf, [figfile], 'png'); %save png version
    close
   
end