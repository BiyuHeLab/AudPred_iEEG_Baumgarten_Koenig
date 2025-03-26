%Aim: Compare neural data timelocked to tone presentation over the course
%of the entire sequence and contrast different pred_p34 values with
%different stat approaches

%% 0) Setup analysis
% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')

NASTD_ECoG_setVars
paths_NASTD_ECoG = NASTD_ECoG_paths(vars.location);

%Add base dir and own script dir
addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));

% 0.2) Determine subject-specific parameters
% sub_list = vars.sub_list(2:end);
sub_list = vars.sub_list;

for i_sub = 1:length(sub_list)
    
    sub = sub_list{i_sub};
    
    timelock_dir = [paths_NASTD_ECoG.ECoGdata_Timelocked '/' sub '/'];
    
    
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    subs_PreProcSettings = NASTD_ECoG_SubPreprocSettings; %load in file with individual preproc infos
    
    %Load in preprocessed neural and behavioral data
    loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];%path to indiv preprocessed ECoG data
    
    % 0.3) Load in preproc data
    tic
    disp(['Loading preprocessed data set for sub: ' sub])
    load(loadfile_ECoGpreprocdata);
    disp(['done loading in ' num2str(toc) ' sec'])
    
    % 0.4) Specify analysis info
    Selected_Channels = setdiff(preprocData_AllTrials.cfg.info_elec.selected.index4EDF, subs_PreProcSettings.(sub).rejectedChan_index);
    ToneDur_text = {'0.2' '0.4'};
    samplefreq = preprocData_AllTrials.fsample;
    
    %% 1) Prepare input data
    %1.1) Select clean/valid trials
    preprocData_CleanTrials = NASTD_ECoG_SelCleanTrials(sub,preprocData_AllTrials);
    
    %1.2) Select only trials with specific Tone Duration (TD)
    ToneDur_text = {'0.2' '0.4'};
    
    for i_TD = 1:length(ToneDur_text)
        preprocData_perTD{i_TD} = NASTD_ECoG_SelTrialsperTD(ToneDur_text{i_TD}, preprocData_CleanTrials);
        disp(['Total trial number for TD ' ToneDur_text{i_TD} 's: ' num2str(length(preprocData_perTD{i_TD}.trial))])
        
    end
    
    %1.3) Read out trial-indices for specific pred_p34
    % label_predP34 = [-1 0 1];% low, medium, high;
    %
    % for i_TD = 1:length(ToneDur_text)
    %     for i_Predp34 = 1:length(label_predP34)
    %         preprocData_perTDperPredp34{i_TD}{i_Predp34} = NASTD_ECoG_SelTrialsperPredp34(label_predP34(i_Predp34), preprocData_perTD{i_TD});
    %         disp(['Total trial number for Predp34: ' label_predP34(i_Predp34) ' in TD ' ToneDur_text{i_TD} 's: ' num2str(length(preprocData_perTDperPredp34{i_TD}{i_Predp34}.trial))])
    %     end
    %
    % end
    
    %1.4) Plot selected Predp34 trials in databrowser
    % cfg             = [];
    % cfg.ylim        = 'maxabs';
    % cfg.viewmode    = 'vertical';
    % cfg.channel     = setdiff(preprocData_perTDperPredp34{1}{1}.cfg.info_elec.selected.index4EDF, subs_PreProcSettings.(sub).rejectedChan_index);
    % ft_databrowser(cfg,preprocData_perTDperPredp34{i_TD}{i_Predp34})
    
%     %1.5) Plot selected TD-specific activity across sequence
%     inputData = {'LF','LogGammaAmp'};
%     plot_STD = 0;
%     save_fig = 0;
%     
%     for i_inputData = 1:length(inputData)
%         for i_chan = 1:length(preprocData_perTD{1}.label) %Selected_Channels'
%             
%             NASTD_ECoG_TLA_PlotTrialAvg_BothTD_perChan...
%                 (sub, inputData{i_inputData}, i_chan, ...
%                 preprocData_perTD, ...
%                 plot_STD, save_fig, paths_NASTD_ECoG)
%             
%             for i_TD = 1:length(ToneDur_text)
%                 NASTD_ECoG_TLA_PlotTrialAvg_perTD_perChan...
%                     (sub, inputData{i_inputData}, ToneDur_text{i_TD}, i_chan, ...
%                     preprocData_perTD{i_TD}, ...
%                     plot_STD, save_fig, paths_NASTD_ECoG)
%             end
%             
%         end
%     end
    
    %% 2) Compute timelocked data per predP34 and TD
    %2.1 Demean, Baseline correct with prestimulus window & LowPass Filter LF signal (all trials per TD)
    BLwin = [-0.5 0]; %window for baseline correction (0 = start of 1st tone)
    LPfreq = 35;
    for i_TD = 1:length(ToneDur_text)
        BLCData_perTD{i_TD} = ...
            NASTD_ECoG_TLA_BLCdataperTD(BLwin, LPfreq, preprocData_perTD{i_TD});
    end
    
    %2.2) Compute timelocked data for specific pred_p34
    for i_TD = 1:length(ToneDur_text)
        TimelockData_perTDPredp34{i_TD} = ...
            NASTD_ECoG_TLA_CompTLdata_perTDPredp34(BLCData_perTD{i_TD});
    end
    
    %Remove other fields to reduce data size
%     test = TimelockData_perTDPredp34{1}
%     for i = 1:3
%     test{i} = rmfield(test{i},'trial')
%     test{i} = rmfield(test{i},'Trial_LogGammaAmp')
%     test{i} = rmfield(test{i},'Avg_LogGammaAmp')
%     test{i} = rmfield(test{i},'STD_LogGammaAmp')
%     end
%     TimelockData_perTDPredp34 = test
    
%     %2.3) Plot GAvg+STD per pred_p34 for each channel
%     inputData = {'LF','LogGammaAmp'};
%     plot_STD = 0;
%     save_fig = 0;
%     
%     for i_inputData = 1:length(inputData)
%         for i_chan = Selected_Channels'
%             
%             %Both TD (subplot)
%             NASTD_ECoG_TLA_PlotPredp34_BothTD_perChan...
%                 (sub, inputData{i_inputData}, i_chan, ...
%                 TimelockData_perTDPredp34, ...
%                 plot_STD, save_fig, paths_NASTD_ECoG)
%             
%             %Per TD (single plot)
%             for i_TD = 1:length(ToneDur_text)
%                 NASTD_ECoG_TLA_PlotPredp34_perTD_perChan...
%                     (sub, inputData{i_inputData}, ToneDur_text{i_TD}, i_chan, ...
%                     TimelockData_perTDPredp34{i_TD}, ...
%                     plot_STD, save_fig, paths_NASTD_ECoG)
%             end
%             
%         end
%     end
    
    %% 3) Perform statistical comparison between pred_p34 conditions (per electrode, for each sample)
    %4.1 Perform cluster-corrected stat. comparison per channel
    inputData = {'LF','LogGammaAmp'};
    ComparedConditions = [1 3]; %which p*34 values should be compared with t-test?
    ToneDur_text = {'0.2' '0.4'};
    
    for i_TD = 1:length(ToneDur_text)
        for i_inputData = 1:length(inputData)
            Predp34DiffElec{i_sub}.Ftest{i_inputData}{i_TD}.Index = [];
%             Predp34DiffElec{i_sub}.ttest{i_inputData}{i_TD}.Index = [];
%             Predp34DiffElec{i_sub}.LinReg{i_inputData}{i_TD}.Index = [];
        end
    end
    
    plot_fig = 0;
    save_fig = 0;
    pval_Predp34comp = 0.05;
    
    for i_inputData = 1:length(inputData)
        for i_TD = 1:length(ToneDur_text)
                            
            counter_ttest = 0;               
            counter_Ftest = 0;
            counter_LinReg = 0;               
             
            for i_chan = Selected_Channels' %Selected channels
                %                  for i_chan = TonalTrackingElec{i_sub}.(inputData{i_inputData}).Index{i_TD} %tonal tracking elecs
                %                  for i_chan = PredictionElec{i_sub}.(inputData{i_inputData}).Pred.Index{i_TD}{i_win} %elecs showing prediction effect
                
                %3.1B F-test (ANOVA; indep samples) comparing 3 conditions (FT cluster-corrected across samples)
                TimelockStat_Predp34_perTD{i_sub}.Ftest{i_chan}{i_TD}{i_inputData} = ...
                    NASTD_ECoG_TLA_StatComp_Fstat...
                    (sub, i_chan, inputData{i_inputData}, pval_Predp34comp, ...
                    TimelockData_perTDPredp34{i_TD}, ...
                    plot_fig, save_fig, paths_NASTD_ECoG);
                
                if TimelockStat_Predp34_perTD{i_sub}.Ftest{i_chan}{i_TD}{i_inputData}.SignDiff == 1
                    counter_Ftest = counter_Ftest + 1;
                    Predp34DiffElec{i_sub}.Ftest{i_inputData}{i_TD}.Index = ...
                        [Predp34DiffElec{i_sub}.Ftest{i_inputData}{i_TD}.Index i_chan];
                    Predp34DiffElec{i_sub}.Ftest{i_inputData}{i_TD}.Label(counter_Ftest) = ...
                        TimelockStat_Predp34_perTD{i_sub}.Ftest{i_chan}{i_TD}{i_inputData}.label;
                    Predp34DiffElec{i_sub}.Ftest{i_inputData}{i_TD}.coords(counter_Ftest,:) = ...
                       preprocData_CleanTrials.hdr.elec.chanpos(i_chan,:);                        
                    Predp34DiffElec{i_sub}.Ftest{i_inputData}{i_TD}.Sum_SignSamples(counter_Ftest) = ...
                        sum(TimelockStat_Predp34_perTD{i_sub}.Ftest{i_chan}{i_TD}{i_inputData}.mask);
                    Predp34DiffElec{i_sub}.Ftest{i_inputData}{i_TD}.SignSamples{counter_Ftest} = ...
                        find(TimelockStat_Predp34_perTD{i_sub}.Ftest{i_chan}{i_TD}{i_inputData}.mask == 1);
                end
                
%                 %3.1A t-test (indep samples) comparing 2 conditions (FT cluster-corrected across samples)
%                 TimelockStat_Predp34_perTD{i_sub}.ttest{i_chan}{i_TD}{i_inputData} = ...
%                     NASTD_ECoG_TLA_StatComp_tstat...
%                     (sub, i_chan, inputData{i_inputData}, ComparedConditions, pval_Predp34comp, ...
%                     TimelockData_perTDPredp34{i_TD}, ...
%                     plot_fig, save_fig, paths_NASTD_ECoG);
% 
%                 if TimelockStat_Predp34_perTD{i_sub}.ttest{i_chan}{i_TD}{i_inputData}.SignDiff == 1
%                     counter_ttest = counter_ttest + 1;
%                     Predp34DiffElec{i_sub}.ttest{i_inputData}{i_TD}.Index = ...
%                         [Predp34DiffElec{i_sub}.ttest{i_inputData}{i_TD}.Index i_chan];
%                     Predp34DiffElec{i_sub}.ttest{i_inputData}{i_TD}.Label(counter_ttest) = ...
%                         TimelockStat_Predp34_perTD{i_sub}.ttest{i_chan}{i_TD}{i_inputData}.label;
%                     Predp34DiffElec{i_sub}.ttest{i_inputData}{i_TD}.coords(counter_ttest,:) = ...
%                        preprocData_CleanTrials.hdr.elec.chanpos(i_chan,:);                    
%                     Predp34DiffElec{i_sub}.ttest{i_inputData}{i_TD}.Sum_SignSamples(counter_ttest) = ...
%                         sum(TimelockStat_Predp34_perTD{i_sub}.ttest{i_chan}{i_TD}{i_inputData}.mask);
%                     Predp34DiffElec{i_sub}.ttest{i_inputData}{i_TD}.SignSamples{counter_ttest} = ...
%                         find(TimelockStat_Predp34_perTD{i_sub}.ttest{i_chan}{i_TD}{i_inputData}.mask == 1);
%                 end
%                 
% 
%                 
%                 %3.1C Linear regression (FT; indep samples) comparing 3 conditions (FT cluster-corrected across samples)
%                 TimelockStat_Predp34_perTD{i_sub}.LinReg{i_chan}{i_TD}{i_inputData} = ...
%                     NASTD_ECoG_TLA_StatComp_LinReg...
%                     (sub, i_chan, inputData{i_inputData}, pval_Predp34comp, ...
%                     TimelockData_perTDPredp34{i_TD}, ...
%                     plot_fig, save_fig, paths_NASTD_ECoG);                
%     
%                 if TimelockStat_Predp34_perTD{i_sub}.LinReg{i_chan}{i_TD}{i_inputData}.SignDiff == 1
%                     counter_LinReg = counter_LinReg + 1;
%                     Predp34DiffElec{i_sub}.LinReg{i_inputData}{i_TD}.Index = ...
%                         [Predp34DiffElec{i_sub}.LinReg{i_inputData}{i_TD}.Index i_chan];
%                     Predp34DiffElec{i_sub}.LinReg{i_inputData}{i_TD}.Label(counter_LinReg) = ...
%                         TimelockStat_Predp34_perTD{i_sub}.LinReg{i_chan}{i_TD}{i_inputData}.label;
%                     Predp34DiffElec{i_sub}.LinReg{i_inputData}{i_TD}.coords(counter_LinReg,:) = ...
%                        preprocData_CleanTrials.hdr.elec.chanpos(i_chan,:);                        
%                     Predp34DiffElec{i_sub}.LinReg{i_inputData}{i_TD}.Sum_SignSamples(counter_LinReg) = ...
%                         sum(TimelockStat_Predp34_perTD{i_sub}.LinReg{i_chan}{i_TD}{i_inputData}.mask);
%                     Predp34DiffElec{i_sub}.LinReg{i_inputData}{i_TD}.SignSamples{counter_LinReg} = ...
%                         find(TimelockStat_Predp34_perTD{i_sub}.LinReg{i_chan}{i_TD}{i_inputData}.mask == 1);
%                 end
                
            end
            
            %3.2 Plot summary plot (subplot with surface showing all
            %electrodes with sign. difference + timelocked p*34 traces with
            %sign. different samples)            
            if ~isempty(Predp34DiffElec{i_sub}.Ftest{i_inputData}{i_TD}.Index)
                inputStat = 'Fstat';
                save_fig = 1;

                NASTD_ECoG_TLA_PlotSignPredp34...
                    (sub, i_TD, inputData{i_inputData}, i_inputData, inputStat, ...
                    TimelockData_perTDPredp34{i_TD}, TimelockStat_Predp34_perTD{i_sub}.Ftest, Predp34DiffElec{i_sub}.Ftest{i_inputData}{i_TD}, ...
                    save_fig, paths_NASTD_ECoG)
            end
            
        end        
    end   
    
    clearvars -except paths_NASTD_ECoG sub_list vars i_sub Predp34DiffElec
    
end
% savefile = [paths_NASTD_ECoG.ECoGdata_Timelocked  'TLA_statPredp34_Allsub.mat'];
% save(savefile, 'Predp34DiffElec' );

%% 4) Plot output
%4.1 Plot brain volume with elecs showing sign. difference from each subject
%Read out all coordinates for all elecs for eachs subject
for i_sub = 1:length(sub_list)
    
    %4.1.1 Load in preprocessed neural and behavioral data to read out coordinate position

    sub = sub_list{i_sub};
    timelock_dir = [paths_NASTD_ECoG.ECoGdata_Timelocked '/' sub '/'];
    
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    subs_PreProcSettings = NASTD_ECoG_SubPreprocSettings; %load in file with individual preproc infos
    loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];%path to indiv preprocessed ECoG data
    tic
    disp(['Loading preprocessed data set for sub: ' sub])
    load(loadfile_ECoGpreprocdata);
    disp(['done loading in ' num2str(toc) ' sec'])    
   
    coords_sub{i_sub} = preprocData_AllTrials.hdr.elec.chanpos;
    
    clear preprocData_AllTrials
end

%Optional: Append sign elecs across TD
temp1 = [];
temp2 = [];

for i_sub = 1:length(sub_list)
    for i_inputData = 1:length(inputData)
        for i_TD = 1:length(ToneDur_text)
            temp1 = [temp1, Predp34DiffElec{i_sub}.ttest{i_inputData}{i_TD}.Index];
            temp2 = [temp2, Predp34DiffElec{i_sub}.Ftest{i_inputData}{i_TD}.Index];            
        end
        Predp34DiffElec{i_sub}.ttest{i_inputData}{3}.Index = unique(temp1);
        Predp34DiffElec{i_sub}.Ftest{i_inputData}{3}.Index = unique(temp2);
        
        temp1 = [];
        temp2 = [];
    end
end

%plot elecs based on indices
inputData = {'LF','LogGammaAmp'};
ToneDur_text = {'0.2' '0.4'};
cmap = 'hot';
view_angle = [270,0]; %(270,0) - LH, (90,0) = RH
StatMeasure = 'Ftest'; %'ttest'
save_Figs = 1;

for i_inputData = 1:length(inputData)
    for i_TD = 1:length(ToneDur_text)
        
        coords = [];
        vals = [];
        chanSize = [];
        
%         %3.2.1A Only sign sensors
%         for i_sub = 1:length(sub_list)            
%             clims = [1,length(sub_list)]; %colormap limits
%             coords = [coords; coords_sub{i_sub}(Predp34DiffElec{i_sub}.(StatMeasure){i_inputData}{i_TD}.Index,:)]; %copy only sign elecs
%             vals = [vals; ones(length(Predp34DiffElec{i_sub}.(StatMeasure){i_inputData}{i_TD}.Index,1))*i_sub]; %used to differentiate subjects
%             chanSize = [chanSize; ones(length(Predp34DiffElec{i_sub}.(StatMeasure){i_inputData}{i_TD}.Index,1))*3];            
%         end
        
        %3.2.1B All sensors but highlight tonal sensors
        for i_sub = 1:length(sub_list)            
            coords = [coords; coords_sub{i_sub}];
            vals = [vals; ones(length(coords_sub{i_sub}),1)*i_sub]; %used to differentiate subjects
            chanSize_sub = ones(length(coords_sub{i_sub}),1);
            chanSize_sub(Predp34DiffElec{i_sub}.(StatMeasure){i_inputData}{i_TD}.Index) = 3; %highlight tonal elecs
%             chanSize_sub(Predp34DiffElec{i_sub}.(StatMeasure){3}{i_inputData}.Index) = 3; %highlight tonal elecs fused across TD
            chanSize = [chanSize; chanSize_sub];
            clims = [1,length(sub_list)+1]; %colormap limits            
        end
        
        PlottedComparison = [StatMeasure ' p*34 (Cluster-corrected, p < 0.05)'];
%         NASTD_ECoG_Plot_PlotEleconSurf...
%             (coords,vals,chanSize,...
%             clims,cmap,view_angle,...
%             inputData{i_inputData}, ToneDur_text{i_TD}) %Single hemisphere, defined by view_angle
        
        NASTD_ECoG_Plot_PlotEleconSurf_SubplotHemis...
            (coords,vals,chanSize,...
            clims,cmap,...
            inputData{i_inputData}, ...
            ToneDur_text{i_TD},PlottedComparison) %Both hemispheres in subplots     
        
        if save_Figs == 1
            path_fig = [paths_NASTD_ECoG.Fig_Timelocked_GroupLvl 'Predp34_Comparison/' StatMeasure '/'];
            mkdir([path_fig]);
            filename     = ['3Dsurf_' StatMeasure 'Predp34_' inputData{i_inputData} '_' ToneDur_text{i_TD} 'msTD_Allsub.png'];
            figfile      = [path_fig filename];
            
            saveas(gcf, [figfile], 'png'); %save png version
            
            close all;
        end        

    end
end
