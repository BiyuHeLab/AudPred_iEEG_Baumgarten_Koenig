%Aim: Compare neural data timelocked to tone presentation over the course
%of the entire sequence and contrast different pred_p34 values

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

%Load in channel indices from tonal tracking
load('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/TLA/TonalTracking_allSub.mat');
%Load in channel indices from sign. prediction effect
% load('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Prediction/PredictionElecIndices_TW10ms_AllSub.mat');
% load('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Prediction/PredictionElecIndices_TW25ms_AllSub.mat');
load('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Prediction/PredictionElecIndices_TW50ms_AllSub.mat');

inputData = {'LF','LogGammaAmp'};
ToneDur_text = {'0.2' '0.4'};
effect = {'Pred'};
pval_plotting = 0.05;
elec_counter = 0;

for i_sub = 1:length(sub_list)
    for i_inputData = 1:length(inputData)
        
        %         sub = sub_list{i_sub};
        %         path_inputdata = [paths_NASTD_ECoG.ECoGdata_Prediction sub '/Exp/'];
        %         NASTD_ECoG_subjectinfo %load subject info file (var: si)
        %         subs_PreProcSettings = NASTD_ECoG_SubPreprocSettings; %load in file with individual preproc infos
        
        temp2 = [];
        
        for i_TD = 1:length(ToneDur_text) %tone duration condition
            
            %             if strcmp(ToneDur_text{i_TD},'0.2')
            %                 tonedur_title = '200';
            %             elseif  strcmp(ToneDur_text{i_TD},'0.4')
            %                 tonedur_title = '400';
            %             end
            %
            %             load([path_inputdata sub '_Predp33_regstatsEXP_' inputData{i_inputData} '_' tonedur_title 'msTD.mat']);
            
            for i_effect = 1:length(effect)
                temp = [];
                %                 for i_win = 1:length(eval([effect{i_effect} 'Effect.pval{1}']))
                
                %                     pval = eval([effect{i_effect} 'Effect.pval{1}{i_win}']);
                %                     [~,~,~,pval_FDRcorrect] = fdr_bh(pval,0.05, 'pdep','no');
                %
                %                     FDRsignElecIndices = find(pval_FDRcorrect < pval_plotting);
                %                     if ~isempty(FDRsignElecIndices)
                %                         temp = [temp, FDRsignElecIndices]; %reads out across TW
                %                     end
                
                
                %                     for i_elec = 1:length(FDRsignElecIndices)
                %                         if ~isempty(FDRsignElecIndices(i_elec))
                %                             elec_counter = elec_counter+1;
                %                             List_FDRsignElecIndices(elec_counter,:) = [i_sub i_inputData i_effect i_TD i_win FDRsignElecIndices(i_elec)];
                %                         end
                %                     end
                for i_win = 1:length(PredictionElec_FDR{i_sub}.(inputData{i_inputData}).Pred.Index{i_TD})
                    
                    if ~isempty(PredictionElec_FDR{i_sub}.(inputData{i_inputData}).Pred.Index{i_TD}{i_win})
                        temp = [temp, PredictionElec_FDR{i_sub}.(inputData{i_inputData}).Pred.Index{i_TD}{i_win}]; %reads out across TW
                    end
                end
                temp = unique(temp);
                PredictionElec_FDR_Allwin{i_sub}.(inputData{i_inputData}).Pred.Index{i_TD} = temp;
            end
            temp2 = [temp2, PredictionElec_FDR_Allwin{i_sub}.(inputData{i_inputData}).Pred.Index{i_TD}];%reads out across TD
            
        end
        PredictionElec_FDR_AllTD{i_sub}.(inputData{i_inputData}).Index = temp2;%reads out across TD
        
    end
end
unique(List_FDRsignElecIndices(:,6))
size(unique(List_FDRsignElecIndices(:,6)))


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
%         inputData = {'LF','LogGammaAmp'};
%         plot_STD = 0;
%         save_fig = 0;
%         
%         for i_inputData = 1:length(inputData)
%             for i_TD = 1:length(ToneDur_text)
%                 
%                 for i_chan = TonalTrackingElec{i_sub}.(inputData{i_inputData}).Index{i_TD}
%                     %                 NASTD_ECoG_TLA_PlotAllTrials_BothTD_perChan(sub, inputData{i_inputData}, i_chan, preprocData_perTD, plot_STD, save_fig, paths_NASTD_ECoG)
%                     
%                     NASTD_ECoG_TLA_PlotAllTrialsperTD_perChan(sub, inputData{i_inputData}, ToneDur_text{i_TD}, i_chan, preprocData_perTD{i_TD}, plot_STD, save_fig, paths_NASTD_ECoG)
%                 end
%                 
%             end
%         end
    %% 2) Compute timelocked data per predP34 and TD
    %2.1 Demean, Baseline correct with prestimulus window & LowPass Filter LF signal (all trials per TD)
    BLwin = [-0.8 0]; %window for baseline correction (0 = start of 1st tone)
    LPfreq = 35;
    for i_TD = 1:length(ToneDur_text)
        BLCData_perTD{i_TD} = NASTD_ECoG_TLA_BLCdataperTD(BLwin, LPfreq, preprocData_perTD{i_TD});
    end
    
    %2.2) Compute timelocked data for specific pred_p34
    for i_TD = 1:length(ToneDur_text)
        TimelockData_perTDPredp34{i_TD} = NASTD_ECoG_TLA_perTDPredp34(BLCData_perTD{i_TD});
    end
%     
%     % %2.3) Plot GAvg+STD per pred_p34 for each channel
%     %     inputData = {'LF','GammaAmp','LogGammaAmp'};
%     %     plot_STD = 0;
%     %     save_fig = 0;
%     %
%     %     for i_inputData = 1:length(inputData)
%     %         for i_chan = Selected_Channels'
%     %             NASTD_ECoG_TLA_PlotPredp34_perChan(sub, inputData{i_inputData}, i_chan, TimelockData_perTDPredp34, plot_STD, save_fig, paths_NASTD_ECoG)
%     %         end
%     %     end
%     
%     
%     % %2.4) Use N1/tonal response as functional localizer
%     %     %Compute activity timelocked to tone presentation to determine electrodes
%     %     %showing increased activity after tone presentation
%     %     BL_win = [-0.5 0]; %window with
%     %     N1_win = [0.09 0.11]; %should focus on N1 at 100s, but take stimulus delay into account
%     %
%     %     plot_fig = 0;
%     %     save_fig = 0;
%     %
%     % %     inputData = {'LF','GammaAmp','LogGammaAmp'};
%     %     inputData = {'LF','LogGammaAmp'};
%     %
%     %     for i_inputData = 1:length(inputData)
%     %
%     %         if i_inputData == 1
%     %             STD_thresh = 5; %subs_PreProcSettings.(sub).N1localizer_STDthresh.LF;
%     %         else
%     %             STD_thresh = 3; %subs_PreProcSettings.(sub).N1localizer_STDthresh.GammaAmp;
%     %         end
%     %
%     %         Index_LocalizedElecs_AllChan.(inputData{i_inputData}) = ...
%     %             NASTD_ECoG_TLA_CompN1localizer(sub, inputData{i_inputData}, ...
%     %             Selected_Channels, BL_win, N1_win, STD_thresh, ...
%     %             preprocData_perTD, plot_fig, save_fig, paths_NASTD_ECoG);
%     %
%     %         for i_TD = 1:length(ToneDur_text)
%     %             if ~isempty(Index_LocalizedElecs_AllChan.(inputData{i_inputData}){i_TD}.Index)
%     %                 TonalTrackingElec{i_sub}.(inputData{i_inputData}).Index{i_TD} = Index_LocalizedElecs_AllChan.(inputData{i_inputData}){i_TD}.Index;
%     %                 TonalTrackingElec{i_sub}.(inputData{i_inputData}).Label{i_TD} = Index_LocalizedElecs_AllChan.(inputData{i_inputData}){i_TD}.Label;
%     %                 if ~isempty(TonalTrackingElec{i_sub}.(inputData{i_inputData}).Index{i_TD})
%     %                     for i_elec = 1:length(TonalTrackingElec{i_sub}.(inputData{i_inputData}).Index{i_TD})
%     %                         TonalTrackingElec{i_sub}.(inputData{i_inputData}).ChanPos{i_TD}(i_elec,:) = ...
%     %                             preprocData_perTD{i_TD}.hdr.elec.chanpos(TonalTrackingElec{i_sub}.(inputData{i_inputData}).Index{i_TD}(i_elec),:);
%     %                     end
%     %                 end
%     %             else
%     %                 TonalTrackingElec{i_sub}.(inputData{i_inputData}).Index{i_TD} = [];
%     %                 TonalTrackingElec{i_sub}.(inputData{i_inputData}).Label{i_TD} = [];
%     %                 TonalTrackingElec{i_sub}.(inputData{i_inputData}).ChanPos{i_TD} = [];
%     %             end
%     %         end
%     %
%     %
%     %     end
%     
%     
%     
    % 3) Perform statistical comparison between pred_p34 conditions (per electrode, for each sample)
    %3.1 Perform cluster-corrected stat. comparison per channel
    %     inputData = {'LF','GammaAmp','LogGammaAmp'};
    inputData = {'LF','LogGammaAmp'};
    ComparedConditions = [1 3]; %which p*34 values should be compared with t-test?
    ToneDur_text = {'0.2' '0.4'};
    
    for i_TD = 1:length(ToneDur_text)
        for i_inputData = 1:length(inputData)
            index_signElec_ttest{i_TD}{i_inputData} = [];
            index_signElec_Ftest{i_TD}{i_inputData} = [];
        end
    end
    
    plot_fig = 1;
    save_fig = 0;
    
    for i_inputData = 1:length(inputData)
        for i_TD = 1:length(ToneDur_text)
            %             for i_chan = TonalTrackingElec{i_sub}.(inputData{i_inputData}).Index{i_TD} %focusses on tonal tracking elecs
            for i_chan = 1:length(PredictionElec_FDR_Allwin{i_sub}.(inputData{i_inputData}).Pred.Index{i_TD})
                
                %                 %3.1A t-test comparing 2 conditions (FT cluster-correction)
                %                 TimelockStat_Predp34_perTD{i_sub}.ttest{i_chan}{i_TD}{i_inputData} = ...
                %                     NASTD_ECoG_TLA_StatComp_tstat(sub, i_chan, inputData{i_inputData}, ComparedConditions, ...
                %                     TimelockData_perTDPredp34{i_TD}, plot_fig, save_fig, paths_NASTD_ECoG);
                %
                %                 if TimelockStat_Predp34_perTD{i_sub}.ttest{i_chan}{i_TD}{i_inputData}.SignDiff == 1
                %                     index_signElec_ttest{i_sub}{i_TD}{i_inputData} = [index_signElec_ttest{i_TD}{i_inputData} i_chan];
                %                 end
                
                %3.1B F-test (ANOVA) comparing 3 conditions (FT cluster-correction)
                TimelockStat_Predp34_perTD{i_sub}.Ftest{i_chan}{i_TD}{i_inputData} = ...
                    NASTD_ECoG_TLA_StatComp_Fstat(sub, PredictionElec_FDR_Allwin{i_sub}.(inputData{i_inputData}).Pred.Index{i_TD}(i_chan), inputData{i_inputData}, ...
                    TimelockData_perTDPredp34{i_TD}, plot_fig, save_fig, paths_NASTD_ECoG);
                
                if TimelockStat_Predp34_perTD{i_sub}.Ftest{i_chan}{i_TD}{i_inputData}.SignDiff == 1
                    index_signElec_Ftest{i_sub}{i_TD}{i_inputData} = [index_signElec_Ftest{i_TD}{i_inputData} i_chan];
                end
            end
        end
    end
        
end
    
    %     clearvars -except paths_NASTD_ECoG sub_list vars i_sub TonalTrackingElec index_signElec_* TimelockStat_Predp34_perTD

%3.2 Plot output
%3.2.1 Plot brain volume with elecs showing tonal tracking from each subject
%Read out all coordinates for all elecs for eachs subject
for i_sub = 1:length(sub_list)
    
    sub = sub_list{i_sub};
    timelock_dir = [paths_NASTD_ECoG.ECoGdata_Timelocked '/' sub '/'];
    
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    subs_PreProcSettings = NASTD_ECoG_SubPreprocSettings; %load in file with individual preproc infos
    %Load in preprocessed neural and behavioral data
    loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];%path to indiv preprocessed ECoG data
    tic
    disp(['Loading preprocessed data set for sub: ' sub])
    load(loadfile_ECoGpreprocdata);
    disp(['done loading in ' num2str(toc) ' sec'])
    
    Selected_Channels = setdiff(preprocData_AllTrials.cfg.info_elec.selected.index4EDF, subs_PreProcSettings.(sub).rejectedChan_index);

    
    coords_sub{i_sub} = preprocData_AllTrials.hdr.elec.chanpos;
    
    for i_inputData = 1:length(inputData)
    PredictionElec_FDR_AllTD{i_sub}.(inputData{i_inputData}).Index = ...
    Selected_Channels(PredictionElec_FDR_AllTD{i_sub}.(inputData{i_inputData}).Index);

    for i_TD = 1:length(ToneDur_text)
    PredictionElec_FDR_Allwin{i_sub}.(inputData{i_inputData}).Pred.Index{i_TD} = ...
        Selected_Channels(PredictionElec_FDR_Allwin{i_sub}.(inputData{i_inputData}).Pred.Index{i_TD});    
    end
    end
    
    clear preprocData_AllTrials
end


%     inputData = {'LF','GammaAmp','LogGammaAmp'};
inputData = {'LF','LogGammaAmp'};
ToneDur_text = {'0.2' '0.4'};
cmap = 'hot';
view_angle = [270,0]; %(270,0) - LH, (90,0) = RH

for i_inputData = 1:length(inputData)
    for i_TD = 1:length(ToneDur_text)
        
        coords = [];
        vals = [];
        chanSize = [];
        
        %             %3.2.1A Only tonal sensors
        %             for i_sub = 1:length(sub_list)
        %
        %                 clims = [1,length(sub_list)]; %colormap limits
        %                 coords = [coords; TonalTrackingElec{i_sub}.(inputData{i_inputData}).ChanPos{i_TD}];
        %                 vals = [vals; ones(length(TonalTrackingElec{i_sub}.(inputData{i_inputData}).Index{i_TD}),1)*i_sub]; %used to differentiate subjects
        %                 chanSize = [chanSize; ones(length(TonalTrackingElec{i_sub}.(inputData{i_inputData}).ChanPos{i_TD}),1)*3];
        %
        %             end
        
        %3.2.1B All sensors but highlight tonal sensors
        for i_sub = 1:length(sub_list)
            
            coords = [coords; coords_sub{i_sub}];
            vals = [vals; ones(length(coords_sub{i_sub}),1)*i_sub]; %used to differentiate subjects
            chanSize_sub = ones(length(coords_sub{i_sub}),1);
%             chanSize_sub(TonalTrackingElec{i_sub}.(inputData{i_inputData}).Index{i_TD}) = 3; %highlight tonal elects
             chanSize_sub(PredictionElec_FDR_Allwin{i_sub}.(inputData{i_inputData}).Pred.Index{i_TD}) = 4; %highlight prediction sensors across TW
%             chanSize_sub(PredictionElec_FDR_AllTD{i_sub}.(inputData{i_inputData}).Index) = 4; %highlight prediction sensors across TD
            chanSize = [chanSize; chanSize_sub];
            clims = [1,length(sub_list)+1]; %colormap limits
            
        end
        
%                      NASTD_ECoG_Plot_PlotEleconSurf(coords,vals,chanSize,clims,cmap,view_angle,inputData{i_inputData}, ToneDur_text{i_TD})
         NASTD_ECoG_Plot_PlotEleconSurf_SubplotHemis(coords,vals,chanSize,clims,cmap,inputData{i_inputData}, ToneDur_text{i_TD})
        
    end
end

%3.2.1 Plot timecourse and highlight sign. differences (per electrode)
%3.2.2 Plot brain+electrode and mark those with sign. difference (in certain TW)