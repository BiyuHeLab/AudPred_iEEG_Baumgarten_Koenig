function SignElecIndices = NASTD_ECoG_Predict_PredEffectTable_lk...
    (validSubjs, sub_list, ...
    FuncInput_DataType, FuncInput_EffectType, FuncInput_ToneDur_text,  ...
    param,...
    plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG)

%Aims:
%Perform FDR correction for electrode-wise prediction effect p-values.
%Create a summary table listing p-values across subjects.
%Plot bar chart comparing number of sign. electrodes across input types and effects.

SignElecIndices = [];

for i_sub = validSubjs
    for i_inputData = 1:length(FuncInput_DataType)
        
        %% 1. Time-windowed data
        %1.1 Load in TW-wise prediction effects and electrode labels
        sub = sub_list{i_sub};
        path_inputdata = [paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
            sub '/Data/'];
        
        for i_TD = 1:length(FuncInput_ToneDur_text)
            
            %TW data
            load([path_inputdata param.Label_TW 'sTW/' sub '_PredEffects_' FuncInput_DataType{i_inputData} '_' FuncInput_ToneDur_text{i_TD} 'sTD.mat'], ...
                'PredEffect','SimplePredErrEffect','ComplexPredErrEffect', 'labels_loadedData');
            
            %1.2 2FDR-corrent across electrodes (per TW)
            for i_effect = 1:length(FuncInput_EffectType)
                for i_win = 1:length(eval([FuncInput_EffectType{i_effect} '.pval{1}']))
                    
                    pval = eval([FuncInput_EffectType{i_effect} '.pval{1}{i_win}']);
                    [~,~,~,pval_FDR] = fdr_bh(pval, param.alpha_FDR, 'pdep','no');
                    %                      pval_FDR = mafdr(pval,'BHFDR', true); %requires bioinformatics toolbox (available for Matlab2017a)
                    
                    tempSignElecIndices = [];
                    tempSignElecIndices.uncorr = find(pval < param.alpha);
                    tempSignElecIndices.FDRcorr = find(pval_FDR < param.alpha_FDR);
                    
                    if ~isempty(tempSignElecIndices.FDRcorr)
                        display([num2str(length(tempSignElecIndices.FDRcorr)) ...
                            ' electrodes showing FDR-corrected sign. TW-effects for sub ' sub_list{i_sub} ...
                            '; InputData: ' FuncInput_DataType{i_inputData} ...
                            '; Effect: ' FuncInput_EffectType{i_effect} ...
                            '; ToneDur: ' FuncInput_ToneDur_text{i_TD}])
                    end
                    
                    %1.3 Store p-vals in common var
                    SignElecIndices.TW.uncorr{i_sub}.(...
                        FuncInput_DataType{i_inputData}).(...
                        FuncInput_EffectType{i_effect}).Index{i_TD}{i_win} = ...
                        tempSignElecIndices.uncorr;
                    SignElecIndices.TW.FDRcorr{i_sub}.(...
                        FuncInput_DataType{i_inputData}).(...
                        FuncInput_EffectType{i_effect}).Index{i_TD}{i_win} = ...
                        tempSignElecIndices.FDRcorr;
                    
                    %1.4 Accumulate in list across subjects {InputData, Effect, ToneDuration}
                    %Create list if not entry not present
                    if ~isfield(SignElecIndices.TW,'List')
                        SignElecIndices.TW.List = [];
                        SignElecIndices.TW.List.FDRcorr = ...
                            cell(length(FuncInput_DataType),...
                            length(FuncInput_EffectType),...
                            length(FuncInput_ToneDur_text));
                        SignElecIndices.TW.List.uncorr = ...
                            cell(length(FuncInput_DataType),...
                            length(FuncInput_EffectType),...
                            length(FuncInput_ToneDur_text));
                    end
                    %Enter FDRcorrceted results in list
                    for i_elec = 1:length(tempSignElecIndices.FDRcorr)
                        if ~isempty(tempSignElecIndices.FDRcorr(i_elec))
                            SignElecIndices.TW.List.FDRcorr{i_inputData,i_effect,i_TD} = ...
                                [SignElecIndices.TW.List.FDRcorr{i_inputData,i_effect,i_TD}; ...
                                sub FuncInput_DataType{i_inputData} FuncInput_EffectType{i_effect} FuncInput_ToneDur_text{i_TD}...
                                'Nan' ...
                                i_win i_win ...
                                tempSignElecIndices.FDRcorr(i_elec) ...
                                labels_loadedData(tempSignElecIndices.FDRcorr(i_elec))...
                                pval_FDR(tempSignElecIndices.FDRcorr(i_elec)) ...
                                pval(tempSignElecIndices.FDRcorr(i_elec))];
                        end
                    end
                    %Enter uncorrected results in list
                    for i_elec = 1:length(tempSignElecIndices.uncorr)
                        if ~isempty(tempSignElecIndices.uncorr(i_elec))
                            SignElecIndices.TW.List.uncorr{i_inputData,i_effect,i_TD} = ...
                                [SignElecIndices.TW.List.uncorr{i_inputData,i_effect,i_TD}; ...
                                sub FuncInput_DataType{i_inputData} FuncInput_EffectType{i_effect} FuncInput_ToneDur_text{i_TD} ...
                                'NaN' ...
                                i_win i_win ...
                                tempSignElecIndices.uncorr(i_elec) ...
                                labels_loadedData(tempSignElecIndices.uncorr(i_elec))...
                                pval_FDR(tempSignElecIndices.uncorr(i_elec)) ...
                                pval(tempSignElecIndices.uncorr(i_elec))];
                        end
                    end
                    %Clean up
                    pval = [];
                end
            end
            %Clean up
            clear PredEffect SimplePredErrEffect ComplexPredErrEffect ...
                labels_loadedData
            
            
            %% 2. Sample-wise data
            %2.1 Load in sample-wise prediction effects and electrode labels
            %Sample-wise data
            load([path_inputdata 'Samplewise/' sub '_PredEffectsCluster_' FuncInput_DataType{i_inputData} '_' FuncInput_ToneDur_text{i_TD} 'sTD_NoClusterCorr.mat'], ...
                'PredEffect', 'SimplePredErrEffect', 'ComplexPredErrEffect', 'labels_loadedData');
            
            %2.2 FDR-correct across electrodes
            %Read out FDR-corrected across time p-values for each electrode
            for i_effect = 1:length(FuncInput_EffectType) 
                currEffect = eval(FuncInput_EffectType{i_effect});
                numElecs = length(currEffect.pval_FDR);

                % Gather minimum FDR p-value per electrode (across time)
                minPvals = nan(1, numElecs);
                for i_elec = 1:numElecs
                    if ~isempty(currEffect.pval_FDR{i_elec})
                        minPvals(i_elec) = min(currEffect.pval_FDR{i_elec});
                    end
                end
                
                % Track uncorrected significant p vals
                tempSignElecIndices = [];
                tempSignElecIndices.uncorr = find(minPvals(1,:) < param.alpha);
                
                if ~isempty(tempSignElecIndices.uncorr)
                    display([num2str(length(tempSignElecIndices.uncorr)) ...
                        ' electrodes showing sign. Clusters for sub ' sub_list{i_sub} ...
                        '; InputData: ' FuncInput_DataType{i_inputData} ...
                        '; Effect: ' FuncInput_EffectType{i_effect} ...
                        '; ToneDur: ' FuncInput_ToneDur_text{i_TD}])
                end
                
                %2.3 Store in common var
                SignElecIndices.Samplewise.uncorr{i_sub}.(...
                    FuncInput_DataType{i_inputData}).(...
                    FuncInput_EffectType{i_effect}).Index{i_TD} = ...
                    tempSignElecIndices.uncorr;
                
                %2.4 Accumulate in list across subjects {InputData, Effect, ToneDuration}
                if ~isfield(SignElecIndices.Samplewise,'List')
                    SignElecIndices.Samplewise.List = [];
                    SignElecIndices.Samplewise.List.uncorr = ...
                        cell(length(FuncInput_DataType),...
                        length(FuncInput_EffectType),...
                        length(FuncInput_ToneDur_text));
                end
                
                % FDR-correct across electrodes
                fdrAdjP = mafdr(minPvals, 'BHFDR', true);  % <- Benjamini-Hochberg
                sigElecs = find(fdrAdjP < param.alpha);
                
                % Initialize FDR list field if it doesn't exist
                if ~isfield(SignElecIndices.Samplewise.List, 'FDR')
                    SignElecIndices.Samplewise.List.FDR = ...
                        cell(length(FuncInput_DataType), ...
                        length(FuncInput_EffectType), ...
                        length(FuncInput_ToneDur_text));
                end

                for idx = 1:length(sigElecs)
                    elecID = sigElecs(idx);
                    SignElecIndices.Samplewise.List.FDR{i_inputData, i_effect, i_TD} = ...
                        [SignElecIndices.Samplewise.List.FDR{i_inputData, i_effect, i_TD}; ...
                        {sub, FuncInput_DataType{i_inputData}, FuncInput_EffectType{i_effect}, FuncInput_ToneDur_text{i_TD}, ...
                         elecID, labels_loadedData{elecID}, fdrAdjP(elecID)}];
                end
                
                %Clean
                minPvals = [];
            end
            %Clean
            clear PredEffect SimplePredErrEffect ComplexPredErrEffect ...
                labels_loadedData
        end
    end
end

%% 3. Transfer array to table
colNames1 = {'sub', 'InputData', 'EffectType', 'ToneDur', ...
    'ClusterNum', 'TW/StartSample', 'TW/StopSample', 'ElecIndex', 'ElecLabel','p_FDR', 'p_uncorr'};
colNames2 = {'sub', 'InputData', 'EffectType', 'ToneDur', ...
    'ElecIndex', 'ElecLabel','p_FDR', 'p_uncorr'};

for i_inputData = 1:length(FuncInput_DataType)
    for i_effect = 1:length(FuncInput_EffectType)
        i_TD = i_inputData;
        if ~isempty(SignElecIndices.TW.List.FDRcorr{i_inputData,i_effect,i_TD})
            SignElecIndices.TW.List.FDRcorr{i_inputData,i_effect,i_TD} = ...
                array2table(SignElecIndices.TW.List.FDRcorr{i_inputData,i_effect,i_TD},...
                'VariableNames',colNames1);
        end
        if ~isempty(SignElecIndices.TW.List.uncorr{i_inputData,i_effect,i_TD})
            SignElecIndices.TW.List.uncorr{i_inputData,i_effect,i_TD} = ...
                array2table(SignElecIndices.TW.List.uncorr{i_inputData,i_effect,i_TD},...
                'VariableNames',colNames1);
        end
        if ~isempty(SignElecIndices.Samplewise.List.uncorr{i_inputData,i_effect,i_TD})
            SignElecIndices.Samplewise.List.uncorr{i_inputData,i_effect,i_TD} = ...
                array2table(SignElecIndices.Samplewise.List.uncorr{i_inputData,i_effect,i_TD},...
                'VariableNames',colNames2);
        end
        if ~isempty(SignElecIndices.Samplewise.List.FDR{i_inputData,i_effect,i_TD})
            SignElecIndices.Samplewise.List.FDR{i_inputData,i_effect,i_TD} = ...
                array2table(SignElecIndices.Samplewise.List.FDR{i_inputData,i_effect,i_TD},...
                'VariableNames',colNames2);
        end
    end
end

%Save table
path_dataoutput = [paths_NASTD_ECoG.ECoGdata_Prediction ...
            '/PredEffects/Allsub_n' num2str(length(validSubjs))...
            '/Data/' param.Label_TW 'sTW/'];
if (~exist(path_dataoutput, 'dir')); mkdir(path_dataoutput); end

savefile = [path_dataoutput 'N' num2str(length(validSubjs)) '_SummaryPredEffects_FTPL.mat'];
save(savefile, 'SignElecIndices');

%% 4. Plot barchart summarizing sign. electrodes per effect & input data type
if plot_poststepFigs == 1
    
    for i_effect = 1:size(SignElecIndices.TW.List.FDRcorr,2)
        for i_inputdata = 1:size(SignElecIndices.TW.List.FDRcorr,1)
            
            SignElecs_TW_FDR(i_effect, i_inputdata) = 0;
            SignElecs_TW_uncorr(i_effect, i_inputdata) = 0;
            SignElecs_Cluster_uncorr(i_effect, i_inputdata) = 0; 
            
            for i_TD = 1:size(SignElecIndices.TW.List.FDRcorr,3)
                if ~isempty(SignElecIndices.TW.List.FDRcorr{i_inputdata,i_effect,i_TD})
                    SignElecs_TW_FDR(i_effect, i_inputdata) = ...
                        SignElecs_TW_FDR(i_effect, i_inputdata) + ...
                        length(SignElecIndices.TW.List.FDRcorr{i_inputdata,i_effect,i_TD}.sub);
                end
                SignElecs_TW_uncorr(i_effect, i_inputdata) = ...
                    SignElecs_TW_uncorr(i_effect, i_inputdata) + ...
                    length(SignElecIndices.TW.List.uncorr{i_inputdata,i_effect,i_TD}.sub);
                SignElecs_Cluster_uncorr(i_effect, i_inputdata) = ...
                    SignElecs_Cluster_uncorr(i_effect, i_inputdata) + ...
                    length(SignElecIndices.Cluster.List.uncorr{i_inputdata,i_effect,i_TD}.sub);
            end
        end
    end
    maxNumSignElecs_TW_FDR = max(max(SignElecs_TW_FDR));
    maxNumSignElecs_TW_uncorr = max(max(SignElecs_TW_uncorr));
    maxNumSignElecs_Cluster_uncorr = max(max(SignElecs_Cluster_uncorr));
    
    h = figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    %Prediction effects

    subplot(3,3,1)
    bar(1:length(SignElecs_TW_uncorr), SignElecs_TW_uncorr(1,:));
    xlabel('Signal Component','FontSize',12)
    ylabel('# sign. TW','FontSize',12)
    XTickLabel = {'LP35Hz', '\alpha', '\beta', '\gamma', 'High \gamma'};
    set(gca,'XTickLabel', XTickLabel,'FontSize',12)
    ylim([0 round(maxNumSignElecs_TW_uncorr*1.1)])
    title({['Prediction Effects'], ['sign. TW (p<' num2str(param.alpha) ' , no MCC)']},'FontSize',12);
    text(1:length(SignElecs_TW_uncorr),SignElecs_TW_uncorr(1,:),num2str(SignElecs_TW_uncorr(1,:)'),'vert','bottom','horiz','center');
    subplot(3,3,2)
    bar(1:length(SignElecs_TW_FDR), SignElecs_TW_FDR(1,:));
    xlabel('Signal Component','FontSize',12)
    ylabel('# sign. TW','FontSize',12)
    XTickLabel = {'LP35Hz', '\alpha', '\beta', '\gamma', 'High \gamma'};
    set(gca,'XTickLabel', XTickLabel,'FontSize',12)
%     ylim([0 round(maxNumSignElecs_TW_FDR*1.1)])
    ylim([0 round(maxNumSignElecs_TW_uncorr*1.1)])
    title({['Prediction Effects'], ['sign. TW (p<' num2str(param.alpha) ' , MCC (FDR) elecs)']},'FontSize',12);
    text(1:length(SignElecs_TW_FDR),SignElecs_TW_FDR(1,:),num2str(SignElecs_TW_FDR(1,:)'),'vert','bottom','horiz','center');
    subplot(3,3,3)
    bar(1:length(SignElecs_Cluster_uncorr), SignElecs_Cluster_uncorr(1,:));
    xlabel('Signal Component','FontSize',12)
    ylabel('# sign. clusters','FontSize',12)
    XTickLabel = {'LP35Hz', '\alpha', '\beta', '\gamma', 'High \gamma'};
    set(gca,'XTickLabel', XTickLabel,'FontSize',12)
    ylim([0 round(maxNumSignElecs_TW_uncorr*1.1)])
    title({['Prediction Effects'], ['sign. clusters (p<' num2str(param.alpha) ', MCC (cluster) samples)']},'FontSize',12);
    text(1:length(SignElecs_Cluster_uncorr),SignElecs_Cluster_uncorr(1,:),num2str(SignElecs_Cluster_uncorr(1,:)'),'vert','bottom','horiz','center');
    
    %Simple Prediction Error effects
    subplot(3,3,4)
    bar(1:length(SignElecs_TW_uncorr), SignElecs_TW_uncorr(2,:));
    xlabel('Signal Component','FontSize',12)
    ylabel('# sign. TW','FontSize',12)
    XTickLabel = {'LP35Hz', '\alpha', '\beta', '\gamma', 'High \gamma'};
    set(gca,'XTickLabel', XTickLabel,'FontSize',12)
    ylim([0 round(maxNumSignElecs_TW_uncorr*1.1)])
    title({['Simple Pred Error Effects'], ['sign. TW (p<' num2str(param.alpha) ', no MCC)']},'FontSize',12);
    text(1:length(SignElecs_TW_uncorr),SignElecs_TW_uncorr(2,:),num2str(SignElecs_TW_uncorr(2,:)'),'vert','bottom','horiz','center');
    subplot(3,3,5)
    bar(1:length(SignElecs_TW_FDR), SignElecs_TW_FDR(2,:));
    xlabel('Signal Component','FontSize',12)
    ylabel('# sign. TW','FontSize',12)
    XTickLabel = {'LP35Hz', '\alpha', '\beta', '\gamma', 'High \gamma'};
    set(gca,'XTickLabel', XTickLabel,'FontSize',12)
    ylim([0 round(maxNumSignElecs_TW_uncorr*1.1)])
%     ylim([0 round(maxNumSignElecs_TW_FDR*1.1)])
    title({['Simple Pred Error Effects'], ['sign. TW (p<' num2str(param.alpha) ', MCC (FDR) elecs)']},'FontSize',12);
    text(1:length(SignElecs_TW_FDR),SignElecs_TW_FDR(2,:),num2str(SignElecs_TW_FDR(2,:)'),'vert','bottom','horiz','center');
    subplot(3,3,6)
    bar(1:length(SignElecs_Cluster_uncorr), SignElecs_Cluster_uncorr(2,:));
    xlabel('Signal Component','FontSize',12)
    ylabel('# sign. clusters','FontSize',12)
    XTickLabel = {'LP35Hz', '\alpha', '\beta', '\gamma', 'High \gamma'};
    set(gca,'XTickLabel', XTickLabel,'FontSize',12)
    ylim([0 round(maxNumSignElecs_TW_uncorr*1.1)])
    title({['Simple Pred Error Effects'], ['sign. cluster (p<' num2str(param.alpha) ', MCC (cluster) samples)']},'FontSize',12);
    text(1:length(SignElecs_Cluster_uncorr),SignElecs_Cluster_uncorr(2,:),num2str(SignElecs_Cluster_uncorr(2,:)'),'vert','bottom','horiz','center');
    
    %Complex Prediction Error effects
    subplot(3,3,7)
    bar(1:length(SignElecs_TW_uncorr), SignElecs_TW_uncorr(3,:));
    xlabel('Signal Component','FontSize',12)
    ylabel('# sign. TW','FontSize',12)
    XTickLabel = {'LP35Hz', '\alpha', '\beta', '\gamma', 'High \gamma'};
    set(gca,'XTickLabel', XTickLabel,'FontSize',12)
    ylim([0 round(maxNumSignElecs_TW_uncorr*1.1)])
    title({['Complex Pred Error  Effects'], ['sign. TW-effects (p<' num2str(param.alpha) ', no MCC)']},'FontSize',12);
    text(1:length(SignElecs_TW_uncorr),SignElecs_TW_uncorr(3,:),num2str(SignElecs_TW_uncorr(3,:)'),'vert','bottom','horiz','center');
    subplot(3,3,8)
    bar(1:length(SignElecs_TW_FDR), SignElecs_TW_FDR(3,:));
    xlabel('Signal Component','FontSize',12)
    ylabel('# sign. TW','FontSize',12)
    XTickLabel = {'LP35Hz', '\alpha', '\beta', '\gamma', 'High \gamma'};
    set(gca,'XTickLabel', XTickLabel,'FontSize',12)
    ylim([0 round(maxNumSignElecs_TW_uncorr*1.1)])
%     ylim([0 round(maxNumSignElecs_TW_FDR*1.1)])
    title({['Complex Pred Error Effects'], ['sign. TW-effects (p<' num2str(param.alpha) ', MCC (FDR) elecs)']},'FontSize',12);
    text(1:length(SignElecs_TW_FDR),SignElecs_TW_FDR(3,:),num2str(SignElecs_TW_FDR(3,:)'),'vert','bottom','horiz','center');
    subplot(3,3,9)
    bar(1:length(SignElecs_Cluster_uncorr), SignElecs_Cluster_uncorr(2,:));
    xlabel('Signal Component','FontSize',12)
    ylabel('# sign. clusters','FontSize',12)
    XTickLabel = {'LP35Hz', '\alpha', '\beta', '\gamma', 'High \gamma'};
    set(gca,'XTickLabel', XTickLabel,'FontSize',12)
    ylim([0 round(maxNumSignElecs_TW_uncorr*1.1)])
    title({['Complex Pred Error Effects'], ['sign. clusters (p<' num2str(param.alpha) ', MCC (cluster) samples)']},'FontSize',12);
    text(1:length(SignElecs_Cluster_uncorr),SignElecs_Cluster_uncorr(3,:),num2str(SignElecs_Cluster_uncorr(3,:)'),'vert','bottom','horiz','center');
    
    sgtitle({['Number of sign. effects (TW or clusters) per effect type & signal component'], ...
        ['(N = ' num2str(length(validSubjs)) ', effects pooled across both TD)']});
    
    %Save summary figure
    if save_poststepFigs == 1
        path_fig = ([paths_NASTD_ECoG.ECoGdata_Prediction ...
            '/PredEffects/Allsub_n' num2str(length(validSubjs))...
            '/Figs/']);
        if (~exist(path_fig, 'dir')); mkdir(path_fig); end
        
        filename     = ['Allsub_n' num2str(length(validSubjs)) ...
            '_SummarySignElecs_perEffect.png'];
        figfile      = [path_fig filename];
        saveas(gcf, figfile, 'png'); %save png version
        close all;
    end
end