function NASTD_ECoG_Connectivity_PlotGC_persub_aggrtone ...
    (subs, ...
    ToneDur_text, ...
    save_poststepFigs, ...
    param, paths_NASTD_ECoG)

%Plot GC results as 1) Overview (GC estimates and p-values per tone group),
%2) temporal GC estimates per tone group, 3) speatcral GC estimates per
%tone group

%% 1) Set paths &  load input data

path_fig = ([paths_NASTD_ECoG.ECoGdata_Connectivity 'GC/Figs/']);
if (~exist(path_fig, 'dir')); mkdir(path_fig); end

path_inputdata = ([paths_NASTD_ECoG.ECoGdata_Connectivity 'GC/']);
label_inputdata = ['GCdata_n' num2str(length(subs)) '_' ...
    param.GC.label_ElecPairSel{1} param.GC.label_ElecPairSel{2} '_' param.GC.InputDataType{1}];

load([path_inputdata label_inputdata], 'GCdata');

%Determine electrode selection labels
if strcmp(param.GC.label_ElecPairSel{1}, 'Pred_Pred')
    title_label1   = 'Pred2Pred';
elseif strcmp(param.GC.label_ElecPairSel{1}, 'PE_PE')
    title_label1   = 'PE2PE';
elseif strcmp(param.GC.label_ElecPairSel{1}, 'Pred_PE')
    title_label1   = 'Pred2PE';
elseif strcmp(param.GC.label_ElecPairSel{1}, 'PE_Pred')
    title_label1   = 'PE2Pred';
end

if strcmp(param.GC.label_ElecPairSel{2}, 'AllRegions')
    title_label2   = 'AllReg';
elseif strcmp(param.GC.label_ElecPairSel{2}, 'Frontal_Temporal')
    title_label2   = 'Front2Temp';
elseif strcmp(param.GC.label_ElecPairSel{2}, 'Temporal_Frontal')
    title_label2   = 'Temp2Front';
end

%% 2) plot GC estimate for time and frequency domain
colormap parula
colormap_parula = colormap;
colormap gray
colormap_gray = colormap;
close

set_visibility = 'off';

for i_sub = 1:length(subs)
    if ~isempty(GCdata{i_sub})
        
        sub = subs{i_sub};        
        path_fig2 = [path_fig 'persub/' sub '/' title_label1 '_' title_label2 '/'];
        if (~exist(path_fig2, 'dir')); mkdir(path_fig2); end
                
        %2.1 Plot Source to Target GC as function of (aggregated) tones over the course of the sequence
        for i_TD = 1:length(ToneDur_text)
            
            %Determine sign. matrix
            sign_matrix_alltonegroups   = [];
            sign_matrix_pertonegroup    = [];
            for i_tonegroup = 1:length(param.GC.tone_aggregation)
                sign_matrix_pertonegroup(:,:,i_tonegroup) = ...
                    GCdata{i_sub}.pval_temporalGC{i_TD}.source2target(:,:,i_tonegroup) < 0.05;
                sign_matrix_pertonegroup(isnan(GCdata{i_sub}.pval_temporalGC{i_TD}.source2target(:,:,i_tonegroup))) = 0;
            end
            sign_matrix_alltonegroups = sum(sign_matrix_pertonegroup,3) > 0;
            
            %Summary GC and p-values across pairs per tone group
            f1 = figure('visible',set_visibility);
            set(f1,'units','normalized','outerposition',[0 0 1 1]) %full screen
            DimSubplot = [length(param.GC.tone_aggregation) 2];
            sgtitle(['Time-domain pw GC estimates per tone group - ' ...
                sub ' - ' ToneDur_text{i_TD} 's TD - ' param.GC.label_ElecPairSel{1} ', ' param.GC.label_ElecPairSel{2} ' - ' param.GC.InputDataType{1}], 'Interpreter', 'none')
            subplot_num = 0;
            
            %Determine Clim
            for i_tonegroup = 1:length(param.GC.tone_aggregation)
                maxC_tempGC{i_TD}(i_tonegroup,1) = max(max(mean(...
                    GCdata{i_sub}.temporalGC{i_TD}.source2target...
                    (:,:,i_tonegroup),3)));
            end
            maxC_tempGC{i_TD} = max(maxC_tempGC{i_TD});
            
            for i_tonegroup = 1:length(param.GC.tone_aggregation)
                %GC values
                subplot_num = subplot_num + 1;
                
                f_sub = subplot(DimSubplot(1),DimSubplot(2),subplot_num);
                imagesc(mean(GCdata{i_sub}.temporalGC{i_TD}.source2target...
                    (:,:,i_tonegroup),3))
                f_sub.XTick =  [1:length(GCdata{i_sub}.ind_pairedelecs_to)];
                f_sub.XTickLabel = GCdata{i_sub}.label_allelecs(GCdata{i_sub}.ind_pairedelecs_to);
                f_sub.XLabel.String = {'Target Elecs'};
                f_sub.YTick =  [1:length(GCdata{i_sub}.ind_pairedelecs_from)];
                f_sub.YTickLabel = GCdata{i_sub}.label_allelecs(GCdata{i_sub}.ind_pairedelecs_from);
                f_sub.YLabel.String = {'Source Elecs'};
                f_sub.CLim = [0 maxC_tempGC{i_TD}];
                f_sub.Colormap = colormap_parula;
                colorbar
                title(['GC estimates - Tone ' ...
                    num2str(param.GC.tone_aggregation{i_tonegroup}(1)) ' to ' num2str(param.GC.tone_aggregation{i_tonegroup}(end))]);
                
                %p-values
                subplot_num = subplot_num + 1;
                
                f_sub = subplot(DimSubplot(1),DimSubplot(2),subplot_num);
                imagesc(GCdata{i_sub}.pval_temporalGC{i_TD}.source2target...
                    (:,:,i_tonegroup))
                f_sub.XTick =  [1:length(GCdata{i_sub}.ind_pairedelecs_to)];
                f_sub.XTickLabel = GCdata{i_sub}.label_allelecs(GCdata{i_sub}.ind_pairedelecs_to);
                f_sub.XLabel.String = {'Target Elecs'};
                f_sub.YTick =  [1:length(GCdata{i_sub}.ind_pairedelecs_from)];
                f_sub.YTickLabel = GCdata{i_sub}.label_allelecs(GCdata{i_sub}.ind_pairedelecs_from);
                f_sub.YLabel.String = {'Source Elecs'};
                f_sub.CLim = [0 0.05];
                f_sub.Colormap = colormap_gray;
                colorbar
                title(['GC p-values (vs. null model, uncorrected) - Tone ' ...
                    num2str(param.GC.tone_aggregation{i_tonegroup}(1)) ' to ' num2str(param.GC.tone_aggregation{i_tonegroup}(end))]);
            end
            
            if istrue(save_poststepFigs)
                filename     = ['Heatmap_' sub '_'...
                    'TemppwGC_TD' num2str(i_TD) '_' ...
                    title_label1 '_' title_label2 '_'  ...
                    param.GC.InputDataType{1} '.png'];
                figfile      = [path_fig2  filename];
                saveas(gcf, figfile, 'png'); %save png version
                close;
            end
            
            %Temporal GC per connection per tone group, 
            %only for connections that are sign. in at least 1 tone epoch            
            DimSubplot = [length(GCdata{i_sub}.ind_pairedelecs_from) length(GCdata{i_sub}.ind_pairedelecs_to)];
            i_subplot = 0;
            f1 = figure('visible',set_visibility);
            set(f1,'units','normalized','outerposition',[0 0 1 1]) %full screen
            %Determine common y-axis limits
            for i_tonegroup = 1:length(param.GC.tone_aggregation)
                maxY_tempGC{i_TD}(i_tonegroup,1) = max(max(GCdata{i_sub}.temporalGC{i_TD}.source2target...
                    (:,:,i_tonegroup)));
            end
            maxY_tempGC{i_TD} = max(maxY_tempGC{i_TD});
            
            for i_sourceelec = 1:length(GCdata{i_sub}.ind_pairedelecs_from)
                for i_targetelec = 1:length(GCdata{i_sub}.ind_pairedelecs_to)
                    
                    i_subplot = i_subplot + 1;
                    f_sub = subplot(DimSubplot(1),DimSubplot(2),i_subplot);
                    
                    if sign_matrix_alltonegroups(i_sourceelec, i_targetelec) > 0
                        
                        meandata_bar = [];
                        for i_tonegroup = 1:length(param.GC.tone_aggregation)
                            meandata_bar(1, i_tonegroup) = ...
                                GCdata{i_sub}.temporalGC{i_TD}.source2target(...
                                i_sourceelec,i_targetelec,i_tonegroup);
                            label_Xtick{i_tonegroup} = ...
                                [num2str(param.GC.tone_aggregation{i_tonegroup}(1)) '-' ...
                                num2str(param.GC.tone_aggregation{i_tonegroup}(end))];
                        end
                        
                        b = bar(meandata_bar, 'facecolor', 'flat');
                        %Color nonsign. tone epcohs differently
                        for i_tonegroup = 1:length(param.GC.tone_aggregation)
                            if sign_matrix_pertonegroup(i_sourceelec, i_targetelec, i_tonegroup) == 0
                                b.CData(i_tonegroup,:) = [0.8 0.8 0.8];
                            else
                                b.CData(i_tonegroup,:) = [0 0 1];
                            end
                        end
                        f_sub.XTickLabel = label_Xtick;
                        f_sub.XTickLabelRotation = 45;
                        f_sub.XLabel.String = {'Tones'};
                        f_sub.YLabel.String = {'GC estimate'};
%                         f_sub.YLim = [0 maxY_tempGC{i_TD}*1.1];
                    else 
                        axis off
                    end
                    title([GCdata{i_sub}.label_allelecs{GCdata{i_sub}.ind_pairedelecs_from(i_sourceelec)} ' -> ' ...
                        GCdata{i_sub}.label_allelecs{GCdata{i_sub}.ind_pairedelecs_to(i_targetelec)}])
                end
            end
            sgtitle(['Time-domain pw GC estimates per tone group - ' ...
                sub ' - ' ToneDur_text{i_TD} 's TD - ' param.GC.label_ElecPairSel{1} ', ' ...
                param.GC.label_ElecPairSel{2} ' - ' param.GC.InputDataType{1}], 'Interpreter', 'none')
            
            if istrue(save_poststepFigs)
                filename     = ['ErrBar_' sub '_'...
                    'TemppwGC_' 'TD' num2str(i_TD) '_' ...
                    title_label1 '_' title_label2 '_'  ...
                    param.GC.InputDataType{1} '.png'];
                figfile      = [path_fig2  filename];
                saveas(gcf, figfile, 'png'); %save png version
                close;
            end
            
            
            %Spectral GC per tone group
            DimSubplot = [length(GCdata{i_sub}.ind_pairedelecs_from) length(GCdata{i_sub}.ind_pairedelecs_to)];
            i_subplot = 0;
            f2 = figure('visible',set_visibility);
            set(f2,'units','normalized','outerposition',[0 0 1 1]) %full screen
            colormap cool
            colormap_cool = colormap;
            %Determine plot settings based on input freq
            if strcmp(param.GC.InputDataType{1}, 'Broadband')
                temp_xlim = [1 150];
                temp_xtick = [0:25:150];
                temp_xticklabel = {'0', '25', '50', '75', '100', '125', '150'};
            elseif strcmp(param.GC.InputDataType{1}, 'HP05toLP30Hz')
                temp_xlim = [1 30];
                temp_xtick = [0:5:30];
                temp_xticklabel = {'0', '5', '10', '15', '20', '25', '30'};
            elseif strcmp(param.GC.InputDataType{1}, 'HighGamma_LogAmp')
                temp_xlim = [70 150];
                temp_xtick = [70:10:150];
                temp_xticklabel = {'70', '80', '90', '100', '110', '120', '130', '140', '150'};
            end
            
            %Determine common y-axis limits
            for i_tonegroup = 1:length(param.GC.tone_aggregation)
                maxY_spectGC{i_TD}(i_tonegroup,1) = ...
                    max(max(max(mean(...
                    GCdata{i_sub}.spectralGC{i_TD}.source2target...
                    (:,:,temp_xlim(1):temp_xlim(2),i_tonegroup),4))));
            end
            maxY_spectGC{i_TD} = max(maxY_spectGC{i_TD});
            
            for i_sourceelec = 1:length(GCdata{i_sub}.ind_pairedelecs_from)
                for i_targetelec = 1:length(GCdata{i_sub}.ind_pairedelecs_to)
                    
                    i_subplot = i_subplot + 1;
                    f_sub = subplot(DimSubplot(1),DimSubplot(2),i_subplot);
                    
                    if sign_matrix_alltonegroups(i_sourceelec, i_targetelec) > 0                        
                        meandata_line = [];
                        for i_tonegroup = 1:length(param.GC.tone_aggregation)
                            meandata_line(i_tonegroup,:) = ...
                                GCdata{i_sub}.spectralGC{i_TD}.source2target(...
                                i_sourceelec,i_targetelec,:,i_tonegroup);
                        end
                        
                        hold on;
                        for i_tonegroup = 1:length(param.GC.tone_aggregation)
                            if sign_matrix_pertonegroup(i_sourceelec, i_targetelec, i_tonegroup) == 0
                                plot(meandata_line(i_tonegroup,:), ...
                                    'Color', colormap_cool(50*i_tonegroup,:), 'LineWidth', 2, 'LineStyle', ':');                                
                            else                            
                                plot(meandata_line(i_tonegroup,:), ...
                                    'Color', colormap_cool(50*i_tonegroup,:), 'LineWidth', 2, 'LineStyle', '-');
                            end
                        end
                        
                        f_sub.XLim = temp_xlim;
                        f_sub.XTick = temp_xtick;
                        f_sub.XTickLabel = temp_xticklabel;
                        
                        f_sub.XLabel.String = {'Frequency [Hz]'};
                        f_sub.YLabel.String = {'GC'};
%                         f_sub.YLim = [0 maxY_spectGC{i_TD}];
                        
                        if i_sourceelec ==1 && i_targetelec == 1
                            legend(label_Xtick)
                        end
                    else 
                        axis off
                    end
                    
                    title([GCdata{i_sub}.label_allelecs{GCdata{i_sub}.ind_pairedelecs_from(i_sourceelec)} ...
                        ' -> ' GCdata{i_sub}.label_allelecs{GCdata{i_sub}.ind_pairedelecs_to(i_targetelec)}])
                    
                end
            end
            sgtitle(['Spectral pw GC estimates per tone group - ' ...
                sub ' - ' ToneDur_text{i_TD} 's TD - ' param.GC.label_ElecPairSel{1} ', ' ...
                param.GC.label_ElecPairSel{2} ' - ' param.GC.InputDataType{1}], 'Interpreter', 'none')
            
            if istrue(save_poststepFigs)
                filename     = ['Line_' sub '_'...
                    'SpectpwGC_' 'TD' num2str(i_TD) '_' ...
                    title_label1 '_' title_label2 '_'  ...
                    param.GC.InputDataType{1} '.png'];
                figfile      = [path_fig2  filename];
                saveas(gcf, figfile, 'png'); %save png version
                close;
            end
            
            
            %5.2 Plot reversed (Target to Source) GC as function of (aggregated) tones over the course of the sequence
            
            %Determine sign. matrix
            sign_matrix_alltonegroups   = [];
            sign_matrix_pertonegroup    = [];
            for i_tonegroup = 1:length(param.GC.tone_aggregation)
                sign_matrix_pertonegroup(:,:,i_tonegroup) = ...
                    GCdata{i_sub}.pval_temporalGC{i_TD}.target2source(:,:,i_tonegroup) < 0.05;
                sign_matrix_pertonegroup(isnan(GCdata{i_sub}.pval_temporalGC{i_TD}.source2target(:,:,i_tonegroup))) = 0;
            end
            sign_matrix_alltonegroups = sum(sign_matrix_pertonegroup,3) > 0;  
            
            %Summary GC and p-values across pairs per tone group
            f1 = figure('visible',set_visibility);
            set(f1,'units','normalized','outerposition',[0 0 1 1]) %full screen
            DimSubplot = [length(param.GC.tone_aggregation) 2];
            sgtitle(['Reversed (target-to-source) Time-domain pw GC estimates per tone group - ' ...
                sub ' - ' ToneDur_text{i_TD} 's TD - ' param.GC.label_ElecPairSel{1} ', ' param.GC.label_ElecPairSel{2} ' - ' param.GC.InputDataType{1}], 'Interpreter', 'none')
            subplot_num = 0;
            
            %Determine Clim
            for i_tonegroup = 1:length(param.GC.tone_aggregation)
                maxC_tempGC{i_TD}(i_tonegroup,1) = max(max(mean(...
                    GCdata{i_sub}.temporalGC{i_TD}.target2source...
                    (:,:,i_tonegroup),3)));
            end
            maxC_tempGC{i_TD} = max(maxC_tempGC{i_TD});
            
            for i_tonegroup = 1:length(param.GC.tone_aggregation)
                %GC values
                subplot_num = subplot_num + 1;
                
                f_sub = subplot(DimSubplot(1),DimSubplot(2),subplot_num);
                imagesc(mean(GCdata{i_sub}.temporalGC{i_TD}.target2source...
                    (:,:,i_tonegroup),3))
                f_sub.XTick =  [1:length(GCdata{i_sub}.ind_pairedelecs_to)];
                f_sub.XTickLabel = GCdata{i_sub}.label_allelecs(GCdata{i_sub}.ind_pairedelecs_to);
                f_sub.XLabel.String = {'Source (revsered Target) Elecs'};
                f_sub.YTick =  [1:length(GCdata{i_sub}.ind_pairedelecs_from)];
                f_sub.YTickLabel = GCdata{i_sub}.label_allelecs(GCdata{i_sub}.ind_pairedelecs_from);
                f_sub.YLabel.String = {'Target (reversed Source) Elecs'};
                f_sub.CLim = [0 maxC_tempGC{i_TD}];
                f_sub.Colormap = colormap_parula;
                colorbar
                title(['GC estimates - Tone ' ...
                    num2str(param.GC.tone_aggregation{i_tonegroup}(1)) ' to ' num2str(param.GC.tone_aggregation{i_tonegroup}(end))]);
                
                %p-values
                subplot_num = subplot_num + 1;
                
                f_sub = subplot(DimSubplot(1),DimSubplot(2),subplot_num);
                imagesc(GCdata{i_sub}.pval_temporalGC{i_TD}.target2source...
                    (:,:,i_tonegroup))
                f_sub.XTick =  [1:length(GCdata{i_sub}.ind_pairedelecs_to)];
                f_sub.XTickLabel = GCdata{i_sub}.label_allelecs(GCdata{i_sub}.ind_pairedelecs_to);
                f_sub.XLabel.String = {'Source (revsered Target) Elecs'};
                f_sub.YTick =  [1:length(GCdata{i_sub}.ind_pairedelecs_from)];
                f_sub.YTickLabel = GCdata{i_sub}.label_allelecs(GCdata{i_sub}.ind_pairedelecs_from);
                f_sub.YLabel.String = {'Target (reversed Source) Elecs'};
                f_sub.CLim = [0 0.05];
                f_sub.Colormap = colormap_gray;
                colorbar
                title(['GC p-values (vs. null model, uncorrected) - Tone ' ...
                    num2str(param.GC.tone_aggregation{i_tonegroup}(1)) ' to ' num2str(param.GC.tone_aggregation{i_tonegroup}(end))]);
            end
            
            if istrue(save_poststepFigs)
                filename     = ['Heatmap_' sub '_'...
                    'TemppwGC_TD' num2str(i_TD) '_reversed' ...
                    title_label1 '_' title_label2 '_'  ...
                    param.GC.InputDataType{1} '.png'];
                if (~exist(path_fig2, 'dir')); mkdir(path_fig2); end
                figfile      = [path_fig2  filename];
                saveas(gcf, figfile, 'png'); %save png version
                close;
            end
            
            %Temporal GC per connection per tone group
            DimSubplot = [length(GCdata{i_sub}.ind_pairedelecs_from) length(GCdata{i_sub}.ind_pairedelecs_to)];
            i_subplot = 0;
            f1 = figure('visible',set_visibility);
            set(f1,'units','normalized','outerposition',[0 0 1 1]) %full screen
            
            %Determine common y-axis limits
            for i_tonegroup = 1:length(param.GC.tone_aggregation)
                maxY_tempGC{i_TD}(i_tonegroup,1) = max(max(GCdata{i_sub}.temporalGC{i_TD}.target2source...
                    (:,:,i_tonegroup)));
            end
            maxY_tempGC{i_TD} = max(maxY_tempGC{i_TD});
            
            for i_sourceelec = 1:length(GCdata{i_sub}.ind_pairedelecs_from)
                for i_targetelec = 1:length(GCdata{i_sub}.ind_pairedelecs_to)
                    
                    i_subplot = i_subplot + 1;
                    f_sub = subplot(DimSubplot(1),DimSubplot(2),i_subplot);
                    
                    if sign_matrix_alltonegroups(i_sourceelec, i_targetelec) > 0
                        meandata_bar = [];
                        for i_tonegroup = 1:length(param.GC.tone_aggregation)
                            meandata_bar(1, i_tonegroup) = ...
                                GCdata{i_sub}.temporalGC{i_TD}.target2source(...
                                i_sourceelec,i_targetelec,i_tonegroup);
                            label_Xtick{i_tonegroup} = ...
                                [num2str(param.GC.tone_aggregation{i_tonegroup}(1)) '-' ...
                                num2str(param.GC.tone_aggregation{i_tonegroup}(end))];
                        end
                        
                        b = bar(meandata_bar, 'facecolor', 'flat');
                        %Color nonsign. tone epcohs differently
                        for i_tonegroup = 1:length(param.GC.tone_aggregation)
                            if sign_matrix_pertonegroup(i_sourceelec, i_targetelec, i_tonegroup) == 0
                                b.CData(i_tonegroup,:) = [0.8 0.8 0.8];
                            else
                                b.CData(i_tonegroup,:) = [0 0 1];
                            end
                        end                        
                        f_sub.XTickLabel = label_Xtick;
                        f_sub.XTickLabelRotation = 45;
                        f_sub.XLabel.String = {'Tones'};
                        f_sub.YLabel.String = {'GC estimate'};
%                         f_sub.YLim = [0 maxY_tempGC{i_TD}*1.1];
                    else 
                        axis off
                    end
                    
                    title([GCdata{i_sub}.label_allelecs{GCdata{i_sub}.ind_pairedelecs_to(i_targetelec)} ...
                        ' -> ' GCdata{i_sub}.label_allelecs{GCdata{i_sub}.ind_pairedelecs_from(i_sourceelec )}])
                end
            end
            sgtitle(['Reversed Time-domain pw GC estimates per tone group - ' ...
                sub ' - ' ToneDur_text{i_TD} 's TD - Reversed' param.GC.label_ElecPairSel{1} ', ' ...
                param.GC.label_ElecPairSel{2} ' - ' param.GC.InputDataType{1}], 'Interpreter', 'none')
            
            if istrue(save_poststepFigs)
                filename     = ['ErrBar_' sub '_'...
                    'TemppwGC_' 'TD' num2str(i_TD) '_reversed' ...
                    title_label1 '_' title_label2 '_'  ...
                    param.GC.InputDataType{1} '.png'];
                figfile      = [path_fig2  filename];
                saveas(gcf, figfile, 'png'); %save png version
                close;
            end
            
            
            %Spectral GC per tone group
            DimSubplot = [length(GCdata{i_sub}.ind_pairedelecs_from) length(GCdata{i_sub}.ind_pairedelecs_to)];
            i_subplot = 0;
            f2 = figure('visible',set_visibility);
            set(f2,'units','normalized','outerposition',[0 0 1 1]) %full screen
            colormap cool
            colormap_cool = colormap;
            
            %Determine common y-axis limits
            for i_tonegroup = 1:length(param.GC.tone_aggregation)
                maxY_spectGC{i_TD}(i_tonegroup,1) = ...
                    max(max(max(mean(...
                    GCdata{i_sub}.spectralGC{i_TD}.target2source...
                    (:,:,temp_xlim(1):temp_xlim(2),i_tonegroup),4))));
            end
            maxY_spectGC{i_TD} = max(maxY_spectGC{i_TD});
            
            for i_sourceelec = 1:length(GCdata{i_sub}.ind_pairedelecs_from)
                for i_targetelec = 1:length(GCdata{i_sub}.ind_pairedelecs_to)
                    
                    i_subplot = i_subplot + 1;
                    f_sub = subplot(DimSubplot(1),DimSubplot(2),i_subplot);
                    
                    if sign_matrix_alltonegroups(i_sourceelec, i_targetelec) > 0
                        
                        meandata_line = [];
                        for i_tonegroup = 1:length(param.GC.tone_aggregation)
                            meandata_line(i_tonegroup,:) = ...
                                GCdata{i_sub}.spectralGC{i_TD}.target2source(...
                                i_sourceelec,i_targetelec,:,i_tonegroup);
                        end
                        
                        hold on;
                        for i_tonegroup = 1:length(param.GC.tone_aggregation)
                            if sign_matrix_pertonegroup(i_sourceelec, i_targetelec, i_tonegroup) == 0
                                plot(meandata_line(i_tonegroup,:), ...
                                    'Color', colormap_cool(50*i_tonegroup,:), 'LineWidth', 2, 'LineStyle', ':');                                
                            else                            
                                plot(meandata_line(i_tonegroup,:), ...
                                    'Color', colormap_cool(50*i_tonegroup,:), 'LineWidth', 2, 'LineStyle', '-');
                            end
                        end

                        f_sub.XLim = temp_xlim;
                        f_sub.XTick = temp_xtick;
                        f_sub.XTickLabel = temp_xticklabel;
                        
                        f_sub.XLabel.String = {'Frequency [Hz]'};
                        f_sub.YLabel.String = {'GC'};
%                         f_sub.YLim = [0 maxY_spectGC{i_TD}];
                        
                        if i_sourceelec ==1 && i_targetelec == 1
                            legend(label_Xtick)
                        end
                    else 
                        axis off
                    end
                    title([GCdata{i_sub}.label_allelecs{GCdata{i_sub}.ind_pairedelecs_to(i_targetelec)} ...
                        ' -> ' GCdata{i_sub}.label_allelecs{GCdata{i_sub}.ind_pairedelecs_from(i_sourceelec)}])
                end
            end
            sgtitle(['Reversed spectral pw GC estimates per tone group - ' ...
                sub ' - ' ToneDur_text{i_TD} 's TD - Reversed' param.GC.label_ElecPairSel{1} ', ' param.GC.label_ElecPairSel{2} ' - ' param.GC.InputDataType{1}], 'Interpreter', 'none')
            
            if istrue(save_poststepFigs)
                filename     = ['Line_' sub '_'...
                    'SpectpwGC_' 'TD' num2str(i_TD) '_reversed' ...
                    title_label1 '_' title_label2 '_'  ...
                    param.GC.InputDataType{1} '.png'];
                figfile      = [path_fig2  filename];
                saveas(gcf, figfile, 'png'); %save png version
                close;
            end
            
        end
    end
end

end