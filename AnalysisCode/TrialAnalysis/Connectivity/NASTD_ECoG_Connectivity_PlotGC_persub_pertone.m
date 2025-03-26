function NASTD_ECoG_Connectivity_PlotGC_persub_pertone ...
    (subs, ...
    ToneDur_text, InputDataType, label_ElecPairSel, ...
    save_poststepFigs, ...
    paths_NASTD_ECoG)

%Plot GC results as 1) Overview (GC estimates and p-values per tone group),
%2) temporal GC estimates per tone group, 3) speatcral GC estimates per
%tone group

%% 1) Set paths &  load input data

path_fig = ([paths_NASTD_ECoG.ECoGdata_Connectivity 'GC/Figs/']);
if (~exist(path_fig, 'dir')); mkdir(path_fig); end

path_inputdata = ([paths_NASTD_ECoG.ECoGdata_Connectivity 'GC/']);
label_inputdata = ['GCdata_n' num2str(length(subs)) '_' ...
    label_ElecPairSel{1} label_ElecPairSel{2} '_' InputDataType{1}];

load([path_inputdata label_inputdata], 'GCdata');

%Determine electrode selection labels
if strcmp(label_ElecPairSel{1}, 'Pred_Pred')
    title_label1   = 'Pred2Pred';
elseif strcmp(label_ElecPairSel{1}, 'PE_PE')
    title_label1   = 'PE2PE';
elseif strcmp(label_ElecPairSel{1}, 'Pred_PE')
    title_label1   = 'Pred2PE';
elseif strcmp(label_ElecPairSel{1}, 'PE_Pred')
    title_label1   = 'PE2Pred';
end

if strcmp(label_ElecPairSel{2}, 'AllRegions')
    title_label2   = 'AllReg';
elseif strcmp(label_ElecPairSel{2}, 'Frontal_Temporal')
    title_label2   = 'Front2Temp';
elseif strcmp(label_ElecPairSel{2}, 'Temporal_Frontal')
    title_label2   = 'Temp2Front';
end

%% 2) plot GC estimate for time and frequency domain
tone_grouping = [1 11; 12 22; 23 33; 34 34];
colormap parula
colormap_parula = colormap;
colormap gray
colormap_gray = colormap;
close

for i_sub = 1:length(subs)
    if ~isempty(GCdata{i_sub})
        sub = subs{i_sub};
        %2.1 Plot Source to Target GC as function of (aggregated) tones over the course of the sequence
        for i_TD = 1:length(ToneDur_text)
            
            %Summary GC and p-values across pairs per tone group
            f1 = figure('visible','off');
            set(f1,'units','normalized','outerposition',[0 0 1 1]) %full screen
            DimSubplot = [length(tone_grouping) 2];
            sgtitle(['Time-domain pw GC estimates per tone group - ' ...
                sub ' - ' ToneDur_text{i_TD} 's TD - ' label_ElecPairSel{1} ', ' label_ElecPairSel{2} ' - ' InputDataType{1}], 'Interpreter', 'none')
            subplot_num = 0;
            
            %Determine Clim
            for i_tonegroup = 1:length(tone_grouping)
                maxC_tempGC{i_TD}(i_tonegroup,1) = max(max(mean(GCdata{i_sub}.temporalGC{i_TD}.source2target...
                    (:,:,tone_grouping(i_tonegroup,1):tone_grouping(i_tonegroup,2)),3)));
            end
            maxC_tempGC{i_TD} = max(maxC_tempGC{i_TD});
            
            for i_tonegroup = 1:length(tone_grouping)
                %GC values
                subplot_num = subplot_num + 1;
                
                f_sub = subplot(DimSubplot(1),DimSubplot(2),subplot_num);
                imagesc(mean(GCdata{i_sub}.temporalGC{i_TD}.source2target...
                    (:,:,tone_grouping(i_tonegroup,1):tone_grouping(i_tonegroup,2)),3))
                f_sub.XTick =  [1:length(GCdata{i_sub}.ind_pairedelecs_to)];
                f_sub.XTickLabel = GCdata{i_sub}.label_allelecs(GCdata{i_sub}.ind_pairedelecs_to);
                f_sub.XLabel.String = {'Target Elecs'};
                f_sub.YTick =  [1:length(GCdata{i_sub}.ind_pairedelecs_from)];
                f_sub.YTickLabel = GCdata{i_sub}.label_allelecs(GCdata{i_sub}.ind_pairedelecs_from);
                f_sub.YLabel.String = {'Source Elecs'};
                f_sub.CLim = [0 maxC_tempGC{i_TD}];
                f_sub.Colormap = colormap_parula;
                colorbar
                title(['GC estimates - Tone ' num2str(tone_grouping(i_tonegroup,1)) ' to ' num2str(tone_grouping(i_tonegroup,2))]);
                %p-values
                subplot_num = subplot_num + 1;
                
                f_sub = subplot(DimSubplot(1),DimSubplot(2),subplot_num);
                imagesc(median(GCdata{i_sub}.pval_temporalGC{i_TD}.source2target...
                    (:,:,tone_grouping(i_tonegroup,1):tone_grouping(i_tonegroup,2)),3))
                f_sub.XTick =  [1:length(GCdata{i_sub}.ind_pairedelecs_to)];
                f_sub.XTickLabel = GCdata{i_sub}.label_allelecs(GCdata{i_sub}.ind_pairedelecs_to);
                f_sub.XLabel.String = {'Target Elecs'};
                f_sub.YTick =  [1:length(GCdata{i_sub}.ind_pairedelecs_from)];
                f_sub.YTickLabel = GCdata{i_sub}.label_allelecs(GCdata{i_sub}.ind_pairedelecs_from);
                f_sub.YLabel.String = {'Source Elecs'};
                f_sub.CLim = [0 0.05];
                f_sub.Colormap = colormap_gray;
                colorbar
                title(['GC p-values (vs. null model, uncorrected) - Tone ' num2str(tone_grouping(i_tonegroup,1)) ' to ' num2str(tone_grouping(i_tonegroup,2))]);
            end
            
            if istrue(save_poststepFigs)
                filename     = ['Heatmap_' sub '_'...
                    'TemppwGC_TD' num2str(i_TD) '_' ...
                    title_label1 '_' title_label2 '_'  ...
                    InputDataType{1} '.png'];
                path_fig2 = [path_fig 'persub/' sub '/'];
                if (~exist(path_fig2, 'dir')); mkdir(path_fig2); end
                figfile      = [path_fig2  filename];
                saveas(gcf, figfile, 'png'); %save png version
                close;
            end
            
            %Temporal GC per connection per tone group
            DimSubplot = [length(GCdata{i_sub}.ind_pairedelecs_from) length(GCdata{i_sub}.ind_pairedelecs_to)];
            i_subplot = 0;
            f1 = figure('visible','off');
            set(f1,'units','normalized','outerposition',[0 0 1 1]) %full screen
            %Determine common y-axis limits
            for i_tonegroup = 1:length(tone_grouping)
                maxY_tempGC{i_TD}(i_tonegroup,1) = max(max(mean(GCdata{i_sub}.temporalGC{i_TD}.source2target...
                    (:,:,tone_grouping(i_tonegroup,1):tone_grouping(i_tonegroup,2)),3)));
            end
            maxY_tempGC{i_TD} = max(maxY_tempGC{i_TD});
            
            for i_sourceelec = 1:length(GCdata{i_sub}.ind_pairedelecs_from)
                for i_targetelec = 1:length(GCdata{i_sub}.ind_pairedelecs_to)
                    i_subplot = i_subplot + 1;
                    f_sub = subplot(DimSubplot(1),DimSubplot(2),i_subplot);
                    bar([mean(GCdata{i_sub}.temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,1:11),3), ...
                        mean(GCdata{i_sub}.temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,12:22),3),  ...
                        mean(GCdata{i_sub}.temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,23:33),3), ...
                        GCdata{i_sub}.temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,34)])
                    hold on
                    errorbar([mean(GCdata{i_sub}.temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,1:11),3), ...
                        mean(GCdata{i_sub}.temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,12:22),3),  ...
                        mean(GCdata{i_sub}.temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,23:33),3), ...
                        GCdata{i_sub}.temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,34)],...
                        [std(GCdata{i_sub}.temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,1:11)), ...
                        std(GCdata{i_sub}.temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,12:22)), ...
                        std(GCdata{i_sub}.temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,23:33)), 0], ...
                        'LineStyle', 'none', 'Color', 'k');
                    
                    f_sub.XTickLabel = {'1-11', '12-22', '23-33', '34'};
                    f_sub.XLabel.String = {'Tones'};
                    f_sub.YLabel.String = {'GC estimate'};
                    f_sub.YLim = [0 maxY_tempGC{i_TD}*1.25];
                    
                    title(['From ' GCdata{i_sub}.label_allelecs{GCdata{i_sub}.ind_pairedelecs_from(i_sourceelec)} ' to ' GCdata{i_sub}.label_allelecs{GCdata{i_sub}.ind_pairedelecs_to(i_targetelec)}])
                end
            end
            sgtitle(['Time-domain pw GC estimates per tone group - ' ...
                sub ' - ' ToneDur_text{i_TD} 's TD - ' label_ElecPairSel{1} ', ' label_ElecPairSel{2} ' - ' InputDataType{1}], 'Interpreter', 'none')
            
            if istrue(save_poststepFigs)
                filename     = ['ErrBar_' sub '_'...
                    'TemppwGC_' 'TD' num2str(i_TD) '_' ...
                    title_label1 '_' title_label2 '_'  ...
                    InputDataType{1} '.png'];
                figfile      = [path_fig2  filename];
                saveas(gcf, figfile, 'png'); %save png version
                close;
            end
            
            
            %Spectral GC per tone group
            DimSubplot = [length(GCdata{i_sub}.ind_pairedelecs_from) length(GCdata{i_sub}.ind_pairedelecs_to)];
            i_subplot = 0;
            f2 = figure('visible','off');
            set(f2,'units','normalized','outerposition',[0 0 1 1]) %full screen
            colormap cool
            colormap_cool = colormap;
            %Determine plot settings based on input freq
            if strcmp(InputDataType{1}, 'Broadband')
                temp_xlim = [1 150];
                temp_xtick = [0:25:150];
                temp_xticklabel = {'0', '25', '50', '75', '100', '125', '150'};
            elseif strcmp(InputDataType{1}, 'HP05toLP30Hz')
                temp_xlim = [0 30];
                temp_xtick = [0:5:30];
                temp_xticklabel = {'0', '5', '10', '15', '20', '25', '30'};
            elseif strcmp(InputDataType{1}, 'HighGamma_LogAmp')
                temp_xlim = [70 150];
                temp_xtick = [70:10:150];
                temp_xticklabel = {'70', '80', '90', '100', '110', '120', '130', '140', '150'};
            end
            
            %Determine common y-axis limits
            for i_tonegroup = 1:length(tone_grouping)
                maxY_spectGC{i_TD}(i_tonegroup,1) = ...
                    max(max(max(mean(...
                    GCdata{i_sub}.spectralGC{i_TD}.source2target...
                    (:,:,temp_xlim(1):temp_xlim(2),tone_grouping(i_tonegroup,1):tone_grouping(i_tonegroup,2)),4))));
            end
            maxY_spectGC{i_TD} = max(maxY_spectGC{i_TD});
            
            for i_sourceelec = 1:length(GCdata{i_sub}.ind_pairedelecs_from)
                for i_targetelec = 1:length(GCdata{i_sub}.ind_pairedelecs_to)
                    i_subplot = i_subplot + 1;
                    f_sub = subplot(DimSubplot(1),DimSubplot(2),i_subplot);
                    plot(squeeze(mean(GCdata{i_sub}.spectralGC{i_TD}.source2target(i_sourceelec,i_targetelec,:,1:11),4)), 'Color', colormap_cool(50,:), 'LineWidth', 1);
                    hold on;
                    plot(squeeze(mean(GCdata{i_sub}.spectralGC{i_TD}.source2target(i_sourceelec,i_targetelec,:,12:22),4)), 'Color', colormap_cool(100,:), 'LineWidth', 1);
                    plot(squeeze(mean(GCdata{i_sub}.spectralGC{i_TD}.source2target(i_sourceelec,i_targetelec,:,23:33),4)), 'Color', colormap_cool(150,:), 'LineWidth', 1);
                    plot(squeeze(mean(GCdata{i_sub}.spectralGC{i_TD}.source2target(i_sourceelec,i_targetelec,:,34),4)), 'Color', colormap_cool(200,:), 'LineWidth', 1);
                    
                    f_sub.XLim = temp_xlim;
                    f_sub.XTick = temp_xtick;
                    f_sub.XTickLabel = temp_xticklabel;
                    
                    f_sub.XLabel.String = {'Frequency [Hz]'};
                    f_sub.YLabel.String = {'GC'};
                    f_sub.YLim = [0 maxY_spectGC{i_TD}];
                    
                    if i_sourceelec ==1 && i_targetelec == 1
                        legend({'Tone 1-11', '12-22', '23-33', '34'})
                    end
                    
                    title(['From ' GCdata{i_sub}.label_allelecs{GCdata{i_sub}.ind_pairedelecs_from(i_sourceelec)} ' to ' GCdata{i_sub}.label_allelecs{GCdata{i_sub}.ind_pairedelecs_to(i_targetelec)}])
                end
            end
            sgtitle(['Spectral pw GC estimates per tone group - ' ...
                sub ' - ' ToneDur_text{i_TD} 's TD - ' label_ElecPairSel{1} ', ' label_ElecPairSel{2} ' - ' InputDataType{1}], 'Interpreter', 'none')
            
            if istrue(save_poststepFigs)
                filename     = ['Line_' sub '_'...
                    'SpectpwGC_' 'TD' num2str(i_TD) '_' ...
                    title_label1 '_' title_label2 '_'  ...
                    InputDataType{1} '.png'];
                figfile      = [path_fig2  filename];
                saveas(gcf, figfile, 'png'); %save png version
                close;
            end
        
        %5.2 Plot reversed (Target to Source) GC as function of (aggregated) tones over the course of the sequence
            
            %Summary GC and p-values across pairs per tone group
            f1 = figure('visible','off');
            set(f1,'units','normalized','outerposition',[0 0 1 1]) %full screen
            DimSubplot = [length(tone_grouping) 2];
            sgtitle(['Reversed (target-to-source) Time-domain pw GC estimates per tone group - ' ...
                sub ' - ' ToneDur_text{i_TD} 's TD - ' label_ElecPairSel{1} ', ' label_ElecPairSel{2} ' - ' InputDataType{1}], 'Interpreter', 'none')
            subplot_num = 0;
            
            for i_tonegroup = 1:length(tone_grouping)
                %GC values
                subplot_num = subplot_num + 1;
                
                f_sub = subplot(DimSubplot(1),DimSubplot(2),subplot_num);
                imagesc(mean(GCdata{i_sub}.temporalGC{i_TD}.target2source...
                    (:,:,tone_grouping(i_tonegroup,1):tone_grouping(i_tonegroup,2)),3))
                f_sub.XTick =  [1:length(GCdata{i_sub}.ind_pairedelecs_to)];
                f_sub.XTickLabel = GCdata{i_sub}.label_allelecs(GCdata{i_sub}.ind_pairedelecs_to);
                f_sub.XLabel.String = {'Source (revsered Target) Elecs'};
                f_sub.YTick =  [1:length(GCdata{i_sub}.ind_pairedelecs_from)];
                f_sub.YTickLabel = GCdata{i_sub}.label_allelecs(GCdata{i_sub}.ind_pairedelecs_from);
                f_sub.YLabel.String = {'Target (reversed Source) Elecs'};
                f_sub.CLim = [0 maxC_tempGC{i_TD}];
                f_sub.Colormap = colormap_parula;
                colorbar
                title(['GC estimates - Tone ' num2str(tone_grouping(i_tonegroup,1)) ' to ' num2str(tone_grouping(i_tonegroup,2))]);
                
                %p-values
                subplot_num = subplot_num + 1;
                
                f_sub = subplot(DimSubplot(1),DimSubplot(2),subplot_num);
                imagesc(median(GCdata{i_sub}.pval_temporalGC{i_TD}.target2source...
                    (:,:,tone_grouping(i_tonegroup,1):tone_grouping(i_tonegroup,2)),3))
                f_sub.XTick =  [1:length(GCdata{i_sub}.ind_pairedelecs_to)];
                f_sub.XTickLabel = GCdata{i_sub}.label_allelecs(GCdata{i_sub}.ind_pairedelecs_to);
                f_sub.XLabel.String = {'Source (revsered Target) Elecs'};
                f_sub.YTick =  [1:length(GCdata{i_sub}.ind_pairedelecs_from)];
                f_sub.YTickLabel = GCdata{i_sub}.label_allelecs(GCdata{i_sub}.ind_pairedelecs_from);
                f_sub.YLabel.String = {'Target (reversed Source) Elecs'};
                f_sub.CLim = [0 0.05];
                f_sub.Colormap = colormap_gray;
                colorbar
                title(['GC p-values (vs. null model, uncorrected) - Tone ' num2str(tone_grouping(i_tonegroup,1)) ' to ' num2str(tone_grouping(i_tonegroup,2))]);
            end
            
            
            if istrue(save_poststepFigs)
                filename     = ['Heatmap_' sub '_'...
                    'TemppwGC_TD' num2str(i_TD) '_reversed' ...
                    title_label1 '_' title_label2 '_'  ...
                    InputDataType{1} '.png'];
                if (~exist(path_fig2, 'dir')); mkdir(path_fig2); end
                figfile      = [path_fig2  filename];
                saveas(gcf, figfile, 'png'); %save png version
                close;
            end
            
            %Temporal GC per connection per tone group
            DimSubplot = [length(GCdata{i_sub}.ind_pairedelecs_from) length(GCdata{i_sub}.ind_pairedelecs_to)];
            i_subplot = 0;
            f1 = figure('visible','off');
            set(f1,'units','normalized','outerposition',[0 0 1 1]) %full screen
            
            for i_sourceelec = 1:length(GCdata{i_sub}.ind_pairedelecs_from)
                for i_targetelec = 1:length(GCdata{i_sub}.ind_pairedelecs_to)
                    i_subplot = i_subplot + 1;
                    f_sub = subplot(DimSubplot(1),DimSubplot(2),i_subplot);
                    bar([mean(GCdata{i_sub}.temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,1:11),3), ...
                        mean(GCdata{i_sub}.temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,12:22),3),  ...
                        mean(GCdata{i_sub}.temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,23:33),3), ...
                        GCdata{i_sub}.temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,34)])
                    hold on
                    errorbar([mean(GCdata{i_sub}.temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,1:11),3), ...
                        mean(GCdata{i_sub}.temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,12:22),3),  ...
                        mean(GCdata{i_sub}.temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,23:33),3), ...
                        GCdata{i_sub}.temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,34)],...
                        [std(GCdata{i_sub}.temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,1:11)), ...
                        std(GCdata{i_sub}.temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,12:22)), ...
                        std(GCdata{i_sub}.temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,23:33)), 0], ...
                        'LineStyle', 'none', 'Color', 'k');
                    
                    f_sub.XTickLabel = {'1-11', '12-22', '23-33', '34'};
                    f_sub.XLabel.String = {'Tones'};
                    f_sub.YLabel.String = {'GC estimate'};
                    f_sub.YLim = [0 maxY_tempGC{i_TD}*1.25];
                   
                    title(['From ' GCdata{i_sub}.label_allelecs{GCdata{i_sub}.ind_pairedelecs_to(i_targetelec)} ' to '  GCdata{i_sub}.label_allelecs{GCdata{i_sub}.ind_pairedelecs_from(i_sourceelec)}])
                end
            end
            sgtitle(['Reversed Time-domain pw GC estimates per tone group - ' ...
                sub ' - ' ToneDur_text{i_TD} 's TD - Reversed' label_ElecPairSel{1} ', ' label_ElecPairSel{2} ' - ' InputDataType{1}], 'Interpreter', 'none')
            
            if istrue(save_poststepFigs)
                filename     = ['ErrBar_' sub '_'...
                    'TemppwGC_' 'TD' num2str(i_TD) '_reversed' ...
                    title_label1 '_' title_label2 '_'  ...
                    InputDataType{1} '.png'];
                figfile      = [path_fig2  filename];
                saveas(gcf, figfile, 'png'); %save png version
                close;
            end
            
            
            %Spectral GC per tone group
            DimSubplot = [length(GCdata{i_sub}.ind_pairedelecs_from) length(GCdata{i_sub}.ind_pairedelecs_to)];
            i_subplot = 0;
            f2 = figure('visible','off');
            set(f2,'units','normalized','outerposition',[0 0 1 1]) %full screen
            colormap cool
            colormap_cool = colormap;
            for i_sourceelec = 1:length(GCdata{i_sub}.ind_pairedelecs_from)
                for i_targetelec = 1:length(GCdata{i_sub}.ind_pairedelecs_to)
                    i_subplot = i_subplot + 1;
                    f_sub = subplot(DimSubplot(1),DimSubplot(2),i_subplot);
                    plot(squeeze(mean(GCdata{i_sub}.spectralGC{i_TD}.target2source(i_sourceelec,i_targetelec,:,1:11),4)), 'Color', colormap_cool(50,:), 'LineWidth', 1);
                    hold on;
                    plot(squeeze(mean(GCdata{i_sub}.spectralGC{i_TD}.target2source(i_sourceelec,i_targetelec,:,12:22),4)), 'Color', colormap_cool(100,:), 'LineWidth', 1);
                    plot(squeeze(mean(GCdata{i_sub}.spectralGC{i_TD}.target2source(i_sourceelec,i_targetelec,:,23:33),4)), 'Color', colormap_cool(150,:), 'LineWidth', 1);
                    plot(squeeze(mean(GCdata{i_sub}.spectralGC{i_TD}.target2source(i_sourceelec,i_targetelec,:,34),4)), 'Color', colormap_cool(200,:), 'LineWidth', 1);
                    
                    if strcmp(InputDataType{1}, 'Broadband')
                        f_sub.XLim = [1 150];
                        f_sub.XTick = [0:25:150];
                        f_sub.XTickLabel = {'0', '25', '50', '75', '100', '125', '150'};
                    elseif strcmp(InputDataType{1}, 'HP05toLP30Hz')
                        f_sub.XLim = [0 30];
                        f_sub.XTick = [0:5:30];
                        f_sub.XTickLabel = {'0', '5', '10', '15', '20', '25', '30'};
                    elseif strcmp(InputDataType{1}, 'HighGamma_LogAmp')
                        f_sub.XLim = [70 150];
                        f_sub.XTick = [70:10:150];
                        f_sub.XTickLabel = {'70', '80', '90', '100', '110', '120', '130', '140', '150'};
                    end
                    f_sub.XLabel.String = {'Frequency [Hz]'};
                    f_sub.YLabel.String = {'GC'};
                    f_sub.YLim = [0 maxY_spectGC{i_TD}];
                    
                    if i_sourceelec ==1 && i_targetelec == 1
                        legend({'Tone 1-11', '12-22', '23-33', '34'})
                    end
                    
                    title(['From ' GCdata{i_sub}.label_allelecs{GCdata{i_sub}.ind_pairedelecs_to(i_targetelec)} ' to '  GCdata{i_sub}.label_allelecs{GCdata{i_sub}.ind_pairedelecs_from(i_sourceelec)}])
                end
            end
            sgtitle(['Reversed spectral pw GC estimates per tone group - ' ...
                sub ' - ' ToneDur_text{i_TD} 's TD - Reversed' label_ElecPairSel{1} ', ' label_ElecPairSel{2} ' - ' InputDataType{1}], 'Interpreter', 'none')
            
            if istrue(save_poststepFigs)
                filename     = ['Line_' sub '_'...
                    'SpectpwGC_' 'TD' num2str(i_TD) '_reversed' ...
                    title_label1 '_' title_label2 '_'  ...
                    InputDataType{1} '.png'];
                figfile      = [path_fig2  filename];
                saveas(gcf, figfile, 'png'); %save png version
                close;
            end
            
        end
    end
end

end