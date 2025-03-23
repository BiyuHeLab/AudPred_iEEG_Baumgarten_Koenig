function NASTD_ECoG_Connectivity_PlotGC_allsub_aggrtone ...
    (subs, ...
    ToneDur_text, SelElecs, AnatReg_CatLabels, ...
    save_poststepFigs, ...
    param, paths_NASTD_ECoG)

%Plot group-level aggregated GC results per anatomical region:
%1) group-average temporal GC estimates per tone group for each region
%2) group-average speatcral GC estimates per tone group

set_visibility = 'off';


%% 1) Set paths &  load input data
path_inputdata = ([paths_NASTD_ECoG.ECoGdata_Connectivity 'GC/']);
% label_inputdata = ['GCdata_n' num2str(length(subs)) '_' ...
%     param.GC.label_ElecPairSel{1} param.GC.label_ElecPairSel{2} '_' param.GC.InputDataType{1} ...
%     param.GC.tone_aggregation_label];
% label_inputdata = ['GCdata_n' num2str(length(subs)) '_' ...
%     param.GC.label_ElecPairSel{1} param.GC.label_ElecPairSel{2} '_' param.GC.InputDataType{1} '_singletones'];
label_inputdata = ['GCdata_n' num2str(length(subs)) '_' ...
    param.GC.label_ElecPairSel{1} param.GC.label_ElecPairSel{2} '_' param.GC.InputDataType{1} param.GC.tone_aggregation_label '_100mswin'];

load([path_inputdata label_inputdata], 'GCdata');

%Determine electrode selection labels
if strcmp(param.GC.label_ElecPairSel{1}, 'Pred_Pred')
    title_label1   = 'Pred2Pred';
elseif strcmp(param.GC.label_ElecPairSel{1}, 'PE_PE')
    title_label1   = 'PE2PE';
elseif strcmp(param.GC.label_ElecPairSel{1}, 'Pred_PE')
    title_label1   = 'Pred2PE';
    toneindex_noreverse = 33;
    toneindex_reverse = 34;    
elseif strcmp(param.GC.label_ElecPairSel{1}, 'PE_Pred')
    title_label1   = 'PE2Pred';
    toneindex_noreverse = 34;
    toneindex_reverse = 33;    
end

if strcmp(param.GC.label_ElecPairSel{2}, 'AllRegions')
    title_label2   = 'AllReg';
elseif strcmp(param.GC.label_ElecPairSel{2}, 'Frontal_Temporal')
    title_label2   = 'Front2Temp';
elseif strcmp(param.GC.label_ElecPairSel{2}, 'Temporal_Frontal')
    title_label2   = 'Temp2Front';
end

subfield_label1 = extractBefore(param.GC.label_ElecPairSel{1}, '_');

% path_fig = ([paths_NASTD_ECoG.ECoGdata_Connectivity 'GC/Figs/allsub/' ...
%     title_label1 '_' title_label2 '/']);
path_fig = ([paths_NASTD_ECoG.ECoGdata_Connectivity 'GC/Figs/allsub/100mswin/' ...
    title_label1 '_' title_label2 '/']);
if (~exist(path_fig, 'dir')); mkdir(path_fig); end
path_fig2 = ([path_fig 'Reversed/']);
if (~exist(path_fig2, 'dir')); mkdir(path_fig2); end


%% 2) Determine and read out anatomical regions

%2.1 Read out anatomical regions for source & target elecs per subject
clear fulllabel* anatreg* anatlabel*
for i_sub = 1:length(subs)
    if ~isempty(GCdata{i_sub})
        clear temp_*
        %Read out anatomical region labels for source elecs per subject
        for i_sourceelec = 1:length(GCdata{i_sub}.ind_pairedelecs_from)
            
            temp_eleclabel_sourceelec = GCdata{i_sub}.label_allelecs{GCdata{i_sub}.ind_pairedelecs_from(i_sourceelec)};
            temp_sublabel = subs{i_sub};
            temp_filt_elec = strncmp(fields(...
                SelElecs.Pairs.(param.GC.label_ElecPairSel{1}).(param.GC.label_ElecPairSel{2}).Pairs_perelec), ...
                [temp_eleclabel_sourceelec '_' temp_sublabel], ...
                length(temp_eleclabel_sourceelec) + length(temp_sublabel) + 1);
            temp_anatlabel_sourceelec = fields(SelElecs.Pairs.(param.GC.label_ElecPairSel{1}).(param.GC.label_ElecPairSel{2}).Pairs_perelec);
            temp_anatlabel_sourceelec = temp_anatlabel_sourceelec{temp_filt_elec};
            
            fulllabel_sourceelec{i_sub}{i_sourceelec, 1} =...
                temp_anatlabel_sourceelec; %Reduced anat label
            temp_anatlabel_sourceelec2 = ...
                cell2mat(SelElecs.([subfield_label1 'Effect']){temp_filt_elec,5});
            anatlabel_sourceelec{i_sub}{i_sourceelec, 1} = ...
                strrep(temp_anatlabel_sourceelec2, ', ', '_'); %Full anat label
            
            %Read out anatomical region label for target elecs (per source elec)
            temp_num_targetelecs = size(SelElecs.Pairs.(param.GC.label_ElecPairSel{1}).(...
                param.GC.label_ElecPairSel{2}).Pairs_perelec.(...
                temp_anatlabel_sourceelec),1);
            for i_targetelec = 1:temp_num_targetelecs
                fulllabel_targetelec{i_sub}{i_targetelec, i_sourceelec} = ...
                    cell2mat(strcat(SelElecs.Pairs.(param.GC.label_ElecPairSel{1}).(...
                    param.GC.label_ElecPairSel{2}).Pairs_perelec.(...
                    temp_anatlabel_sourceelec){i_targetelec,1}, '_', ...
                    SelElecs.Pairs.(param.GC.label_ElecPairSel{1}).(...
                    param.GC.label_ElecPairSel{2}).Pairs_perelec.(...
                    temp_anatlabel_sourceelec){i_targetelec,2}, '_', ...
                    SelElecs.Pairs.(param.GC.label_ElecPairSel{1}).(...
                    param.GC.label_ElecPairSel{2}).Pairs_perelec.(...
                    temp_anatlabel_sourceelec){i_targetelec,5}));
                fulllabel_targetelec{i_sub}{i_targetelec, i_sourceelec} = ...
                    strrep(fulllabel_targetelec{i_sub}{i_targetelec, i_sourceelec}, ', ', '_');
                anatlabel_targetelec{i_sub}{i_targetelec, i_sourceelec} = ...
                    cell2mat(SelElecs.Pairs.(param.GC.label_ElecPairSel{1}).(...
                    param.GC.label_ElecPairSel{2}).Pairs_perelec.(...
                    temp_anatlabel_sourceelec){i_targetelec,5});
                anatlabel_targetelec{i_sub}{i_targetelec, i_sourceelec} = ...
                    strrep(anatlabel_targetelec{i_sub}{i_targetelec, i_sourceelec}, ', ', '_');
            end
        end
    end
end

%2.2 Read out all present anatomical locations across all subjects
anatlabel_sourceelec_allsub = [];
anatlabel_targetelec_allsub = [];

for i_sub = 1:length(subs)
    if ~isempty(anatlabel_sourceelec{i_sub})
        
        anatlabel_sourceelec_allsub = ...
            [anatlabel_sourceelec_allsub; anatlabel_sourceelec{i_sub}];
        
        for i_sourceelec = 1:length(anatlabel_sourceelec{i_sub})
            for i_targetelec = 1:size(anatlabel_targetelec{i_sub},1)
                if ~isempty(anatlabel_targetelec{i_sub}{i_targetelec, i_sourceelec})
                    anatlabel_targetelec_allsub = ...
                        [anatlabel_targetelec_allsub; ...
                        anatlabel_targetelec{i_sub}(i_targetelec, i_sourceelec)];
                end
            end
        end
    end
end

anatreg.label_unique = unique([anatlabel_sourceelec_allsub; anatlabel_targetelec_allsub]);

%Create filter limiting analyzed anatomical regions to provided label list
for i_catlabels = 1:length(anatreg.label_unique)
    anatreg.filter(i_catlabels,1) = ...
    any(strcmp(AnatReg_CatLabels, anatreg.label_unique{i_catlabels}));
end
anatreg.filter_label = [num2str(sum(anatreg.filter)) 'AnatReg'];
  
%Flexibly determine if and how to group anat regions
%For each unique label, create index and index of all sub-labels
anatreg.labelindex_aggregated = [];
for i_catlabels = 1:length(anatreg.label_unique)
    temp_index = ...
        find(~cellfun('isempty', strfind(...
        anatreg.label_unique, anatreg.label_unique{i_catlabels})));    
    if ~isempty(temp_index)
        anatreg.labelindex_aggregated(i_catlabels,:) = temp_index;
    else
        anatreg.labelindex_aggregated(i_catlabels,:) = NaN;
    end
end
%Adjust anat cat labels
anatreg.label_aggregated = [];
for i_catlabels = 1:length(anatreg.label_unique)
    if strcmp(anatreg.label_unique{i_catlabels}, 'AntPFC')
        anatreg.label_aggregated{i_catlabels,1} = 'Frontal';
    elseif strcmp(anatreg.label_unique{i_catlabels}, 'AntPFC_IFG')
        anatreg.label_aggregated{i_catlabels,1} = 'IFG';
    elseif strcmp(anatreg.label_unique{i_catlabels}, 'AntPFC_PrecentralG')
        anatreg.label_aggregated{i_catlabels,1} = 'PrecG';
    elseif strcmp(anatreg.label_unique{i_catlabels}, 'VentralT')
        anatreg.label_aggregated{i_catlabels,1} = 'Temporal';
    elseif strcmp(anatreg.label_unique{i_catlabels}, 'VentralT_STG')
        anatreg.label_aggregated{i_catlabels,1} = 'STG';
    elseif strcmp(anatreg.label_unique{i_catlabels}, 'VentralT_MTG')
        anatreg.label_aggregated{i_catlabels,1} = 'MTG';
    elseif strcmp(anatreg.label_unique{i_catlabels}, 'SupParLob')
        anatreg.label_aggregated{i_catlabels,1} = 'Parietal';
    elseif strcmp(anatreg.label_unique{i_catlabels}, 'SupParLob_PostcentralG')
        anatreg.label_aggregated{i_catlabels,1} = 'PostcG';
    elseif strcmp(anatreg.label_unique{i_catlabels}, 'SupParLob_SupramarginalG')
        anatreg.label_aggregated{i_catlabels,1} = 'SupmarG';
    elseif strcmp(anatreg.label_unique{i_catlabels}, 'OccipitalL')
        anatreg.label_aggregated{i_catlabels,1} = 'OcciL';
    end
end


%2.3 Average GC estimates across subjects per selected anatomical region
%Define anatreg index and label for each source and target elec per subject
for i_sub = 1:length(subs)
    if ~isempty(GCdata{i_sub})
        
        for i_sourceelec = 1:length(GCdata{i_sub}.ind_pairedelecs_from)
            for i_anatlabels = 1:length(anatreg.label_aggregated)
                if contains(anatlabel_sourceelec{i_sub}{i_sourceelec}, ...
                        anatreg.label_unique{i_anatlabels})
                    anatreg.sourceelec_index{i_sub}(i_sourceelec,1) = i_anatlabels;
                    anatreg.sourceelec_label{i_sub}{i_sourceelec,1} = anatreg.label_aggregated{i_anatlabels};
                end
            end
            %Ensure to use target elecs of a source elec where no
            %target elec is missing (possible in case of double elecs)
            if ~any(cellfun(@isempty,anatlabel_targetelec{i_sub}(:,i_sourceelec)))
                sel_sourceelec = i_sourceelec;
            end            
        end       
        for i_targetelec = 1:length(GCdata{i_sub}.ind_pairedelecs_to)
            for i_anatlabels = 1:length(anatreg.label_aggregated)
                if ~isempty(anatlabel_targetelec{i_sub}(:,sel_sourceelec))
                    if contains(anatlabel_targetelec{i_sub}{i_targetelec,sel_sourceelec}, ...
                            anatreg.label_unique{i_anatlabels})
                        anatreg.targetelec_index{i_sub}(i_targetelec,1) = ...
                            i_anatlabels;
                        anatreg.targetelec_label{i_sub}{i_targetelec,1} = ...
                            anatreg.label_aggregated{i_anatlabels};
                    end
                end
            end
        end
        
    end
end

%Average temporal and spectral GC values across subjects for each anatomical region
%Create empty structs
Gavg_tempGC_peranatreg.source2target = ...
    nan(length(anatreg.label_aggregated), ...
    length(anatreg.label_aggregated), ...
    size(GCdata{i_sub}.temporalGC{1}.source2target,3), ...
    length(subs), ...
    length(ToneDur_text)); %Dim: AnatLabelSource, AnatLabelTarget, Epochs/Tones, Subs, TD
Gavg_tempGC_peranatreg.target2source = Gavg_tempGC_peranatreg.source2target;

Gavg_spectGC_peranatreg.source2target = ...
    nan(length(anatreg.label_aggregated), ...
    length(anatreg.label_aggregated), ...
    size(GCdata{i_sub}.spectralGC{1}.source2target,3), ...
    size(GCdata{i_sub}.spectralGC{1}.source2target,4), ...
    length(subs), ...
    length(ToneDur_text)); %Dim: AnatLabelSource, AnatLabelTarget, Freqs, Epochs/Tones, Subs, TD
Gavg_spectGC_peranatreg.target2source = Gavg_spectGC_peranatreg.source2target;

NumGC_peranatreg = ...
    zeros(length(anatreg.label_aggregated), ...
    length(anatreg.label_aggregated), ...
    length(ToneDur_text)); %Dim: AnatLabelSource, AnatLabelTarget, TD

%Aggregate GC estimates for similar anat reg in common struct
for i_sub = 1:length(subs)
    if ~isempty(GCdata{i_sub})        
        for i_TD = 1:length(ToneDur_text)
            
            for i_anatlabels_source = 1:length(anatreg.label_aggregated)
                for i_sourceelec = 1:length(GCdata{i_sub}.ind_pairedelecs_from)
                    if any(anatreg.sourceelec_index{i_sub}(i_sourceelec) == anatreg.labelindex_aggregated(i_anatlabels_source,:))
                        
                        for i_anatlabels_target = 1:length(anatreg.label_aggregated)
                            temp_tempGC.source2target = []; %Proxy vars for averaging
                            temp_tempGC.target2source = [];                            
                            temp_spectGC.target2source = [];
                            temp_spectGC.source2target = [];
                            for i_targetelec = 1:length(GCdata{i_sub}.ind_pairedelecs_to)
                                if any(anatreg.targetelec_index{i_sub}(i_targetelec) == anatreg.labelindex_aggregated(i_anatlabels_target,:))
                                    %Aggregate temporal GC estimates across elecs from identical source-target anat regs                                   
                                    temp_tempGC.source2target = ...
                                        [temp_tempGC.source2target; ... 
                                        squeeze(GCdata{i_sub}.temporalGC{i_TD}.source2target(i_sourceelec, i_targetelec, :))'];
                                    temp_tempGC.target2source = ...
                                        [temp_tempGC.target2source; ... 
                                        squeeze(GCdata{i_sub}.temporalGC{i_TD}.target2source(i_sourceelec, i_targetelec, :))'];                                    
                                    temp_spectGC.source2target(:,:,i_targetelec) = ...
                                        squeeze(GCdata{i_sub}.spectralGC{i_TD}.source2target(i_sourceelec, i_targetelec, :, :));
                                    temp_spectGC.target2source(:,:,i_targetelec) = ...
                                        squeeze(GCdata{i_sub}.spectralGC{i_TD}.target2source(i_sourceelec, i_targetelec, :, :));                                                                        
                                    %Count pairings for current anat region combination
                                    NumGC_peranatreg(i_anatlabels_source, i_anatlabels_target, i_TD) = ...
                                        NumGC_peranatreg(i_anatlabels_source, i_anatlabels_target, i_TD) + 1;
                                end
                            end
                            
                            %Average GC across elecs from identical source-target anat regs
                            Gavg_tempGC_peranatreg.source2target(i_anatlabels_source, i_anatlabels_target, :, i_sub, i_TD) = ...
                                nanmean(temp_tempGC.source2target,1);
                            Gavg_tempGC_peranatreg.target2source(i_anatlabels_target, i_anatlabels_source, :, i_sub, i_TD) = ...
                                nanmean(temp_tempGC.target2source,1);                            
                            Gavg_spectGC_peranatreg.source2target(i_anatlabels_source, i_anatlabels_target, :, :, i_sub, i_TD) = ...
                                nanmean(temp_spectGC.source2target,3);
                            Gavg_spectGC_peranatreg.target2source(i_anatlabels_target, i_anatlabels_source, :, :, i_sub, i_TD) = ...                                
                                nanmean(temp_spectGC.target2source,3);
                        end
                    end
                end
            end
        end
    end
end

%Optional: Perform 1 sample KS test to test GC values across subs (per pairing, time point, TD)
%for normality (0 = null hypothesis of normal dist, 1 = rejection of null)

for i_TD = 1:length(ToneDur_text)
    KS_mat{i_TD} = ...
        nan(length(anatreg.label_aggregated), ...
        length(anatreg.label_aggregated), ...
        size(GCdata{i_sub}.temporalGC{1}.source2target,3));

    for i_sourceelec = 1:length(anatreg.label_aggregated)
        for i_targetelec = 1:length(anatreg.label_aggregated)
            for i_toneepoch = 1:length(param.GC.tone_aggregation)
                if any(~isnan(squeeze(Gavg_tempGC_peranatreg.source2target(i_sourceelec, i_targetelec, i_toneepoch, :, i_TD))))
                    KS_mat{i_TD}(i_sourceelec, i_targetelec, i_toneepoch) = ...
                        kstest(squeeze(...
                        Gavg_tempGC_peranatreg.source2target(...
                        i_sourceelec, i_targetelec, i_toneepoch, :, i_TD)), ...
                        'Alpha', 0.05);
                end
            end
        end
    end
    
    KStestdata.sum_nonan{i_TD} = length(find(~isnan(KS_mat{i_TD}(:,:,:))));
    KStestdata.sum_sign{i_TD} = length(find(KS_mat{i_TD}(:,:,:) == 1));
    disp([num2str(KStestdata.sum_sign{i_TD}) ' / ' ...
        num2str(KStestdata.sum_nonan{i_TD}) ...
        ' entries not normal distributed'])
end


%Average GC estimates across subs
%Mean
Gavg_tempGC_peranatreg.source2target = ...
    squeeze(nanmean(Gavg_tempGC_peranatreg.source2target,4));
Gavg_tempGC_peranatreg.target2source = ...
    squeeze(nanmean(Gavg_tempGC_peranatreg.target2source,4));%Dim: AnatLabelSource, AnatLabelTarget, Epochs/Tones, TD
Gavg_spectGC_peranatreg.source2target = ...
    squeeze(nanmean(Gavg_spectGC_peranatreg.source2target,5));
Gavg_spectGC_peranatreg.target2source = ...
    squeeze(nanmean(Gavg_spectGC_peranatreg.target2source,5));%Dim: AnatLabelSource, AnatLabelTarget, Freqs, Epochs/Tones, TD

%Median
% Gavg_tempGC_peranatreg.source2target = ...
%     squeeze(nanmedian(Gavg_tempGC_peranatreg.source2target,4));
% Gavg_tempGC_peranatreg.target2source = ...
%     squeeze(nanmedian(Gavg_tempGC_peranatreg.target2source,4));%Dim: AnatLabelSource, AnatLabelTarget, Epochs/Tones, TD
% Gavg_spectGC_peranatreg.source2target = ...
%     squeeze(nanmedian(Gavg_spectGC_peranatreg.source2target,5));
% Gavg_spectGC_peranatreg.target2source = ...
%     squeeze(nanmedian(Gavg_spectGC_peranatreg.target2source,5));%Dim: AnatLabelSource, AnatLabelTarget, Freqs, Epochs/Tones, TD


%% 3. Plot GC estimate for time and frequency domain
colormap_parula = [0.7 0.7 0.7; parula(10000)]; %Ensures that NaN values are greyed out

%3.1 Plot Source to Target GC as function of (aggregated) tones over the course of the sequence
for i_TD = 1:length(ToneDur_text)
    
    %Summary GC and number of pairings across anat reg per tone group
    f1 = figure('visible',set_visibility);
    set(f1,'units','normalized','outerposition',[0 0 1 1]) %full screen
    DimSubplot = [length(param.GC.tone_aggregation) 1];
    sgtitle(['Time-domain Gavg (n = ' num2str(length(subs)) ...
        ') pw GC estimates per tone group across anat regions - ' ...
        ToneDur_text{i_TD} 's TD - ' ...
        param.GC.label_ElecPairSel{1} ', ' param.GC.label_ElecPairSel{2} ' - ' ...
        param.GC.InputDataType{1}], 'Interpreter', 'none')
    subplot_num = 0;
    
    %Determine Clim
    for i_tonegroup = 1:length(param.GC.tone_aggregation)
        maxC_tempGC{i_TD}(i_tonegroup,1) = max(max(nanmean(...
            Gavg_tempGC_peranatreg.source2target...
            (anatreg.filter,anatreg.filter,i_tonegroup, i_TD))));
    end
    maxC_tempGC{i_TD} = max(maxC_tempGC{i_TD});
    
    for i_tonegroup = 1:length(param.GC.tone_aggregation)
        %GC values
        subplot_num = subplot_num + 1;
        
        f_sub = subplot(DimSubplot(1),DimSubplot(2),subplot_num);
        imagesc(Gavg_tempGC_peranatreg.source2target...
            (anatreg.filter,anatreg.filter,i_tonegroup,i_TD))
        f_sub.XTick =  [1:length(anatreg.label_aggregated(anatreg.filter))];
        f_sub.XTickLabel = anatreg.label_aggregated(anatreg.filter);
        a = get(gca, 'XTickLabel');
        set(gca,'XTickLabel', a, 'fontsize', 8);
        f_sub.XLabel.String = {'Target Elecs'};
        f_sub.XLabel.FontSize = 12;
        f_sub.YTick =  [1:length(anatreg.label_aggregated(anatreg.filter))];
        f_sub.YTickLabel = anatreg.label_aggregated(anatreg.filter);
        f_sub.YLabel.String = {'Source Elecs'};
        f_sub.YLabel.FontSize = 12;
        f_sub.CLim = [0 maxC_tempGC{i_TD}];
        f_sub.Colormap = colormap_parula;
        set(gca,'TickLabelInterpreter', 'none');
        
        colorbar
        title(['GC estimates - Tone ' ...
            num2str(param.GC.tone_aggregation{i_tonegroup}(1)) ...
            ' to ' num2str(param.GC.tone_aggregation{i_tonegroup}(end))], ...
            'fontsize', 12);
        
        %Add text specifying number of pairings to each subfield
        NumGC_peranatreg_filter = NumGC_peranatreg(anatreg.filter, anatreg.filter, i_TD);
        for i_imagefield_x = 1:size(Gavg_tempGC_peranatreg.source2target...
                (anatreg.filter,anatreg.filter,i_tonegroup,i_TD),1)
            for i_imagefield_y = 1:size(Gavg_tempGC_peranatreg.source2target...
                    (anatreg.filter,anatreg.filter,i_tonegroup,i_TD),2)
                text(i_imagefield_y,i_imagefield_x, ...
                    num2str(NumGC_peranatreg_filter(i_imagefield_x, i_imagefield_y)), ...
                    'fontsize', 6)
            end
        end
    end
    
    if istrue(save_poststepFigs)
        filename     = ['Heatmap_AllsubN' num2str(length(subs)) '_'...
            'TemppwGC_' anatreg.filter_label '_TD' num2str(i_TD) '_' ...
            title_label1 '_' title_label2 '_'  ...
            param.GC.InputDataType{1} param.GC.tone_aggregation_label '.png'];
        figfile      = [path_fig  filename];
        saveas(gcf, figfile, 'png'); %save png version
        close;
    end
    
    %Summary GC and number of pairings across anat reg for toneindex_noreverse tone group
    if ismember(toneindex_noreverse, cellfun(@(c)[c],param.GC.tone_aggregation))
        f1 = figure('visible',set_visibility);
        set(f1,'units','normalized','outerposition',[0 0 1 1]) %full screen
        DimSubplot = [1 1];
        sgtitle(['Time-domain Gavg (n = ' num2str(length(subs)) ...
            ') pw GC estimates during p' num2str(toneindex_noreverse) ' across anat regions - ' ...
            ToneDur_text{i_TD} 's TD - ' ...
            param.GC.label_ElecPairSel{1} ', ' param.GC.label_ElecPairSel{2} ' - ' ...
            param.GC.InputDataType{1}], 'Interpreter', 'none')
        subplot_num = 0;
        
        %Determine Clim
        for i_tonegroup = 1:length(param.GC.tone_aggregation)
            if param.GC.tone_aggregation{i_tonegroup} == toneindex_noreverse
                index_noreverse = i_tonegroup;
                maxC_tempGC{i_TD} = max(max(nanmean(...
                    Gavg_tempGC_peranatreg.source2target...
                    (anatreg.filter,anatreg.filter,i_tonegroup, i_TD))));
            end
        end
        
        for i_tonegroup = index_noreverse
            %GC values
            subplot_num = subplot_num + 1;
            
            f_sub = subplot(DimSubplot(1),DimSubplot(2),subplot_num);
            imagesc(Gavg_tempGC_peranatreg.source2target...
                (anatreg.filter,anatreg.filter,i_tonegroup,i_TD))
            f_sub.XTick =  [1:length(anatreg.label_aggregated(anatreg.filter))];
            f_sub.XTickLabel = anatreg.label_aggregated(anatreg.filter);
            a = get(gca, 'XTickLabel');
            set(gca,'XTickLabel', a, 'fontsize', 8);
            f_sub.XLabel.String = {'Target Elecs'};
            f_sub.XLabel.FontSize = 12;
            f_sub.YTick =  [1:length(anatreg.label_aggregated(anatreg.filter))];
            f_sub.YTickLabel = anatreg.label_aggregated(anatreg.filter);
            f_sub.YLabel.String = {'Source Elecs'};
            f_sub.YLabel.FontSize = 12;
            f_sub.CLim = [0 maxC_tempGC{i_TD}];
            f_sub.Colormap = colormap_parula;
            set(gca,'TickLabelInterpreter', 'none');
            
            colorbar
            title(['Time-domain GAvg (n = ' num2str(length(subs)) ') GC estimates - Tone ' ...
                num2str(param.GC.tone_aggregation{i_tonegroup}(1)) ' to ' num2str(param.GC.tone_aggregation{i_tonegroup}(end))]);
            
            %Add text specifying number of pairings to each subfield
            NumGC_peranatreg_filter = NumGC_peranatreg(anatreg.filter, anatreg.filter, i_TD);
            for i_imagefield_x = 1:size(Gavg_tempGC_peranatreg.source2target...
                    (anatreg.filter,anatreg.filter,i_tonegroup,i_TD),1)
                for i_imagefield_y = 1:size(Gavg_tempGC_peranatreg.source2target...
                        (anatreg.filter,anatreg.filter,i_tonegroup,i_TD),2)
                    text(i_imagefield_y,i_imagefield_x, ...
                        num2str(NumGC_peranatreg_filter(i_imagefield_x, i_imagefield_y)))
                end
            end
        end
        
        if istrue(save_poststepFigs)
            filename     = ['Heatmap_AllsubN' num2str(length(subs)) '_'...
                'TemppwGC_per' anatreg.filter_label '_p' num2str(toneindex_noreverse) '_TD' num2str(i_TD) '_' ...
                title_label1 '_' title_label2 '_'  ...
                param.GC.InputDataType{1} param.GC.tone_aggregation_label '.png'];
            figfile      = [path_fig  filename];
            saveas(gcf, figfile, 'png'); %save png version
            close;
        end
    end
    
    %3.2 Plot Source to Target temporal GC as function of (aggregated) tones over the course of the sequence
    DimSubplot = [length(anatreg.label_aggregated(anatreg.filter)) ...
        length(anatreg.label_aggregated(anatreg.filter))];
    i_subplot = 0;
    f1 = figure('visible',set_visibility);
    set(f1,'units','normalized','outerposition',[0 0 1 1]) %full screen
    %Determine common y-axis limits
    for i_tonegroup = 1:length(param.GC.tone_aggregation)
        maxY_tempGC{i_TD}(i_tonegroup,1) = max(max(max(...
            Gavg_tempGC_peranatreg.source2target...
            (anatreg.filter,anatreg.filter,i_tonegroup, i_TD))));
    end
    maxY_tempGC{i_TD} = max(maxY_tempGC{i_TD});
    
    %Plot barcharts across tone epochs for each anat reg connection
    for i_sourceelec = 1:length(anatreg.label_aggregated(anatreg.filter))
        for i_targetelec = 1:length(anatreg.label_aggregated(anatreg.filter))
            
            i_subplot = i_subplot + 1;
            f_sub = subplot(DimSubplot(1),DimSubplot(2),i_subplot);
            
            Gavg_tempGC_peranatreg_filter.source2target = ...
                Gavg_tempGC_peranatreg.source2target(anatreg.filter,anatreg.filter,:,i_TD);
            if any(~isnan(Gavg_tempGC_peranatreg_filter.source2target(i_sourceelec, i_targetelec,:)))
                
                meandata_bar = [];
                for i_tonegroup = 1:length(param.GC.tone_aggregation)
                    meandata_bar(1, i_tonegroup) = ...
                        Gavg_tempGC_peranatreg_filter.source2target...
                        (i_sourceelec,i_targetelec,i_tonegroup);
                    label_Xtick{i_tonegroup} = ...
                        [num2str(param.GC.tone_aggregation{i_tonegroup}(1)) '-' ...
                        num2str(param.GC.tone_aggregation{i_tonegroup}(end))];
                end
                
                b = bar(meandata_bar, 'facecolor', [0.7 0.7 0.7]);
                
                %Optional: test for linear trend over tone groups
                x = 1:length(meandata_bar);
                linregmodel = fitlm(x, meandata_bar, 'linear');
                linCoeffs = polyfit(x, meandata_bar, 1) ;              
                f = polyval(linCoeffs, x);                
                if linregmodel.Coefficients{2,4} < 0.05
                    b.FaceColor = [0 0 1];
                end                
                hold on
                plot(x,meandata_bar,'.',x,f,'-','LineWidth',2, 'Color', 'k')
                %To do: add p-value to plot
                
                %Adjust axes                
                f_sub.XTickLabel = label_Xtick;
                f_sub.XTickLabelRotation = 45;
                a = get(gca, 'XTickLabel');
                set(gca,'XTickLabel', a, 'fontsize', 6);
                f_sub.XLabel.String = {'Tones'};
                f_sub.XLabel.FontSize = 6;
                f_sub.YLabel.String = {'GC estimate'};
                f_sub.YLabel.FontSize = 6;
%                 f_sub.YLim = [0 maxY_tempGC{i_TD}*1.1];
            else
                axis off
            end
            anatreg.label_aggregated_filter = ...
                anatreg.label_aggregated(anatreg.filter);
            NumGC_peranatreg_filter = ...
                NumGC_peranatreg(anatreg.filter, anatreg.filter, i_TD);            
            title({[anatreg.label_aggregated_filter{i_sourceelec} ' -> '], ...
                [anatreg.label_aggregated_filter{i_targetelec} ...
                ' (' num2str(NumGC_peranatreg_filter(i_sourceelec, i_targetelec)) ')']}, ...
                'Interpreter', 'none', 'FontSize', 8)
        end
    end
    sgtitle(['Time-domain Gavg (n = ' num2str(length(subs)) ') pw GC estimates per tone group per anat region - ' ...
        ToneDur_text{i_TD} 's TD - ' param.GC.label_ElecPairSel{1} ', ' ...
        param.GC.label_ElecPairSel{2} ' - ' param.GC.InputDataType{1}], ...
        'Interpreter', 'none')
    
    if istrue(save_poststepFigs)
        filename     = ['ErrBar_AllsubN' num2str(length(subs)) ...
            '_TemppwGC_' anatreg.filter_label '_' 'TD' num2str(i_TD) '_' ...
            title_label1 '_' title_label2 '_'  ...
            param.GC.InputDataType{1} param.GC.tone_aggregation_label '.png'];
        figfile      = [path_fig  filename];
        saveas(gcf, figfile, 'png'); %save png version
        close;
    end
    
    
    %3.3 Plot Source to Target spectral GC as function of (aggregated) tones over the course of the sequence
    DimSubplot = [length(anatreg.label_aggregated(anatreg.filter)) ...
        length(anatreg.label_aggregated(anatreg.filter))];
    i_subplot = 0;
    f2 = figure('visible',set_visibility);
    set(f2,'units','normalized','outerposition',[0 0 1 1]) %full screen
    colormap_cool = cool(length(param.GC.tone_aggregation));
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
            max(max(max(nanmean(...
            Gavg_spectGC_peranatreg.source2target...
            (anatreg.filter,anatreg.filter,temp_xlim(1):temp_xlim(2),i_tonegroup, i_TD),4))));
    end
    maxY_spectGC{i_TD} = max(maxY_spectGC{i_TD});
    
    %Plot line plots across tone epochs for each anat reg connection
    for i_sourceelec = 1:length(anatreg.label_aggregated(anatreg.filter))
        for i_targetelec = 1:length(anatreg.label_aggregated(anatreg.filter))
            
            i_subplot = i_subplot + 1;
            f_sub = subplot(DimSubplot(1),DimSubplot(2),i_subplot);
            
            Gavg_spectGC_peranatreg_filter.source2target = ...
                Gavg_spectGC_peranatreg.source2target(anatreg.filter,anatreg.filter,:,:,i_TD);
            
            if any(any(~isnan(Gavg_spectGC_peranatreg_filter.source2target(i_sourceelec, i_targetelec,:, :))))
                
                meandata_line = [];
                for i_tonegroup = 1:length(param.GC.tone_aggregation)
                    meandata_line(i_tonegroup,:) = ...
                        Gavg_spectGC_peranatreg_filter.source2target(...
                        i_sourceelec,i_targetelec,:,i_tonegroup);
                end
                
                hold on;
                for i_tonegroup = 1:length(param.GC.tone_aggregation)
                    plot(meandata_line(i_tonegroup,:), ...
                        'Color', colormap_cool(i_tonegroup,:), 'LineWidth', 2, 'LineStyle', '-');
                end
                
                f_sub.XLim = temp_xlim;
                f_sub.XTick = temp_xtick;
                f_sub.XTickLabel = temp_xticklabel;
                a = get(gca, 'XTickLabel');
                set(gca,'XTickLabel', a, 'fontsize', 4);
                f_sub.XLabel.String = {'Frequency [Hz]'};
                f_sub.YLabel.String = {'GC'};
                f_sub.XLabel.FontSize = 6;
                f_sub.YLabel.FontSize = 6;
                %                         f_sub.YLim = [0 maxY_spectGC{i_TD}];
                
                if i_sourceelec ==1 && i_targetelec == 1
                    legend(label_Xtick)
                end
            else
                axis off
            end
            
            anatreg.label_aggregated_filter = ...
                anatreg.label_aggregated(anatreg.filter);
            NumGC_peranatreg_filter = ...
                NumGC_peranatreg(anatreg.filter, anatreg.filter, i_TD);            
            title({[anatreg.label_aggregated_filter{i_sourceelec} ' -> '], ...
                [anatreg.label_aggregated_filter{i_targetelec} ...
                ' (' num2str(NumGC_peranatreg_filter(i_sourceelec, i_targetelec)) ')']}, ...
                'Interpreter', 'none', 'FontSize', 6)            
        end
    end
    sgtitle(['Spectral Gavg (n = ' num2str(length(subs)) ') pw GC estimates per tone group per anat region - ' ...
        ToneDur_text{i_TD} 's TD - ' param.GC.label_ElecPairSel{1} ', ' ...
        param.GC.label_ElecPairSel{2} ' - ' param.GC.InputDataType{1}], ...
        'Interpreter', 'none')
    
    if istrue(save_poststepFigs)
        filename     = ['Line_AllsubN' num2str(length(subs)) ...
            '_SpectpwGC_' anatreg.filter_label '_' 'TD' num2str(i_TD) '_' ...
            title_label1 '_' title_label2 '_'  ...
            param.GC.InputDataType{1} param.GC.tone_aggregation_label '.png'];
        figfile      = [path_fig  filename];
        saveas(gcf, figfile, 'png'); %save png version
        close;
    end
    
    
    %% 4. Plot reversed (Target to Source) GC as function of (aggregated) tones over the course of the sequence
    %4.1 Summary GC and number of pairings across anat reg per tone group
    f1 = figure('visible',set_visibility);
    set(f1,'units','normalized','outerposition',[0 0 1 1]) %full screen
    DimSubplot = [length(param.GC.tone_aggregation) 1];
    sgtitle(['Reversed (target-to-source) time-domain Gavg (n = ' num2str(length(subs)) ...
        ') pw GC estimates per tone group across anat regions - ' ...
        ToneDur_text{i_TD} 's TD - ' ...
        param.GC.label_ElecPairSel{1} ', ' param.GC.label_ElecPairSel{2} ' - ' ...
        param.GC.InputDataType{1}], 'Interpreter', 'none')
    subplot_num = 0;
    
    %Determine Clim
    for i_tonegroup = 1:length(param.GC.tone_aggregation)
        maxC_tempGC{i_TD}(i_tonegroup,1) = max(max(nanmean(...
            Gavg_tempGC_peranatreg.target2source...
            (anatreg.filter,anatreg.filter,i_tonegroup, i_TD))));
    end
    maxC_tempGC{i_TD} = max(maxC_tempGC{i_TD});
    
    for i_tonegroup = 1:length(param.GC.tone_aggregation)
        %GC values
        subplot_num = subplot_num + 1;
        
        f_sub = subplot(DimSubplot(1),DimSubplot(2),subplot_num);
        imagesc(Gavg_tempGC_peranatreg.target2source...
            (anatreg.filter,anatreg.filter,i_tonegroup,i_TD))
        f_sub.XTick =  [1:length(anatreg.label_aggregated(anatreg.filter))];
        f_sub.XTickLabel = anatreg.label_aggregated(anatreg.filter);
        a = get(gca, 'XTickLabel');
        set(gca,'XTickLabel', a, 'fontsize', 8);
        f_sub.XLabel.String = {'Source (rev. Target) Elecs'};
        f_sub.XLabel.FontSize = 10;
        f_sub.YTick =  [1:length(anatreg.label_aggregated(anatreg.filter))];
        f_sub.YTickLabel = anatreg.label_aggregated(anatreg.filter);
        f_sub.YLabel.String = {'Target (rev. Source) Elecs'};
        f_sub.YLabel.FontSize = 10;
        f_sub.CLim = [0 maxC_tempGC{i_TD}];
        f_sub.Colormap = colormap_parula;
        set(gca,'TickLabelInterpreter', 'none');
        
        colorbar
        title(['GC estimates - Tone ' ...
            num2str(param.GC.tone_aggregation{i_tonegroup}(1)) ...
            ' to ' num2str(param.GC.tone_aggregation{i_tonegroup}(end))], ...
            'fontsize', 12);
        
        %Add text specifying number of pairings to each subfield
        NumGC_peranatreg_filter = NumGC_peranatreg(anatreg.filter, anatreg.filter, i_TD);
        for i_imagefield_x = 1:size(Gavg_tempGC_peranatreg.target2source...
                (anatreg.filter,anatreg.filter,i_tonegroup,i_TD),1)
            for i_imagefield_y = 1:size(Gavg_tempGC_peranatreg.target2source...
                    (anatreg.filter,anatreg.filter,i_tonegroup,i_TD),2)
                text(i_imagefield_y,i_imagefield_x, ...
                    num2str(NumGC_peranatreg_filter(i_imagefield_x, i_imagefield_y)), ...
                    'fontsize', 6)
            end
        end
    end
    
    if istrue(save_poststepFigs)
        filename     = ['Heatmap_AllsubN' num2str(length(subs)) '_'...
            'TemppwGC_' anatreg.filter_label '_TD' num2str(i_TD) '_reversed_' ...
            title_label1 '_' title_label2 '_'  ...
            param.GC.InputDataType{1} param.GC.tone_aggregation_label '.png'];
        figfile      = [path_fig2  filename];
        saveas(gcf, figfile, 'png'); %save png version
        close;
    end
    
    %Summary GC and number of pairings across anat reg for toneindex_reverse tone group
    if ismember(toneindex_reverse, cellfun(@(c)[c],param.GC.tone_aggregation))
        f1 = figure('visible',set_visibility);
        set(f1,'units','normalized','outerposition',[0 0 1 1]) %full screen
        DimSubplot = [1 1];
        sgtitle(['Reversed (target-to-source) time-domain Gavg (n = ' num2str(length(subs)) ...
            ') pw GC estimates during p' num2str(toneindex_reverse) ' across anat regions - ' ...
            ToneDur_text{i_TD} 's TD - ' ...
            param.GC.label_ElecPairSel{1} ', ' param.GC.label_ElecPairSel{2} ' - ' ...
            param.GC.InputDataType{1}], 'Interpreter', 'none')
        subplot_num = 0;
        
        %Determine Clim
        for i_tonegroup = 1:length(param.GC.tone_aggregation)
            if param.GC.tone_aggregation{i_tonegroup} == toneindex_reverse
                index_reverse = i_tonegroup;
                maxC_tempGC{i_TD} = max(max(nanmean(...
                    Gavg_tempGC_peranatreg.target2source...
                    (anatreg.filter,anatreg.filter,i_tonegroup, i_TD))));
            end
        end
        
        for i_tonegroup = index_reverse
            %GC values
            subplot_num = subplot_num + 1;
            
            f_sub = subplot(DimSubplot(1),DimSubplot(2),subplot_num);
            imagesc(Gavg_tempGC_peranatreg.target2source...
                (anatreg.filter,anatreg.filter,i_tonegroup,i_TD))
            f_sub.XTick =  [1:length(anatreg.label_aggregated(anatreg.filter))];
            f_sub.XTickLabel = anatreg.label_aggregated(anatreg.filter);
            a = get(gca, 'XTickLabel');
            set(gca,'XTickLabel', a, 'fontsize', 8);
            f_sub.XLabel.String = {'Source (revsered Target) Elecs'};
            f_sub.XLabel.FontSize = 12;
            f_sub.YTick =  [1:length(anatreg.label_aggregated(anatreg.filter))];
            f_sub.YTickLabel = anatreg.label_aggregated(anatreg.filter);
            f_sub.YLabel.String = {'Target (reversed Source) Elecs'};
            f_sub.YLabel.FontSize = 12;
            f_sub.CLim = [0 maxC_tempGC{i_TD}];
            f_sub.Colormap = colormap_parula;
            set(gca,'TickLabelInterpreter', 'none');
            
            colorbar
            title(['Time-domain GAvg (n = ' num2str(length(subs)) ') GC estimates - Tone ' ...
                num2str(param.GC.tone_aggregation{i_tonegroup}(1)) ' to ' num2str(param.GC.tone_aggregation{i_tonegroup}(end))]);
            
            %Add text specifying number of pairings to each subfield
            NumGC_peranatreg_filter = NumGC_peranatreg(anatreg.filter, anatreg.filter, i_TD);
            for i_imagefield_x = 1:size(Gavg_tempGC_peranatreg.target2source...
                    (anatreg.filter,anatreg.filter,i_tonegroup,i_TD),1)
                for i_imagefield_y = 1:size(Gavg_tempGC_peranatreg.target2source...
                        (anatreg.filter,anatreg.filter,i_tonegroup,i_TD),2)
                    text(i_imagefield_y,i_imagefield_x, ...
                        num2str(NumGC_peranatreg_filter(i_imagefield_x, i_imagefield_y)))
                end
            end
        end
        
        if istrue(save_poststepFigs)
            filename     = ['Heatmap_AllsubN' num2str(length(subs)) '_'...
                'TemppwGC_' anatreg.filter_label '_p' num2str(toneindex_reverse) '_TD' num2str(i_TD) '_reversed_' ...
                title_label1 '_' title_label2 '_'  ...
                param.GC.InputDataType{1} param.GC.tone_aggregation_label '.png'];
            figfile      = [path_fig2  filename];
            saveas(gcf, figfile, 'png'); %save png version
            close;
        end
    end
    
    %4.2 Plot reversed Source to Target temporal GC as function of (aggregated) tones over the course of the sequence
    DimSubplot = [length(anatreg.label_aggregated(anatreg.filter)) ...
        length(anatreg.label_aggregated(anatreg.filter))];
    i_subplot = 0;
    f1 = figure('visible',set_visibility);
    set(f1,'units','normalized','outerposition',[0 0 1 1]) %full screen
    %Determine common y-axis limits
    for i_tonegroup = 1:length(param.GC.tone_aggregation)
        maxY_tempGC{i_TD}(i_tonegroup,1) = max(max(max(...
            Gavg_tempGC_peranatreg.target2source...
            (anatreg.filter,anatreg.filter,i_tonegroup, i_TD))));
    end
    maxY_tempGC{i_TD} = max(maxY_tempGC{i_TD});
    
    %Plot barcharts across tone epochs for each anat reg connection
    for i_sourceelec = 1:length(anatreg.label_aggregated(anatreg.filter))
        for i_targetelec = 1:length(anatreg.label_aggregated(anatreg.filter))
            
            i_subplot = i_subplot + 1;
            f_sub = subplot(DimSubplot(1),DimSubplot(2),i_subplot);

            Gavg_tempGC_peranatreg_filter.target2source = ...
                Gavg_tempGC_peranatreg.target2source(anatreg.filter,anatreg.filter,:,i_TD);            
            if any(~isnan(Gavg_tempGC_peranatreg_filter.target2source(i_sourceelec, i_targetelec,:)))
                
                meandata_bar = [];
                for i_tonegroup = 1:length(param.GC.tone_aggregation)
                    meandata_bar(1, i_tonegroup) = ...
                        Gavg_tempGC_peranatreg_filter.target2source...
                        (i_sourceelec,i_targetelec,i_tonegroup);
                    label_Xtick{i_tonegroup} = ...
                        [num2str(param.GC.tone_aggregation{i_tonegroup}(1)) '-' ...
                        num2str(param.GC.tone_aggregation{i_tonegroup}(end))];
                end
                
                b = bar(meandata_bar, 'facecolor', 'flat');
                
                f_sub.XTickLabel = label_Xtick;
                f_sub.XTickLabelRotation = 45;
                a = get(gca, 'XTickLabel');
                set(gca,'XTickLabel', a, 'fontsize', 4);
                f_sub.XLabel.String = {'Tones'};
                f_sub.XLabel.FontSize = 6;
                f_sub.YLabel.String = {'GC estimate'};
                f_sub.YLabel.FontSize = 6;
                %                         f_sub.YLim = [0 maxY_tempGC{i_TD}*1.1];
            else
                axis off
            end
            anatreg.label_aggregated_filter = ...
                anatreg.label_aggregated(anatreg.filter);
            NumGC_peranatreg_filter = ...
                NumGC_peranatreg(anatreg.filter, anatreg.filter, i_TD);            
            title({[anatreg.label_aggregated_filter{i_sourceelec} ' -> '], ...
                [anatreg.label_aggregated_filter{i_targetelec} ...
                ' (' num2str(NumGC_peranatreg_filter(i_sourceelec, i_targetelec)) ')']}, ...
                'Interpreter', 'none', 'FontSize', 6)
        end
    end
    sgtitle(['Reversed time-domain Gavg (n = ' num2str(length(subs)) ...
        ') pw GC estimates per tone group per anat region - ' ...
        ToneDur_text{i_TD} 's TD - ' param.GC.label_ElecPairSel{1} ', ' ...
        param.GC.label_ElecPairSel{2} ' - ' param.GC.InputDataType{1}], ...
        'Interpreter', 'none')
    
    if istrue(save_poststepFigs)
        filename     = ['ErrBar_AllsubN' num2str(length(subs)) ...
            '_TemppwGC_' anatreg.filter_label '_' 'TD' num2str(i_TD) '_reversed_' ...
            title_label1 '_' title_label2 '_'  ...
            param.GC.InputDataType{1} param.GC.tone_aggregation_label '.png'];
        figfile      = [path_fig2  filename];
        saveas(gcf, figfile, 'png'); %save png version
        close;
    end    
    
    
    %4.3 Plot Source to Target spectral GC as function of (aggregated) tones over the course of the sequence
    DimSubplot = [length(anatreg.label_aggregated(anatreg.filter)) ...
        length(anatreg.label_aggregated(anatreg.filter))];
    i_subplot = 0;
    f2 = figure('visible',set_visibility);
    set(f2,'units','normalized','outerposition',[0 0 1 1]) %full screen
    colormap_cool = cool(length(param.GC.tone_aggregation));
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
            max(max(max(nanmean(...
            Gavg_spectGC_peranatreg.target2source...
            (anatreg.filter,anatreg.filter,temp_xlim(1):temp_xlim(2),i_tonegroup, i_TD),4))));
    end
    maxY_spectGC{i_TD} = max(maxY_spectGC{i_TD});
    
    %Plot line plots across tone epochs for each anat reg connection
    for i_sourceelec = 1:length(anatreg.label_aggregated(anatreg.filter))
        for i_targetelec = 1:length(anatreg.label_aggregated(anatreg.filter))
            
            i_subplot = i_subplot + 1;
            f_sub = subplot(DimSubplot(1),DimSubplot(2),i_subplot);
            
            Gavg_spectGC_peranatreg_filter.target2source = ...
                Gavg_spectGC_peranatreg.target2source(anatreg.filter,anatreg.filter,:,:,i_TD);            
            if any(any(~isnan(Gavg_spectGC_peranatreg_filter.target2source(i_sourceelec, i_targetelec,:, :))))
                
                meandata_line = [];
                for i_tonegroup = 1:length(param.GC.tone_aggregation)
                    meandata_line(i_tonegroup,:) = ...
                        Gavg_spectGC_peranatreg_filter.target2source(...
                        i_sourceelec,i_targetelec,:,i_tonegroup);
                end
                
                hold on;
                for i_tonegroup = 1:length(param.GC.tone_aggregation)
                    plot(meandata_line(i_tonegroup,:), ...
                        'Color', colormap_cool(i_tonegroup,:), 'LineWidth', 2, 'LineStyle', '-');
                end
                
                f_sub.XLim = temp_xlim;
                f_sub.XTick = temp_xtick;
                f_sub.XTickLabel = temp_xticklabel;
                a = get(gca, 'XTickLabel');
                set(gca,'XTickLabel', a, 'fontsize', 4);
                f_sub.XLabel.String = {'Frequency [Hz]'};
                f_sub.YLabel.String = {'GC'};
                f_sub.XLabel.FontSize = 6;
                f_sub.YLabel.FontSize = 6;
                %                         f_sub.YLim = [0 maxY_spectGC{i_TD}];
                
                if i_sourceelec ==1 && i_targetelec == 1
                    legend(label_Xtick)
                end
            else
                axis off
            end
            
            anatreg.label_aggregated_filter = ...
                anatreg.label_aggregated(anatreg.filter);
            NumGC_peranatreg_filter = ...
                NumGC_peranatreg(anatreg.filter, anatreg.filter, i_TD);            
            title({[anatreg.label_aggregated_filter{i_sourceelec} ' -> '], ...
                [anatreg.label_aggregated_filter{i_targetelec} ...
                ' (' num2str(NumGC_peranatreg_filter(i_sourceelec, i_targetelec)) ')']}, ...
                'Interpreter', 'none', 'FontSize', 6)   
        end
    end
    sgtitle(['Reversed spectral Gavg (n = ' num2str(length(subs)) ...
        ') pw GC estimates per tone group per anat region - ' ...
        ToneDur_text{i_TD} 's TD - ' param.GC.label_ElecPairSel{1} ', ' ...
        param.GC.label_ElecPairSel{2} ' - ' param.GC.InputDataType{1}], ...
        'Interpreter', 'none')
    
    if istrue(save_poststepFigs)
        filename     = ['Line_AllsubN' num2str(length(subs)) ...
            '_SpectpwGC_' anatreg.filter_label '_' 'TD' num2str(i_TD) '_reversed_' ...
            title_label1 '_' title_label2 '_'  ...
            param.GC.InputDataType{1} param.GC.tone_aggregation_label '.png'];
        figfile      = [path_fig2  filename];
        saveas(gcf, figfile, 'png'); %save png version
        close;
    end    
    
end

end