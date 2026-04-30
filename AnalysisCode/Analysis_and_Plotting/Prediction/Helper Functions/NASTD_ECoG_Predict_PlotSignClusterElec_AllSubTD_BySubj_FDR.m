function NASTD_ECoG_Predict_PlotSignClusterElec_AllSubTD_BySubj_FDR...
    (subs, ...
    FuncInput_EffectType, FuncInput_DataType, FuncInput_ToneDur_text, ...
    param, ...
    save_poststepFigs, paths_NASTD_ECoG)

%Aim: Plot electrodes showing sign. cluster-corrected prediction effects
%per input data & TD condition aggregated over subjects and tone durations
%(TD) and projected on 1 Hemisphere.
%Based on sample-wise analysis, thus no time windows
%but instead entire analysis window.
%Plot Figure with 5 subfigures:
% - clusterstat
% - p-value derivative
% - onset
% - offset
% - duration

% set(0, 'DefaultFigureVisible', 'on');

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
%Add base dir and own script dir
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

if param.FDRcorrect == 1
    FDR_label = 'FDRcorr';
else
    FDR_label = 'uncorr';
end

path_fig = ([paths_NASTD_ECoG.ECoGdata_Prediction ...
    '/PredEffects/Allsub_n' num2str(length(subs)) ...
    '/Figs/PredEffects_Surf/' FuncInput_EffectType ...
    '/ClusterCorr/' param.ElecSelect '/' FDR_label '/']);
if (~exist(path_fig, 'dir')); mkdir(path_fig); end


%% 1) Load prediction effect data and aggregate relevant info across subjects
clear Param2plot
coords_allsub                       = [];
labels_allsub                       = [];
anatlabels_allsub                   = [];
Param2plot.all_subs.sub_index       = [];
Param2plot.all_subs.TD_index        = [];
usedElecs_chanposIndex              = [];

for i_sub = 1:length(subs)
    
    tic
    sub = subs{i_sub};
    disp(['-- Loading data for sub: ' sub ' --'])
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;
    
    for i_TD = 1:length(FuncInput_ToneDur_text)
        %Load prediction data and select current effect
        path_inputdata = [paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
            sub '/Data/Samplewise/'];
        load([path_inputdata sub '_PredEffectsCluster_' ...
            FuncInput_DataType '_' ...
            FuncInput_ToneDur_text{i_TD} 'sTD.mat'], ...
            FuncInput_EffectType , 'labels_loadedData');
        CurrentEffect = eval(FuncInput_EffectType);
        clear PredEffect SimplePredErrEffect ComplexPredErrEffect
        
        %Also load ECoG preproc data for channel labels and position
        loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];
        load(loadfile_ECoGpreprocdata);
        
        SampleFreq              = DataClean_AllTrials.fsample;
        nSensors_all            = size(CurrentEffect.stats,2);
        nSamples(i_TD)          = size(CurrentEffect.clusterstat{1}.cluster_timecourse,2);
        
        %Load stimulus correlation data and select current effect
        if strcmp(param.ElecSelect, 'StimCorr')
            path_inputdata = [paths_NASTD_ECoG.ECoGdata_StimCorr '/' sub '/Data/'];
            load([path_inputdata sub '_StimCorr_' FuncInput_DataType '_' ...
                FuncInput_ToneDur_text{i_TD} 'sTD.mat'], ...
                'corr_ttest', 'SensorLabels');
            filt_signelecs_StimCorr = corr_ttest.p < 0.05;%only sign. stim corr elecs
        else
            filt_signelecs_StimCorr  = true(length(labels_loadedData),1);  %all elecs
        end
        ind_signelecs_StimCorr  = find(filt_signelecs_StimCorr);
        nSensors_sel            = sum(filt_signelecs_StimCorr);
        
        %Read electrode labels, coordinates, and anatomical labels for all
        %sign. StimCorr electrodes and aggregate them across subjects
        for i_elec = 1:length(ind_signelecs_StimCorr)
            usedElecs_chanposIndex(i_elec,1) = ...
                find(strcmp(labels_loadedData{ind_signelecs_StimCorr(i_elec)}, ...
                DataClean_AllTrials.elec.label));
        end
        coords_sub{i_sub, i_TD}   = ...
            DataClean_AllTrials.elec.chanpos(usedElecs_chanposIndex,:);
        coords_allsub       = [coords_allsub; coords_sub{i_sub, i_TD}];
        anatlabels_allsub   = ...
            [anatlabels_allsub; DataClean_AllTrials.elec.T1AnatLabel(usedElecs_chanposIndex)];
        for i_elec = 1:length(ind_signelecs_StimCorr)
            labels_allsub{end+1,1} = [labels_loadedData{ind_signelecs_StimCorr(i_elec)} ' ' sub];
        end
        
        %Determine max. number of sign. clusters (uncorrected) and restrict matrix to this range
        maxnum_cluster = 10; %Across subject max cluster number estimate
        %     maxnum_cluster = 0; %Individual determination doesn't work since different across subjects
        %     for i_elec = 1:nSensors_sel
        %         temp_signcluster = [];
        %         temp_signcluster = ...
        %             sum(CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_pval ...
        %             < param.pval_plotting);
        %         if temp_signcluster > maxnum_cluster
        %            maxnum_cluster = temp_signcluster;
        %         end
        %     end
        
        %Store cluster information from all
        Param2plot.per_sub.clusterstat{i_sub, i_TD}           = nan(nSensors_sel, maxnum_cluster);
        Param2plot.per_sub.pval{i_sub, i_TD}                  = nan(nSensors_sel, maxnum_cluster);
        Param2plot.per_sub.pval_derivative{i_sub, i_TD}       = nan(nSensors_sel, maxnum_cluster);
        Param2plot.per_sub.effect_onset{i_sub, i_TD}          = nan(nSensors_sel, maxnum_cluster);
        Param2plot.per_sub.effect_offset{i_sub, i_TD}         = nan(nSensors_sel, maxnum_cluster);
        Param2plot.per_sub.effect_duration{i_sub, i_TD}       = nan(nSensors_sel, maxnum_cluster);
        
        for i_elec = 1:nSensors_sel
            %Determine clusterorder based on minimal p-value
            clusterorder_minp = [];
            [~, clusterorder_minp] = ...
                sort(CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_pval);
            %For each cluster, read out corresponding information and order it according to minpval
            index_cluster_placement = 0;
            for i_cluster = clusterorder_minp
                index_cluster_placement = index_cluster_placement + 1;
                Param2plot.per_sub.clusterstat{i_sub, i_TD}(i_elec, index_cluster_placement) = ... %Clusterstat
                    CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_statSum(i_cluster);
                Param2plot.per_sub.pval{i_sub, i_TD}(i_elec, index_cluster_placement) = ... %pval
                    CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_pval(i_cluster);
                %If there is a valid cluster,
                if ~isnan(CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_pval(i_cluster))
                    %read out cluster timing info
                    Param2plot.per_sub.effect_onset{i_sub, i_TD}(i_elec,index_cluster_placement ) = ... %First sample of cluster to compute relative onset
                        CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_samples{i_cluster}(1) / nSamples(i_TD);
                    Param2plot.per_sub.effect_offset{i_sub, i_TD}(i_elec, index_cluster_placement) = ...%Last sample of cluster to compute relative offset
                        CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_samples{i_cluster}(end) / nSamples(i_TD);
                    Param2plot.per_sub.effect_duration{i_sub, i_TD}(i_elec, index_cluster_placement) = ...%Relative cluster duration in samples
                        length(CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_samples{i_cluster}) / nSamples(i_TD);
                end
            end
        end
        
        %Restrict cluster data to current electrode selection (StimCorr or All)
        for i_elec = 1:nSensors_sel
            if filt_signelecs_StimCorr(ind_signelecs_StimCorr(i_elec)) == 0 %Restrict p-values to selected elecs (StimCorr or All)
                Param2plot.per_sub.pval{i_sub, i_TD}(i_elec, :) = NaN;
            end
        end
        
        %Compute derivative of p-value to estimate strength of effect
        Param2plot.per_sub.pval_derivative{i_sub, i_TD} = ...
            -(log10(Param2plot.per_sub.pval{i_sub, i_TD}));
        
        %Perform FDR_correction on cluster-p-values
        %NOTE: FDR-correction is performed across electrodes for each cluster
        %seperately, since it is much too strict if we treat different clusters
        %as different electrodes/measurements
        cluster_pvalFDR_perelec = [];
        for i_cluster = 1:size(Param2plot.per_sub.pval{i_sub, i_TD},2)
            [~, cluster_critpFDR,~ , cluster_pvalFDR_perelec(:,i_cluster)] = ...
                fdr_bh(Param2plot.per_sub.pval{i_sub, i_TD}(:,i_cluster), param.pval_FDR, 'pdep','no');
        end
        Param2plot.per_sub.pval_FDR{i_sub, i_TD} = cluster_pvalFDR_perelec;
        if param.FDRcorrect == 1
            FDR_label = 'FDRcorr';
            %Replace uncorrected p-value matrix with FDR-corrected p-values
            Param2plot.per_sub.pval{i_sub, i_TD} = cluster_pvalFDR_perelec;
        end
        
        %Read out cluster timecourse for sign. clusters only (to later easily
        %use it for plotting sign. clusters)
        for i_elec = 1:length(Param2plot.per_sub.pval{i_sub, i_TD})
            %Set empty time course as proxy
            Param2plot.per_sub.cluster_timecourse{i_sub, i_TD}(i_elec, :) = ...
                zeros(1,length(CurrentEffect.clusterstat...
                {ind_signelecs_StimCorr(i_elec)}.cluster_timecourse));
            
            for i_cluster = 1:length(Param2plot.per_sub.pval{i_sub, i_TD}(i_elec,:))
                %If a cluster is sign, enter its timecourse in the blank proxy
                %nd denote it with the cluster index
                if Param2plot.per_sub.pval{i_sub, i_TD}(i_elec,i_cluster) < param.pval_plotting
                    Param2plot.per_sub.cluster_timecourse{i_sub, i_TD}...
                        (i_elec, find(CurrentEffect.clusterstat...
                        {ind_signelecs_StimCorr(i_elec)}.cluster_timecourse ...
                        == i_cluster)) = i_cluster;
                end
            end
        end
        
        %Data aggregation over subjects
        %Array to differentiate subject entries
        Param2plot.all_subs.sub_index = ...
            [Param2plot.all_subs.sub_index; ...
            ones(length(coords_sub{i_sub, i_TD}),1)*i_sub];        
        %Array to differentiate TD entries for each electrode
        for i_elec = 1:nSensors_all
            Param2plot.all_subs.TD_index    = ...
                [Param2plot.all_subs.TD_index; ...
                FuncInput_ToneDur_text{i_TD} 's'];
        end
        
        %Clusterinfo aggregation over subjects (Dim: elec * cluster)
        if i_sub == 1 && i_TD == 1
            Param2plot.all_subs.clusterstat     = ...
                Param2plot.per_sub.clusterstat{i_sub, i_TD};
            Param2plot.all_subs.pval            = ...
                Param2plot.per_sub.pval{i_sub, i_TD};
            Param2plot.all_subs.pval_FDR        = ...
                Param2plot.per_sub.pval_FDR{i_sub, i_TD};            
            Param2plot.all_subs.pval_derivative = ...
                Param2plot.per_sub.pval_derivative{i_sub, i_TD};
            Param2plot.all_subs.effect_onset    = ...
                Param2plot.per_sub.effect_onset{i_sub, i_TD};
            Param2plot.all_subs.effect_offset   = ...
                Param2plot.per_sub.effect_offset{i_sub, i_TD};
            Param2plot.all_subs.effect_duration = ...
                Param2plot.per_sub.effect_duration{i_sub, i_TD};
        else
            Param2plot.all_subs.clusterstat     = ...
                [Param2plot.all_subs.clusterstat; ...
                Param2plot.per_sub.clusterstat{i_sub, i_TD}];
            Param2plot.all_subs.pval            = ...
                [Param2plot.all_subs.pval; ...
                Param2plot.per_sub.pval{i_sub, i_TD}];
            Param2plot.all_subs.pval_FDR        = ...
                [Param2plot.all_subs.pval_FDR; ...
                Param2plot.per_sub.pval_FDR{i_sub, i_TD}];
            Param2plot.all_subs.pval_derivative = ...
                [Param2plot.all_subs.pval_derivative; ...
                Param2plot.per_sub.pval_derivative{i_sub, i_TD}];
            Param2plot.all_subs.effect_onset    = ...
                [Param2plot.all_subs.effect_onset; ...
                Param2plot.per_sub.effect_onset{i_sub, i_TD}];
            Param2plot.all_subs.effect_offset   = ...
                [Param2plot.all_subs.effect_offset; ...
                Param2plot.per_sub.effect_offset{i_sub, i_TD}];
            Param2plot.all_subs.effect_duration = ...
                [Param2plot.all_subs.effect_duration; ...
                Param2plot.per_sub.effect_duration{i_sub, i_TD}];
        end
        
        %Cleanup
        usedElecs_chanposIndex = [];
        clear CurrentEffect labels_loadedData corr_ttest DataClean_AllTrials temp*
    end
    disp([' -- done loading in ' num2str(toc) ' sec --'])
end

%% 2) Determine sign. electrodes
SignElecs.array         = any(Param2plot.all_subs.pval_FDR < param.pval_plotting,2); %1D filter denoting sign. elecs (independent of number of clusters)
SignElecs.index         = find(SignElecs.array);
SignElecs.num_elecs     = length(SignElecs.index);
SignElecs.num_cluster   = sum(sum(Param2plot.all_subs.pval_FDR < param.pval_plotting));

for i_elec = 1:length(SignElecs.index)
    SignElecs.onset_minpcluster(i_elec,1) = ...
        Param2plot.all_subs.effect_onset(SignElecs.index(i_elec),1);
    SignElecs.offset_minpcluster(i_elec,1) = ...
        Param2plot.all_subs.effect_offset(SignElecs.index(i_elec),1);
    SignElecs.duration_minpcluster(i_elec,1) = ...
        Param2plot.all_subs.effect_duration(SignElecs.index(i_elec),1);
end

%Determine labels of sign. elecs and add number to subtitle
SignElecs.labels        = labels_allsub(SignElecs.index);
SignElecs.anatlabels    = anatlabels_allsub(SignElecs.index);
SignElecs.fulllabels    = [];
for i_elec = 1:length(SignElecs.labels)
    SignElecs.fulllabels{i_elec,1} = ...
        [SignElecs.labels{i_elec}, ' ', ...
        Param2plot.all_subs.TD_index(SignElecs.index(i_elec),:), ' ', ...
        SignElecs.anatlabels{i_elec}];
end

%Create electrode labels for plotting
for i_elec = 1:length(labels_allsub)%No electrode labels
    labels_plotting_empty{i_elec} = '';
end
counter_elecs           = 0;
labels_plotting_number  = labels_plotting_empty;
% Remove electrode labels
for i_elec = SignElecs.index'
    counter_elecs = counter_elecs +1;
    labels_plotting_number{i_elec} = '';
end

%% 3) Prepare & Plot Figure 1: Clusterstat & p-value derivative
 
%3.1 Prepare figure
h = figure('visible','on'); %ensures figure doesn't pop up during plotting
set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
set(gcf,'renderer','opengl');

DimSubplot = [1, 2];    
CounterSubplot = 1;
CounterSurfplot = 1;
%3.2 Adjust header/title
%3.2 Adjust header/title
switch FuncInput_EffectType
    case 'PredEffect'
        effect_text = ['Effect: Prediction (explain p33 activity by p*34); p < ' num2str(param.pval_plotting) ' ' FDR_label];
    case 'SimplePredErrEffect'
        effect_text = ['Effect: Simple Prediction error (explain p34 activity by absolute p33-p34 difference); p < ' num2str(param.pval_plotting) ' ' FDR_label];
    case 'ComplexPredErrEffect'
        effect_text = ['Effect: Complex Prediction error (explain p34 activity by absolute p*34-p34 difference); p < ' num2str(param.pval_plotting) ' ' FDR_label];
end

Fig_title = {['Group level (n = ' num2str(length(subs)) ')'] ...
    [effect_text] ...
    ['Input data: ' FuncInput_DataType ', Pooled TD']} ;
sgtitle(Fig_title,'FontSize',20,'Interpreter','none')

% 3.3 Adapted plotting for p-value derivative with FDR overlay
for i_elec = 1:length(Param2plot.all_subs.pval_derivative)
    PlotInput(i_elec,1) = max(Param2plot.all_subs.pval_derivative(i_elec,:));
end

% Only plot significant clusters
PlotInput(PlotInput <= 1.3) = NaN;

% Colorlim
clims = [0 3];
Label_Colorbar = '- log10 (cluster p-value) ';

% prepare struct
dat.dimord      = 'chan_time';
dat.time        = 0;
dat.label       = labels_allsub;
dat.avg         = PlotInput;
dat.sign_elecs  = SignElecs.array;

% Sizes
SizeFactor_default = 2;  % small spheres for cluster-significant not FDR
SizeFactor_FDR = 4;     % large spheres for FDR significant

% Prepare sizes vector: default size for all
chanSize = ones(1, length(dat.avg)) * SizeFactor_default;

% Set large size for FDR-significant electrodes
chanSize(SignElecs.array) = SizeFactor_FDR;

% Project all electrodes on left hemisphere (flip x)
coords_plot = coords_allsub;
coords_plot(:,1) = abs(coords_plot(:,1)) * -1;

% Prepare colormap and text color
cmap = parula;
dat.textcolor_rgb = [1 0 0];

% Plot surface with small spheres colored by pval_derivative for cluster-significant electrodes
% BUT we want the coloring only for cluster significant but NOT FDR electrodes
% So mask color values for FDR electrodes to NaN (so they are not colored by cmap)
PlotInput_forcolor = PlotInput;
PlotInput_forcolor(SignElecs.array) = NaN;

SubplotPosition = [0 -0.1 0 0];
    
% Use your plotting function with cluster-only electrodes colored
sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH(...
    coords_plot, labels_plotting_number, PlotInput_forcolor, ~SignElecs.array,...
    chanSize, clims, cmap, dat.textcolor_rgb, ...
    DimSubplot, CounterSubplot, SubplotPosition, [], []);

hold on;

% Overlay large yellow spheres for FDR significant electrodes (if any)
% Plot by subject color
last5 = cellfun(@(x) x(end-4:end), SignElecs.labels, 'UniformOutput', false); % extract last five characters of labels which indicates subject ID
[uniqueTags, ~, tagIndices] = unique(last5, 'stable');
subjID_sigelecs      = nan(size(PlotInput));
subjID_sigelecs(SignElecs.index) = tagIndices;
nSubs = max(tagIndices);            % number of unique subjects
cmap_subjects = lines(nSubs);       % distinct colors per subject

main_ax = gca;  % main axes handle

% Create overlay axes
overlay_ax = axes('Position', get(main_ax, 'Position'), ...
                  'Color', 'none', ...
                  'XLim', main_ax.XLim, ...
                  'YLim', main_ax.YLim, ...
                  'ZLim', main_ax.ZLim, ...
                  'View', main_ax.View, ...
                  'DataAspectRatio', main_ax.DataAspectRatio, ...
                  'Box', 'off', ...
                  'Visible', 'off');

hold(overlay_ax, 'on');

drawnow;  % ensure small spheres are committed

if any(SignElecs.array)
    [sphereX, sphereY, sphereZ] = sphere(30); % increase resolution for better 3D look
    for i = find(SignElecs.array)'
        cx = coords_plot(i,1);
        cy = coords_plot(i,2);
        cz = coords_plot(i,3) + 5;  % vertical offset
        radius = SizeFactor_FDR;
        subjColor = cmap_subjects(subjID_sigelecs(i), :);

        h_sphere = surf(radius*sphereX + cx, ...
                        radius*sphereY + cy, ...
                        radius*sphereZ + cz, ...
                        'FaceColor', subjColor, ...
                        'EdgeColor', 'none', ...      
                        'FaceLighting', 'gouraud', ...
                        'SpecularStrength', 0.5, ...
                        'DiffuseStrength', 0.7, ...
                        'AmbientStrength', 0.4);

        set(h_sphere, 'Clipping', 'off');
        uistack(h_sphere, 'top'); % ensure it’s last
    end
end

sp_handle_surf{CounterSurfplot} = sp_handle_surf_temp.L;
CounterSubplot                  = CounterSubplot + 1;

% Colorbar settings
h = colorbar;
h.Position(4) = h.Position(4) * 0.5;              % Half height
h.Position(2) = h.Position(2) + h.Position(4)/2;  % Move up to center it
h.Position(1) = h.Position(1) + 0.02;             % Move slightly to the right
h.Label.String = Label_Colorbar;
h.FontSize = 16;
caxis(clims);

% --- Adjust existing colorbar (for small spheres) ---
h.Color = [0 0 0]; % black text for clarity
h.Label.String = 'Small spheres: -log10(cluster p-value)';
h.FontSize = 16;
caxis(clims);

% --- Add second colorbar for subject ID colors (large spheres) ---
% Get position of first colorbar for reference
cb_pos = h.Position;

% Define second colorbar size and position (to the right, enough space)
cb_width = 0.1;   % width of new colorbar axis
cb_height = 0.35; % height, enough to hold patches and labels
cb_x = cb_pos(1) + cb_pos(3) + 0.06; % right of first colorbar + some padding
cb_y = cb_pos(2) + (cb_pos(4) - cb_height)/2; % vertically centered relative to first

% Create new axes for the subject ID colorbar
cb_ax = axes('Position', [cb_x, cb_y, cb_width, cb_height], ...
             'Visible', 'off', 'XLim', [0 1], 'YLim', [0 nSubs]);

hold(cb_ax, 'on');

% For each subject, draw a colored rectangle stacked vertically
w = 0.2;          % width of each rectangle (narrower)
h = 1;            % full height so rectangles stack without gaps
x_rect = 0.3;     % x-position of rectangles (move right)
x_text = x_rect + w + 0.05;  % text slightly to the right of rectangles

for i = 1:nSubs
    % Draw rectangles stacked without vertical gaps, shifted right
    rectangle('Position', [x_rect, i - 1, w, h], ...
              'FaceColor', cmap_subjects(i,:), 'EdgeColor', 'none', ...
              'Parent', cb_ax);
    % Draw labels next to rectangles
    text(x_text, i - 0.5, num2str(i), ...
        'Parent', cb_ax, ...
        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', ...
        'FontSize', 12, 'FontWeight', 'bold');
end

% Set Y ticks at center of each patch
yticks(cb_ax, 0.5:1:(nSubs-0.5));
yticklabels(cb_ax, arrayfun(@num2str, 1:nSubs, 'UniformOutput', false));

% Remove X ticks and labels
cb_ax.XTick = [];
cb_ax.XColor = 'none';

% Label colorbar title on top
text(cb_ax, 0.8, nSubs/2 + 0.6, 'Large spheres: Subject IDs', ...
    'HorizontalAlignment', 'center', ...
    'FontWeight', 'normal', ...     % no bold
    'FontSize', 16, ...             % same font size
    'Rotation', 90);                % vertical text

% Reverse Y axis so 1 is at top (optional)
set(cb_ax, 'YDir', 'reverse');

%% Add table with info about electrodes
% subplot(DimSubplot(1), DimSubplot(2), 2)
% textbox_info1 = {...
%     [num2str(length(subs)) ' subjects'] ...
%     [num2str(SignElecs.num_elecs) ' / ' num2str(length(PlotInput)) ' FDR-corr. sign. elec / all elec'] ...
%     [num2str(length(PlotInput(~isnan(PlotInput)))) ' sig. clusters total'] []};
% textbox_info2 = {};
% if length(SignElecs.index) < 30 %one row
%     for i_elec = 1:length(SignElecs.index)
%         sel_elec = SignElecs.index(i_elec);
%         textbox_info2{i_elec} = ...
%             [labels_plotting_number{sel_elec} ' = ' SignElecs.fulllabels{i_elec}] ;
%     end
%     textbox = [textbox_info1 textbox_info2];
%     t = text(0, 0.5, 0, textbox, 'FontSize',12,'Interpreter','none');
% else
%     for i_elec = 1:25 %two rows
%         sel_elec = SignElecs.index(i_elec);
%         textbox_info2{i_elec} = ...
%             [labels_plotting_number{sel_elec} ' = ' SignElecs.fulllabels{i_elec}] ;
%     end
%     for i_elec = 26:length(SignElecs.index)
%         sel_elec = SignElecs.index(i_elec);
%         textbox_info3{i_elec} = ...
%             [labels_plotting_number{sel_elec} ' = ' SignElecs.fulllabels{i_elec}] ;
%     end
%     textbox = [textbox_info1 textbox_info2];
%     t = text(-0.3, 0.5, 0, textbox, 'FontSize',12,'Interpreter','none');
%     t2 = text(0.6, 1, 0, textbox_info3, 'FontSize',12,'Interpreter','none');
% end
% set(gca,'visible','off')
% CounterSubplot = CounterSubplot + 1;
% 
% 
% %% Plot on both hemispheres
% h = figure('visible','on'); %ensures figure doesn't pop up during plotting
% set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
% set(gcf,'renderer','opengl');
% 
% DimSubplot = [1, 2];    
% CounterSubplot = 1;
% CounterSurfplot = 1;
% %3.2 Adjust header/title
% %3.2 Adjust header/title
% switch FuncInput_EffectType
%     case 'PredEffect'
%         effect_text = ['Effect: Prediction (explain p33 activity by p*34); p < ' num2str(param.pval_plotting) ' ' FDR_label];
%     case 'SimplePredErrEffect'
%         effect_text = ['Effect: Simple Prediction error (explain p34 activity by absolute p33-p34 difference); p < ' num2str(param.pval_plotting) ' ' FDR_label];
%     case 'ComplexPredErrEffect'
%         effect_text = ['Effect: Complex Prediction error (explain p34 activity by absolute p*34-p34 difference); p < ' num2str(param.pval_plotting) ' ' FDR_label];
% end
% 
% Fig_title = {['Group level (n = ' num2str(length(subs)) ')'] ...
%     [effect_text] ...
%     ['Input data: ' FuncInput_DataType ', Pooled TD']} ;
% sgtitle(Fig_title,'FontSize',20,'Interpreter','none')
% 
% % 3.3 Adapted plotting for p-value derivative with FDR overlay
% for i_elec = 1:length(Param2plot.all_subs.pval_derivative)
%     PlotInput(i_elec,1) = max(Param2plot.all_subs.pval_derivative(i_elec,:));
% end
% 
% % Only plot significant clusters
% PlotInput(PlotInput <= 1.3) = NaN;
% 
% % Colorlim
% clims = [0 3];
% Label_Colorbar = '- log10 (cluster p-value) ';
% 
% % prepare struct
% dat.dimord      = 'chan_time';
% dat.time        = 0;
% dat.label       = labels_allsub;
% dat.avg         = PlotInput;
% dat.sign_elecs  = SignElecs.array;
% 
% % Sizes
% SizeFactor_default = 2;  % small spheres for cluster-significant not FDR
% SizeFactor_FDR = 4;     % large spheres for FDR significant
% 
% % Prepare sizes vector: default size for all
% chanSize = ones(1, length(dat.avg)) * SizeFactor_default;
% 
% % Set large size for FDR-significant electrodes
% chanSize(SignElecs.array) = SizeFactor_FDR;
% 
% % Project all electrodes on left hemisphere (flip x)
% coords_plot = coords_allsub;
% 
% % Prepare colormap and text color
% cmap = parula;
% dat.textcolor_rgb = [1 0 0];
% 
% % Plot surface with small spheres colored by pval_derivative for cluster-significant electrodes
% % BUT we want the coloring only for cluster significant but NOT FDR electrodes
% % So mask color values for FDR electrodes to NaN (so they are not colored by cmap)
% PlotInput_forcolor = PlotInput;
% PlotInput_forcolor(SignElecs.array) = NaN;
% 
% SubplotPosition = [0 -0.1 0 0];
%     
% % Use your plotting function with cluster-only electrodes colored
% % Plot per hemisphere
% sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_SubplotHemis(...
%     coords_plot, labels_plotting_number, PlotInput_forcolor, ~SignElecs.array,...
%     chanSize, clims, cmap, ...
%     DimSubplot, CounterSubplot, SubplotPosition, [], []);
% 
% hold on;
% 
% % Overlay large yellow spheres for FDR significant electrodes (if any)
% if any(SignElecs.array)
%     % Get left and right hemisphere indices (based on x-coordinate)
%     left_inds = find(coords_plot(:,1) < 0 & SignElecs.array);
%     right_inds = find(coords_plot(:,1) > 0 & SignElecs.array);
%     
%     % Sphere parameters
%     [sphereX, sphereY, sphereZ] = sphere(20);
%     radius = SizeFactor_FDR;
%     z_shift = 5; % your strong z-shift
%     
%     % ----- Plot spheres on LEFT hemisphere -----
%     axes(sp_handle_surf_temp.L);  % Switch to left subplot
%     hold on;
%     for i = left_inds'
%         cx = coords_plot(i,1);
%         cy = coords_plot(i,2);
%         cz = coords_plot(i,3) + z_shift;
%         surf(radius*sphereX + cx, ...
%              radius*sphereY + cy, ...
%              radius*sphereZ + cz, ...
%              'FaceColor', 'yellow', ...
%              'EdgeColor', 'none', ...
%              'FaceAlpha', 1, ...
%              'FaceLighting', 'gouraud', ...
%              'SpecularStrength', 0.3, ...
%              'DiffuseStrength', 0.6, ...
%              'AmbientStrength', 0.4, ...
%              'Clipping', 'off');
%     end
% 
%     % ----- Plot spheres on RIGHT hemisphere -----
%     axes(sp_handle_surf_temp.R);  % Switch to right subplot
%     hold on;
%     for i = right_inds'
%         cx = coords_plot(i,1);
%         cy = coords_plot(i,2);
%         cz = coords_plot(i,3) + z_shift;
%         surf(radius*sphereX + cx, ...
%              radius*sphereY + cy, ...
%              radius*sphereZ + cz, ...
%              'FaceColor', 'yellow', ...
%              'EdgeColor', 'none', ...
%              'FaceAlpha', 1, ...
%              'FaceLighting', 'gouraud', ...
%              'SpecularStrength', 0.3, ...
%              'DiffuseStrength', 0.6, ...
%              'AmbientStrength', 0.4, ...
%              'Clipping', 'off');
%     end
% end
% 
% %sp_handle_surf{CounterSurfplot} = sp_handle_surf_temp.L;
% CounterSubplot                  = CounterSubplot + 1;
% 
% % Colorbar settings
% h = colorbar;
% h.Position(4) = h.Position(4) * 0.5;              % Half height
% h.Position(2) = h.Position(2) + h.Position(4)/2;  % Move up to center it
% h.Position(1) = h.Position(1) + 0.05;             % Move slightly to the right
% h.Label.String = Label_Colorbar;
% h.FontSize = 16;
% caxis(clims);



%% 4 Save figure
if save_poststepFigs == 1
    filename     = ['Surf1HSignElec' param.ElecSelect ...
        '_Allsubn' num2str(length(subs)) '_' ...
        FuncInput_EffectType '_' FuncInput_DataType ...
        '_AllTD_Stat' FDR_label '_FDRBySubjColor.png'];
    figfile      = [path_fig filename];
    saveas(gcf, figfile, 'png'); %save png version
    close all;
end

end