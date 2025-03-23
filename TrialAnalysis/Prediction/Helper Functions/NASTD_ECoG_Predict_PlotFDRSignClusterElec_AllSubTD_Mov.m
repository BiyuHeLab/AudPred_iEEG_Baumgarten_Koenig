function NASTD_ECoG_Predict_PlotFDRSignClusterElec_AllSubTD_Mov...
    (subs, ...
    FuncInput_EffectType, FuncInput_DataType, FuncInput_ToneDur_text, ...
    param, ...
    save_poststepFigs, paths_NASTD_ECoG)

%Aim: Plot electrodes showing sign. cluster-corrected prediction effects
%per input data & TD condition aggregated over subjects and tone durations 
%(TD) and projected on 1 Hemisphere. 
%Use only electrodes that show sign. Stim-Corr effect.
%Plot effects with temporal resolution, save single-figure output, and
%create movie from single images

%Based on sample-wise analysis, thus no time windows in analysis.
%Time windows added for movie plotting

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

path_vid = ([paths_NASTD_ECoG.ECoGdata_Prediction ...
    '/PredEffects/Allsub_n' num2str(length(subs)) ...
    '/Figs/PredEffects_Surf/' FuncInput_EffectType ...
    '/Vid/' param.ElecSelect '/' FDR_label '/']);
if (~exist(path_vid, 'dir')); mkdir(path_vid); end
path_frames = ([path_vid 'Frames/']);
if (~exist(path_frames, 'dir')); mkdir(path_frames); end

nFrames = 100;


%% 1) Load prediction effect data and aggregate relevant info across subjects
clear Param2plot
coords_allsub                       = [];
labels_allsub                       = [];
anatlabels_allsub                   = [];
Param2plot.all_subs.sub_index       = [];
Param2plot.all_subs.TD_index        = [];
Param2plot.all_subs.filt_StimCorr   = [];
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
        if param.FDRcorrect == 1
            FDR_label = 'FDRcorr';
            cluster_pvalFDR_perelec = [];
            for i_cluster = 1:size(Param2plot.per_sub.pval{i_sub, i_TD},2)
                [~, cluster_critpFDR,~ , cluster_pvalFDR_perelec(:,i_cluster)] = ...
                    fdr_bh(Param2plot.per_sub.pval{i_sub, i_TD}(:,i_cluster), param.pval_FDR, 'pdep','no');
            end
            %Replace uncorrected p-value matrix with FDR-corrected p-values
            Param2plot.per_sub.pval{i_sub, i_TD} = cluster_pvalFDR_perelec;
        else
            FDR_label = 'uncorr';
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
        
        %Array to differentiate subject and TD entries for each electrode
        Param2plot.all_subs.sub_index   = ...
            [Param2plot.all_subs.sub_index; ...
            ones(length(coords_sub{i_sub}),1)*i_sub];
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
        clear CurrentEffect labels_loadedData corr_ttest DataClean_AllTrials
    end
    disp([' -- done loading in ' num2str(toc) ' sec --'])      
end

%% 2) Determine sign. electrodes
SignElecs.array         = any(Param2plot.all_subs.pval < param.pval_plotting,2); %1D filter denoting sign. elecs (independent of number of clusters)
SignElecs.index         = find(SignElecs.array);
SignElecs.num_elecs     = length(SignElecs.index);
SignElecs.num_cluster   = sum(sum(Param2plot.all_subs.pval < param.pval_plotting));

% %Option 1: Determine when minpval cluster is sign. 
% for i_elec = 1:length(SignElecs.index)
%     SignElecs.onset_minpcluster(i_elec,1) = ...
%         Param2plot.all_subs.effect_onset(SignElecs.index(i_elec),1);
%     SignElecs.offset_minpcluster(i_elec,1) = ...
%         Param2plot.all_subs.effect_offset(SignElecs.index(i_elec),1);
%     SignElecs.duration_minpcluster(i_elec,1) = ...
%         Param2plot.all_subs.effect_duration(SignElecs.index(i_elec),1);
% end

%Option 2: Determine when all clusters are sign.
%Create elec*samplepoints (relative - common across TD) matrix
%Mark sign. clusters per electrode
TimingMat_signcluster       = false(length(Param2plot.all_subs.pval), nFrames);
SignElecs.array_allcluster  = find(Param2plot.all_subs.pval < param.pval_plotting);
[SignElecs.row_allcluster, SignElecs.column_allcluster] = ...
    find(Param2plot.all_subs.pval < param.pval_plotting);

for i_signcluster = 1:SignElecs.num_cluster
        %Find electrode of current cluster
        curr_elec = SignElecs.row_allcluster(i_signcluster);
%         labels_allsub(curr_elec)
        %Find cluster number of current cluster
        curr_cluster = SignElecs.column_allcluster(i_signcluster);        
        %Find nearest entry in timing matrix that corresponds toonset/offset
        curr_onsetsample = ...
            ceil(round(Param2plot.all_subs.effect_onset(curr_elec, curr_cluster),2) * nFrames);
            if curr_onsetsample == 0
                curr_onsetsample = 1;
            end
        curr_offsetsample = ...
            ceil(round(Param2plot.all_subs.effect_offset(curr_elec, curr_cluster),2) * nFrames);
        %Fill corresponding entries in timing matrix
        TimingMat_signcluster(curr_elec,curr_onsetsample:curr_offsetsample) = true;        
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
for i_elec = SignElecs.index'
    counter_elecs = counter_elecs +1;
    labels_plotting_number{i_elec} = num2str(counter_elecs);
end

%% 3) Prepare & Plot frames for each TW
if ~isempty(SignElecs.index)
    
    for i_frame = 1:nFrames
        
        %Set up frame
        h = figure('visible','off'); %ensures figure doesn't pop up during plotting
        % set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
        set(gcf,'renderer','opengl');
        DimSubplot      = [1 1];
        CounterSubplot  = 1;
        CounterSurfplot = 1;
        SizeFactor      = 3;
        SubplotPosition = [0.05 -0.02 0 0];
        
        %Determine to-be-plotted electrodes based on cluster-timing
%         %Option 1: Only minpval sign. cluster
%         dat.sign_elecs = false(length(Param2plot.all_subs.pval_derivative),1);
%         dat.sign_elecs(SignElecs.index(...
%             find(SignElecs.onset_minpcluster <= (1/nFrames * i_frame) & ...
%             SignElecs.offset_minpcluster >= (1/nFrames * (i_frame + 1))))) = true;
        
        %Option 2: all sign. cluster
        dat.sign_elecs = TimingMat_signcluster(:,i_frame);
        
        %Determine plotting parameters
        for i_elec = 1:length(Param2plot.all_subs.pval_derivative)
            PlotInput(i_elec,1) = Param2plot.all_subs.pval_derivative(i_elec,1);
            %select custer based on maximal clusterstat for each elctrode
        end
        
        %Colorlim
        dv_max                  = max(max(Param2plot.all_subs.pval_derivative));
        dv_min                  = min(min(Param2plot.all_subs.pval_derivative));
        dv_absmax               = max([abs(dv_max) abs(dv_min)]);
        %     clims_maxclusterstat    = [-dv_absmax dv_absmax];
        clims_maxclusterstat    = [0 3];
        clims                   = clims_maxclusterstat;
        
        %Data struct
        dat.dimord          = 'chan_time';
        dat.time            = 0;
        dat.label           = labels_allsub;
        dat.avg             = PlotInput;
        dat.textcolor_rgb   = [1 0 0];
        
        chanSize    = ones(1,length(dat.avg))*SizeFactor; %electrode size (arbitrary)
        cmap        = 'parula';
        
        %Project all electrodes on one hemisphere,
        coords_allsub(:,1) = abs(coords_allsub(:,1)) * -1;
        
        sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH...
            (coords_allsub, labels_plotting_empty, dat.avg, dat.sign_elecs,...
            chanSize, clims, cmap, dat.textcolor_rgb, ...
            DimSubplot, CounterSubplot, SubplotPosition, [], []);
        
        sp_handle_surf{CounterSurfplot} = sp_handle_surf_temp.L;
        
        %Set title
        t = title({['Time relative to tone onset:'],...
            [num2str(1/nFrames * i_frame) ' - ' ...
            num2str(1/nFrames * (i_frame+1))]});
        t.FontSize = 12;
        t.FontWeight = 'normal';
        %     t.Position(3) = t.Position(3) + (t.Position(3)/25);
        
        %Set Colorbar
        Label_Colorbar          = '- log10 (cluster p-value)';
        ColorbarPosition    = [0.1 0.2 0 0.5]; %Left of surface
        h                   = colorbar;
        h.Position(1)       = ColorbarPosition(1); %sets colorbar to the right
        h.Position(2)       = ColorbarPosition(2); %sets colorbar higher
        h.Position(4)       = h.Position(4)*0.75; %makes colorbar shorter
        h.Label.String      = Label_Colorbar;
        h.FontSize          = 12;
        caxis(clims)
        
        disp(['Processing Frame ' num2str(i_frame)]);
        
        %Save frame
        filename     = ['Frame' num2str(i_frame) '.png'];
        figfile      = [path_frames filename];
        saveas(gcf, figfile, 'png'); %save png version
        close all;
    end
    
    %% Aggregate frames to movie
    cd(path_vid);
    writerObj = VideoWriter(...
        ['Video_' FuncInput_EffectType '_' ...
        FuncInput_DataType '_PooledTD_' ...
        param.ElecSelect '_' ...
        FDR_label '.avi']);
    
    writerObj.FrameRate = 10;
    open(writerObj)
    
    for i_frame = 1:(nFrames-1)
        frame = sprintf([path_frames 'Frame' num2str(i_frame) '.png']);
        input = imread(frame);
        
        writeVideo(writerObj, input);
    end
    close(writerObj);
    
end

end