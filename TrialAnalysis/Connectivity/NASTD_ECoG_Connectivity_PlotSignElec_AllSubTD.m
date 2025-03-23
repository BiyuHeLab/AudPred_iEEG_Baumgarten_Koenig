function NASTD_ECoG_Connectivity_PlotSignElec_AllSubTD...
    (subs, ...
    FuncInput_EffectType, FuncInput_DataType, FuncInput_ToneDur_text, ...
    param, ...
    save_poststepFigs, paths_NASTD_ECoG)

%Aim: 

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

path_fig = ([paths_NASTD_ECoG.ECoGdata_Connectivity ...
    'ElecSelect/Allsub_n' num2str(length(subs)) ...
    '/Figs/' FuncInput_EffectType '/']);
if (~exist(path_fig, 'dir')); mkdir(path_fig); end


%% 1) Load prediction effect data and aggregate relevant info across subjects
clear Param2plot
coords_allsub                       = [];
labels_allsub                       = [];
anatlabels_allsub                   = [];
Param2plot.all_subs.sub_index       = [];
Param2plot.all_subs.TD_index        = [];
Param2plot.all_subs.label_AnatCat   = [];
Param2plot.all_subs.index_AnatCat   = [];
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
        
        %Categorize electrodes according to anatomical regions
        AnatReg_allSubs{i_sub, i_TD} = ...
            NASTD_ECoG_AssignAnatRegions(DataClean_AllTrials, labels_loadedData);
        
        %Aggregate category labels and indices across subjects
        Param2plot.all_subs.label_AnatCat  = ...
            [Param2plot.all_subs.label_AnatCat; ...
            AnatReg_allSubs{i_sub, i_TD}.Info_perelec(:,3)];
        Param2plot.all_subs.index_AnatCat  = ...
            [Param2plot.all_subs.index_AnatCat; ...
            AnatReg_allSubs{i_sub, i_TD}.CatIndex];
        
        %Determine max. number of sign. clusters (uncorrected) and restrict matrix to this range
        maxnum_cluster = 10; %Across subject max cluster number estimate
        
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
        clear CurrentEffect labels_loadedData corr_ttest DataClean_AllTrials temp*
    end
    disp([' -- done loading in ' num2str(toc) ' sec --'])
end

%% 2) Determine sign. electrodes
SignElecs.array         = any(Param2plot.all_subs.pval < param.pval_plotting,2); %1D filter denoting sign. elecs (independent of number of clusters)
SignElecs.index         = find(SignElecs.array);
SignElecs.num_elecs     = length(SignElecs.index);
SignElecs.num_cluster   = sum(sum(Param2plot.all_subs.pval < param.pval_plotting));

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
for i_elec = SignElecs.index'
    counter_elecs = counter_elecs +1;
    labels_plotting_number{i_elec} = num2str(counter_elecs);
end

%% 3) Prepare & Plot Figure 1: Clusterstat & p-value derivative

%3.1 Prepare figure
h = figure('visible','on'); %ensures figure doesn't pop up during plotting
set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
set(gcf,'renderer','opengl');

DimSubplot      = [2,2];
SizeFactor      = 3;

%3.2 Adjust header/title
%3.2 Adjust header/title
switch FuncInput_EffectType
    case 'PredEffect'
        effect_text = ['Effect: Prediction (explain p33 activity by p*34); p < ' num2str(param.pval_plotting) ' ' FDR_label];
    case 'ComplexPredErrEffect'
        effect_text = ['Effect: Complex Prediction error (explain p34 activity by absolute p*34-p34 difference); p < ' num2str(param.pval_plotting) ' ' FDR_label];
end

Fig_title = {['Group level (n = ' num2str(length(subs)) ') - Electrodes for Connectivity Analysis'] ...
    [effect_text] ...
    ['Input data: ' FuncInput_DataType ', Pooled TD']} ;
sgtitle(Fig_title,'FontSize',20,'Interpreter','none')

%Project all electrodes on one hemisphere,
coords_allsub(:,1) = abs(coords_allsub(:,1)) * -1;

%Color coding: Subject-Number
for i_elec = 1:length(Param2plot.all_subs.pval_derivative)
    PlotInput(i_elec,1) = max(Param2plot.all_subs.sub_index(i_elec,:));
end

%Colorlim
clims = [1 9];
Label_Colorbar = 'Subject Number';

%prepare struct
dat.dimord          = 'chan_time';
dat.time            = 0;
dat.label           = labels_allsub;
dat.avg             = PlotInput;
dat.sign_elecs      = SignElecs.array;
dat.textcolor_rgb   = [1 0 0];

chanSize    = ones(1,length(dat.avg))*SizeFactor; %electrode size (arbitrary)
cmap        = distinguishable_colors(9);

%Plot surface without connectivity markers
sp_handle_surf_temp1 = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH...
    (coords_allsub, labels_plotting_empty, dat.avg, SignElecs.array, ...
    chanSize, clims, cmap, dat.textcolor_rgb, ...
    DimSubplot, 1, [0 -0.1 0 0], [], []);
%     %Enlarge subplot
%     sp_handle_surf_temp1.L.Position(3) =  sp_handle_surf_temp1.L.Position(3) * 1.1;
%     sp_handle_surf_temp1.L.Position(4) =  sp_handle_surf_temp1.L.Position(4) * 1.1;

%Colorbar
ColorbarPosition    = [0.05 0.4 0 0.6]; %Left of surface
h = colorbar;
h.Position(1)       = ColorbarPosition(1); %sets colorbar to the right
h.Position(2)       = ColorbarPosition(2); %sets colorbar higher
h.Position(4)       = ColorbarPosition(4)*0.8; %shortens colorbar
h.Label.String      = Label_Colorbar;
h.FontSize          = 16;
caxis(clims)

%Plot surface with connectivity markers (lines between electrodes)
sp_handle_surf_temp2 = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH_Connect...
    (coords_allsub, labels_plotting_empty, dat.avg, SignElecs.array, ...
    Param2plot.all_subs.sub_index, Param2plot.all_subs.index_AnatCat, ...
    chanSize, clims, cmap, dat.textcolor_rgb, ...
    DimSubplot, 2, [0 -0.1 0 0], [], []);
%     %Enlarge subplot
%     sp_handle_surf_temp2.Position(3) =  sp_handle_surf_temp2.Position(3) * 1.1;
%     sp_handle_surf_temp2.Position(4) =  sp_handle_surf_temp2.Position(4) * 1.1;


%Compute number of electrode-connections per anat. parcellation per subject
%Determine anat parcel for all sign. elecs per subject
temp_anatparcel_signelecs = zeros(length(subs),10);
for i_sub = 1:length(subs)
    temp_index_signelecs = [];    
    %Determine sign. elecs for current subject
    temp_index_signelecs = find(Param2plot.all_subs.sub_index(SignElecs.array) == i_sub);
    Param2plot.all_subs.label_AnatCat(SignElecs.index(temp_index_signelecs))
    temp_index_anatcat = ...
    min(Param2plot.all_subs.index_AnatCat(SignElecs.index(temp_index_signelecs),:),[],2)';
%Note: min means that we read out the coarsest parcellation possible    
    for i_anatcat = unique(temp_index_anatcat)
        temp_anatparcel_signelecs(i_sub,i_anatcat) = ...
            sum(min(Param2plot.all_subs.index_AnatCat...
            (SignElecs.index(temp_index_signelecs),:),[],2)' == i_anatcat);
    end   
end
SignElecs_perAnatParcel = table({'NY688';'NY704';'NY708';'NY723';'NY742';'NY751';'NY787';'NY794';'NY798'}, ...
    temp_anatparcel_signelecs(:,1),temp_anatparcel_signelecs(:,2),temp_anatparcel_signelecs(:,3), ...
    temp_anatparcel_signelecs(:,4),temp_anatparcel_signelecs(:,5),temp_anatparcel_signelecs(:,6), ...
    temp_anatparcel_signelecs(:,7),temp_anatparcel_signelecs(:,8),temp_anatparcel_signelecs(:,9), ...
    temp_anatparcel_signelecs(:,10));
SignElecs_perAnatParcel.Properties.VariableNames = ...
    {'Subject', AnatReg_allSubs{1,1}.CatLabels{1:end}};

%Determine number of cross-anat parcel connections per anat_parcel and subject 
for i_sub = 1:length(subs)
    %Count possible connections per anat parcel to electrodes in other parcels
    temp_numconnections_peranatcat = ...
        zeros(length(AnatReg_allSubs{1,1}.CatLabels),length(AnatReg_allSubs{1,1}.CatLabels));    
    for i_curranatcat = 1:length(AnatReg_allSubs{1,1}.CatLabels)
        for i_connectanatcat = 1:length(AnatReg_allSubs{1,1}.CatLabels)
            temp_numconnections_peranatcat(i_curranatcat,i_connectanatcat) = ...
            temp_anatparcel_signelecs(i_sub,i_curranatcat) * ...
            temp_anatparcel_signelecs(i_sub,i_connectanatcat);
        end
    end
    %NaN diagonals (i.e., connectivity within anat parcels) and redundant matrix half 
    temp_numconnections_peranatcat(1:1+size(temp_numconnections_peranatcat,1):end) ...
        = NaN;
    for i_curranatcat = 1:length(temp_numconnections_peranatcat)
        temp_numconnections_peranatcat(...
            (i_curranatcat + 1):length(temp_numconnections_peranatcat),i_curranatcat) ...
            = NaN;
    end
    
    %Place results in table
    Connections_perAnatParcel{i_sub} = table({AnatReg_allSubs{1,1}.CatLabels{1:end}}', ...
        temp_numconnections_peranatcat(:,1),temp_numconnections_peranatcat(:,2),temp_numconnections_peranatcat(:,3), ...
        temp_numconnections_peranatcat(:,4),temp_numconnections_peranatcat(:,5),temp_numconnections_peranatcat(:,6), ...
        temp_numconnections_peranatcat(:,7),temp_numconnections_peranatcat(:,8),temp_numconnections_peranatcat(:,9), ...
        temp_numconnections_peranatcat(:,10));       
    Connections_perAnatParcel{i_sub}.Properties.VariableNames = ...
        {' ', AnatReg_allSubs{1,1}.CatLabels{1:end}};  
end
clear temp*


%Text listing electrode connections per subject
subplot(DimSubplot(1), DimSubplot(2), 3)
textbox_info1 = {... %Summary across subjects
    [num2str(length(subs)) ' subjects'] ...
    [num2str(SignElecs.num_elecs) ' / ' num2str(length(SignElecs.array)) ' sign. elec / all elec'] ...
    [num2str(SignElecs.num_cluster) ' sign. clusters total']};
t = text(-0.3, 0.6, 0, textbox_info1, 'FontSize',10,'Interpreter','none','FontWeight', 'bold');

%Electrode connections
ypos = [1.1 0.5 1.1 0.5 1.1 0.5 1.1 0.5 1.1];
for i_sub = 1:length(subs) %1 column per subject
    textbox_info_sub{i_sub} = {};
    sel_elecs = find(Param2plot.all_subs.sub_index(SignElecs.index) == i_sub)'; 
    text_connect = {};
    
    text_anatcat = [];
    sel_elecs_anatcat = ...
        min(Param2plot.all_subs.index_AnatCat(SignElecs.index(sel_elecs),:),[],2);    
    for i_anatcat = unique(sel_elecs_anatcat)'
        sum(sel_elecs_anatcat == i_anatcat);
        text_anatcat = ...
            [text_anatcat ' ' num2str(sum(sel_elecs_anatcat == i_anatcat)) ...
            ' ' AnatReg_allSubs{1,1}.CatLabels{i_anatcat} ';'];
    end 
    
    i_entry = 0;
    for i_curranatcat = 1:length(AnatReg_allSubs{1,1}.CatLabels)
        for i_connectanatcat = 1:length(AnatReg_allSubs{1,1}.CatLabels)
            if Connections_perAnatParcel{i_sub}{i_curranatcat, i_connectanatcat+1} > 0
                i_entry = i_entry + 1;
                text_connect{i_entry} = [
                   AnatReg_allSubs{1,1}.CatLabels{i_curranatcat} '-' ...
                   AnatReg_allSubs{1,1}.CatLabels{i_connectanatcat} ': ' ...
                   num2str(Connections_perAnatParcel{i_sub}{i_curranatcat, i_connectanatcat+1}) ...
                   ' connections'];                    
            end
        end
    end 
    
    textbox_info_sub{i_sub} = ...
        {[subs{i_sub}] [num2str(length(sel_elecs)) ' sign. elecs'] ...
        text_anatcat [] text_connect{:}};
    t = text((i_sub*0.28)-0.3, ypos(i_sub), 0, textbox_info_sub{i_sub}, ...
        'FontSize',8,'Interpreter','none', 'FontWeight', 'bold',...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
end
% %Single electrodes
% ypos = [1.1 0.5 1.1 0.5 1.1 0.5 1.1 0.5 1.1];
% for i_sub = 1:length(subs) %1 column per subject
%     textbox_info_sub{i_sub} = {};
%     sel_elecs = find(Param2plot.all_subs.sub_index(SignElecs.index) == i_sub)';
%     
%     sel_elecs_anatcat = min(Param2plot.all_subs.index_AnatCat(SignElecs.index(sel_elecs),:),[],2);
%     text_anatcat = [];
%     for i_anatcat = unique(sel_elecs_anatcat)'
%         sum(sel_elecs_anatcat == i_anatcat);
%         text_anatcat = ...
%             [text_anatcat ' ' num2str(sum(sel_elecs_anatcat == i_anatcat)) ...
%             ' ' AnatReg_allSubs{1,1}.CatLabels{i_anatcat} ';'];
%     end
%     
%     elec_label = [];
%     for i_elec = 1:length(sel_elecs)
%         elec_label{i_elec} = ...
%             [SignElecs.labels{sel_elecs(i_elec)} ' ' ...
%             Param2plot.all_subs.label_AnatCat{SignElecs.index(sel_elecs(i_elec))}];
%     end
%     
%     textbox_info_sub{i_sub} = ...
%         {[subs{i_sub} ' - ' num2str(length(sel_elecs)) ' sign. elecs'] ...
%         text_anatcat []...
%         elec_label{:}};
%     t = text((i_sub*0.28)-0.3, ypos(i_sub), 0, textbox_info_sub{i_sub}, ...
%         'FontSize',8,'Interpreter','none', 'FontWeight', 'bold',...
%         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
% end
set(gca,'visible','off')


%% 4 Save figure
if save_poststepFigs == 1
    filename     = ['Surf1H_ConnectSignElec_' param.ElecSelect ...
        '_Allsubn' num2str(length(subs)) '_' ...
        FuncInput_EffectType '_' FuncInput_DataType ...
        '_AllTD_p' num2str(param.pval_plotting) '_' FDR_label '.png'];
    figfile      = [path_fig filename];
    saveas(gcf, figfile, 'png'); %save png version
    close all;
end

end