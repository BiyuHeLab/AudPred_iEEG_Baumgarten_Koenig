addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
%Add project base path
NASTD_ECoG_setVars
paths_NASTD_ECoG = NASTD_ECoG_paths;
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
%Add project base and script dir

%Determine subjects
sub_list = vars.sub_list;
subs = sub_list(vars.validSubjs);

%Load in file with individual preproc infos
subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;
ToneDur_text = {'0.2' '0.4'};

plot_poststepFigs = 0;
save_poststepFigs = 1;
addpath('/isilon/LFMI/VMdrive/Thomas/toolboxes/mvgc_v1.3/');
load('/isilon/LFMI/VMdrive/Thomas//NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Connectivity/ElecSelect/Allsub_n9/SelElecs_p0.05uncorr.mat') %p < 0.05 thresh,uncorr elec selection    

param.GC            = [];
param.GC.fs         = 512;
param.GC.downsample = 0; %Cave: Downsampling lead to problems with tone sample selection
    if param.GC.downsample == 0
        param.GC.newfs  = param.GC.fs;
    else
        param.GC.newfs  = param.GC.fs/param.GC.downsample;
    end
param.GC.nvars      = 2;
param.GC.regmode    = 'OLS';   % VAR model estimation regression mode ('OLS', 'LWR' or empty for default
param.GC.maxmorder  = 50;   % maximum model order for model order estimation, rule of thumb = number of samples per input data snippet, but not > 100
% param.GC.morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)

param.GC.tstat      = 'F'; % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
param.GC.alpha      = 0.05;   % significance level for significance test
param.GC.mhtc       = 'FDRD'; % multiple hypothesis test correction (see routine 'significance')

param.GC.ElecPairEffect   = {'Pred_Pred','PE_PE','Pred_PE','PE_Pred'};
param.GC.ElecPairRegion   = 'AllRegions';
param.GC.InputDataType    = {'Broadband'}; %{'Broadband', 'HP05toLP30Hz', 'HighGamma_LogAmp'}; 

param.GC.tone_aggregation = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};

param.GC.tone_aggregation_label = cellfun(@(c)[c],param.GC.tone_aggregation);
param.GC.tone_aggregation_label = num2str((param.GC.tone_aggregation_label),'_%d');
param.GC.tone_aggregation_label = erase(param.GC.tone_aggregation_label, ' ');

param.GC.epochdur_ms = 'full'; %'full', 200, 100
if strcmp(param.GC.epochdur_ms, 'full')
    param.GC.epochdur_ms_label = 'fullTW';
else
   param.GC.epochdur_ms_label = [num2str(param.GC.epochdur_ms) 'msTW'];
end

AnatReg_CatLabels = ... %Select which anatomical regions should be plotted
    {'AntPFC'; ...
    'VentralT'; ...
    'SupParLob'};

%Plot group-level aggregated GC results for each anatomical parcel pairing
%for each effect type, resoectively. Focused on comparisons between p1 to
%p31. %Test slope of linear fit vs. 0

set_visibility  = 'on';
Yaxis_scaling   = 'true';
digits(6) %show 6 digits after decimal point

%% 1) Set paths &  load input data
% Load data for all tones from 1 to 32
path_inputdata = ([paths_NASTD_ECoG.ECoGdata_Connectivity 'GC/']);
label_inputdata = ['GCdata_n9_' ...
    param.GC.ElecPairRegion '_' ...
    param.GC.InputDataType{1} param.GC.tone_aggregation_label '_' ...
    param.GC.epochdur_ms_label '_ensemnorm'];
GCdata = load([path_inputdata label_inputdata], 'GCdata');
GCdata = GCdata.GCdata;

%% 2) Determine and read out anatomical regions
for i_effect = 1:length(param.GC.ElecPairEffect)

    %2.0 Determine electrode selection labels
    if strcmp(param.GC.ElecPairEffect{i_effect}, 'Pred_Pred')
        title_label1   = 'Pred2Pred';
    elseif strcmp(param.GC.ElecPairEffect{i_effect}, 'PE_PE')
        title_label1   = 'PE2PE';
    elseif strcmp(param.GC.ElecPairEffect{i_effect}, 'Pred_PE')
        title_label1   = 'Pred2PE';
    elseif strcmp(param.GC.ElecPairEffect{i_effect}, 'PE_Pred')
        title_label1   = 'PE2Pred';
    end

    if strcmp(param.GC.ElecPairRegion, 'AllRegions')
        title_label2   = 'AllReg';
    elseif strcmp(param.GC.ElecPairRegion, 'Frontal_Temporal')
        title_label2   = 'Front2Temp';
    elseif strcmp(param.GC.ElecPairRegion, 'Temporal_Frontal')
        title_label2   = 'Temp2Front';
    end

    subfield_label1 = extractBefore(param.GC.ElecPairEffect{i_effect}, '_');    

    %2.1 Read out detailed electrode label 
    %and correspoding anatomical region label
    %for each source & target elecs per subject
    clear fulllabel* anatlabel*
    for i_sub = 1:length(subs)
        if ~isempty(GCdata{i_sub}.temporalGC{i_effect})
            clear temp_*
            %Read out anatomical region labels for source elecs per subject
            for i_sourceelec = 1:length(GCdata{i_sub}.ind_pairedelecs_from{i_effect})

                temp_eleclabel_sourceelec = ...
                    GCdata{i_sub}.label_allelecs{GCdata{i_sub}.ind_pairedelecs_from{i_effect}(i_sourceelec)};
                temp_sublabel = subs{i_sub};
                temp_filt_elec = strncmp(fields(...
                    SelElecs.Pairs.(param.GC.ElecPairEffect{i_effect}).(param.GC.ElecPairRegion).Pairs_perelec), ...
                    [temp_eleclabel_sourceelec '_' temp_sublabel], ...
                    length(temp_eleclabel_sourceelec) + length(temp_sublabel) + 1);
                temp_anatlabel_sourceelec = ...
                    fields(SelElecs.Pairs.(param.GC.ElecPairEffect{i_effect}).(param.GC.ElecPairRegion).Pairs_perelec);
                temp_anatlabel_sourceelec = temp_anatlabel_sourceelec{temp_filt_elec};

                fulllabel_sourceelec{i_sub}{i_sourceelec, 1} =...
                    temp_anatlabel_sourceelec; %Reduced anat label
                temp_anatlabel_sourceelec2 = ...
                    cell2mat(SelElecs.([subfield_label1 'Effect']){temp_filt_elec,5});
                anatlabel_sourceelec{i_sub}{i_sourceelec, 1} = ...
                    strrep(temp_anatlabel_sourceelec2, ', ', '_'); %Full anat label

                %Read out anatomical region label for target elecs (per source elec)

                temp_num_targetelecs = size(SelElecs.Pairs.(... %Number of target elecs for current source elec
                    param.GC.ElecPairEffect{i_effect}).(...
                    param.GC.ElecPairRegion).Pairs_perelec.(...
                    temp_anatlabel_sourceelec),1);

                for i_targetelec = 1:temp_num_targetelecs
                    fulllabel_targetelec{i_sub}{i_targetelec, i_sourceelec} = ...
                        cell2mat(strcat(...
                        SelElecs.Pairs.(param.GC.ElecPairEffect{i_effect}).(...
                        param.GC.ElecPairRegion).Pairs_perelec.(...
                        temp_anatlabel_sourceelec){i_targetelec,1}, '_', ... %Elec label
                        SelElecs.Pairs.(param.GC.ElecPairEffect{i_effect}).(...
                        param.GC.ElecPairRegion).Pairs_perelec.(...
                        temp_anatlabel_sourceelec){i_targetelec,2}, '_', ... %Sub label
                        SelElecs.Pairs.(param.GC.ElecPairEffect{i_effect}).(...
                        param.GC.ElecPairRegion).Pairs_perelec.(...
                        temp_anatlabel_sourceelec){i_targetelec,5})); %Anat label
                    fulllabel_targetelec{i_sub}{i_targetelec, i_sourceelec} = ...
                        strrep(...
                        fulllabel_targetelec{i_sub}{i_targetelec, i_sourceelec}, ...
                        ', ', '_');
                    anatlabel_targetelec{i_sub}{i_targetelec, i_sourceelec} = ...
                        cell2mat(...
                        SelElecs.Pairs.(param.GC.ElecPairEffect{i_effect}).(...
                        param.GC.ElecPairRegion).Pairs_perelec.(...
                        temp_anatlabel_sourceelec){i_targetelec,5});
                    anatlabel_targetelec{i_sub}{i_targetelec, i_sourceelec} = ...
                        strrep(...
                        anatlabel_targetelec{i_sub}{i_targetelec, i_sourceelec}, ...
                        ', ', '_');
                end
            end
            %In cases where source and target elec are identical (e.g.,
            %for Pred_Pred/PE_PE effects) add an empty entry in place of
            %the identical elec along the diagonal, to keep the matrix 
            %similar to the GC output matrix
            if strcmp(param.GC.ElecPairEffect{i_effect},  'Pred_Pred') || ...
                    strcmp(param.GC.ElecPairEffect{i_effect},  'PE_PE')

                temp_fulllabel_targetelec = ...
                    cell(length(GCdata{i_sub}.ind_pairedelecs_to{i_effect}), ...
                    length(GCdata{i_sub}.ind_pairedelecs_from{i_effect}));
                temp_anatlabel_targetelec = ...
                    cell(length(GCdata{i_sub}.ind_pairedelecs_to{i_effect}), ...
                    length(GCdata{i_sub}.ind_pairedelecs_from{i_effect}));

                for i_sourceelec = 1:size(temp_fulllabel_targetelec,2)
                    outputposchange_targetelec = 0;
                    for i_targetelec = 1:size(fulllabel_targetelec{i_sub},1)
                        if i_sourceelec == i_targetelec
                            outputposchange_targetelec = outputposchange_targetelec + 1;
                        end
                        temp_fulllabel_targetelec{i_targetelec + outputposchange_targetelec, i_sourceelec} = ...
                            fulllabel_targetelec{i_sub}{i_targetelec, i_sourceelec};
                        temp_anatlabel_targetelec{i_targetelec + outputposchange_targetelec, i_sourceelec} = ...
                            anatlabel_targetelec{i_sub}{i_targetelec, i_sourceelec};                        
                    end
                end
                fulllabel_targetelec{i_sub} = temp_fulllabel_targetelec;
                anatlabel_targetelec{i_sub} = temp_anatlabel_targetelec;
            end  
            %Note: For Pred_PE and PE_Pred, in cases where there is the
            %same electrode as source and target, also a NaN value is
            %introduced. However, here the NaN value is placed at the end 
            %of the (source-elec-specific) column and not at the respective
            %position of the target elec. Since later read-out is
            %individual across target elecs for each source elec, this does
            %not interfere withthe selection of anatomical regions
        end
    end

    %2.2 Read out all present anatomical locations across all subjects
    anatlabel_sourceelec_allsub = [];
    anatlabel_targetelec_allsub = [];
    anatreg{i_effect}.label_unique = ...
        {'AntPFC'; 'AntPFC_IFG'; 'AntPFC_PrecentralG'; ...
        'VentralT'; 'VentralT_STG'; 'VentralT_MTG'; ...
        'SupParLob'; 'SupParLob_PostcentralG'; 'SupParLob_SupramarginalG'; ...
        'OccipitalL'};

    %Create filter limiting to-be-analyzed anatomical regions to the ones selected
    for i_catlabels = 1:length(anatreg{i_effect}.label_unique)
        anatreg{i_effect}.filter(i_catlabels,1) = ...
            any(strcmp(AnatReg_CatLabels, anatreg{i_effect}.label_unique{i_catlabels}));
    end
    anatreg{i_effect}.filter_label = [num2str(sum(anatreg{i_effect}.filter)) 'AnatReg'];

    %Flexibly determine if and how to suballocate anat regions
    %For each unique label, create index and index of all sub-labels
    anatreg{i_effect}.labelindex_aggregated = [];
    for i_catlabels = 1:length(anatreg{i_effect}.label_unique)
        temp_index = ...
            find(~cellfun('isempty', strfind(...
            anatreg{i_effect}.label_unique, anatreg{i_effect}.label_unique{i_catlabels})));
        if ~isempty(temp_index)
            anatreg{i_effect}.labelindex_aggregated(i_catlabels,:) = temp_index;
        else
            anatreg{i_effect}.labelindex_aggregated(i_catlabels,:) = NaN;
        end
    end

    %Adjust anat cat labels to more fitting anatomical labels
    anatreg{i_effect}.label_aggregated = [];
    for i_catlabels = 1:length(anatreg{i_effect}.label_unique)
        if strcmp(anatreg{i_effect}.label_unique{i_catlabels}, 'AntPFC')
            anatreg{i_effect}.label_aggregated{i_catlabels,1} = 'Frontal';
        elseif strcmp(anatreg{i_effect}.label_unique{i_catlabels}, 'AntPFC_IFG')
            anatreg{i_effect}.label_aggregated{i_catlabels,1} = 'IFG';
        elseif strcmp(anatreg{i_effect}.label_unique{i_catlabels}, 'AntPFC_PrecentralG')
            anatreg{i_effect}.label_aggregated{i_catlabels,1} = 'PrecG';
        elseif strcmp(anatreg{i_effect}.label_unique{i_catlabels}, 'VentralT')
            anatreg{i_effect}.label_aggregated{i_catlabels,1} = 'Temporal';
        elseif strcmp(anatreg{i_effect}.label_unique{i_catlabels}, 'VentralT_STG')
            anatreg{i_effect}.label_aggregated{i_catlabels,1} = 'STG';
        elseif strcmp(anatreg{i_effect}.label_unique{i_catlabels}, 'VentralT_MTG')
            anatreg{i_effect}.label_aggregated{i_catlabels,1} = 'MTG';
        elseif strcmp(anatreg{i_effect}.label_unique{i_catlabels}, 'SupParLob')
            anatreg{i_effect}.label_aggregated{i_catlabels,1} = 'Parietal';
        elseif strcmp(anatreg{i_effect}.label_unique{i_catlabels}, 'SupParLob_PostcentralG')
            anatreg{i_effect}.label_aggregated{i_catlabels,1} = 'PostcG';
        elseif strcmp(anatreg{i_effect}.label_unique{i_catlabels}, 'SupParLob_SupramarginalG')
            anatreg{i_effect}.label_aggregated{i_catlabels,1} = 'SupmarG';
        elseif strcmp(anatreg{i_effect}.label_unique{i_catlabels}, 'OccipitalL')
            anatreg{i_effect}.label_aggregated{i_catlabels,1} = 'OcciL';
        end
    end


    %2.3 Assign anatreg index and label for each source and target elec per subject
    for i_sub = 1:length(subs)
        if ~isempty(GCdata{i_sub}.temporalGC{i_effect})
            for i_sourceelec = 1:length(GCdata{i_sub}.ind_pairedelecs_from{i_effect})
                for i_anatlabels = 1:length(anatreg{i_effect}.label_aggregated)
                    if contains(anatlabel_sourceelec{i_sub}{i_sourceelec}, ...
                            anatreg{i_effect}.label_unique{i_anatlabels})
                        anatreg{i_effect}.sourceelec_index{i_sub}(i_sourceelec,1) = ...
                            i_anatlabels;
                        anatreg{i_effect}.sourceelec_label{i_sub}{i_sourceelec,1} = ...
                            anatreg{i_effect}.label_aggregated{i_anatlabels};
                    end
                end
                %Labels of target elecs vary depending on their source elecs
                for i_targetelec = 1:size(anatlabel_targetelec{i_sub},1)
                    for i_anatlabels = 1:length(anatreg{i_effect}.label_aggregated)
                        if ~isempty(anatlabel_targetelec{i_sub}{i_targetelec,i_sourceelec})
                            if contains(anatlabel_targetelec{i_sub}{i_targetelec,i_sourceelec}, ...
                                    anatreg{i_effect}.label_unique{i_anatlabels})
                                anatreg{i_effect}.targetelec_index{i_sub}(i_targetelec,i_sourceelec) = ...
                                    i_anatlabels;
                                anatreg{i_effect}.targetelec_label{i_sub}{i_targetelec,i_sourceelec} = ...
                                    anatreg{i_effect}.label_aggregated{i_anatlabels};
                            end
                        else
                            anatreg{i_effect}.targetelec_index{i_sub}(i_targetelec,i_sourceelec) = ...
                                nan;
                            anatreg{i_effect}.targetelec_label{i_sub}{i_targetelec,i_sourceelec} = ...
                                nan;
                        end
                    end
                end
            end
        end
    end


    %2.4 Pool elec pair-wise GC estimates across subjects per selected anatomical region
    %Create empty structs
    for i_anatlabels_source = 1:length(anatreg{i_effect}.label_aggregated)
        for i_anatlabels_target = 1:length(anatreg{i_effect}.label_aggregated)
            for i_TD = 1:length(ToneDur_text)                
                Gavg_tempGC_peranatreg.source2target{i_effect}{...
                    i_anatlabels_source,i_anatlabels_target,i_TD} = [];
                Gavg_tempGC_peranatreg.target2source{i_effect}{...
                    i_anatlabels_source,i_anatlabels_target,i_TD} = [];

                Gavg_spectGC_peranatreg.source2target{i_effect}{...
                    i_anatlabels_source,i_anatlabels_target,i_TD} = [];
                Gavg_spectGC_peranatreg.target2source{i_effect}{...
                    i_anatlabels_source,i_anatlabels_target,i_TD} = [];

                Gavg_tempGC_peranatreg.subject_index{i_effect}{...
                    i_anatlabels_source,i_anatlabels_target,i_TD} = [];
            end
        end
    end
    NumGC_peranatreg{i_effect} = ...
        zeros(length(anatreg{i_effect}.label_aggregated), ...
        length(anatreg{i_effect}.label_aggregated), ...
        length(ToneDur_text)); %Dim: AnatLabelSource, AnatLabelTarget, TD

    %Pool GC estimates for similar anat reg across electrode pairs    
    for i_TD = 1:length(ToneDur_text)        
        for i_anatlabels_source = 1:length(anatreg{i_effect}.label_aggregated)
            for i_anatlabels_target = 1:length(anatreg{i_effect}.label_aggregated)

                %Empty structs generated for each source/target anatomical region, 
                %since we pool over this variable
                temp_tempGC.source2target   = []; 
                temp_tempGC.target2source   = [];
                temp_spectGC.target2source  = [];
                temp_spectGC.source2target  = [];
                temp_tempGC.subject_index   = [];
                temp_spectGC_repcount       = 0;

                for i_sub = 1:length(subs)
                    if ~isempty(GCdata{i_sub}.temporalGC{i_effect}) && ...
                            isreal(GCdata{i_sub}.temporalGC{i_effect}{i_TD}.source2target)                        
                        for i_sourceelec = 1:length(GCdata{i_sub}.ind_pairedelecs_from{i_effect})
                            if any(anatreg{i_effect}.sourceelec_index{i_sub}(i_sourceelec) == ... 
                                    anatreg{i_effect}.labelindex_aggregated(i_anatlabels_source,:))
                                    %Filter: Select GC value only if corresponing source elec corresponds to current anat reg

                                for i_targetelec = 1:size(anatreg{i_effect}.targetelec_index{i_sub},1) 
                                    if any(anatreg{i_effect}.targetelec_index{i_sub}(i_targetelec,i_sourceelec) == ...
                                            anatreg{i_effect}.labelindex_aggregated(i_anatlabels_target,:))
                                            %Filter: Select GC value only if corresponing target elec corresponds to current anat reg

                                        label_targetelec = strtok(fulllabel_targetelec{i_sub}{i_targetelec,i_sourceelec},'_');
                                        disp(['Found GC value for source loc [' anatreg{i_effect}.label_aggregated{i_anatlabels_source} ...
                                            '(' GCdata{i_sub}.label_allelecs{GCdata{i_sub}.ind_pairedelecs_from{i_effect}(i_sourceelec)} ')] ' ...
                                            'to target loc [' anatreg{i_effect}.label_aggregated{i_anatlabels_target}  '('...
                                            label_targetelec ')]' ...
                                            ' for sub: ' subs{i_sub} ...
                                            ', TD: ' num2str(i_TD)]);

                                        %Pool GC across electrode pairs from same anatomical region
                                        temp_tempGC.source2target = ...
                                            [temp_tempGC.source2target; ...
                                            squeeze(GCdata{i_sub}.temporalGC{i_effect}{i_TD...
                                            }.source2target(i_sourceelec, i_targetelec, :))'];
                                        temp_tempGC.target2source = ...
                                            [temp_tempGC.target2source; ...
                                            squeeze(GCdata{i_sub}.temporalGC{i_effect}{i_TD...
                                            }.target2source(i_sourceelec, i_targetelec, :))'];

                                        temp_spectGC_repcount = temp_spectGC_repcount + 1;
                                        temp_spectGC.source2target(:,:,temp_spectGC_repcount) = ...
                                            squeeze(GCdata{i_sub}.spectralGC{i_effect}{i_TD...
                                            }.source2target(i_sourceelec, i_targetelec, :, :));
                                        temp_spectGC.target2source(:,:,temp_spectGC_repcount) = ...
                                            squeeze(GCdata{i_sub}.spectralGC{i_effect}{i_TD...
                                            }.target2source(i_sourceelec, i_targetelec, :, :));

                                        temp_tempGC.subject_index = ...
                                            [temp_tempGC.subject_index; i_sub];
                                    end
                                end
                            end
                        end
                    end
                end
                %Count total number of valid electrode pairs for current anat region
                if ~isempty(temp_tempGC.source2target)
                    NumGC_peranatreg{i_effect}(i_anatlabels_source, i_anatlabels_target, i_TD) = ...
                        sum(~isnan(temp_tempGC.source2target(:,1)));
                else
                    NumGC_peranatreg{i_effect}(i_anatlabels_source, i_anatlabels_target, i_TD) = ...
                        0;
                end

                %Store pooled data in common struct
                Gavg_tempGC_peranatreg.source2target{i_effect}{...
                    i_anatlabels_source,i_anatlabels_target,i_TD} = ...
                    temp_tempGC.source2target;
                Gavg_tempGC_peranatreg.target2source{i_effect}{...
                    i_anatlabels_source,i_anatlabels_target,i_TD} = ...
                    temp_tempGC.target2source;
                Gavg_spectGC_peranatreg.source2target{i_effect}{...
                    i_anatlabels_source,i_anatlabels_target,i_TD} = ...
                     temp_spectGC.source2target;
                Gavg_spectGC_peranatreg.target2source{i_effect}{...
                    i_anatlabels_source,i_anatlabels_target,i_TD} = ...
                    temp_spectGC.target2source;

                Gavg_tempGC_peranatreg.subject_index{i_effect}{...
                    i_anatlabels_source,i_anatlabels_target,i_TD} = ...
                    temp_tempGC.subject_index;
            end
        end
    end

    %Optional: Perform 1 sample KS test to test GC values across elec pairs (per time point and TD)
    %for normality (0 = null hypothesis of normal dist, 1 = rejection of null)    
%         for i_TD = 1:length(ToneDur_text)
%             KS_mat{i_TD} = ...
%                 nan(length(anatreg{i_effect}.label_aggregated), ...
%                 length(anatreg{i_effect}.label_aggregated), ...
%                 size(GCdata{i_sub}.temporalGC{i_effect}{1}.source2target,3));
%             for i_sourceelec = 1:length(anatreg{i_effect}.label_aggregated)
%                 for i_targetelec = 1:length(anatreg{i_effect}.label_aggregated)
%                     for i_toneepoch = 1:length(param.GC.tone_aggregation)
%                         if any((squeeze(Gavg_tempGC_peranatreg.source2target{i_effect}...
%                                 {i_sourceelec, i_targetelec, i_TD})))
%                             KS_mat{i_TD}(i_sourceelec, i_targetelec, i_toneepoch) = ...
%                                 kstest(...
%                                 Gavg_tempGC_peranatreg.source2target{i_effect}{...
%                                 i_sourceelec, i_targetelec, i_TD}(:,i_toneepoch),...
%                                 'Alpha', 0.05);
%                         end
%                     end
%                 end
%             end
%             KStestdata.sum_nonan{i_effect}{i_TD}  = length(find(~isnan(KS_mat{i_TD}(:,:,:))));
%             KStestdata.sum_sign{i_effect}{i_TD}   = length(find(KS_mat{i_TD}(:,:,:) == 1));
%             disp([num2str(KStestdata.sum_sign{i_effect}{i_TD}) ' / ' ...
%                 num2str(KStestdata.sum_nonan{i_effect}{i_TD}) ...
%                 ' entries not normally distributed for ' title_label1 ' - ' title_label2])
%         end
end

% In order to plot across both tone sequence conditions and across tone
% durations:
Gavg_tempGC_peranatreg_avg = Gavg_tempGC_peranatreg;  % Copy the structure to preserve its shape

% To plot for each tone sequence separately
ToneSeq_GC = Gavg_tempGC_peranatreg;

% Loop through the cells of source2target in both Gavg_tempGC_peranatreg1 and Gavg_tempGC_peranatreg2
for i = 1:4
    % Initialize a 10x10 cell array for the averaged data
    averaged_cell = cell(10, 10);
    avg_cell = cell(10,10); %to plot GC across tone sequence 1
    
    for j = 1:10
        for k = 1:10
            % Extract the relevant data from both Gavg_tempGC_peranatreg1 and Gavg_tempGC_peranatreg2
            data1_1 = Gavg_tempGC_peranatreg.source2target{i}{j, k, 1}; %data for tone duration 1
            data1_2 = Gavg_tempGC_peranatreg.source2target{i}{j, k, 2}; %data for tone duration 2
            
            % Average the data across both structures and both time points (1 and 2)
            avg_data = (data1_1 + data1_2) / 2;
            
            % Store the result in the averaged_cell array
            averaged_cell{j, k} = avg_data;
            avg_cell{j,k} = avg_data;
        end
    end

    % Store the averaged data back into the output structure
    Gavg_tempGC_peranatreg_avg.source2target{i} = averaged_cell;
    ToneSeq_GC.source2target{i} = avg_cell;
end
Gavg_tempGC_peranatreg_avg = rmfield(Gavg_tempGC_peranatreg_avg, {'target2source', 'subject_index'});
ToneSeq_GC = rmfield(ToneSeq_GC, {'target2source', 'subject_index'});

%% 3. Plot GC estimate for time and frequency domain within (ie., Pred_pred/PE_PE) effects
%Plot median and IQR of GC estimates for each anatomical region as a
%shaded line plot. Different anatomical pairs are plotted in subplots.
%Plotting is separetly done for within (ie., Pred_pred/PE_PE) and between
%(Pred_PE, PE_Pred) effects.

% The anatomical pairs we want to save plots for
index_pairs = { ...
    [1, 2], ...
    [2, 1], ...
    [1, 3], ...
    [3, 1], ...
    [2, 3], ...
    [3, 2] ...
};
color_list = {'r', 'g', 'b', 'm', 'c', 'y'};  % Red, Green, Blue, Magenta, Cyan, Yellow
figures_dir = '/isilon/LFMI/VMdrive/Lua/NASTD/Figures/';

% Loop over each pair of indices
for i_effect = 1:2
    temp_currfilt = find(anatreg{i_effect}.filter == 1);
    for pair_idx = 1:(length(index_pairs)-1)
        % Extract the current pair of indices
        temp_currfilt_pair = index_pairs{pair_idx};

        % Create a new figure for each plot
        f = figure('visible', set_visibility);
        set(f, 'units', 'normalized', 'outerposition', [0 0 1 1]) % full screen

        % Create an axes object for the plot
        f_sub1 = axes(f);  % Use axes within the figure

        % Get the data for the current pair of indices
        data_to_plot = nanmean(Gavg_tempGC_peranatreg_avg.source2target{i_effect}{temp_currfilt(temp_currfilt_pair(1)), temp_currfilt(temp_currfilt_pair(2))});
        plot_color = color_list{pair_idx};

        % Create the plot (Line plot with no error bars)
        hold(f_sub1, 'on');
        plot(f_sub1, (1:length(param.GC.tone_aggregation)) - 0.05, ...
            data_to_plot, 'Color', plot_color, 'Marker', 'd', 'MarkerFaceColor', plot_color, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerSize', 10);
        hold(f_sub1, 'off');

        % Set axis properties
        Ylim_temp = f_sub1.YLim;
        f_sub1.XLim = [1 - 0.3, length(param.GC.tone_aggregation) + 0.3];

        if istrue(Yaxis_scaling)
            f_sub1.YLim = [Ylim_temp(1), Ylim_temp(2)];
            f_sub1.YLim = [5e-4, 15e-3];
        end

        f_sub1.XTick = 1:length(param.GC.tone_aggregation);
        f_sub1.XTickLabel = param.GC.tone_aggregation;
        a = get(gca, 'XTickLabel');
        set(gca, 'XTickLabel', a, 'fontsize', 8);

        % Set labels
        f_sub1.XLabel.String = {'Tones'};
        f_sub1.XLabel.FontSize = 12;
        f_sub1.YLabel.String = {'GC estimate'};
        f_sub1.YLabel.FontSize = 12;
        
        for i_elec = 1:size(Gavg_tempGC_peranatreg_avg.source2target{i_effect}{temp_currfilt(temp_currfilt_pair(1)), temp_currfilt(temp_currfilt_pair(2))},1)
            coeff_linfit(i_elec,:) = ...
                polyfit...
                (1:length(param.GC.tone_aggregation), ...
                (Gavg_tempGC_peranatreg_avg.source2target{i_effect}{temp_currfilt(temp_currfilt_pair(1)), temp_currfilt(temp_currfilt_pair(2))}(i_elec,:)), 1);
        end
        
        % Compute predicted values for each electrode
        x_vals = (1:length(param.GC.tone_aggregation)) - 0.05;
        predicted_values = zeros(size(coeff_linfit, 1), length(x_vals));
        for i_elec = 1:size(coeff_linfit, 1)
            predicted_values(i_elec, :) = coeff_linfit(i_elec, 1) * x_vals + coeff_linfit(i_elec, 2);
        end

        % Compute the mean trendline by averaging predicted values across electrodes
        mean_trendline = mean(predicted_values, 1);

        % Add the mean trendline to the plot
        hold(f_sub1, 'on');
        plot(f_sub1, x_vals, mean_trendline - 0.001, '--', 'Color', plot_color, 'LineWidth', 2);  % Black dashed line for the trendline


        % Add title and legend
        legend(f_sub1, [anatreg{i_effect}.label_aggregated{temp_currfilt(temp_currfilt_pair(1))} '->' ...
                anatreg{i_effect}.label_aggregated{temp_currfilt(temp_currfilt_pair(2))} ...
                ' (' num2str(NumGC_peranatreg{i_effect}(temp_currfilt(temp_currfilt_pair(1)), temp_currfilt(temp_currfilt_pair(2)))) ' elec pairs)'], ...
                'location', 'northeast');
        legend(f_sub1, 'FontSize', 14);  

        set(f_sub1, 'TickLabelInterpreter', 'none');
        title(f_sub1, ['Pred -> Pred | ' ...
            anatreg{i_effect}.label_aggregated{temp_currfilt(temp_currfilt_pair(1))} ' - ' ...
            anatreg{i_effect}.label_aggregated{temp_currfilt(temp_currfilt_pair(2))}], ...
            'fontsize', 12);

        % Save the figure to a file with transparent background
        savePath = [figures_dir, 'GC_1to31ALL_', param.GC.ElecPairEffect{i_effect}, '_', anatreg{i_effect}.label_aggregated{temp_currfilt(temp_currfilt_pair(1))}, anatreg{i_effect}.label_aggregated{temp_currfilt(temp_currfilt_pair(2))}, '_withTrendline.eps'];
        set(f_sub1, 'Color', 'none');  % Set figure background to transparent

        % Save the figure as EPS with transparent background
        print(f, savePath, '-depsc', '-r300', '-painters'); 
        
        % Close the figure after saving
        close(f);
    end
end

% Figures for both tone sequences for parietal to frontal
for i_effect = 1:2
    temp_currfilt = find(anatreg{i_effect}.filter == 1);
    pair_idx = 4
        
    % Extract the current pair of indices
    temp_currfilt_pair = index_pairs{pair_idx};

    % Create a new figure for each plot
    f = figure('visible', set_visibility);
    set(f, 'units', 'normalized', 'outerposition', [0 0 1 1]) % full screen

    % Create an axes object for the plot
    f_sub1 = axes(f);  % Use axes within the figure

    % Get the data for the current pair of indices
    data_to_plot1 = nanmean(ToneSeq1_GC.source2target{i_effect}{temp_currfilt(temp_currfilt_pair(1)), temp_currfilt(temp_currfilt_pair(2))});
    data_to_plot2 = nanmean(ToneSeq2_GC.source2target{i_effect}{temp_currfilt(temp_currfilt_pair(1)), temp_currfilt(temp_currfilt_pair(2))});
    plot_color1 = color_list{pair_idx};
    plot_color2 = color_list{pair_idx + 1};

    % Create the plot (Line plot with no error bars)
    hold(f_sub1, 'on');
    plot(f_sub1, (1:length(param.GC.tone_aggregation)) - 0.05, ...
        data_to_plot1, 'Color', plot_color1, 'Marker', 'd', 'MarkerFaceColor', plot_color1, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerSize', 10);
    plot(f_sub1, (1:length(param.GC.tone_aggregation)) + 0.05, ...
        data_to_plot2, 'Color', plot_color2, 'Marker', 'd', 'MarkerFaceColor', plot_color2, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerSize', 10);
    hold(f_sub1, 'off');

    % Set axis properties
    Ylim_temp = f_sub1.YLim;
    f_sub1.XLim = [1 - 0.3, length(param.GC.tone_aggregation) + 0.3];

    if istrue(Yaxis_scaling)
        f_sub1.YLim = [Ylim_temp(1) - 0.1, Ylim_temp(2)];
        f_sub1.YLim = [5e-4, 15e-3];
    end

    f_sub1.XTick = 1:length(param.GC.tone_aggregation1);
    f_sub1.XTickLabel = param.GC.tone_aggregation1;
    a = get(gca, 'XTickLabel');
    set(gca, 'XTickLabel', a, 'fontsize', 8);

    % Set labels
    f_sub1.XLabel.String = {'Tones'};
    f_sub1.XLabel.FontSize = 12;
    f_sub1.YLabel.String = {'GC estimate'};
    f_sub1.YLabel.FontSize = 12;

    for i_elec = 1:size(ToneSeq1_GC.source2target{i_effect}{temp_currfilt(temp_currfilt_pair(1)), temp_currfilt(temp_currfilt_pair(2))},1)
        coeff_linfit1(i_elec,:) = ...
            polyfit...
            (1:length(param.GC.tone_aggregation1), ...
            (ToneSeq1_GC.source2target{i_effect}{temp_currfilt(temp_currfilt_pair(1)), temp_currfilt(temp_currfilt_pair(2))}(i_elec,:)), 1);
    end
    
    for i_elec = 1:size(ToneSeq2_GC.source2target{i_effect}{temp_currfilt(temp_currfilt_pair(1)), temp_currfilt(temp_currfilt_pair(2))},1)
        coeff_linfit2(i_elec,:) = ...
            polyfit...
            (1:length(param.GC.tone_aggregation1), ...
            (ToneSeq2_GC.source2target{i_effect}{temp_currfilt(temp_currfilt_pair(1)), temp_currfilt(temp_currfilt_pair(2))}(i_elec,:)), 1);
    end

    % Compute predicted values for each electrode
    x_vals = (1:length(param.GC.tone_aggregation1)) - 0.05;
    predicted_values1 = zeros(size(coeff_linfit1, 1), length(x_vals));
    predicted_values2 = zeros(size(coeff_linfit2, 1), length(x_vals));

    for i_elec = 1:size(coeff_linfit1, 1)
        predicted_values1(i_elec, :) = coeff_linfit1(i_elec, 1) * x_vals + coeff_linfit1(i_elec, 2);
    end
    
    for i_elec = 1:size(coeff_linfit2, 1)
         predicted_values2(i_elec, :) = coeff_linfit2(i_elec, 1) * x_vals + coeff_linfit2(i_elec, 2);
    end

    % Compute the mean trendline by averaging predicted values across electrodes
    mean_trendline1 = mean(predicted_values1, 1);
    mean_trendline2 = mean(predicted_values2, 1);

    % Add the mean trendline to the plot
    hold(f_sub1, 'on');
    plot(f_sub1, x_vals, mean_trendline1 - 0.001, '--', 'Color', plot_color1, 'LineWidth', 2);  % Black dashed line for the trendline
    plot(f_sub1, x_vals + 0.1, mean_trendline2 - 0.001, '--', 'Color', plot_color2, 'LineWidth', 2);  % Black dashed line for the trendline
    
    all_data = [data_to_plot1, data_to_plot2, mean_trendline1, mean_trendline2];
    new_ylim = [min(all_data) - 0.005, max(all_data) + 0.005];  % Add a small buffer for better visualization
    f_sub1.YLim = new_ylim;

    % Add title and legend
    legend(f_sub1, ...
    {[anatreg{i_effect}.label_aggregated{temp_currfilt(temp_currfilt_pair(1))} '->' ...
      anatreg{i_effect}.label_aggregated{temp_currfilt(temp_currfilt_pair(2))} ...
      ' (Tone Seq. 1)'], ...
     [anatreg{i_effect}.label_aggregated{temp_currfilt(temp_currfilt_pair(1))} '->' ...
      anatreg{i_effect}.label_aggregated{temp_currfilt(temp_currfilt_pair(2))} ...
      ' (Tone Seq. 2)']}, ...
    'location', 'northwest', 'FontSize', 14);
    
    set(f_sub1, 'TickLabelInterpreter', 'none');
    title(f_sub1, ['Pred -> Pred | ' ...
        anatreg{i_effect}.label_aggregated{temp_currfilt(temp_currfilt_pair(1))} ' - ' ...
        anatreg{i_effect}.label_aggregated{temp_currfilt(temp_currfilt_pair(2))}], ...
        'fontsize', 12);

    % Save the figure to a file with transparent background
    savePath = [figures_dir, 'GC_1to31_', param.GC.ElecPairEffect{i_effect}, '_', anatreg{i_effect}.label_aggregated{temp_currfilt(temp_currfilt_pair(1))}, anatreg{i_effect}.label_aggregated{temp_currfilt(temp_currfilt_pair(2))}, '_withTrendline_bothToneSeqs.eps'];
    set(f_sub1, 'Color', 'none');  % Set figure background to transparent

    % Save the figure as EPS with transparent background
    print(f, savePath, '-depsc', '-r300', '-painters'); 

    % Close the figure after saving
    close(f);
end