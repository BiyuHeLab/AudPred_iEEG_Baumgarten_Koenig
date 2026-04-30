function SelElecs = NASTD_ECoG_Connectivity_ReadOutAllElecs_AllSubTD...
    (subs, ...
    plot_poststepFigs, ...
    paths_NASTD_ECoG)

%Aim: Read out all existing electrodes.
%These electrodes will form the basis for GC analysis.

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
%Add base dir and own script dir
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

path_fig = ([paths_NASTD_ECoG.ECoGdata_Connectivity ...
    'ElecSelect/Allsub_n' num2str(length(subs)) ...
    '/Figs/ElecPairs/AllElecs/']);
if (~exist(path_fig, 'dir')); mkdir(path_fig); end


%% 1) Load prediction effect data and aggregate relevant info across subjects
tic

clear AllElecs
allElecs_chanposIndex               = [];
AllElecs.all_subs.sub_index         = [];
AllElecs.all_subs.label_AnatCat     = [];
AllElecs.all_subs.index_AnatCat     = [];

AllElecs.all_subs.label_elec        = [];
AllElecs.all_subs.label_anat        = [];
AllElecs.all_subs.coords_elec       = [];

i_TD = 1; %TD doesn't matter since all elecs used
for i_sub = 1:length(subs)
    
    sub = subs{i_sub};
    disp(['-- Loading data for sub: ' sub ' --'])
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;
    
    
    %Load ECoG preproc data for channel labels and position
    loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];
    load(loadfile_ECoGpreprocdata);
    
    %Read out electrode number, label, coordinates, and anatomical labels and aggregate them across subjects
    SampleFreq                          = ...
        DataClean_AllTrials.fsample;
    nSensors_all                        = ...
        size(DataClean_AllTrials.cfg.info_elec.selected.Label,1);
    AllElecs.per_sub.elec_labels{i_sub} = ...
        DataClean_AllTrials.cfg.info_elec.selected.Label;
    
    for i_elec = 1:nSensors_all %Find index of preprocessed elecs in all elecs
        allElecs_chanposIndex(i_elec,1) = ...
            find(strcmp(DataClean_AllTrials.cfg.info_elec.selected.Label{i_elec}, ...
            DataClean_AllTrials.elec.label));
    end
    
    coords_sub{i_sub}   = ...
        DataClean_AllTrials.elec.chanpos(allElecs_chanposIndex,:);
    AllElecs.all_subs.coords_elec = ...
        [AllElecs.all_subs.coords_elec; coords_sub{i_sub}];
    AllElecs.all_subs.label_anat = ...
        [AllElecs.all_subs.label_anat; ...
        DataClean_AllTrials.elec.T1AnatLabel(allElecs_chanposIndex)];
    
    %Add additional elec + sub label
    for i_elec = 1:nSensors_all
        AllElecs.all_subs.label_elec{end+1,1} = ...OK
            [DataClean_AllTrials.cfg.info_elec.selected.Label{i_elec} ' ' sub];
    end
    
    %Categorize electrodes according to anatomical regions
    AnatReg_allSubs{i_sub} = ...
        NASTD_ECoG_AssignAnatRegions...
        (DataClean_AllTrials, DataClean_AllTrials.cfg.info_elec.selected.Label);
    %Aggregate category labels and indices across subjects
    AllElecs.all_subs.label_AnatCat  = ...
        [AllElecs.all_subs.label_AnatCat; ...
        AnatReg_allSubs{i_sub}.Info_perelec(:,3)];
    AllElecs.all_subs.index_AnatCat  = ...
        [AllElecs.all_subs.index_AnatCat; ...
        AnatReg_allSubs{i_sub}.CatIndex];
    
    %Data aggregation over subjects
    %Array to differentiate subject entries
    AllElecs.all_subs.sub_index = ...
        [AllElecs.all_subs.sub_index; ...
        ones(length(coords_sub{i_sub}),1)*i_sub];
    
    %Cleanup
    allElecs_chanposIndex = [];
    clear DataClean_AllTrials temp* coords_sub
    
    disp([' -- All electrodes processed for ' sub ' --']) 
end
disp(['-- All electrodes read out for all subs after ' num2str(round(toc/60),2) ' min --'])


%% 2) Summarize electrodes in table
SelElecs = struct;

index_elecall = [];
index_elecall = [1:length(AllElecs.all_subs.label_elec)]';
index_elecall = num2str(index_elecall);

index_elecpersub = [];
for i_sub = 1:length(subs)
    temp = [];
    temp = [1:length(AllElecs.per_sub.elec_labels{i_sub})]';
    index_elecpersub = [index_elecpersub; temp];
end
index_elecpersub = num2str(index_elecpersub);

varnames = ...
    {'Electrode Label', 'Subject Label', 'Sub index',...
    'AnatCat Label', 'Detailed Label', 'Elec Coords', ...
    'Elec Index (all elecs)', 'Elec Index (Ssub elecs)'};
SelElecs.AllElecs = table(...
    extractBefore(AllElecs.all_subs.label_elec,' '), ...
    extractAfter(AllElecs.all_subs.label_elec,' '), ...
    AllElecs.all_subs.sub_index, ...
    AllElecs.all_subs.label_AnatCat, ...
    AllElecs.all_subs.label_anat, ...
    AllElecs.all_subs.coords_elec, ...
    index_elecall, ...
    index_elecpersub, ...
    'VariableNames', varnames);
SelElecs.Properties.Description = ...
    ['Information for all preprocessed electrodes'];


%% 3) Determine possible pairings between sign. electrodes and plot single subject pairings
%Determine pairings between electrodes within the same subject 
%and specific spatial (e.g., frontal - temporal) combinations.

SelElecs.Pairs = struct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_poststepFigs == 1
    
    %% 3.1 All regions
    SelElecs.Pairs.All_All.AllRegions.Pairs_perelec = struct;
    SelElecs.Pairs.All_All.AllRegions.Ind_UniquePairings = [];

    for i_elec = 1:height(SelElecs.AllElecs)
        %Determine current elec information
        temp_label  = SelElecs.AllElecs{i_elec,1};
        temp_sub    = SelElecs.AllElecs{i_elec,2};
        temp_anatcat = SelElecs.AllElecs{i_elec,4};
        %Find other elecs from same subject
        temp_index_otherelecs = ...
            find(strcmp(temp_sub, SelElecs.AllElecs{:,2}) & ...
            ~strcmp(temp_label, SelElecs.AllElecs{:,1}));
        %Create label for subfield (elec label + sub + anatcat)
        if ~isempty(strfind(temp_anatcat{1}, ','))
            temp_fieldname = ...
                [temp_label{1} '_' temp_sub{1} '_' extractBefore(temp_anatcat{1}, ',')];
        else
            temp_fieldname = ...
                [temp_label{1} '_' temp_sub{1} '_' temp_anatcat{1}];
        end
        %Store paired elec information in new struct
        if ~isempty(temp_index_otherelecs)
            SelElecs.Pairs.All_All.AllRegions.Pairs_perelec.(temp_fieldname) = ...
                SelElecs.AllElecs(temp_index_otherelecs,:);
        else
            SelElecs.Pairs.All_All.AllRegions.Pairs_perelec.(temp_fieldname) = [];
        end
        clear temp*
    end
    
    %Determine number of possible pairings for each subject
    for i_sub = 1:length(subs)
        sub = subs(i_sub);
        index_sign_elecs_currsub = find(strcmp(sub, SelElecs.AllElecs{:,2}));
        
        %Determine number of possible pairings for each subject
        [p,q] = meshgrid(index_sign_elecs_currsub, index_sign_elecs_currsub);
        temp_pairs1 = [p(:) q(:)];
        if size(temp_pairs1,1) > 1
            %Keep only electrode pairings between non-identical electrodes
            temp_pairs2 = [];
            for i_pairs = 1:length(temp_pairs1)
                if temp_pairs1(i_pairs,1) ~= temp_pairs1(i_pairs,2)
                    temp_pairs2 = [temp_pairs2; temp_pairs1(i_pairs,:)];
                end
            end
            %delete double pairings (order irrelevant)
            i_doubles = false(length(temp_pairs2),1);
            i_compared = false(length(temp_pairs2),1);
            for i_pairs = 1:length(temp_pairs2)
                temp_doubles = find(all(temp_pairs2(i_pairs,:) == fliplr(temp_pairs2),2));
                i_compared(i_pairs) = true;
                if i_compared(temp_doubles) == false
                    i_doubles(temp_doubles) = true;
                end
            end
            %Read out number and indices of pairings
            SelElecs.Pairs.All_All.AllRegions.Num_UniquePairings_persub(1,i_sub) = ...
                size(temp_pairs2(~i_doubles,:),1);
            SelElecs.Pairs.All_All.AllRegions.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs2(~i_doubles,:);
        elseif size(temp_pairs1,1) == 1 && temp_pairs1(1) ~= temp_pairs1(2)
            SelElecs.Pairs.All_All.AllRegions.Num_UniquePairings_persub(1,i_sub) = ...
                1;
            SelElecs.Pairs.All_All.AllRegions.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs1;
        else
            SelElecs.Pairs.All_All.AllRegions.Num_UniquePairings_persub(1,i_sub) = ...
                0;
            SelElecs.Pairs.All_All.AllRegions.Ind_UniquePairings_persub{i_sub} = ...
                nan;
        end
        
        %aggregate pairings across subjects
        SelElecs.Pairs.All_All.AllRegions.Ind_UniquePairings = ...
            [SelElecs.Pairs.All_All.AllRegions.Ind_UniquePairings; ...
            SelElecs.Pairs.All_All.AllRegions.Ind_UniquePairings_persub{i_sub}];        
    end
    clear temp*
    
    %Plot figure showing Number and location of each pair for each subject
    figure
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    DimSubplot = [3 3];
    sgtitle(['Possible Electrode Pairs - All to All - All regions'])
    
    for i_sub = 1:length(subs)        
        if any(strcmp(subs(i_sub),SelElecs.AllElecs{:,2}))
            if ~isnan(SelElecs.Pairs.All_All.AllRegions.Ind_UniquePairings_persub{i_sub})
                curr_elec_ind = unique(SelElecs.Pairs.All_All.AllRegions.Ind_UniquePairings_persub{i_sub});
            else
                curr_elec_ind = find(strcmp(subs(i_sub),SelElecs.AllElecs{:,2}));
            end
            elec_coords = SelElecs.AllElecs{curr_elec_ind,6};
            elec_coords(:,1) = abs(elec_coords(:,1)) * -1;%Project all electrodes on one hemisphere
            elec_labels = {};
            for i_elec = 1:size(elec_coords,1) %No electrode labels
                elec_labels{i_elec} = cell2mat(SelElecs.AllElecs{curr_elec_ind(i_elec),1});
            end
            elec_data = ones(size(elec_coords,1),1)*i_sub;
            sub_index = ones(size(elec_coords,1),1)*i_sub;
            elec_signindex = true(size(elec_coords,1),1);
            elec_typeindex = [ones(length(curr_elec_ind),1) * 3]; %1 = source, 2 = target, 3 = both
            elec_pairings = SelElecs.Pairs.All_All.AllRegions.Ind_UniquePairings_persub{i_sub};
            effect_comp = 'All_All';
            chanSize = ones(1,length(elec_data))*3; %electrode size (arbitrary)
            clims = [1 9];
            cmap  = distinguishable_colors(9); %different colors per sub
            cmap  = ones(length(subs), 3).*[0.6207    0.3103    0.2759]; %all dark red
            
            textcolor_rgb   = [0 0 0];
            sp_title = {[subs{i_sub} ' - ' num2str(size(elec_coords,1)) ' All elecs - ' ...
                num2str(SelElecs.Pairs.All_All.AllRegions.Num_UniquePairings_persub(i_sub)) ' elec pairs']};
            
            NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH_Pairs...
                (elec_coords, elec_labels, elec_data, elec_signindex, elec_typeindex, ...
                sub_index, elec_pairings, effect_comp, ...
                chanSize, clims, cmap, textcolor_rgb, ...
                DimSubplot, i_sub, [0 -0.02 0 0], [], sp_title);
        end
    end
    filename     = ['Surf1H_' ...
        'Allsubn' num2str(length(subs)) '_' ...
        'PairsAllElec_AllReg.png'];
    figfile      = [path_fig filename];
    saveas(gcf, figfile, 'png'); %save png version
    close;
    clear elec*
    
    %% 3.2 Frontal - temporal
    SelElecs.Pairs.All_All.Frontal_Temporal.Pairs_perelec = struct;
    SelElecs.Pairs.All_All.Frontal_Temporal.Ind_UniquePairings = [];

    for i_elec = 1:height(SelElecs.AllElecs)
        %Determine current elec information
        temp_label  = SelElecs.AllElecs{i_elec,1};
        temp_sub    = SelElecs.AllElecs{i_elec,2};
        temp_anatcat = SelElecs.AllElecs{i_elec,4};
        if strfind(temp_anatcat{1}, 'AntPFC') %If current electrode is frontal
            
            temp_index_otherelecs = ...%Find other sign. elecs from same subject in temporal areas
                find(strcmp(temp_sub, SelElecs.AllElecs{:,2}) & ...
                ~strcmp(temp_label, SelElecs.AllElecs{:,1}) & ...
                ~cellfun('isempty',strfind(SelElecs.AllElecs{:,4}, 'VentralT')));
            if ~isempty(strfind(temp_anatcat{1}, ',')) %Create label for subfield
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' extractBefore(temp_anatcat{1}, ',')];
            else
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' temp_anatcat{1}];
            end
            
            if ~isempty(temp_index_otherelecs) %Store pair information in new struct
                SelElecs.Pairs.All_All.Frontal_Temporal.Pairs_perelec.(temp_fieldname) = ...
                    SelElecs.AllElecs(temp_index_otherelecs,:);
            else
                SelElecs.Pairs.All_All.Frontal_Temporal.Pairs_perelec.(temp_fieldname) = [];
            end
        end
    end
    %Determine number of possible pairings for each subject
    for i_sub = 1:length(subs)
        sub = subs(i_sub);
        sign_elecs_currsub1 = ...
            find(strcmp(sub, SelElecs.AllElecs{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.AllElecs{:,4}, 'AntPFC')));
        sign_elecs_currsub2 = ...
            find(strcmp(sub, SelElecs.AllElecs{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.AllElecs{:,4}, 'VentralT')));
        %Determine number of possible pairings for each subject
        [p,q] = meshgrid(sign_elecs_currsub1, sign_elecs_currsub2);
        temp_pairs1 = [p(:) q(:)];
        if size(temp_pairs1,1) > 1
            %delete same-same pairings
            temp_pairs2 = [];
            for i_pairs = 1:length(temp_pairs1)
                if temp_pairs1(i_pairs,1) ~= temp_pairs1(i_pairs,2)
                    temp_pairs2 = [temp_pairs2; temp_pairs1(i_pairs,:)];
                end
            end
            %delete doubles
            i_doubles = false(length(temp_pairs2),1);
            i_compared = false(length(temp_pairs2),1);
            for i_pairs = 1:length(temp_pairs2)
                temp_doubles = find(all(temp_pairs2(i_pairs,:) == fliplr(temp_pairs2),2));
                i_compared(i_pairs) = true;
                
                if i_compared(temp_doubles) == false
                    i_doubles(temp_doubles) = true;
                end
            end
            SelElecs.Pairs.All_All.Frontal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                size(temp_pairs2(~i_doubles,:),1);
            SelElecs.Pairs.All_All.Frontal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs2(~i_doubles,:);
        elseif size(temp_pairs1,1) == 1 && temp_pairs1(1) ~= temp_pairs1(2)
            SelElecs.Pairs.All_All.Frontal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                1;
            SelElecs.Pairs.All_All.Frontal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs1;
        else
            SelElecs.Pairs.All_All.Frontal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                0;
            SelElecs.Pairs.All_All.Frontal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                nan;
        end
        %aggregate pairings across subjects
        SelElecs.Pairs.All_All.Frontal_Temporal.Ind_UniquePairings = ...
            [SelElecs.Pairs.All_All.Frontal_Temporal.Ind_UniquePairings; ...
            SelElecs.Pairs.All_All.Frontal_Temporal.Ind_UniquePairings_persub{i_sub}];        
    end
    
    %Plot figure showing Number and location of each pair for each subject
    figure
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    DimSubplot = [3 3];
    sgtitle(['Possible Electrode Pairs - All to All - {\color{blue}Frontal} to {\color{red}Temporal}'])
    for i_sub = 1:length(subs)
        
        if any(strcmp(subs(i_sub),SelElecs.AllElecs{:,2}))
            if ~isnan(SelElecs.Pairs.All_All.Frontal_Temporal.Ind_UniquePairings_persub{i_sub})
                curr_elec_ind = unique(SelElecs.Pairs.All_All.Frontal_Temporal.Ind_UniquePairings_persub{i_sub});
            else
                curr_elec_ind = find(strcmp(subs(i_sub),SelElecs.AllElecs{:,2}));
            end
            elec_coords = SelElecs.AllElecs{curr_elec_ind,6};
            elec_coords(:,1) = abs(elec_coords(:,1)) * -1;%Project all electrodes on one hemisphere
            elec_labels = {};
            for i_elec = 1:size(elec_coords,1) %No electrode labels
                elec_labels{i_elec} = cell2mat(SelElecs.AllElecs{curr_elec_ind(i_elec),1});
            end
            elec_data = ones(size(elec_coords,1),1)*i_sub;
            sub_index = ones(size(elec_coords,1),1)*i_sub;
            elec_signindex = true(size(elec_coords,1),1);
            elec_typeindex = ones(size(elec_coords,1),1)*3;
            elec_typeindex(~cellfun('isempty',strfind(SelElecs.AllElecs{curr_elec_ind,4}, 'AntPFC'))) = 1; %1 = source, 2 = target
            elec_typeindex(~cellfun('isempty',strfind(SelElecs.AllElecs{curr_elec_ind,4}, 'VentralT'))) = 2; %1 = source, 2 = target
            elec_pairings = SelElecs.Pairs.All_All.Frontal_Temporal.Ind_UniquePairings_persub{i_sub};
            effect_comp = 'All_All';
            chanSize = ones(1,length(elec_data))*3; %electrode size (arbitrary)
            clims = [1 9];
            cmap  = ones(length(subs), 3).*[0.6207    0.3103    0.2759];
            textcolor_rgb   = [0 0 0];
            sp_title = {[subs{i_sub} ' - ' num2str(size(elec_coords,1)) ' elecs - ' ...
                num2str(SelElecs.Pairs.All_All.Frontal_Temporal.Num_UniquePairings_persub(i_sub)) ' elec pairs']};
            
            NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH_Pairs...
                (elec_coords, elec_labels, elec_data, elec_signindex, elec_typeindex, ...
                sub_index, elec_pairings, effect_comp, ...
                chanSize, clims, cmap, textcolor_rgb, ...
                DimSubplot, i_sub, [0 -0.02 0 0], [], sp_title);
        end
    end
    filename     = ['Surf1H_' ...
        'Allsubn' num2str(length(subs)) '_' ...
        'PairsAllElec_Front2Temp.png'];
    figfile      = [path_fig filename];
    saveas(gcf, figfile, 'png'); %save png version
    close;
    
    %% 3.3 Frontal - Parietal
    SelElecs.Pairs.All_All.Frontal_Parietal.Pairs_perelec = struct;
    SelElecs.Pairs.All_All.Frontal_Parietal.Ind_UniquePairings = [];
    
    for i_elec = 1:height(SelElecs.AllElecs)
        %Determine current elec information
        temp_label  = SelElecs.AllElecs{i_elec,1};
        temp_sub    = SelElecs.AllElecs{i_elec,2};
        temp_anatcat = SelElecs.AllElecs{i_elec,4};
        if strfind(temp_anatcat{1}, 'AntPFC') %If current electrode is frontal
            
            temp_index_otherelecs = ...%Find other sign. elecs from same subject in temporal areas
                find(strcmp(temp_sub, SelElecs.AllElecs{:,2}) & ...
                ~strcmp(temp_label, SelElecs.AllElecs{:,1}) & ...
                ~cellfun('isempty',strfind(SelElecs.AllElecs{:,4}, 'SupParLob')));
            if ~isempty(strfind(temp_anatcat{1}, ',')) %Create label for subfield
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' extractBefore(temp_anatcat{1}, ',')];
            else
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' temp_anatcat{1}];
            end
            
            if ~isempty(temp_index_otherelecs) %Store pair information in new struct
                SelElecs.Pairs.All_All.Frontal_Parietal.Pairs_perelec.(temp_fieldname) = ...
                    SelElecs.AllElecs(temp_index_otherelecs,:);
            else
                SelElecs.Pairs.All_All.Frontal_Parietal.Pairs_perelec.(temp_fieldname) = [];
            end
        end
    end
    %Determine number of possible pairings for each subject
    for i_sub = 1:length(subs)
        sub = subs(i_sub);
        sign_elecs_currsub1 = ...
            find(strcmp(sub, SelElecs.AllElecs{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.AllElecs{:,4}, 'AntPFC')));
        sign_elecs_currsub2 = ...
            find(strcmp(sub, SelElecs.AllElecs{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.AllElecs{:,4}, 'SupParLob')));
        %Determine number of possible pairings for each subject
        [p,q] = meshgrid(sign_elecs_currsub1, sign_elecs_currsub2);
        temp_pairs1 = [p(:) q(:)];
        if size(temp_pairs1,1) > 1
            %delete same-same pairings
            temp_pairs2 = [];
            for i_pairs = 1:length(temp_pairs1)
                if temp_pairs1(i_pairs,1) ~= temp_pairs1(i_pairs,2)
                    temp_pairs2 = [temp_pairs2; temp_pairs1(i_pairs,:)];
                end
            end
            %delete doubles
            i_doubles = false(length(temp_pairs2),1);
            i_compared = false(length(temp_pairs2),1);
            for i_pairs = 1:length(temp_pairs2)
                temp_doubles = find(all(temp_pairs2(i_pairs,:) == fliplr(temp_pairs2),2));
                i_compared(i_pairs) = true;
                
                if i_compared(temp_doubles) == false
                    i_doubles(temp_doubles) = true;
                end
            end
            SelElecs.Pairs.All_All.Frontal_Parietal.Num_UniquePairings_persub(1,i_sub) = ...
                size(temp_pairs2(~i_doubles,:),1);
            SelElecs.Pairs.All_All.Frontal_Parietal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs2(~i_doubles,:);
        elseif size(temp_pairs1,1) == 1 && temp_pairs1(1) ~= temp_pairs1(2)
            SelElecs.Pairs.All_All.Frontal_Parietal.Num_UniquePairings_persub(1,i_sub) = ...
                1;
            SelElecs.Pairs.All_All.Frontal_Parietal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs1;
        else
            SelElecs.Pairs.All_All.Frontal_Parietal.Num_UniquePairings_persub(1,i_sub) = ...
                0;
            SelElecs.Pairs.All_All.Frontal_Parietal.Ind_UniquePairings_persub{i_sub} = ...
                nan;
        end
        %aggregate pairings across subjects
        SelElecs.Pairs.All_All.Frontal_Parietal.Ind_UniquePairings = ...
            [SelElecs.Pairs.All_All.Frontal_Parietal.Ind_UniquePairings; ...
            SelElecs.Pairs.All_All.Frontal_Parietal.Ind_UniquePairings_persub{i_sub}];        
    end
    
    %Plot figure showing Number and location of each pair for each subject
    figure
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    DimSubplot = [3 3];
    sgtitle(['Possible Electrode Pairs - All to All - {\color{blue}Frontal} to {\color{green}Parietal}'])
    for i_sub = 1:length(subs)
        
        if any(strcmp(subs(i_sub),SelElecs.AllElecs{:,2}))
            if ~isnan(SelElecs.Pairs.All_All.Frontal_Parietal.Ind_UniquePairings_persub{i_sub})
                curr_elec_ind = unique(SelElecs.Pairs.All_All.Frontal_Parietal.Ind_UniquePairings_persub{i_sub});
            else
                curr_elec_ind = find(strcmp(subs(i_sub),SelElecs.AllElecs{:,2}));
            end
            elec_coords = SelElecs.AllElecs{curr_elec_ind,6};
            elec_coords(:,1) = abs(elec_coords(:,1)) * -1;%Project all electrodes on one hemisphere
            elec_labels = {};
            for i_elec = 1:size(elec_coords,1) %No electrode labels
                elec_labels{i_elec} = cell2mat(SelElecs.AllElecs{curr_elec_ind(i_elec),1});
            end
            elec_data = ones(size(elec_coords,1),1)*i_sub;
            sub_index = ones(size(elec_coords,1),1)*i_sub;
            elec_signindex = true(size(elec_coords,1),1);
            elec_typeindex = ones(size(elec_coords,1),1)*3;
            elec_typeindex(~cellfun('isempty',strfind(SelElecs.AllElecs{curr_elec_ind,4}, 'AntPFC'))) = 1; %1 = source, 2 = target
            elec_typeindex(~cellfun('isempty',strfind(SelElecs.AllElecs{curr_elec_ind,4}, 'SupParLob'))) = 4; %1 = source, 2 = target
            elec_pairings = SelElecs.Pairs.All_All.Frontal_Parietal.Ind_UniquePairings_persub{i_sub};
            effect_comp = 'All_All';
            chanSize = ones(1,length(elec_data))*3; %electrode size (arbitrary)
            clims = [1 9];
            cmap  = ones(length(subs), 3).*[0.6207    0.3103    0.2759];
            textcolor_rgb   = [0 0 0];
            sp_title = {[subs{i_sub} ' - ' num2str(size(elec_coords,1)) ' elecs - ' ...
                num2str(SelElecs.Pairs.All_All.Frontal_Parietal.Num_UniquePairings_persub(i_sub)) ' elec pairs']};
            
            NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH_Pairs...
                (elec_coords, elec_labels, elec_data, elec_signindex, elec_typeindex, ...
                sub_index, elec_pairings, effect_comp, ...
                chanSize, clims, cmap, textcolor_rgb, ...
                DimSubplot, i_sub, [0 -0.02 0 0], [], sp_title);
        end
    end
    filename     = ['Surf1H_' ...
        'Allsubn' num2str(length(subs)) '_' ...
        'PairsAllElec_Front2Par.png'];
    figfile      = [path_fig filename];
    saveas(gcf, figfile, 'png'); %save png version
    close;
 
    %% 3.4 Parietal - Temporal
    SelElecs.Pairs.All_All.Parietal_Temporal.Pairs_perelec = struct;
    SelElecs.Pairs.All_All.Parietal_Temporal.Ind_UniquePairings = [];

    for i_elec = 1:height(SelElecs.AllElecs)
        %Determine current elec information
        temp_label  = SelElecs.AllElecs{i_elec,1};
        temp_sub    = SelElecs.AllElecs{i_elec,2};
        temp_anatcat = SelElecs.AllElecs{i_elec,4};
        if strfind(temp_anatcat{1}, 'SupParLob') %If current electrode is frontal
            
            temp_index_otherelecs = ...%Find other sign. elecs from same subject in temporal areas
                find(strcmp(temp_sub, SelElecs.AllElecs{:,2}) & ...
                ~strcmp(temp_label, SelElecs.AllElecs{:,1}) & ...
                ~cellfun('isempty',strfind(SelElecs.AllElecs{:,4}, 'VentralT')));
            if ~isempty(strfind(temp_anatcat{1}, ',')) %Create label for subfield
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' extractBefore(temp_anatcat{1}, ',')];
            else
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' temp_anatcat{1}];
            end
            
            if ~isempty(temp_index_otherelecs) %Store pair information in new struct
                SelElecs.Pairs.All_All.Parietal_Temporal.Pairs_perelec.(temp_fieldname) = ...
                    SelElecs.AllElecs(temp_index_otherelecs,:);
            else
                SelElecs.Pairs.All_All.Parietal_Temporal.Pairs_perelec.(temp_fieldname) = [];
            end
        end
    end
    %Determine number of possible pairings for each subject
    for i_sub = 1:length(subs)
        sub = subs(i_sub);
        sign_elecs_currsub1 = ...
            find(strcmp(sub, SelElecs.AllElecs{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.AllElecs{:,4}, 'SupParLob')));
        sign_elecs_currsub2 = ...
            find(strcmp(sub, SelElecs.AllElecs{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.AllElecs{:,4}, 'VentralT')));
        %Determine number of possible pairings for each subject
        [p,q] = meshgrid(sign_elecs_currsub1, sign_elecs_currsub2);
        temp_pairs1 = [p(:) q(:)];
        if size(temp_pairs1,1) > 1
            %delete same-same pairings
            temp_pairs2 = [];
            for i_pairs = 1:length(temp_pairs1)
                if temp_pairs1(i_pairs,1) ~= temp_pairs1(i_pairs,2)
                    temp_pairs2 = [temp_pairs2; temp_pairs1(i_pairs,:)];
                end
            end
            %delete doubles
            i_doubles = false(length(temp_pairs2),1);
            i_compared = false(length(temp_pairs2),1);
            for i_pairs = 1:length(temp_pairs2)
                temp_doubles = find(all(temp_pairs2(i_pairs,:) == fliplr(temp_pairs2),2));
                i_compared(i_pairs) = true;
                
                if i_compared(temp_doubles) == false
                    i_doubles(temp_doubles) = true;
                end
            end
            SelElecs.Pairs.All_All.Parietal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                size(temp_pairs2(~i_doubles,:),1);
            SelElecs.Pairs.All_All.Parietal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs2(~i_doubles,:);
        elseif size(temp_pairs1,1) == 1 && temp_pairs1(1) ~= temp_pairs1(2)
            SelElecs.Pairs.All_All.Parietal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                1;
            SelElecs.Pairs.All_All.Parietal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs1;
        else
            SelElecs.Pairs.All_All.Parietal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                0;
            SelElecs.Pairs.All_All.Parietal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                nan;
        end
        %aggregate pairings across subjects
        SelElecs.Pairs.All_All.Parietal_Temporal.Ind_UniquePairings = ...
            [SelElecs.Pairs.All_All.Parietal_Temporal.Ind_UniquePairings; ...
            SelElecs.Pairs.All_All.Parietal_Temporal.Ind_UniquePairings_persub{i_sub}];        
    end
    
    %Plot figure showing Number and location of each pair for each subject
    figure
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    DimSubplot = [3 3];
    sgtitle(['Possible Electrode Pairs - All to All - {\color{green}Parietal} to {\color{red}Temporal}'])
    for i_sub = 1:length(subs)
        
        if any(strcmp(subs(i_sub),SelElecs.AllElecs{:,2}))
            if ~isnan(SelElecs.Pairs.All_All.Parietal_Temporal.Ind_UniquePairings_persub{i_sub})
                curr_elec_ind = unique(SelElecs.Pairs.All_All.Parietal_Temporal.Ind_UniquePairings_persub{i_sub});
            else
                curr_elec_ind = find(strcmp(subs(i_sub),SelElecs.AllElecs{:,2}));
            end
            elec_coords = SelElecs.AllElecs{curr_elec_ind,6};
            elec_coords(:,1) = abs(elec_coords(:,1)) * -1;%Project all electrodes on one hemisphere
            elec_labels = {};
            for i_elec = 1:size(elec_coords,1) %No electrode labels
                elec_labels{i_elec} = cell2mat(SelElecs.AllElecs{curr_elec_ind(i_elec),1});
            end
            elec_data = ones(size(elec_coords,1),1)*i_sub;
            sub_index = ones(size(elec_coords,1),1)*i_sub;
            elec_signindex = true(size(elec_coords,1),1);
            elec_typeindex = ones(size(elec_coords,1),1)*3;
            elec_typeindex(~cellfun('isempty',strfind(SelElecs.AllElecs{curr_elec_ind,4}, 'SupParLob'))) = 4; %1 = source, 2 = target
            elec_typeindex(~cellfun('isempty',strfind(SelElecs.AllElecs{curr_elec_ind,4}, 'VentralT'))) = 2; %1 = source, 2 = target
            elec_pairings = SelElecs.Pairs.All_All.Parietal_Temporal.Ind_UniquePairings_persub{i_sub};
            effect_comp = 'All_All';
            chanSize = ones(1,length(elec_data))*3; %electrode size (arbitrary)
            clims = [1 9];
            cmap  = ones(length(subs), 3).*[0.6207    0.3103    0.2759];
            textcolor_rgb   = [0 0 0];
            sp_title = {[subs{i_sub} ' - ' num2str(size(elec_coords,1)) ' elecs - ' ...
                num2str(SelElecs.Pairs.All_All.Parietal_Temporal.Num_UniquePairings_persub(i_sub)) ' elec pairs']};
            
            NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH_Pairs...
                (elec_coords, elec_labels, elec_data, elec_signindex, elec_typeindex, ...
                sub_index, elec_pairings, effect_comp, ...
                chanSize, clims, cmap, textcolor_rgb, ...
                DimSubplot, i_sub, [0 -0.02 0 0], [], sp_title);
        end
    end
    filename     = ['Surf1H_' ...
        'Allsubn' num2str(length(subs)) '_' ...
        'PairsAllElec_Par2Temp.png'];
    figfile      = [path_fig filename];
    saveas(gcf, figfile, 'png'); %save png version
    close;    
    
end


%% 4. Plot group-level pairings per anat selection
    figure
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    DimSubplot = [1 3];
    sgtitle(['All subs (n = ' num2str(length(subs)) ') - All to All'])
    
    effect_comp = 'All_All';
    cmap  = distinguishable_colors(3);
    clims = [1 3];
    textcolor_rgb   = [0 0 0];
    
    %Frontal - Temporal
    curr_elec_ind = [];
    for i_sub = 1:length(subs)
        if any(strcmp(subs(i_sub),SelElecs.AllElecs{:,2}))
            if ~isnan(SelElecs.Pairs.All_All.Frontal_Temporal.Ind_UniquePairings_persub{i_sub})
                curr_elec_ind = [curr_elec_ind; ...
                    unique(SelElecs.Pairs.All_All.Frontal_Temporal.Ind_UniquePairings_persub{i_sub})];
            end
        end
    end
    elec_coords = SelElecs.AllElecs{curr_elec_ind,6};
    elec_coords(:,1) = abs(elec_coords(:,1)) * -1;%Project all electrodes on one hemisphere
    elec_labels = {};
    for i_elec = 1:size(elec_coords,1) %No electrode labels
        elec_labels{i_elec} = cell2mat(SelElecs.AllElecs{curr_elec_ind(i_elec),1});
    end
    elec_data = ones(size(elec_coords,1),1);
    elec_signindex = true(size(elec_coords,1),1);
    elec_typeindex = ones(size(elec_coords,1),1)*3;
    elec_typeindex(~cellfun('isempty',strfind(SelElecs.AllElecs{curr_elec_ind,4}, 'AntPFC'))) = 1;
    elec_typeindex(~cellfun('isempty',strfind(SelElecs.AllElecs{curr_elec_ind,4}, 'VentralT'))) = 2;
    elec_pairings = SelElecs.Pairs.All_All.Frontal_Temporal.Ind_UniquePairings;
    elec_pairings = elec_pairings(~isnan(elec_pairings(:,1)),:);
    
    chanSize = ones(1,length(elec_data))*3; %electrode size (arbitrary)
    sp_title = {['{\color{blue}Frontal} to {\color{red}Temporal}  {\color{black}electrodes} '], ...
        [ num2str(size(elec_coords,1)) ' All elecs - ' ...
        num2str(sum(SelElecs.Pairs.All_All.Frontal_Temporal.Num_UniquePairings_persub)) ...
        ' elec pairs']};
    
    NASTD_ECoG_Plot_SubplotSignElecsSurf_ColorAnat_LH_Pairs...
        (elec_coords, elec_labels, elec_typeindex, elec_signindex, ...
        elec_pairings, effect_comp, ...
        chanSize, clims, cmap, textcolor_rgb, ...
        DimSubplot, 1, [0 -0.02 0 0], [], sp_title);
    clear elec*
    
    %Frontal - Parietal
    curr_elec_ind = [];
    for i_sub = 1:length(subs)
        if any(strcmp(subs(i_sub),SelElecs.AllElecs{:,2}))
            if ~isnan(SelElecs.Pairs.All_All.Frontal_Parietal.Ind_UniquePairings_persub{i_sub})
                curr_elec_ind = [curr_elec_ind; ...
                    unique(SelElecs.Pairs.All_All.Frontal_Parietal.Ind_UniquePairings_persub{i_sub})];
            end
        end
    end
    elec_coords = SelElecs.AllElecs{curr_elec_ind,6};
    elec_coords(:,1) = abs(elec_coords(:,1)) * -1;%Project all electrodes on one hemisphere
    elec_labels = {};
    for i_elec = 1:size(elec_coords,1) %No electrode labels
        elec_labels{i_elec} = cell2mat(SelElecs.AllElecs{curr_elec_ind(i_elec),1});
    end
    elec_data = ones(size(elec_coords,1),1);
    elec_signindex = true(size(elec_coords,1),1);
    elec_typeindex = ones(size(elec_coords,1),1)*3;
    elec_typeindex(~cellfun('isempty',strfind(SelElecs.AllElecs{curr_elec_ind,4}, 'AntPFC'))) = 1;
    elec_typeindex(~cellfun('isempty',strfind(SelElecs.AllElecs{curr_elec_ind,4}, 'SupParLob'))) = 3;
    elec_pairings = SelElecs.Pairs.All_All.Frontal_Parietal.Ind_UniquePairings;
    elec_pairings = elec_pairings(~isnan(elec_pairings(:,1)),:);
    
    chanSize = ones(1,length(elec_data))*3; %electrode size (arbitrary)    
    sp_title = {['{\color{blue}Frontal} to {\color{green}Parietal}  {\color{black}electrodes} '], ...
        [ num2str(size(elec_coords,1)) ' All elecs - ' ...
        num2str(sum(SelElecs.Pairs.All_All.Frontal_Parietal.Num_UniquePairings_persub)) ...
        ' elec pairs']};
    
    NASTD_ECoG_Plot_SubplotSignElecsSurf_ColorAnat_LH_Pairs...
        (elec_coords, elec_labels, elec_typeindex, elec_signindex, ...
        elec_pairings, effect_comp, ...
        chanSize, clims, cmap, textcolor_rgb, ...
        DimSubplot, 2, [0 -0.02 0 0], [], sp_title);
    clear elec*
    
    %Parietal - Temporal
    curr_elec_ind = [];
    for i_sub = 1:length(subs)
        if any(strcmp(subs(i_sub),SelElecs.AllElecs{:,2}))
            if ~isnan(SelElecs.Pairs.All_All.Parietal_Temporal.Ind_UniquePairings_persub{i_sub})
                curr_elec_ind = [curr_elec_ind; ...
                    unique(SelElecs.Pairs.All_All.Parietal_Temporal.Ind_UniquePairings_persub{i_sub})];
            end
        end
    end
    elec_coords = SelElecs.AllElecs{curr_elec_ind,6};
    elec_coords(:,1) = abs(elec_coords(:,1)) * -1;%Project all electrodes on one hemisphere
    elec_labels = {};
    for i_elec = 1:size(elec_coords,1) %No electrode labels
        elec_labels{i_elec} = cell2mat(SelElecs.AllElecs{curr_elec_ind(i_elec),1});
    end
    elec_data = ones(size(elec_coords,1),1);
    elec_signindex = true(size(elec_coords,1),1);
    elec_typeindex = ones(size(elec_coords,1),1)*3;
    elec_typeindex(~cellfun('isempty',strfind(SelElecs.AllElecs{curr_elec_ind,4}, 'SupParLob'))) = 3;
    elec_typeindex(~cellfun('isempty',strfind(SelElecs.AllElecs{curr_elec_ind,4}, 'VentralT'))) = 2;
    elec_pairings = SelElecs.Pairs.All_All.Parietal_Temporal.Ind_UniquePairings;
    elec_pairings = elec_pairings(~isnan(elec_pairings(:,1)),:);
    
    chanSize = ones(1,length(elec_data))*3; %electrode size (arbitrary)  
    sp_title = {['{\color{green}Parietal} to {\color{red}Temporal}  {\color{black}electrodes} '], ...
        [ num2str(size(elec_coords,1)) ' All elecs - ' ...
        num2str(sum(SelElecs.Pairs.All_All.Parietal_Temporal.Num_UniquePairings_persub)) ...
        ' elec pairs']};
    
    NASTD_ECoG_Plot_SubplotSignElecsSurf_ColorAnat_LH_Pairs...
        (elec_coords, elec_labels, elec_typeindex, elec_signindex, ...
        elec_pairings, effect_comp, ...
        chanSize, clims, cmap, textcolor_rgb, ...
        DimSubplot, 3, [0 -0.02 0 0], [], sp_title);
    clear elec*
    

filename     = ['Surf1H_' ...
    'Allsubn' num2str(length(subs)) '_' ...
    'PairsAllElec.png'];
figfile      = [path_fig filename];
saveas(gcf, figfile, 'png'); %save png version
close;

end