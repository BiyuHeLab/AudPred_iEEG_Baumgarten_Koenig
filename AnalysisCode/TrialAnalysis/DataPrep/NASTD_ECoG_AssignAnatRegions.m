function AnatReg = NASTD_ECoG_AssignAnatRegions...
    (DataClean_AllTrials, labels_loadedData)

%Aim: To categorize T1 electrode labels into higher-level anatomical categories

%1) Determine .elec subfield index for current electrodes
ind_selelecs_MNI = [];
for i_elec = 1:length(labels_loadedData)
    ind_selelecs_MNI = ...
        [ind_selelecs_MNI; ...
        find(strcmp(DataClean_AllTrials.elec.label, ...
        labels_loadedData{i_elec}))];
end

%2) Select respective labels and corresponding anat region description
AnatReg.Info_perelec = cell(length(ind_selelecs_MNI), 5);
for i_elec = 1:length(ind_selelecs_MNI)
    AnatReg.Info_perelec{i_elec,1} = ...
        DataClean_AllTrials.elec.label{ind_selelecs_MNI(i_elec)};
    AnatReg.Info_perelec{i_elec,2} = ...
        DataClean_AllTrials.elec.T1AnatLabel{ind_selelecs_MNI(i_elec)};
end

%3) Categorize electrodes in anatomical regions
AnatReg.Search.selLabels = ... %Labels and indices used to identify electrode locations
    {'frontal', 1; 'precentral', 2; 'opercularis', 3; 'triangularis', 4; 'orbitalis', 5; ... %Frontal
    'STG', 6; 'MTG', 7; 'temporalpole', 8; 'temporal', 9; 'fusiform', 10; 'entorhinal' , 11; 'bankssts', 12; 'parahippocampal', 13; ... %Temporal
    'parietal', 14; 'postcentral', 15; 'supramarginal', 16; 'cuneus', 17; ... %Parietal
    'occipital', 18; 'calcarine', 19}; %Occipital
AnatReg.Search.CatLabels = ... %Anatomical categories we use to differentiate regions
    {'AntPFC', [1 2 3 4 5]; 'PrecentralG', 2; 'IFG', [3 4 5]; ...
    'VentralT', [6 7 8 9 10 11 12 13]; 'STG', [6, 12]; 'MTG', 7; ...
    'SupParLob', [14 15 16 17]; 'PostcentralG', 15; 'SupramarginalG', 16; ...
    'OccipitalL', [18 19]};

%4) Scan each electrode and check if its name contains an anatomical label
for i_elec = 1:length(AnatReg.Info_perelec)
    temp_currlabels  = [];
    temp_currindx    = [];
    for i_anatreg = 1:length(AnatReg.Search.selLabels)
        if contains(...
                AnatReg.Info_perelec{i_elec,2}, ...
                AnatReg.Search.selLabels{i_anatreg,1})
            temp_currlabels = ...
                [temp_currlabels, ' ', AnatReg.Search.selLabels{i_anatreg,1}];
            temp_currindx = ...
                [temp_currindx, AnatReg.Search.selLabels{i_anatreg,2}];
        end
        AnatReg.Info_perelec{i_elec,4} = temp_currlabels;
        AnatReg.Info_perelec{i_elec,5} = temp_currindx;
    end
    if isempty(temp_currlabels)
        disp([' -- No fitting anatomical label found for: ' ...
            AnatReg.Info_perelec{i_elec,1} ' - ' ...
            AnatReg.Info_perelec{i_elec,2} ' -- ']);
    end
end

%5) Depending on the found label, place electrode in respective category
for i_elec = 1:length(AnatReg.Info_perelec)
    temp_category = cell(1, length(AnatReg.Search.CatLabels));
    for i_category = 1:length(AnatReg.Search.CatLabels)
        for i_entry = 1:length(AnatReg.Search.CatLabels{i_category,2})
            if any(AnatReg.Info_perelec{i_elec,5} == ...
                    AnatReg.Search.CatLabels{i_category,2}(i_entry))
                temp_category{i_category} = ...
                    AnatReg.Search.CatLabels{i_category,1};
            end
        end
    end
    temp_category2 = [];
    for i_entry = 1:length(temp_category)
        if ~isempty(temp_category{i_entry})
            if isempty(temp_category2)
                temp_category2 = ...
                    [temp_category{i_entry}];
            else
                temp_category2 = ...
                    [temp_category2 ', ' temp_category{i_entry}];
            end
        end
    end
    if isempty(temp_category2)
        disp([' -- No fitting anatomical category found for: ' ...
            AnatReg.Info_perelec{i_elec,1} ' - ' ...
            AnatReg.Info_perelec{i_elec,2} ' -- ']);
    end
    AnatReg.Info_perelec{i_elec,3} = temp_category2;
end

%6) Provide summary output
temp_table1                  = cell2table(AnatReg.Info_perelec);
AnatReg.CatLabels            = cell(size(AnatReg.Search.CatLabels,1),1);
for i_catlabel = 1:length(AnatReg.CatLabels)
    AnatReg.CatLabels{i_catlabel} = AnatReg.Search.CatLabels{i_catlabel,1};
end
% AnatReg.CatLabels            = table2cell(temp_table2);
AnatReg.Num_SelElecs         = length(AnatReg.Info_perelec);
AnatReg.Num_Cat  = size(AnatReg.CatLabels, 1);
AnatReg.Num_ElecsperCat = nan(size(AnatReg.CatLabels, 1),1);
AnatReg.CatIndex             = nan(length(AnatReg.Info_perelec),5);

%For every category, check how many electrodes are in there
for i_category = 1:AnatReg.Num_Cat
    
    AnatReg.Num_ElecsperCat(i_category) = ... %Count label-matching electrodes
        [sum(contains(...
        temp_table1{:,3}, AnatReg.CatLabels{i_category}))];
end

for i_elec = 1:length(AnatReg.Info_perelec)
    placing_counter = 1;
    for i_category = 1:AnatReg.Num_Cat
        
        if contains(temp_table1{i_elec,3}, AnatReg.CatLabels{i_category})
            AnatReg.CatIndex(i_elec, placing_counter) = ... %copy fitting category indices
                i_category;
            placing_counter = placing_counter + 1;
        end
    end
end

% %Optional: Modify .CatIndex to ensure color plot grouping of similar categories
% for i_row = 1:size(AnatReg.CatIndex,1)
%     for i_column = 1:size(AnatReg.CatIndex,2)
%         if AnatReg.CatIndex(i_row, i_column) == 2 
%             AnatReg.CatIndex(i_row, i_column) = 1.3;
%         elseif AnatReg.CatIndex(i_row, i_column) == 3
%             AnatReg.CatIndex(i_row, i_column) = 1.6;
%         elseif AnatReg.CatIndex(i_row, i_column) == 4
%             AnatReg.CatIndex(i_row, i_column) = 3;
%         elseif AnatReg.CatIndex(i_row, i_column) == 5
%             AnatReg.CatIndex(i_row, i_column) = 3.3;            
%         elseif AnatReg.CatIndex(i_row, i_column) == 6
%             AnatReg.CatIndex(i_row, i_column) = 3.6;            
%         elseif AnatReg.CatIndex(i_row, i_column) == 7
%             AnatReg.CatIndex(i_row, i_column) = 5;            
%         elseif AnatReg.CatIndex(i_row, i_column) == 8
%             AnatReg.CatIndex(i_row, i_column) = 5.3;            
%         elseif AnatReg.CatIndex(i_row, i_column) == 9
%             AnatReg.CatIndex(i_row, i_column) = 5.6;
%         elseif AnatReg.CatIndex(i_row, i_column) == 10
%             AnatReg.CatIndex(i_row, i_column) = 7;
%         end
%     end
% end
            
            
temp_table2(:,1)    = cell2table(AnatReg.CatLabels);
temp_table2(:,2)    = table(AnatReg.Num_ElecsperCat);
temp_table2.Properties.VariableNames = {'CategoryLabel', 'NumElectrodes'};
AnatReg.Cat         = temp_table2;

% disp([' -- Anatomical categrization performed for: ' ...
%     num2str(AnatReg.Num_SelElecs) ' electrodes and ' num2str(AnatReg.Num_Cat) ' categories -- '])
% disp(AnatReg.Cat)

AnatReg = orderfields(AnatReg);

clear temp*