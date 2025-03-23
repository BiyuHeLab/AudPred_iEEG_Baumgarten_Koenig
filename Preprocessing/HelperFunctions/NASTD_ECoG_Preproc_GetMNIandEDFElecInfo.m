function [Elec_info, SelectElec_info]= ...
    NASTD_ECoG_Preproc_GetMNIandEDFElecInfo ...
    (sub, MNIelec_table, T1ElecAnat_table, LabelEDF_all)

%Reads info from MNI(location), T1 (location & anatomical label),
%and EDF(recording) files. String-compares and matches the electrode labels.
%From matching labels, strip and grid electrodes are selected for further processing. 
%Output is a summary file with electrode label, type, number info for 
%1) all electrodes and 2) selected electrodes

Elec_info = struct;

%% 1.1) Prepare MNI data
%1.1.1 Transfer table info to respective subfields
Elec_info.NumMNI_all = size(MNIelec_table,1); 
%determine number of localized electrodes based on .txt file
Elec_info.LabelMNI_all = MNIelec_table{:,1}; 
%determine corresponding electrode label
Elec_info.LabelMNI_adjusted_all = MNIelec_table{:,1};
Elec_info.Type_all = MNIelec_table{:,5}; 
%determine electrode type (G = grid, S = strip, D = depth)

for i_elec = 1:size(T1ElecAnat_table,1)
    if strcmp(Elec_info.LabelMNI_all{i_elec}, T1ElecAnat_table{i_elec,1})
        Elec_info.T1AnatLabel_MNIorder = T1ElecAnat_table{:,5};
    else
            disp(['Mismatch in MNI and T1 channel labeling for ' ...
                Elec_info.LabelMNI_all{i_elec} '(MNI) and '...
                T1ElecAnat_table{i_elec,1} '(T1)'])
            pause
    end
end
%determine corresponding anatomical label


%Adjust MNI label names
%(MNI labels have an additional 0 before digit numbers. Delete only those 
%0 to enable EDF-MNI electrode label matching.)
for i_MNILabel = 1:Elec_info.NumMNI_all
    pos_zero = strfind(Elec_info.LabelMNI_adjusted_all{i_MNILabel,1},'0'); 
    %find location of all zeros in label name
    
    pos_nondecimal_zero = []; %restrict zeros to nondecimals
    for i_zeros = 1:length(pos_zero)
        if isempty(regexp(Elec_info.LabelMNI_adjusted_all{i_MNILabel,1}...
                (pos_zero(i_zeros)-1),'[1-9]','match')) 
            %if there is no non-zero digit before the zero (e.g. 10)
            if pos_zero(i_zeros) ~= length(Elec_info.LabelMNI_adjusted_all{i_MNILabel,1}) 
                %if the current 0 is not at the end of the label (e.g., 100)
                pos_nondecimal_zero = [pos_nondecimal_zero, pos_zero(i_zeros)];
            end
        end
    end
    
    Elec_info.LabelMNI_adjusted_all{i_MNILabel,1}(pos_nondecimal_zero) = []; 
    %remove nondecimal zeros
    
end
clear pos*

%     %Get quick overview to figure out which labels mismatch between MNI and EDF file
%     [Elec_info.LabelMNI_adjusted_all, LabelEDF_all(1:size(Elec_info.LabelMNI_adjusted_all(:,1),1))]

%Individual label adjustments
%Changes MNI labels which differ from EDF labels to the corresponding EDF form.
Elec_info = NASTD_ECoG_Preproc_IndivMNI2EDFLabelAdjust(sub,Elec_info);

%Count who many electrodes of each type are present
Elec_info.NumStrip_all = 0;
Elec_info.NumGrid_all = 0;
Elec_info.NumDepth_all = 0;
Elec_info.NumUndef_all = 0;

for i_elec = 1:length(Elec_info.Type_all)
    if strcmp(Elec_info.Type_all(i_elec),'S') %Strip
        index_selectStripGridElec_info(i_elec) = 1;
        Elec_info.NumStrip_all = Elec_info.NumStrip_all+1;
    elseif strcmp(Elec_info.Type_all(i_elec),'G') %Grid
        index_selectStripGridElec_info(i_elec) = 1;
        Elec_info.NumGrid_all = Elec_info.NumGrid_all+1;
    elseif strcmp(Elec_info.Type_all(i_elec),'EG') %Grid
        index_selectStripGridElec_info(i_elec) = 1;
        Elec_info.NumGrid_all = Elec_info.NumGrid_all+1;
    elseif strcmp(Elec_info.Type_all(i_elec),'D') %Depth
        index_selectStripGridElec_info(i_elec) = 0;
        Elec_info.NumDepth_all = Elec_info.NumDepth_all+1;
    else
        index_selectStripGridElec_info(i_elec) = NaN;
        Elec_info.NumUndef_all = Elec_info.NumUndef_all+1;
    end
end

%% 1.2) EDF
Elec_info.LabelEDF_all = LabelEDF_all;

%% 2) Match electrode labels between .EDF and .txt files 
%(i.e., determine which channels are present in both files)

%2.1) Find MNI indices corresponding to EDF elecs and corresponding labels (Original EDF to adjusted MNI)
[Elec_info.EDF2MNI_ElecLabelmatch, Elec_info.EDF2MNI_NOElecLabelmatch] = ...
    NASTD_ECoG_Preproc_MatchElecLabels...
    (sub, Elec_info.LabelEDF_all, Elec_info.LabelMNI_adjusted_all);

Elec_info.Num_EDFMNImatch = sum((~isnan(Elec_info.EDF2MNI_ElecLabelmatch)));
Elec_info.Num_EDFMNINOmatch = sum((isnan(Elec_info.EDF2MNI_ElecLabelmatch)));

%Display matching and non-matching electrodes
%TASK: Manually check if new labels agree
disp('TASK 3 (manual): Check manually if listed EDF-to-MNI labels match.');
pause
for i_channel = 1:length(Elec_info.EDF2MNI_ElecLabelmatch)
    if isnan(Elec_info.EDF2MNI_ElecLabelmatch(i_channel))
        disp(['EDFChannel ' num2str(i_channel) ': ' ...
            Elec_info.LabelEDF_all{i_channel} ' - No MNI coordinates found']);
    else
        disp(['EDFChannel ' num2str(i_channel) ': ' ...
            Elec_info.LabelEDF_all{i_channel} ' - original MNILabel ' ...
            Elec_info.LabelMNI_all{Elec_info.EDF2MNI_ElecLabelmatch(i_channel)}]);
    end
end
pause 

disp('TASK 4 (manual): Check EDF channels without corresponding MNI-match.');
pause
%     disp([Elec_info.LabelEDF_all(end-(subs_PreProcSettings.(sub).number_nonEEGchan-1):end)']);
disp([Elec_info.EDF2MNI_NOElecLabelmatch']);
pause

%2.2) Determine indices of matching electrodes for both EDF and MNI files
Elec_info.index4EDF_matchingElec_info = ...
    find(~isnan(Elec_info.EDF2MNI_ElecLabelmatch)); 
%for each EDF elec, read out the corresponding MNI-channel index

Elec_info.index4MNI_matchingElec_info = ...
    Elec_info.EDF2MNI_ElecLabelmatch(Elec_info.index4EDF_matchingElec_info); 
%for each MNI elec, read out the corresponding EDF-channel index

%2.3) Determine EDF and MNI labels for matching electrodes
Elec_info.LabelEDF_match = ...
    Elec_info.LabelEDF_all(...
    Elec_info.index4EDF_matchingElec_info); 
%gives EDF-labels for matching electrodes

Elec_info.LabelMNI_adjusted_match = ...
    Elec_info.LabelMNI_adjusted_all(...
    Elec_info.EDF2MNI_ElecLabelmatch(...
    ~isnan(Elec_info.EDF2MNI_ElecLabelmatch))); 
%gives MNI-labels for matching electrodes

%2.4) Ensure same order for both labels (align MNI order to EDF order)
[~, order] = ismember(Elec_info.LabelEDF_match, Elec_info.LabelMNI_adjusted_match);
disp('TASK 5 (manual): Check if EDF and MNI labels are identically ordered.');
pause
[Elec_info.LabelEDF_match, Elec_info.LabelMNI_adjusted_match] %displays parallel ordered arrays - should match
if ~isequal(Elec_info.LabelEDF_match, Elec_info.LabelMNI_adjusted_match) %checks is ordered arrays are the same
    disp('EDF - MNI label mismatch !!');
end
pause

%2.5) Align order for elec typ as well (align MNI type info to EDF order)
Elec_info.TypeMNI_match = Elec_info.Type_all(Elec_info.index4MNI_matchingElec_info); 

%2.6) Align order of anatomical labels as well (align MNI-type order to EDF order)
Elec_info.T1AnatLabel_EDForder = Elec_info.T1AnatLabel_MNIorder(Elec_info.index4MNI_matchingElec_info); 

%2.6) separately save labels of non-matching electrodes
Elec_info.LabelEDF_NOmatch = Elec_info.EDF2MNI_NOElecLabelmatch;

%% 3) Select strip & grid electrodes (and get rid of depth electrodes)
%3.1) Determine index for electrode subtypes and count respective subtypes
Elec_info.NumStrip_match = 0;
Elec_info.NumGrid_match = 0;
Elec_info.NumDepth_match = 0;
Elec_info.NumUndef_match = 0;

for i_elec = 1:length(Elec_info.TypeMNI_match)
    if strcmp(Elec_info.TypeMNI_match(i_elec),'S') %Strip
        Elec_info.index4EDF_selectStripGridElec_info(i_elec) = 1;
        Elec_info.NumStrip_match = Elec_info.NumStrip_match + 1;
    elseif strcmp(Elec_info.TypeMNI_match(i_elec),'G') %Grid
        Elec_info.index4EDF_selectStripGridElec_info(i_elec) = 1;
        Elec_info.NumGrid_match = Elec_info.NumGrid_match + 1;
    elseif strcmp(Elec_info.TypeMNI_match(i_elec),'EG') %Smaller Grid
        Elec_info.index4EDF_selectStripGridElec_info(i_elec) = 1;
        Elec_info.NumGrid_match = Elec_info.NumGrid_match + 1;
    elseif strcmp(Elec_info.TypeMNI_match(i_elec),'D') %Depth
        Elec_info.index4EDF_selectStripGridElec_info(i_elec) = 0;
        Elec_info.NumDepth_match = Elec_info.NumDepth_match + 1;
    else
        Elec_info.index4EDF_selectStripGridElec_info(i_elec) = NaN; %Not determinable
        Elec_info.NumUndef_match = Elec_info.NumUndef_match +1 ;
    end
end
Elec_info.index4EDF_selectStripGridElec_info = ...
    logical(Elec_info.index4EDF_selectStripGridElec_info);

%% 4) Copy selected electrodes (matching strip+grid) info into separate subfield
SelectElec_info = [];
SelectElec_info.Label = Elec_info.LabelEDF_match(Elec_info.index4EDF_selectStripGridElec_info);
SelectElec_info.Type = Elec_info.TypeMNI_match(Elec_info.index4EDF_selectStripGridElec_info);
SelectElec_info.index4EDF = Elec_info.index4EDF_matchingElec_info(Elec_info.index4EDF_selectStripGridElec_info);
SelectElec_info.index4MNI = Elec_info.index4MNI_matchingElec_info(Elec_info.index4EDF_selectStripGridElec_info);
SelectElec_info.T1AnatLabel_EDForder = Elec_info.T1AnatLabel_EDForder(Elec_info.index4EDF_selectStripGridElec_info);

%% 5) Determine trigger and ECG channels
%5.1) Find Trigger channel
if strcmp(sub,'NY798')
  Elec_info.index_triggerchan = ...
      find(strcmp(Elec_info.LabelEDF_all,'DC2')); %for NY798 its DC2    
else
  Elec_info.index_triggerchan = ...
      find(strcmp(Elec_info.LabelEDF_all,'DC1')); %find trigger channel (usually DC1)  
end

%5.2) Find Cardic/ECG channels (bipolar recording)
Elec_info.index_ECGchan = [];
%cardiac channels are differently named, either C*** or EKG * - check for both options
if sum(contains(Elec_info.LabelEDF_all, 'EKG')) > 0 %check if channels labeled 'EKG' are present
    Elec_info.index_ECGchan = find(contains(Elec_info.LabelEDF_all, 'EKG'));
    Elec_info.index_ECGchan = [Elec_info.index_ECGchan(:)'];
elseif sum(contains(Elec_info.LabelEDF_all, 'ECG')) > 0 %check if channels labeled 'EKG' are present
    Elec_info.index_ECGchan = find(contains(Elec_info.LabelEDF_all, 'ECG'));
    Elec_info.index_ECGchan = [Elec_info.index_ECGchan(:)'];
elseif sum(contains(Elec_info.LabelEDF_all,'C')) > 0 %check if channels labeled with 'C' are present
    chan_C = find(contains(Elec_info.LabelEDF_all,'C')); %find cardiac channel (usually C***)
    for i_channel = 1:length(chan_C)
        pos_C{i_channel} = isstrprop(Elec_info.LabelEDF_all{chan_C(i_channel),1},'alpha');  %find location of all letters in label name
        if pos_C{i_channel}(1) == 1 && sum(pos_C{i_channel}) == 1
            Elec_info.index_ECGchan = [Elec_info.index_ECGchan; chan_C(i_channel)];
        end
    end
    Elec_info.index_ECGchan = [Elec_info.index_ECGchan(:)'];
    clear chan_C pos_C
    if isempty(Elec_info.index_ECGchan)
        disp('!! no ECG channels found !!');
    end
else
    disp('!! no ECG channels found !!');
end

%% 6) Order subfields
subfield_order = {'LabelEDF_all','LabelEDF_match','LabelEDF_NOmatch',...
    'LabelMNI_all','LabelMNI_adjusted_all', 'LabelMNI_adjusted_match',...
    'Type_all','TypeMNI_match',...
    'T1AnatLabel_MNIorder', 'T1AnatLabel_EDForder',...
    'NumMNI_all','Num_EDFMNImatch','Num_EDFMNINOmatch',...
    'NumStrip_all','NumGrid_all','NumDepth_all','NumUndef_all',...
    'NumStrip_match','NumGrid_match','NumDepth_match','NumUndef_match',...
    'EDF2MNI_ElecLabelmatch','EDF2MNI_NOElecLabelmatch',...
    'index4MNI_matchingElec_info','index4EDF_matchingElec_info',...
    'index4EDF_selectStripGridElec_info',...
    'index_triggerchan','index_ECGchan'};
Elec_info = orderfields(Elec_info, subfield_order);

end