function [Elec_info, SelectElec_info]= sfa_expt4_getMNIandEDFElecInfo(sub, MNIelec_table, LabelEDF_all)
%Reads info from both MNI(location) and EDF(recording) files. Compares and
%matches the electrode labels. From matching labels, strip and grid
%electrodes are selected for further processing. Output is a summary file
%with electrode label, type, number info for 1) all electrodes and 2)
%selected electrodes
    subs_PreProcSettings = sfa_expt4_subjectPreProcSettings; %load in file with individual preproc infos

    Elec_info = struct;

%1.1) MNI
%1.1.1 Transfer table infor to respective subfields
    Elec_info.NumMNI_all = size(MNIelec_table,1); %determine number of electrodes based on .txt file
    Elec_info.LabelMNI_all = MNIelec_table{:,1};    
    Elec_info.LabelMNI_adjusted_all = MNIelec_table{:,1};
    Elec_info.Type_all = MNIelec_table{:,5}; %G = grid, S = strip, D = depth
   
%1.1.2 Adjust MNI label names
    %2.1.2.1 MNI labels have an additional 0 before 1 digit numbers. 
    %Delete those 0 to enable EDF-MNI electrode label matching
for i_MNILabel = 1:Elec_info.NumMNI_all
    pos_zero = strfind(Elec_info.LabelMNI_adjusted_all{i_MNILabel,1},'0'); %find location of zero in label name
    pos_letters = isstrprop(Elec_info.LabelMNI_adjusted_all{i_MNILabel,1},'alpha');  %find location of all letters in label name
    pos_lastletter = find(pos_letters,1,'last');  %find location of all letters in label name
   
    if ~isempty(pos_zero) %if there is a 0 in the current label
        if pos_zero == pos_lastletter+1 %and if there is a letter in front of the 0
            Elec_info.LabelMNI_adjusted_all{i_MNILabel,1}(regexp(Elec_info.LabelMNI_adjusted_all{i_MNILabel,1},'0')) = []; %remove 0
        end
    end   
end
clear pos*
    
    %1.1.2.2 Individual label adjustments
    Elec_info = sfa_expt4_IndivMNI2EDFLabelAdjust(sub,Elec_info); 

%1.1.3 Count who many electrodes of each type are present
    Elec_info.NumStrip_all = 0;
    Elec_info.NumGrid_all = 0;
    Elec_info.NumDepth_all = 0;
    Elec_info.NumUndef_all = 0;

for i_elec = 1:length(Elec_info.Type_all)
    if strcmp(Elec_info.Type_all(i_elec),'S')
        index_selectStripGridElec_info(i_elec) = 1;
        Elec_info.NumStrip_all = Elec_info.NumStrip_all+1;
    elseif strcmp(Elec_info.Type_all(i_elec),'G')
        index_selectStripGridElec_info(i_elec) = 1;           
        Elec_info.NumGrid_all = Elec_info.NumGrid_all+1;
   elseif strcmp(Elec_info.Type_all(i_elec),'D')
        index_selectStripGridElec_info(i_elec) = 0; 
        Elec_info.NumDepth_all = Elec_info.NumDepth_all+1;
    else
        index_selectStripGridElec_info(i_elec) = NaN;
        Elec_info.NumUndef_all = Elec_info.NumUndef_all+1;
    end
end

%1.2) EDF
    Elec_info.LabelEDF_all = LabelEDF_all;

%% 2) Match electrode labels between .EDF and .txt files (i.e., determine which channels are present in both)
    %2.1) Find corresponding labels (Original EDF to adjusted MNI)
    [Elec_info.EDF2MNI_ElecLabelmatch, Elec_info.EDF2MNI_NOElecLabelmatch] = ...
        sfa_expt4_matchElecLabels(sub, Elec_info.LabelEDF_all, Elec_info.LabelMNI_adjusted_all);
    Elec_info.Num_EDFMNImatch = sum((~isnan(Elec_info.EDF2MNI_ElecLabelmatch)));
    Elec_info.Num_EDFMNINOmatch = sum((isnan(Elec_info.EDF2MNI_ElecLabelmatch)));

    %Display matching and non-matching electrodes
    for i_channel = 1:Elec_info.Num_EDFMNImatch
        if isnan(Elec_info.EDF2MNI_ElecLabelmatch(i_channel))
            disp(['EDFChannel ' num2str(i_channel) ': ' ...
                Elec_info.LabelEDF_all{i_channel} '. Not Found']);
        else
            disp(['EDFChannel ' num2str(i_channel) ': ' ...
                Elec_info.LabelEDF_all{i_channel} '. MNILabel ' ...
                Elec_info.LabelMNI_all{Elec_info.EDF2MNI_ElecLabelmatch(i_channel)}]);
        end   
    end
    disp('Non-compared EDF channels (non-iEEG channels):');
    disp([Elec_info.LabelEDF_all(end-(subs_PreProcSettings.(sub).number_nonEEGchan-1):end)']);

    %2.2) Determine position/index of matching electrodes for both EDF and MNI files 
    Elec_info.index4EDF_matchingElec_info = find(~isnan(Elec_info.EDF2MNI_ElecLabelmatch)); %read out index of matching channels for EDF file
    Elec_info.index4MNI_matchingElec_info = Elec_info.EDF2MNI_ElecLabelmatch(Elec_info.index4EDF_matchingElec_info);
    %2.3) Determine label in EDF and MNI style of matching electrodes
    Elec_info.LabelEDF_match = Elec_info.LabelEDF_all(Elec_info.index4EDF_matchingElec_info); %gives EDF-labels for matching electrodes
    Elec_info.LabelMNI_adjusted_match = Elec_info.LabelMNI_adjusted_all(sort(Elec_info.EDF2MNI_ElecLabelmatch(~isnan(Elec_info.EDF2MNI_ElecLabelmatch)))); %gives MNI-labels for matching electrodes
    %2.4) Align order for both label versions (align MNI order to EDF order)
    [~, order] = ismember(Elec_info.LabelEDF_match, Elec_info.LabelMNI_adjusted_match);
    % [Elec_info.LabelEDF_match, Elec_info.LabelMNI_adjusted_match(order,1)] %displays parallel ordered varrays
    if ~isequal(Elec_info.LabelEDF_match, Elec_info.LabelMNI_adjusted_match(order,1)) %checks is ordered arrays are the same
        disp('EDF - MNI label mismatch !!');
    end
    Elec_info.LabelMNI_adjusted_match_EDFordered = Elec_info.LabelMNI_adjusted_match(order,1); %reorder MNI labels in the same order as EDF labels
    %2.5) Align order for elec typ as well (align MNI type info to EDF order)
    Elec_info.TypeMNI_match = Elec_info.Type_all(Elec_info.index4MNI_matchingElec_info); %gives MNI-labels for matching electrodes
    Elec_info.TypeMNI_match_EDFordered = Elec_info.TypeMNI_match(order,1);
    %2.6) separately save labels of non-matching electrodes
    Elec_info.LabelEDF_NOmatch = Elec_info.EDF2MNI_NOElecLabelmatch;
    
%% 3) Select strip & grid electrodes (and get rid of depth electrodes)
    %3.1) Determine index for electrode subtypes and count respective subtypes
    Elec_info.NumStrip_match = 0;
    Elec_info.NumGrid_match = 0;
    Elec_info.NumDepth_match = 0;
    Elec_info.NumUndef_match = 0;

    for i_elec = 1:length(Elec_info.TypeMNI_match_EDFordered)
        if strcmp(Elec_info.TypeMNI_match_EDFordered(i_elec),'S') %Strip
            Elec_info.index4EDF_selectStripGridElec_info(i_elec) = 1;
            Elec_info.NumStrip_match = Elec_info.NumStrip_match+1;
        elseif strcmp(Elec_info.TypeMNI_match_EDFordered(i_elec),'G') %Grip
            Elec_info.index4EDF_selectStripGridElec_info(i_elec) = 1;           
            Elec_info.NumGrid_match = Elec_info.NumGrid_match+1;
       elseif strcmp(Elec_info.TypeMNI_match_EDFordered(i_elec),'D') %Depth
            Elec_info.index4EDF_selectStripGridElec_info(i_elec) = 0; 
            Elec_info.NumDepth_match = Elec_info.NumDepth_match+1;
        else
            Elec_info.index4EDF_selectStripGridElec_info(i_elec) = NaN; %Not determinable
            Elec_info.NumUndef_match = Elec_info.NumUndef_match+1;
        end
    end
    Elec_info.index4EDF_selectStripGridElec_info = logical(Elec_info.index4EDF_selectStripGridElec_info);

%% 4) Copy selected electrodes (matching strip+grid) info into separate subfield
SelectElec_info = [];
SelectElec_info.Label = Elec_info.LabelEDF_match(Elec_info.index4EDF_selectStripGridElec_info);
SelectElec_info.Type = Elec_info.TypeMNI_match_EDFordered(Elec_info.index4EDF_selectStripGridElec_info);
SelectElec_info.index4EDF = Elec_info.index4EDF_matchingElec_info(Elec_info.index4EDF_selectStripGridElec_info);
SelectElec_info.index4MNI = Elec_info.index4MNI_matchingElec_info(Elec_info.index4EDF_selectStripGridElec_info);

%% 5) Determine trigger and ECG channels
    %5.1) Find Trigger channel
    Elec_info.index_triggerchan = find(strcmp(Elec_info.LabelEDF_all,'DC1')); %find trigger channel (usually DC1)
    
    %5.2) Find Cardic/ECG channels (bipolar recording)    
    Elec_info.index_ECGchan = [];
    %cardiac channels are differently named, either C*** or EKG * - check for both options
    if sum(contains(Elec_info.LabelEDF_all, 'EKG')) > 0 %check if channels labeled 'EKG' are present
        Elec_info.index_ECGchan = find(contains(Elec_info.LabelEDF_all, 'EKG'));
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
    'LabelMNI_all','LabelMNI_adjusted_all',...
    'LabelMNI_adjusted_match','LabelMNI_adjusted_match_EDFordered',...
    'Type_all','TypeMNI_match','TypeMNI_match_EDFordered',...
    'NumMNI_all','Num_EDFMNImatch','Num_EDFMNINOmatch',...
    'NumStrip_all','NumGrid_all','NumDepth_all','NumUndef_all',...
    'NumStrip_match','NumGrid_match','NumDepth_match','NumUndef_match',...
    'EDF2MNI_ElecLabelmatch','EDF2MNI_NOElecLabelmatch',...
    'index4MNI_matchingElec_info','index4EDF_matchingElec_info',...
    'index4EDF_selectStripGridElec_info',...
    'index_triggerchan','index_ECGchan'};
    Elec_info = orderfields(Elec_info,subfield_order);
    
end