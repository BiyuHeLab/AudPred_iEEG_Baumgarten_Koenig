function [EDF2MNI_ElecLabelmatch, EDF2MNI_NOElecLabelmatch] = ...
    NASTD_ECoG_Preproc_MatchElecLabels ...
    (sub, EDFLabels, MNILabels)

%Finds corresponding labels for EDF(recording) and MNI(location) data.
%Important: Uses adjusted MNI labels (i.e., 0 before single-digits removed
%and individual labels adjusted). 

%Output lists corresponding indices of MNI labels for each EDF label
subs_PreProcSettings = ...
    NASTD_ECoG_Preproc_SubPreprocSettings; %load in file with individual preproc infos
number_nonEEGchan = ...
    subs_PreProcSettings.(sub).number_nonEEGchan; %read out number of non-iEEG chan (e.g., DC, Pleth)

newEDFLabels = {};
for i_chan = 1:(length(EDFLabels) - number_nonEEGchan) %skip non-iEEG chan
    newEDFLabels{i_chan} = EDFLabels{i_chan};
    %newEDFDataLabels{i_chan} = strrep(newEDFDataLabels{i_chan},'_',''); %Change string
end

EDF2MNI_ElecLabelmatch = nan(length(newEDFLabels),1);  %LH, RH
   
%Match indices & labels        
for i_MNILabel = 1:size(MNILabels,1) %for each MNI label
    for i_newEDFLabel = 1:length(newEDFLabels) %check for each EDF Label if it matches
        if strcmp(MNILabels{i_MNILabel,1}, newEDFLabels{i_newEDFLabel})%exact match
            EDF2MNI_ElecLabelmatch(i_newEDFLabel,1) = i_MNILabel;
            %if it matches, place the MNI index at the index-position of the EDF elec
        end
    end
end

%Also read out those labels (EDF) where no match has been found
EDF2MNI_NOElecLabelmatch = newEDFLabels(isnan(EDF2MNI_ElecLabelmatch))';

end
