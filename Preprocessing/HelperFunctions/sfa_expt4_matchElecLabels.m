function [EDF2MNI_ElecLabelmatch, EDF2MNI_NOElecLabelmatch]= sfa_expt4_matchElecLabels(sub, EDFLabels, MNILabels)
%Finds corresponding labels for EDF(recording) and MNI(location) data.
%Important: Uses adjusted MNI labels (i.e., 0 before single-digits removed
%and individual labels adjusted). %Output lists corresponding indices of
%MNI labels for each EDF label
subs_PreProcSettings = sfa_expt4_subjectPreProcSettings; %load in file with individual preproc infos
number_nonEEGchan = subs_PreProcSettings.(sub).number_nonEEGchan; %read out number of non EEG chan (e.g., DC, Pleth)


newEDFLabels = {};
for i_chan = 1:(length(EDFLabels)-number_nonEEGchan) %skip non-EEG chan
    newEDFLabels{i_chan} = EDFLabels{i_chan};
    %newEDFDataLabels{i_chan} = strrep(newEDFDataLabels{i_chan},'_',''); %Change string
end

EDF2MNI_ElecLabelmatch = nan(length(newEDFLabels),1);  %LH, RH

% %MNI labels have an additional 0 before 1 digit numbers. Delete those 0
% for i_MNILabel = 1:size(MNILabels,1)
%     pos_zero = strfind(MNILabels{i_MNILabel,1},'0'); %find location of zero in label name
%     pos_letters = isstrprop(MNILabels{i_MNILabel,1},'alpha');  %find location of all letters in label name
%     pos_lastletter = find(pos_letters,1,'last');  %find location of all letters in label name
%    
%     if ~isempty(pos_zero) %if there is a 0 in the current label
%         if pos_zero == pos_lastletter+1 %and if there is a letter in front of the 0
%             MNILabels{i_MNILabel,1}(regexp(MNILabels{i_MNILabel,1},'0')) = []; %remove 0
%         end
%     end
% end
   
%Match labels        
for i_MNILabel = 1:size(MNILabels,1) %for each MNI label
    for i_newEDFLabel = 1:length(newEDFLabels) %check for each EDF Label if it macthes
        if strcmp(MNILabels{i_MNILabel,1},newEDFLabels{i_newEDFLabel}) == 1 %exact match
            EDF2MNI_ElecLabelmatch(i_newEDFLabel,1) = i_MNILabel;
            %output list corresponding indices of MNI labels for each EDF label            
        end
    end
end

%Also read out those labels (EDF) where no match has been found
EDF2MNI_NOElecLabelmatch = newEDFLabels(isnan(EDF2MNI_ElecLabelmatch))';

