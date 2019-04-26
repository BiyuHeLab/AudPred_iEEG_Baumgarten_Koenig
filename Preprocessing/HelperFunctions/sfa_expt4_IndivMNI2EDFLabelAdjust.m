function Elec_info = sfa_expt4_IndivMNI2EDFLabelAdjust(sub,Elec_info)
%Changes MNI labels which differ from EDF labels to the corresponding EDF
%form. Info which channels and which new labels are selected is taken from
%sfa_expt4_subjectPreProcSettings file.

subs_PreProcSettings = sfa_expt4_subjectPreProcSettings; %load in file with individual preproc infos

index_mismatchchannels = subs_PreProcSettings.(sub).MNI2EDFmismatches.index_mismatchchannels; %read out index of label-mismatching channels
for i_MNILabel = 1:length(index_mismatchchannels)
    Elec_info.LabelMNI_adjusted_all{index_mismatchchannels(i_MNILabel),1} = ...
        subs_PreProcSettings.(sub).MNI2EDFmismatches.correctedlabel_mismatchchannels{i_MNILabel}; %replace old label with new label
disp(['Original MNI label ' subs_PreProcSettings.(sub).MNI2EDFmismatches.label_mismatchchannels{i_MNILabel} ... 
    ' is changed to EDF-corresponding label ' subs_PreProcSettings.(sub).MNI2EDFmismatches.correctedlabel_mismatchchannels{i_MNILabel} ]);

end
