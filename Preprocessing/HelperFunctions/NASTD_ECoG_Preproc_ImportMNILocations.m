function elec_mni_frv_new = NASTD_ECoG_Preproc_ImportMNILocations...
    (sub, MNIfilename, numElectrodes)

%Aim: Import MNI coordinates for all electrodes
elec_mni_frv_new = struct;
if strcmp(sub, 'NY798') %Problems with text file in this sub, thus using a converted excel sheet
    temp_table ...
        = readtable(MNIfilename);
    MNIlabels = temp_table{:,1};
    MNIx = temp_table{:,2};
    MNIy = temp_table{:,3};
    MNIz = temp_table{:,4};
else
    [MNIlabels, MNIx, MNIy, MNIz] ...
        = importElecsfile(MNIfilename,1, numElectrodes);
end

elec_mni_frv_new.chanpos = zeros(numElectrodes,3);
elec_mni_frv_new.chantype = cell(numElectrodes,1);
elec_mni_frv_new.chanunit = cell(numElectrodes,1);
elec_mni_frv_new.coordsys = 'mni';
elec_mni_frv_new.elecpos = zeros(numElectrodes,3);
elec_mni_frv_new.label = cell(numElectrodes,1);
elec_mni_frv_new.tra = zeros(numElectrodes,numElectrodes);
elec_mni_frv_new.tra = diag(ones(numElectrodes,1));
elec_mni_frv_new.unit = 'mm';
elec_mni_frv_new.cfg = struct;

for i_electrode = 1:length(MNIx)
    elec_mni_frv_new.chanpos(i_electrode,1) = MNIx(i_electrode);
    elec_mni_frv_new.chanpos(i_electrode,2) = MNIy(i_electrode);
    elec_mni_frv_new.chanpos(i_electrode,3) = MNIz(i_electrode);
    elec_mni_frv_new.chanunit{i_electrode} = 'V';
    elec_mni_frv_new.chantype{i_electrode} = 'ieeg';
    elec_mni_frv_new.label{i_electrode} = MNIlabels{i_electrode};
end

elec_mni_frv_new.elecpos = elec_mni_frv_new.chanpos;
