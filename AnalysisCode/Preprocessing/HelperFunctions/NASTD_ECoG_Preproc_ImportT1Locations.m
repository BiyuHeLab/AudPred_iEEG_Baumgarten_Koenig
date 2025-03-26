function elec_T1_frv_new = NASTD_ECoG_Preproc_ImportT1Locations...
    (T1filename,numElectrodes)

%Aim: Import MNI coordinates for all electrodes

elec_T1_frv_new = struct;
[T1labels,T1x,T1y,T1z] = importElecsfile(T1filename,1, numElectrodes);

elec_T1_frv_new.chanpos = zeros(numElectrodes,3);
elec_T1_frv_new.chantype = cell(numElectrodes,1);
elec_T1_frv_new.chanunit = cell(numElectrodes,1);
elec_T1_frv_new.coordsys = 'T1';
elec_T1_frv_new.elecpos = zeros(numElectrodes,3);
elec_T1_frv_new.label = cell(numElectrodes,1);
elec_T1_frv_new.tra = zeros(numElectrodes,numElectrodes);
elec_T1_frv_new.tra = diag(ones(numElectrodes,1));
elec_T1_frv_new.unit = 'mm';
elec_T1_frv_new.cfg = struct;


for i_electrode = 1:length(T1x)
    elec_T1_frv_new.chanpos(i_electrode,1) = T1x(i_electrode);
    elec_T1_frv_new.chanpos(i_electrode,2) = T1y(i_electrode);
    elec_T1_frv_new.chanpos(i_electrode,3) = T1z(i_electrode);
    elec_T1_frv_new.chanunit{i_electrode} = 'V';
    elec_T1_frv_new.chantype{i_electrode} = 'ieeg';
    elec_T1_frv_new.label{i_electrode} = T1labels{i_electrode};
end
elec_T1_frv_new.elecpos = elec_T1_frv_new.chanpos;
