function elec = NASTD_ECoG_Preproc_CreateElecStruct...
    (elec_mni_frv_new, ...
    EDFlabel, ...
    compareLabels, T1AnatLabel, ...
    EKGchannels, DCchannels)

% Aim: Create FT-like electrode structure containing all relevant info

numEDFelectrodes = length(compareLabels);

elec = struct;

elec.chanpos = zeros(numEDFelectrodes,3);
elec.chantype = cell(numEDFelectrodes,1);
elec.chanunit = cell(numEDFelectrodes,1);
elec.T1AnatLabel = cell(numEDFelectrodes,1);
elec.coordsys = elec_mni_frv_new.coordsys;
elec.elecpos = zeros(numEDFelectrodes,3);
elec.label = cell(numEDFelectrodes,1);
elec.tra = zeros(numEDFelectrodes,numEDFelectrodes);
elec.tra = diag(ones(numEDFelectrodes,1));
elec.unit = 'mm';
elec.cfg = struct;

for i_electrode = 1:numEDFelectrodes
    if isnan(compareLabels(i_electrode))
        elec.chantype{i_electrode} = 'missingLocation';
        elec.chanunit{i_electrode} = 'uV';
        elec.T1AnatLabel{i_electrode} = 'missingLocation';
        
    else
        elec.chanpos(i_electrode,:) = elec_mni_frv_new.chanpos(compareLabels(i_electrode),:);
        elec.chanunit(i_electrode)  = elec_mni_frv_new.chanunit(compareLabels(i_electrode));
        elec.chantype(i_electrode)  = elec_mni_frv_new.chantype(compareLabels(i_electrode));
        elec.T1AnatLabel{i_electrode} = T1AnatLabel{compareLabels(i_electrode)};
        
    end
    elec.label{i_electrode} = EDFlabel{i_electrode};
end

for i_EKGchannels = EKGchannels
    elec.chanunit{i_EKGchannels} = 'uV';
    elec.chantype{i_EKGchannels} = 'ecg';
end
    
for i_DCchannels = DCchannels
    elec.chanunit{i_DCchannels} = 'uV';
    elec.chantype{i_DCchannels} = 'trigger';
end

