%TJB: Plot all electrodes of all subjects to show overall coverage
%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')

NASTD_ECoG_setVars
paths_NASTD_ECoG = NASTD_ECoG_paths;

% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer));
%Add project base and script dir

%Determine subjects
sub_list = vars.sub_list;

%Load in file with individual preproc infos
subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;

ToneDur_text = {'0.2' '0.4'};

plot_poststepFigs = 1;
save_poststepFigs = 0;

%% 0.2) Determine subject-specific parameters (whole-recording)
% i_sub = 1; %proxy to load subs_PreProcSettings
% sub_list = vars.sub_list; %patients
% sub = sub_list{i_sub};
% 
% NASTD_subjectinfo %load subject info file (var: si)
% 
% subs_PreProcSettings = NASTD_SubjectPreProcSettings; %load in file with individual preproc infos

%subs = si.sub_list;
saveplot = 1;
% 
% path_fig = [paths_NASTD.Fig_ECoGdataraw_GroupAvg '/ElectrodeCoverage/'];
% mkdir([path_fig]);

%Fill out SOZ electrodes
SOZ_Elecs_persub = struct();
SOZ_Elecs_persub.NY688 = {'G35','G36','G43','G44'};
SOZ_Elecs_persub.NY704 = {'G41','G42','AT1','AT2','AT3','MT1', 'MT2', 'MT3', 'DPMT1', 'DPMT2', 'DPMT3', 'DPMT4', 'DPMT5'};
SOZ_Elecs_persub.NY708 = {'G11','G12','G13','G14'};
SOZ_Elecs_persub.NY723 = {'RDMT2','RDMT3','LDMT3','LDMT4'};
SOZ_Elecs_persub.NY742 = {};
SOZ_Elecs_persub.NY751 = {'G59','G60','G61','G62', 'G63', 'G64', 'PT3', 'PT4'};
SOZ_Elecs_persub.NY787 = {'G3','G4','G25','OF3','OF4'};
SOZ_Elecs_persub.NY794 = {'G45','G46','LAT1','LAT2','LAT3','LAT4','LPT1','LPT2','LPT3','LPT4','LPF1','LOF3','LOF4'};
SOZ_Elecs_persub.NY798 = {};

RejElecs = NASTD_ECoG_Preproc_SubPreprocSettings; %field rejecedtedChan_label

%% 1) Load ECoG preproc data for channel labels and position
for i_sub = vars.validSubjs
    sub = sub_list{i_sub};
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];%path to indiv preprocessed ECoG data
    
    % 0.3) Load in preproc data
    tic
    disp(['Loading preprocessed data set for sub: ' sub])
    load(loadfile_ECoGpreprocdata);
    disp(['done loading in ' num2str(toc) ' sec'])
    
    % Get electrode labels and positions for current subject
    labels = DataClean_AllTrials.label(DataClean_AllTrials.cfg.info_elec.selected.index4EDF);
    coords = DataClean_AllTrials.elec.chanpos(DataClean_AllTrials.cfg.info_elec.selected.index4EDF,:);
    nElec  = length(labels);
    
    %% 2) Plot all electrodes (nly location, not functional parameter)
    
    %2Plot electrodes as spheres on MNI brain with each sphere color-coded
    %according to k-prime val
    %         coords = data_ECoGfiltref_trials.elec.chanpos(1:subs_PreProcSettings.(sub).number_ECoGchan,;); %all ECoG chan
    nElec  = size(coords,1);
    
    % --- Build a "vals" vector: 1 = non‑SOZ, 2 = SOZ, 3 = Rejected ---
    vals = ones(nElec,1);  % default = non‑SOZ
    chanSize = ones(nElec,1) * 1.5;
    if isfield(SOZ_Elecs_persub, sub) %Mark SOZ Electrodes
        sozLab = SOZ_Elecs_persub.(sub);
        if ~isempty(sozLab)
            isSOZ_here = ismember(labels, sozLab);
            vals(isSOZ_here) = 2;  % mark SOZ as 2
        end
    end
    if isfield(RejElecs, sub)
        RejLab = RejElecs.(sub).rejectedChan_label;
        if ~isempty(RejLab)
            isRej_here = ismember(labels, RejLab);
            chanSize(isRej_here) = 3;  % mark Rejected electrode as 3
        end
    end
    
    clims    = [1 2];                                          
    cmap     = [0 0 0; 1 0 0];       
    view_angle = [270,0]; %(270,0) - LH, (90,0) = RH
    
    %Option 1: %both hemispheres
    c = NASTD_ECoG_Plot_PlotElecsSurf_SubplotHemis(coords,vals,chanSize,clims,cmap);
    
    %Option 2: Project all electrodes on left hemisphere,
    % coords(:,1) = abs(coords(:,1)) * -1;
    % c = NASTD_ECoG_Plot_PlotElecsSurf(coords,vals,chanSize,clims,cmap, view_angle, [], []);
    set(gcf, 'Position', [100 100 1000 500]);  % [left bottom width height]
    colorbar('off');  % remove colorbar
    
    Figtitle = {['All Electrodes from subject #' num2str(i_sub)]};
    t = sgtitle(Figtitle);
    t.FontSize = 16;
    
    if saveplot
        filename     = ['AllElecs_MNIsurf_' sub '_withSOZandRejected.png'];
        figfile      = ['/isilon/LFMI/VMdrive/Lua/NASTD/Figs/' filename];
        
        saveas(gcf, [figfile], 'png'); %save png version
    end
    close all
end