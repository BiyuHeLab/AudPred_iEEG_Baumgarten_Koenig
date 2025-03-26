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
% subs_PreProcSettings = NASTD_SubjectPreProcSettings; %load in file with individual preproc infos
% 
% subs = si.sub_list;
% saveplot = 1;
% 
% path_fig = [paths_NASTD.Fig_ECoGdataraw_GroupAvg '/ElectrodeCoverage/'];
% mkdir([path_fig]);


%% 1) Load ECoG preproc data for channel labels and position
    data_AllElec = struct;
    data_AllElec.label = [];
    data_AllElec.sub = [];
    data_AllElec.elec.chanpos = [];
    elecs_persub = [];
for i_sub = vars.validSubjs
    
    sub = sub_list{i_sub};
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];%path to indiv preprocessed ECoG data
    
    % 0.3) Load in preproc data
    tic
    disp(['Loading preprocessed data set for sub: ' sub])
    load(loadfile_ECoGpreprocdata);
    disp(['done loading in ' num2str(toc) ' sec'])
    
    %append electrodes to one common summary struct
    data_AllElec.label = [data_AllElec.label; DataClean_AllTrials.label(DataClean_AllTrials.cfg.info_elec.selected.index4EDF)];
    data_AllElec.sub = [data_AllElec.sub; ones(length(DataClean_AllTrials.cfg.info_elec.selected.index4EDF),1)*i_sub];
    elecs_persub{i_sub} = find(data_AllElec.sub == i_sub);    
    data_AllElec.elec.chanpos = [data_AllElec.elec.chanpos; DataClean_AllTrials.elec.chanpos(DataClean_AllTrials.cfg.info_elec.selected.index4EDF,:)];
end  

%% 2) Plot all electrodes (nly location, not functional parameter)

        %2Plot electrodes as spheres on MNI brain with each sphere color-coded
        %according to k-prime val     
%         coords = data_ECoGfiltref_trials.elec.chanpos(1:subs_PreProcSettings.(sub).number_ECoGchan,;); %all ECoG chan
        coords = data_AllElec.elec.chanpos; %MNI coordinates for selected electrodes
        vals = ones(length(data_AllElec.elec.chanpos),1); %parameter of interest for resp. electrodes
            for i_sub = vars.validSubjs
                vals(elecs_persub{i_sub}) = vals(elecs_persub{i_sub}).*(i_sub*10);
            end
                    
        chanSize = (vals./vals)*1.5; %electrode size (arbitrary)
        clims = [0,90]; %colormap limits
%         clims = [1,max(vals)]; %colormap limits
        cmap = 'jet';       
        view_angle = [270,0]; %(270,0) - LH, (90,0) = RH
        
        %Option 1: %both hemispheres
        c = NASTD_ECoG_Plot_PlotElecsSurf_SubplotHemis(coords,vals,chanSize,clims,cmap); 
                
        %Option 2: Project all electrodes on left hemisphere,
        coords(:,1) = abs(coords(:,1)) * -1;
        c = NASTD_ECoG_Plot_PlotElecsSurf(coords,vals,chanSize,clims,cmap, view_angle, [], []); 
           
        
        Figtitle = {['All Electrodes from ' num2str(length(subs)) ' subjects']};
        title(Figtitle)

        
        if saveplot
            filename     = ['AllElecs_MNIsurf_' num2str(length(vars.validSubjs)) 'subs.png'];
            figfile      = [path_fig filename];                

            saveas(gcf, [figfile], 'png'); %save png version
        end        
