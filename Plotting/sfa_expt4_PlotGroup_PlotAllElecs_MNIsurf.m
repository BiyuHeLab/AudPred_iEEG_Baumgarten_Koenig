%TJB: Plot all electrodes of all subjects to show overall coverage
%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/data/gogodisk4/thomas/AuditoryPrediction/iEEG/')
% addpath('//gogo.sb.nyumc.org/data/gogodisk4/thomas/AuditoryPrediction/iEEG/')

sfa_expt4_setVars
paths_sfa_expt4 = sfa_expt4_paths(vars.location);

%Add base dir and own script dir
addpath(genpath(paths_sfa_expt4.BaseDir));
addpath(genpath(paths_sfa_expt4.ScriptsDir));
addpath(genpath(paths_sfa_expt4.Freesurfer)); %path to freesurfer where read_surf function is

% addpath(genpath('/data/gogodisk4/thomas/AuditoryPrediction/iEEG/Scripts/FromRichard'));

%% 0.2) Determine subject-specific parameters (whole-recording)
i_sub = 1; %proxy to load subs_PreProcSettings
sub_list = vars.sub_list; %patients
sub = sub_list{i_sub};

sfa_expt4_subjectinfo %load subject info file (var: si)
subs_PreProcSettings = sfa_expt4_subjectPreProcSettings; %load in file with individual preproc infos

subs = si.sub_list;
saveplot = 1;

path_fig = [paths_sfa_expt4.Fig_ECoGdataraw_GroupAvg '/ElectrodeCoverage/'];
mkdir([path_fig]);


%% 1) Load ECoG preproc data for channel labels and position
    data_AllElec = struct;
    data_AllElec.label = [];
    data_AllElec.sub = [];
    data_AllElec.elec.chanpos = [];
    elecs_persub = [];
for i_sub = 1:length(subs)
    sub = sub_list{i_sub};
    load([paths_sfa_expt4.Preproc_ECoGdata sub '/' sub '_ECoGdata_trials_refLNfiltdetrend.mat']);%path to indiv preprocessed ECoG data

    %append electrodes to one common summary struct
    data_AllElec.label = [data_AllElec.label; data_ECoGfiltref_trials.label(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF)];
    data_AllElec.sub = [data_AllElec.sub; ones(length(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF),1)*i_sub];
    elecs_persub{i_sub} = find(data_AllElec.sub == i_sub);    
    data_AllElec.elec.chanpos = [data_AllElec.elec.chanpos; data_ECoGfiltref_trials.elec.chanpos(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF,:)];
end  

%% 2) Plot all electrodes (nly location, not functional parameter)

        %2Plot electrodes as spheres on MNI brain with each sphere color-coded
        %according to k-prime val     
%         coords = data_ECoGfiltref_trials.elec.chanpos(1:subs_PreProcSettings.(sub).number_ECoGchan,;); %all ECoG chan
        coords = data_AllElec.elec.chanpos; %MNI coordinates for selected electrodes
        vals = ones(length(data_AllElec.elec.chanpos),1); %parameter of interest for resp. electrodes
            for i_sub = 1:length(subs)
                vals(elecs_persub{i_sub}) = vals(elecs_persub{i_sub}).*i_sub;
            end
                    
        chanSize = (vals./vals)*2; %electrode size (arbitrary)
        clims = [0,3]; %colormap limits
%         clims = [1,max(vals)]; %colormap limits
        cmap = 'winter';       
        view_angle = [270,0]; %(270,0) - LH, (90,0) = RH
        c = sfa_expt4_Plot_PlotEleconSurf(coords,vals,chanSize,clims,cmap,view_angle);
        
        Figtitle = {['All Electrodes from ' num2str(length(subs)) ' subjects']};
        title(Figtitle)

        
        if saveplot
            filename     = ['AllElecs_MNIsurf_' num2str(length(subs)) 'subs.png'];
            figfile      = [path_fig filename];                

            saveas(gcf, [figfile], 'png'); %save png version
        end        
