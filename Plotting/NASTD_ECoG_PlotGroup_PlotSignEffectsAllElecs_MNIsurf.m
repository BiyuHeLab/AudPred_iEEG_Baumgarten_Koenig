%TJB: Plot allsign effects (HisTRack, Prediction) for all electrodes of all subjects to show overall coverage
%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
% addpath('/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
% addpath('//gogo.sb.nyumc.org/data/gogodisk4/thomas/AuditoryPrediction/iEEG/')
NASTD_setVars
paths_NASTD = NASTD_paths(vars.location);

%Add base dir and own script dir
addpath(genpath(paths_NASTD.BaseDir));
addpath(genpath(paths_NASTD.ScriptsDir));
addpath(genpath(paths_NASTD.Freesurfer)); %path to freesurfer where read_surf function is

%% 0.2) Determine subject-specific parameters (whole-recording)
i_sub = 1; %proxy to load subs_PreProcSettings
sub_list = vars.sub_list; %patients
sub = sub_list{i_sub};

NASTD_subjectinfo %load subject info file (var: si)
subs_PreProcSettings = NASTD_SubjectPreProcSettings; %load in file with individual preproc infos

subs = si.sub_list;
saveplot = 1;

path_fig = [paths_NASTD.Fig_ECoGdataraw_GroupAvg '/ElectrodeCoverage/'];
mkdir([path_fig]);

tonedur_text = {'0.2' '0.4'};

nReps = 100;
pval = 0.05;
ShuffleSets = 'averaged'; %'appended'

baseline_correct = [0];
BC_text = 'NBC';

effects = {'HisTrack_ERP' 'HisTrack_Gamma' 'Prediction_ERP' 'Prediction_Gamma' 'PredictionError_ERP' 'PredictionError_Gamma', 'SimplePredictionError_ERP'};
% effect2plot = 'HisTrack_ERP';
% effect2plot = 'HisTrack_Gamma';
% effect2plot = 'Prediction_ERP';
% effect2plot = 'Prediciton_Gamma';
% effect2plot = 'PredictionError_ERP';
% effect2plot = 'PredicitonError_Gamma';
% effect2plot = 'SimplePredictionError_ERP';

%% 1) Load single subject data and construct summary struct

%1.3 Load in Effect data
for i_effect = 1:length(effects)    
%     effect2plot = effects{i_effect}  
    
    %1.1) Create empyt proxy
    data_AllElec = struct;
    data_AllElec.label = [];
    data_AllElec.sub = [];
    data_AllElec.elec.chanpos = [];
    elecs_persub = cell(length(subs),1);
    data_AllElec.effect = cell(1,2);
        data_AllElec.effect{1} = cell(1,4);
        data_AllElec.effect{2} = cell(1,8);
    data_AllElec.sign_elecs = cell(1,2);
        data_AllElec.sign_elecs{1} = cell(1,4);
        data_AllElec.sign_elecs{2} = cell(1,8);
  
for i_sub = 1:length(subs) 
    
%1.2) Load in preproc data for electrode labels and positions  
    
    sub = sub_list{i_sub};
    load([paths_NASTD.Preproc_ECoGdata sub '/' sub '_data_preproc.mat']);%path to indiv preprocessed ECoG data
    %append electrodes to one common summary struct
    data_AllElec.label = [data_AllElec.label; data_ECoGfiltref_trials.label(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF)];
    data_AllElec.sub = [data_AllElec.sub; ones(length(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF),1)*i_sub];
    elecs_persub{i_sub} = find(data_AllElec.sub == i_sub);    
    data_AllElec.elec.chanpos = [data_AllElec.elec.chanpos; data_ECoGfiltref_trials.elec.chanpos(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF,:)];
          
    for i_tonedur = 1:length(tonedur_text) %tone duration condition
        
            % 1.2.1) Define analysis parameters
            fsample = data_ECoGfiltref_trials.fsample;
            toneDur_inSecs  = str2num(tonedur_text{i_tonedur});
            nSamplesPerTone = toneDur_inSecs * fsample;

            win_size    = 25;  %25 data points with fs 512 = ~49ms as mentioned in manuscript
            s_per_win = (1/fsample)*win_size;
            win_overlap = 0; %as mentioned in manuscript

            windows = [1 win_size];
            while windows(end,end) < nSamplesPerTone
                windows = [windows; windows(end,:) + (win_size - win_overlap)];
            end

            if windows(end,end) > nSamplesPerTone
                windows(end,:) = [];
            end

            
            %1.2.2) Load Subject effect data   
            NASTD_subjectinfo %load subject info file (var: si)
            subs_PreProcSettings = NASTD_SubjectPreProcSettings; %load in file with individual preproc infos

            %Tone duration condition for data load-in
            if strcmp(tonedur_text{i_tonedur},'0.2')
                tonedur_title = '200';
            elseif  strcmp(tonedur_text{i_tonedur},'0.4')
                tonedur_title = '400';
            end

            %History Tracking - ERP range%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(effect2plot,'HisTrack_ERP')  
                
                path_load = [paths_NASTD.ECoGdata_HisTrack sub '/not_baseline_corrected/EXPvsSHUFFLE/'];
                path_fig = [paths_NASTD.Fig_HisTrack_SingleSubs 'all/not_baseline_corrected/EXPvsSHUFFLE/ERP/'];
                mkdir([path_fig]);
                
                Scale_Limits = [1,16]; %Plotting Scale
                Fig_title = ['sign k (p = ' num2str(pval) ') with minimum test error for linear regression on tone history (low Freqs)'];

                if strcmp(ShuffleSets, 'appended')
                    load([path_load sub '_HisTrackSLIDE16_regstatsSHUFFLED' num2str(nReps) 'Reps_' BC_text '_' tonedur_title 'msToneDur_kModelComp_appSets.mat']); %appended shuffle-sets
                elseif strcmp(ShuffleSets, 'averaged')
                    load([path_load sub '_HisTrackSLIDE16_regstatsSHUFFLED' num2str(nReps) 'Reps_' BC_text '_' tonedur_title 'msToneDur_kModelComp_avgSets.mat']); %averaged shuffle-sets
                end    %Parameter of interest: sig2_test_k_pval{i_win}(i_sensor) 

                for i_win = 1:size(windows,1)

                    data_AllElec.effect{i_tonedur}{i_win} = [data_AllElec.effect{i_tonedur}{i_win}; sig2_test_k{i_win}(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF)'];
                    data_AllElec.sign_elecs{i_tonedur}{i_win} = [data_AllElec.sign_elecs{i_tonedur}{i_win}; (sig2_test_k_pval{i_win}(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF') < pval)'];

                end
                
            %History Tracking - Gamma range%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif strcmp(effect2plot,'HisTrack_Gamma')  
                
                path_load = [paths_NASTD.ECoGdata_HisTrack sub '/not_baseline_corrected/EXPvsSHUFFLE/'];
                path_fig = [paths_NASTD.Fig_HisTrack_SingleSubs 'all/not_baseline_corrected/EXPvsSHUFFLE/Gamma/'];
                mkdir([path_fig]);
                
                Scale_Limits = [1,16]; %Plotting Scale
                Fig_title = ['sign k (p = ' num2str(pval) ') with minimum test error for linear regression on tone history (Gamma Freqs)'];

                if strcmp(ShuffleSets, 'appended')
                    load([path_load sub '_HisTrackGAMMAEnvSLIDE16_regstatsSHUFFLED' num2str(nReps) 'Reps_' BC_text '_' tonedur_title 'msToneDur_kModelComp_appSets.mat']); %appended shuffle-sets
                elseif strcmp(ShuffleSets, 'averaged')
                    load([path_load sub '_HisTrackGAMMAEnvSLIDE16_regstatsSHUFFLED' num2str(nReps) 'Reps_' BC_text '_' tonedur_title 'msToneDur_kModelComp_avgSets.mat']); %averaged shuffle-sets
                end    %Parameter of interest: sig2_test_k_pval{i_win}(i_sensor) 
             
                for i_win = 1:size(windows,1)

                    data_AllElec.effect{i_tonedur}{i_win} = [data_AllElec.effect{i_tonedur}{i_win}; sig2_test_k{i_win}(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF)'];
                    data_AllElec.sign_elecs{i_tonedur}{i_win} = [data_AllElec.sign_elecs{i_tonedur}{i_win}; (sig2_test_k_pval{i_win}(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF') < pval)'];

                end

            %Prediction Effect - ERP range%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif strcmp(effect2plot,'Prediction_ERP')
                
                path_load = [paths_NASTD.ECoGdata_Prediction sub '/not_baseline_corrected/Exp/'];
                path_fig = [paths_NASTD.Fig_Prediction_SingleSubs 'all/not_baseline_corrected/ERP/'];
                mkdir([path_fig]);
                
                Scale_Limits = [-3, 3]; %Plotting Scale
                Fig_title = ['Prediction - Linear regression (Neural activity during p33 - p*34) t-stat (Low Freqs)'];
               
                load([path_load sub '_PredictionERF_regstatsEXP_' BC_text '_' tonedur_title 'msToneDur.mat']);

                for i_win = 1:size(windows,1)

                    data_AllElec.effect{i_tonedur}{i_win} = [data_AllElec.effect{i_tonedur}{i_win}; pred_t1{1}{i_win}'];
                    data_AllElec.sign_elecs{i_tonedur}{i_win} = [data_AllElec.sign_elecs{i_tonedur}{i_win}; (pred_t1_p{1}{i_win} < pval)'];

                end                

            %Prediction Effect - Gamma range%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif strcmp(effect2plot,'Prediction_Gamma')
                
                path_load = [paths_NASTD.ECoGdata_Prediction sub '/not_baseline_corrected/Exp/'];
                path_fig = [paths_NASTD.Fig_Prediction_SingleSubs 'all/not_baseline_corrected/Gamma/'];
                mkdir([path_fig]);
                
                Scale_Limits = [-3, 3]; %Plotting Scale
                Fig_title = ['Prediction - Linear regression (Neural activity during p33 - p*34) t-stat (Gamma Freqs)'];
                
                load([path_load sub '_PredictionGamma_regstatsEXP_' BC_text '_' tonedur_title 'msToneDur.mat']);

                for i_win = 1:size(windows,1)

                    data_AllElec.effect{i_tonedur}{i_win} = [data_AllElec.effect{i_tonedur}{i_win}; pred_t1{1}{i_win}'];
                    data_AllElec.sign_elecs{i_tonedur}{i_win} = [data_AllElec.sign_elecs{i_tonedur}{i_win}; (pred_t1_p{1}{i_win} < pval)'];

                end     
                
            %Prediction Error Effect - ERP range%%%%%%%%%%%%%%%%%%%%%%%%%%%                
            elseif strcmp(effect2plot,'PredictionError_ERP')
                
                path_load = [paths_NASTD.ECoGdata_Prediction sub '/not_baseline_corrected/Exp/'];
                path_fig = [paths_NASTD.Fig_Prediction_SingleSubs 'all/not_baseline_corrected/ERP/'];
                mkdir([path_fig]);
                
                Scale_Limits = [-3, 3]; %Plotting Scale
                Fig_title = ['Prediction Error - Linear regression (Neural activity during p34 - difference (p*34-p34) - p34) t-stat (Low Freqs)'];
               
                load([path_load sub '_PredictionERF_regstatsEXP_' BC_text '_' tonedur_title 'msToneDur.mat']);

                for i_win = 1:size(windows,1)

                    data_AllElec.effect{i_tonedur}{i_win} = [data_AllElec.effect{i_tonedur}{i_win}; pred_err_t1{1}{i_win}'];
                    data_AllElec.sign_elecs{i_tonedur}{i_win} = [data_AllElec.sign_elecs{i_tonedur}{i_win}; (pred_err_t1_p{1}{i_win} < pval)'];

                end                

            %Prediction Error Effect - Gamma range%%%%%%%%%%%%%%%%%%%%%%%%%                
            elseif strcmp(effect2plot,'PredictionError_Gamma')
                
                path_load = [paths_NASTD.ECoGdata_Prediction sub '/not_baseline_corrected/Exp/'];
                path_fig = [paths_NASTD.Fig_Prediction_SingleSubs 'all/not_baseline_corrected/Gamma/'];
                mkdir([path_fig]);
                Scale_Limits = [-3, 3]; %Plotting Scale
                Fig_title = ['Prediction Error - Linear regression (Neural activity during p34 - difference (p*34-p34) - p34) t-stat (Gamma Freqs)'];
                
                load([path_load sub '_PredictionGamma_regstatsEXP_' BC_text '_' tonedur_title 'msToneDur.mat']);

                for i_win = 1:size(windows,1)

                    data_AllElec.effect{i_tonedur}{i_win} = [data_AllElec.effect{i_tonedur}{i_win}; pred_err_t1{1}{i_win}'];
                    data_AllElec.sign_elecs{i_tonedur}{i_win} = [data_AllElec.sign_elecs{i_tonedur}{i_win}; (pred_err_t1_p{1}{i_win} < pval)'];

                end   
                
            %SimplePrediction Error Effect - ERP range%%%%%%%%%%%%%%%%%%%%%%%%%%%                
            elseif strcmp(effect2plot,'SimplePredictionError_ERP')
                
                path_load = [paths_NASTD.ECoGdata_Prediction sub '/not_baseline_corrected/Exp/'];
                path_fig = [paths_NASTD.Fig_Prediction_SingleSubs 'all/not_baseline_corrected/ERP/'];
                mkdir([path_fig]);
                
                Scale_Limits = [-3, 3]; %Plotting Scale
                Fig_title = ['Simple Prediction Error - Linear regression (Neural activity during p34 - difference (p*34-p33)) t-stat (Low Freqs)'];
               
                load([path_load sub '_NEWPredictionERF_regstatsEXP_' BC_text '_' tonedur_title 'msToneDur.mat']);

                for i_win = 1:size(windows,1)

                    data_AllElec.effect{i_tonedur}{i_win} = [data_AllElec.effect{i_tonedur}{i_win}; simplepred_err_t1{1}{i_win}'];
                    data_AllElec.sign_elecs{i_tonedur}{i_win} = [data_AllElec.sign_elecs{i_tonedur}{i_win}; (simplepred_err_t1_p{1}{i_win} < pval)'];

                end 
                
    end  
    end
end

%% 2) Plot effect for all electrodes across all subs
for i_tonedur = 1:length(tonedur_text) %tone duration condition
            % 1.2.1) Define analysis parameters
            fsample = data_ECoGfiltref_trials.fsample;
            toneDur_inSecs  = str2num(tonedur_text{i_tonedur});
            nSamplesPerTone = toneDur_inSecs * fsample;

            win_size    = 25;  %25 data points with fs 512 = ~49ms as mentioned in manuscript
            s_per_win = (1/fsample)*win_size;
            win_overlap = 0; %as mentioned in manuscript

            windows = [1 win_size];
            while windows(end,end) < nSamplesPerTone
                windows = [windows; windows(end,:) + (win_size - win_overlap)];
            end

            if windows(end,end) > nSamplesPerTone
                windows(end,:) = [];
            end
            
    for i_win = 1:size(windows,1)

        %determine start & end points for respective window
        w1 = num2str(1000 * windows(i_win, 1) / fsample , 3 ); 
        w2 = num2str(1000 * windows(i_win, 2) / fsample , 3 );


        win_title = ['time window = [' w1 ' ms - ' w2 ' ms]'];
          
        clear dat
        dat.dimord = 'chan_time';
        dat.time   = 0;

        dat.label  = data_AllElec.label; %seleced ECoG chan       
        
        dat.avg = data_AllElec.effect{i_tonedur}{i_win}; %read out k values per sensor        
        dat.sign_elecs = logical(data_AllElec.sign_elecs{i_tonedur}{i_win});
        
        
%         %2.1)Plot 2D line plot with k-prime values for meaningfull channel order (e.g., anterior to posterior)           
%         x = (1:length(dat.avg)); %proxy x axis array
%         y = dat.avg';%proxy y axis array
%         c = [1, 16]; %colorbar
%         
%         h = figure;
%         set(gcf,'units','normalized','outerposition',[0 0 1 1])
%         hold on
%         subFigtitle = {['All subs - ' model_comp_title] ['ToneDur:' tonedur_text{i_tonedur} 'ms - ERF ' win_title ' - p' num2str(pval)]};
%         suptitle(subFigtitle)       
%            
%             plot(x,y)    
%             xlim([1, length(x)])
%             ylim([0, 16])
%             grid on
%             set(gca,'xtick',1:length(x))
%             ax = gca;
%             ax.XMinorGrid = 'on';
%             ax.XAxis.TickLabel = dat.label;
%             set(gca,'FontSize',10,'XTickLabelRotation',90) 
%             hold on
%             plot(x(dat.sign_elecs),y(dat.sign_elecs'),'.r','Markersize',30) %highlight sign. channels
        
        
        %2.2)Plot 3D surfe plot with k-prime values for meaningfull channel order (e.g., anterior to posterior)      
        coords = data_AllElec.elec.chanpos; %MNI coordinates for selected electrodes
        vals = dat.avg ; %parameter of interest for resp. electrodes
        sign_index = dat.sign_elecs;
        chanSize = (vals./vals)*3.5; %electrode size for sign elecs (arbitrary)
        clims = Scale_Limits; %colormap limits
%         clims = [1,max(vals)]; %colormap limits
        cmap = 'jet';       
        view_angle = [270,0]; %(270,0) - LH, (90,0) = RH        
        c = NASTD_Plot_PlotSignElecsSurf(coords,vals,sign_index,chanSize,clims,cmap,view_angle);        
        Figtitle = {['All subs (n = ' num2str(length(sub_list)) ') - ' Fig_title],['ToneDur:' tonedur_text{i_tonedur} 'ms - ERF ' win_title ' - p' num2str(pval)]};
        title(Figtitle)
        
        if saveplot
            filename     = [effect2plot '_AllSubs_SignElecs2MNIsurf_' w1 '-' w2 'msTW_' tonedur_text{i_tonedur} 'msTD_LH.png'];
            figfile      = [path_fig filename];                
            saveas(gcf, [figfile], 'png'); %save png version
        end
        
        
        
        view_angle = [90,0]; %(270,0) - LH, (90,0) = RH
        c = NASTD_Plot_PlotSignElecsSurf(coords,vals,sign_index,chanSize,clims,cmap,view_angle);
        Figtitle = {['All subs (n = ' num2str(length(sub_list)) ') - ' Fig_title],['ToneDur:' tonedur_text{i_tonedur} 'ms - ERF ' win_title ' - p' num2str(pval)]};
        title(Figtitle)
        
        if saveplot
            filename     = [effect2plot '_AllSubs_SignElecs2MNIsurf_' w1 '-' w2 'msTW_' tonedur_text{i_tonedur} 'msTD_RH.png'];
            figfile      = [path_fig filename];               
            saveas(gcf, [figfile], 'png'); %save png version
        end


    end
    
    close all;
end

end
