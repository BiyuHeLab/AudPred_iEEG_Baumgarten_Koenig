function sfa_expt4_Predict_PlotPredict_SingleSub(sub, baseline_correct, tonedur_text, ks, toneIndex, save_poststepFigs)
%TJB: Based on original script 'continuous_prediction_jenn.m'
%Aim: Plot prediciton effect and prediction-error per channel for single subjects
  
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

set(0, 'DefaultFigureVisible', 'on');

%% 0.2) Determine subject- and condition- specific parameters (whole-recording)
sfa_expt4_subjectinfo %load subject info file (var: si)
subs_PreProcSettings = sfa_expt4_subjectPreProcSettings; %load in file with individual preproc infos

%Tone duration condition for data load-in
if strcmp(tonedur_text,'0.2')
    tonedur_title = '200';
elseif  strcmp(tonedur_text,'0.4')
    tonedur_title = '400';
end

if baseline_correct == 1
    path_load = [paths_sfa_expt4.ECoGdata_Prediction sub '/baseline_corrected/Exp/'];
    path_fig = [paths_sfa_expt4.Fig_Prediction_SingleSubs sub '/baseline_corrected/'];
    BC_text = 'BC';

elseif baseline_correct == 0
    path_load = [paths_sfa_expt4.ECoGdata_Prediction sub '/not_baseline_corrected/Exp/'];
    path_fig = [paths_sfa_expt4.Fig_Prediction_SingleSubs sub '/not_baseline_corrected/'];
    BC_text = 'NBC';
end

%% 1) Load prediction effect data and raw data (for electrode info)
load([path_load sub '_PredictionERF_regstatsEXP_' BC_text '_' tonedur_title 'msToneDur.mat']);

%Also load ECoG preproc data for channel labels and position
loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];%path to indiv preprocessed ECoG data
load(loadfile_ECoGpreprocdata);

%% 2) Define analysis parameters
%1.2 Define analysis parameters
fsample = data_ECoGfiltref_trials.fsample;
toneDur_inSecs  = str2num(tonedur_text);
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

windows_inms = (windows / fsample) * 1000;

%% 3) Plot measures of association between ERF window and tone pitch

%3.1) Decide which measure of association should be plotted 
%pred_t1 = t-values for linear regression between ERF at window X of 
    %tone 33 and predicted tone pitch for tone 34
%pred_t2 = t-values fr quadratic regression between ERF at window X of 
    %tone 33, predicted tone pitch for tone 34, and tone pitch at tone 33
%pred_err_t1 = prediction error = t-value for linear regression ERF at
    %window X of tone 34, prediction error (difference between actually
    %presented final tone (log(p34)) and computationally determined final
    %tone), and tone pitch of final tone (p34)

%association = {'pred_t1', 'pred_t2', 'pred_err_t1'}; %incl. quadratic version
association = {'pred_t1', 'pred_err_t1'};

for i = 1:length(association) %loop across measures of association
    
    for i_k = 1:length(ks) %loop across tones on which prediction is focused on
    
        if ks(i_k) == 1 && i <= 2 %if only first tone is in focus, do nothing, since there is no prediciton possible
                    
        else
        
           %3.2 find zlimits for constant plotting across time windows
           dv_absmax = 0;
           for i_win = 1:size(windows,1)
                dv_win_absmax = eval(['max( abs( ' association{i} '{i_k}{i_win} ) )']);
                if dv_win_absmax > dv_absmax
                    dv_absmax = dv_win_absmax;
                end
           end

            %3.3 Plot data for current time window
            for i_win = 1:size(windows,1)

                %3.3.1 Set up title information
                w1 = num2str( 1000 * windows(i_win, 1) / fsample , 3 );
                w2 = num2str( 1000 * windows(i_win, 2) / fsample , 3 );

                switch association{i}
                    case 'pred_t1'
                        association_title = 'prediction linear regression t-stat';
                        pval = pred_t1_p{i_k}{i_win};
                        title_text = ['evaluated at tone #' num2str(toneIndex) ', prediction using k = ' num2str(ks(i_k)) ' tones'];

                    case 'pred_t2'
                        association_title = 'prediction quadratic regression t-stat';
                        pval = pred_t2_p{i_k}{i_win};
                        title_text = ['evaluated at tone #' num2str(toneIndex) ', prediction using k = ' num2str(ks(i_k)) ' tones'];

                    case 'pred_err_t1'
                        association_title = 'abs(prediction error) linear regression t-stat';
                        pval = pred_err_t1_p{i_k}{i_win};
                        title_text = ['evaluated at tone #' num2str(toneIndex+1) ', prediction using k = ' num2str(ks(i_k)) ' tones'];

    %                 case 'pred_err_t2'
    %                     assoc_title = 'prediction error quadratic regression t-stat';
    %                     pval = pred_err_t2_p{i_k}{i_win};
    %                     title_text = ['evaluated at tone #' num2str(toneIndex+1) ', prediction using k = ' num2str(ks(i_k)) ' tones'];                    


                end

                switch baseline_correct
                    case 1
                        title_text = [title_text ', baseline corrected by first 20 ms of ERF'];
                    case 0
                        title_text = [title_text ', not baseline corrected'];
                    case -1
                        title_text = [title_text ', baseline corrected by first 20 ms of sequence'];
                end

                win_title = ['time window = [' w1 ' ms - ' w2 ' ms]'];

                %3.3.2 Set up data struct for plotting
                clear dat
                dat.dimord = 'chan_time';
                dat.time   = 0;        
        %         dat.label  = data_ECoGfiltref_trials.label(1:subs_PreProcSettings.(sub).number_ECoGchan); %all ECoG chan
                dat.label  = data_ECoGfiltref_trials.label(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF); %seleced ECoG chan
        
                dv = eval([association{i} '{i_k}{i_win}'';']);
                dv(pval > .05) = 0; 
                
                dat.avg = dv; %all ECoG chan     
                
                %1)Plot matrix with association measure for meaningfull channel
                %order (e.g., anterior to posterior)           
                x = (1:length(dat.avg)); %proxy x axis array
                y = dat.avg';%proxy y axis array
                c = [-(max([dv_absmax, 1.1])), max([dv_absmax, 1.1])]; %colorbar

                h = figure;
                set(gcf,'units','normalized','outerposition',[0 0 1 1])
                hold on
                subFigtitle = {[sub ' - ' association_title] ['ToneDur:' tonedur_text 'ms - ERF ' win_title ' - ' title_text]};
                suptitle(subFigtitle)

                for i_subfig = 1:2
                    range = (length(dat.avg)/2*(i_subfig-1))+1:(length(dat.avg)/2*(i_subfig));

                subplot(2,1,i_subfig);      
                    plot(x(range),y(range))    

                    xlim([range(1), range(end)])
                    ylim([-(max([dv_absmax, 1.1])), max([dv_absmax, 1.1])])
                    grid on
                    set(gca,'xtick',1:length(x))
                    ax = gca;
                    ax.XMinorGrid = 'on';
                    ax.XAxis.TickLabel = dat.label(range);
                    set(gca,'FontSize',10,'XTickLabelRotation',90)            
                end    

                if save_poststepFigs == 1
                    path_fig1 = [path_fig '2DLineplot/'];
                    mkdir([path_fig1]);
                    filename     = [sub '_' association{i} '_' w1 '-' w2 'msTimeWin_' tonedur_text 'msToneDur.png'];
                    figfile      = [path_fig1 filename];                

                    saveas(gcf, [figfile], 'png'); %save png version

                    delete(h);            
                end                
                
                
                %2) Plot electrodes as spheres on MNI brain with each sphere color-coded
                %according to k-prime val     
        %         coords = data_ECoGfiltref_trials.elec.chanpos(1:subs_PreProcSettings.(sub).number_ECoGchan,;); %all ECoG chan
                coords = data_ECoGfiltref_trials.elec.chanpos(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF,:); %MNI coordinates for selected electrodes
                vals = dat.avg; %parameter of interest for resp. electrodes
                chanSize = (vals./vals)*2.5; %electrode size (arbitrary)
                clims = [-(max([dv_absmax, 1.1])), max([dv_absmax, 1.1])]; %colorbar
        %         clims = [1,max(vals)]; %colormap limits
                cmap = 'jet';       
                view_angle = [270,0]; %(270,0) - LH, (90,0) = RH
                c = sfa_expt4_Plot_PlotEleconSurf(coords,vals,chanSize,clims,cmap,view_angle);

                Figtitle = {[sub ' - ' association_title] ['ToneDur:' tonedur_text 'ms - ERF ' win_title ' - ' title_text]};
                title(Figtitle)

                if save_poststepFigs == 1
                    path_fig2 = [path_fig '3DSurfplot/'];
                    mkdir([path_fig2]);
                    filename     = [sub '_' association{i} '_' w1 '-' w2 'msTimeWin_' tonedur_text 'msToneDur.png'];
                    figfile      = [path_fig2 filename];                

                    saveas(gcf, [figfile], 'png'); %save png version

                    delete(h); 
                    close all
                    
                end 
        
        if strcmp(sub,'NY723') %For subject 4 also plot RH
            view_angle = [90,0]; %(270,0) - LH, (90,0) = RH
            c = sfa_expt4_Plot_PlotEleconSurf(coords,vals,chanSize,clims,cmap,view_angle);

            Figtitle = {[sub ' - ' association_title] ['ToneDur:' tonedur_text 'ms - ERF ' win_title ' - ' title_text]};
            title(Figtitle)
            
            if save_poststepFigs == 1
                path_fig2 = [path_fig '3DSurfplot/'];
                mkdir([path_fig2]);
                filename     = [sub '_' association{i} '_' w1 '-' w2 'msTimeWin_' tonedur_text 'msToneDur_RH.png'];
                figfile      = [path_fig2 filename];                

                saveas(gcf, [figfile], 'png'); %save png version

                delete(h);            
            end         
        
        end
        
                
            end
        end
    end
    
end

end

