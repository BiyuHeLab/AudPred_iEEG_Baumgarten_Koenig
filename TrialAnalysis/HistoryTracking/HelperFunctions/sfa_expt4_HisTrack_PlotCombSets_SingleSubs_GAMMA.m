function sfa_expt4_HisTrack_PlotCombSets_SingleSubs_GAMMA(sub, tonedur_text, baseline_correct, saveplot)
%TJB: New version of tone_history_CV_plot.m, plotting part
%Aim: Plot topoplots of optimal k'values fom combined sets 
%for each sbject, time window, for a specific tone duration condition

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
    path_load = [paths_sfa_expt4.ECoGdata_HisTrack sub '/baseline_corrected/Exp/'];
    path_fig = [paths_sfa_expt4.Fig_HisTrack_SingleSubs sub '/baseline_corrected/'];
    BC_text = 'BC';

elseif baseline_correct == 0
    path_load = [paths_sfa_expt4.ECoGdata_HisTrack sub '/not_baseline_corrected/Exp/'];
    path_fig = [paths_sfa_expt4.Fig_HisTrack_SingleSubs sub '/not_baseline_corrected/'];
    BC_text = 'NBC';
end


%% 1) Load combined k-value data
load([path_load sub '_HisTrackGAMMAEnvSLIDE16_regstatsEXP_' BC_text '_' tonedur_title 'msToneDur_kModelComp_CombSets.mat']);

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

%% 2) Plot measures of association of ERF window with pitch
% model_comp = {'linear_AIC', 'linear_BIC', 'sig2_test', 'R2_test'};
model_comp = {'sig2_test'}; %define which version of model selection used
%sig2_test = sum of squared residuals divided by length of model order
%sig2_test_k = k-prime-value (i.e., preferred model order/Number of tones back best explaining MEG data)

for i_model = 1:length(model_comp)

    %2.1 find zlimits for constant plotting across time windows
    dv_absmax = 0;
    for i_win = 1:size(windows,1)
        dv_win_absmax = eval(['max( abs( ' model_comp{i_model} '_k{i_win} ) );']); %read out k'-value to determine colour bar scaling
        if dv_win_absmax > dv_absmax
            dv_absmax = dv_win_absmax;
        end
    end    
    
    %2.2 Plot data for current time window
    for i_win = 1:size(windows,1)

        %determine start & end points for respective window
        w1 = num2str(1000 * windows(i_win, 1) / fsample , 3 ); 
        w2 = num2str(1000 * windows(i_win, 2) / fsample , 3 );

        switch model_comp{i_model}
            case 'linear_AIC'
                model_comp_title = 'k with best AIC weight for linear regression on tone history';

            case 'linear_BIC'
                model_comp_title = 'k with best BIC weight for linear regression on tone history';

            case 'sig2_test'
                model_comp_title = 'k with minimum test error for linear regression on tone history';

        end

        switch baseline_correct
            case 1
                title_text = ['baseline corrected by first 20 ms of ERF'];
            case 0
                title_text = ['not baseline corrected'];
            case -1
                title_text = ['baseline corrected by first 20 ms of sequence'];
        end

        win_title = ['time window = [' w1 ' ms - ' w2 ' ms]'];

    
           
        clear dat
        dat.dimord = 'chan_time';
        dat.time   = 0;

        %Determine to be plotted channels (selected Chan - rejected Chan)
        index_selectedElecs = data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF; %Choose selected Elecs
        for i_chan = 1:length(subs_PreProcSettings.(sub).rejectedChan_index) %Delete rejected Elecs
        index_selectedElecs(find(index_selectedElecs == subs_PreProcSettings.(sub).rejectedChan_index(i_chan))) = [];
        end

        
        dat.label  = data_ECoGfiltref_trials.label(index_selectedElecs); %seleced ECoG chan    
        
        dv = eval([model_comp{i_model} '_k{i_win}'';']); %read out k values per sensor
%         dat.avg = dv; %all ECoG chan
        dat.avg = dv(index_selectedElecs); %seleced ECoG chan
        
        %1)Plot matrix with k-prime values for meaningfull channel
        %order (e.g., anterior to posterior)           
 
        x = (1:length(dat.avg)); %proxy x axis array
        y = dat.avg';%proxy y axis array
        c = [1, max([dv_absmax, 1.1])]; %colorbar
        
        h = figure;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        hold on
        subFigtitle = {[sub ' - ' model_comp_title] ['ToneDur:' tonedur_text 'ms - ERF ' win_title ' - ' title_text]};
        suptitle(subFigtitle)       
           
            plot(x,y)    
            xlim([1, length(x)])
            ylim([0, max([dv_absmax, 1.1])])
            grid on
            set(gca,'xtick',1:length(x))
            ax = gca;
            ax.XMinorGrid = 'on';
            ax.XAxis.TickLabel = dat.label;
            set(gca,'FontSize',10,'XTickLabelRotation',90)  
     
        if saveplot
            path_fig1 = [path_fig '2DLineplot/Gamma/'];
            mkdir([path_fig1]);
            filename     = [sub '_K_GAMMAEnv_2DLine_' w1 '-' w2 'msTW_' tonedur_text 'msToneDur.png'];

            figfile      = [path_fig1 filename];                

            saveas(gcf, [figfile], 'png'); %save png version

            delete(h);
            close all 
        end

        %2) Plot electrodes as spheres on MNI brain with each sphere color-coded
        %according to k-prime val     
%         coords = data_ECoGfiltref_trials.elec.chanpos(1:subs_PreProcSettings.(sub).number_ECoGchan,;); %all ECoG chan
        coords = data_ECoGfiltref_trials.elec.chanpos(index_selectedElecs,:); %MNI coordinates for selected electrodes
        vals = dat.avg; %parameter of interest for resp. electrodes
        chanSize = (vals./vals)*2.5; %electrode size (arbitrary)
        clims = [1,16]; %colormap limits
%         clims = [1,max(vals)]; %colormap limits
        cmap = 'jet';       
        view_angle = [270,0]; %(270,0) - LH, (90,0) = RH
        c = sfa_expt4_Plot_PlotEleconSurf(coords,vals,chanSize,clims,cmap,view_angle);
        
        Figtitle = {[sub ' - ' model_comp_title] ['ToneDur:' tonedur_text 'ms - ERF ' win_title ' - ' title_text]};
        title(Figtitle)

        
        if saveplot
            path_fig2 = [path_fig '3DSurfplot/Gamma/'];
            mkdir([path_fig2]);
            filename     = [sub '_K_GAMMAEnv_3DSurf_' w1 '-' w2 'msTimeWin_' tonedur_text 'msToneDur.png'];
            figfile      = [path_fig2 filename];                

            saveas(gcf, [figfile], 'png'); %save png version

            delete(h);
            close all
        end  
        
         if strcmp(sub,'NY723') %For subject 4 also plot RH
            view_angle = [90,0]; %(270,0) - LH, (90,0) = RH
            c = sfa_expt4_Plot_PlotEleconSurf(coords,vals,chanSize,clims,cmap,view_angle);

            Figtitle = {[sub ' - ' model_comp_title] ['ToneDur:' tonedur_text 'ms - ERF ' win_title ' - ' title_text]};
            title(Figtitle)
            
            if saveplot
                path_fig2 = [path_fig '3DSurfplot/Gamma/'];
                mkdir([path_fig2]);
                filename     = [sub '_K_GAMMAEnv_3DSurfRH_' w1 '-' w2 'msTimeWin_' tonedur_text 'msToneDur.png'];
                figfile      = [path_fig2 filename];                

                saveas(gcf, [figfile], 'png'); %save png version

                delete(h); 
                close all               
            end         
        
        end
       
        
        
    end
    
%      close all
end
