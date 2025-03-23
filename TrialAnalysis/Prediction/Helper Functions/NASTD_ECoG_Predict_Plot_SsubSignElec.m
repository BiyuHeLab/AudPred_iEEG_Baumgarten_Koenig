function NASTD_ECoG_Predict_Plot_SsubSignElec...
    (sub, ToneDur_text, inputData, ...
    predictive_sequencerange, effect, SamplesTW, Label_TW, BL_label, ...
    FDR_correct, pval_plotting, ...
    save_poststepFigs, paths_NASTD_ECoG, vars)
%Aim: Plot electrodes showing uncorrected prediciton effect on MNI volume for single subjects

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')

%Add base dir and own script dir
addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

path_inputdata = [paths_NASTD_ECoG.ECoGdata_Prediction sub '/Exp/' Label_TW 'sTW/'];
path_fig = [paths_NASTD_ECoG.Fig_Prediction_SingleSubs sub '/Exp/' Label_TW 'sTW/'];

%% 0.2) Determine subject-specific parameters (whole-recording)
NASTD_ECoG_subjectinfo %load subject info file (var: si)
subs_PreProcSettings = NASTD_ECoG_SubPreprocSettings; %load in file with individual preproc infos

% set(0, 'DefaultFigureVisible', 'on');

%Tone duration condition for data load-in
if strcmp(ToneDur_text,'0.2')
    tonedur_title = '200';
elseif  strcmp(ToneDur_text,'0.4')
    tonedur_title = '400';
elseif  strcmp(ToneDur_text,'all')
    tonedur_title = 'all';
end


%% 1) Load prediction effect data and preproc data (for electrode info)
load([path_inputdata sub '_Predp33_regstatsEXP_' inputData '_' tonedur_title 'msTD_' BL_label '.mat']);

%Also load ECoG preproc data for channel labels and position
loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];%path to indiv preprocessed ECoG data
tic
disp(['Loading preprocessed data set for sub: ' sub])
load(loadfile_ECoGpreprocdata);
disp(['done loading in ' num2str(toc) ' sec'])

%% 2) Define analysis parameters
%1.2 Define analysis parameters
fsample = preprocData_AllTrials.fsample;
if  ~strcmp(ToneDur_text,'all')
    toneDur_inSecs  = str2num(ToneDur_text);
elseif strcmp(ToneDur_text,'all')
    toneDur_inSecs  = 0.2;
end
nSamplesPerTone = toneDur_inSecs * fsample;
i_SelectChan = setdiff(preprocData_AllTrials.cfg.info_elec.selected.index4EDF, subs_PreProcSettings.(sub).rejectedChan_index); %selected and clean ECoG chan

win_size    = SamplesTW;  %25 data points with fs 512 = ~49ms as mentioned in manuscript
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
disp(['-Plotted effect = ' effect])
%3.1) Decide which effect should be plotted
%Pred = t-values for linear regression between ERF at window X of
%tone 33 and predicted tone pitch for tone 34
%PredError = t-value for linear regression ERF at
%window X of tone 34, prediction error (difference between actually
%presented final tone (log(p34)) and computationally determined final
%tone), and tone pitch of final tone (p34)
%SimplePredError = t-value for linear regression ERF at
%window X of tone 34, simple prediction error (difference between actually
%presented final tone log(p34) and log(p33)

for i_k = 1:length(predictive_sequencerange) %loop across tones on which prediction is focused on
    %3.2 find zlimits from t-values for constant plotting across time windows
    dv_absmax = 0;
    for i_win = 1:size(windows,1)
        dv_win_absmax = eval(['max( abs( ' effect 'Effect.tval{i_k}{i_win} ) )']);
        if dv_win_absmax > dv_absmax
            dv_absmax = dv_win_absmax;
        end
    end
    
    %3.3 Plot data for current time window
    for i_win = 1:size(windows,1)
        
        %3.3.1 Set up title information
        w1 = num2str( 1000 * windows(i_win, 1) / fsample , 3 );
        w2 = num2str( 1000 * windows(i_win, 2) / fsample , 3 );
        
        switch effect
            case 'Pred'
                effect_title = ['Prediction (predict p*34 based on p33 activity [lin reg t-stat]) - ' inputData ' - ' BL_label];
                pval = PredEffect.pval{i_k}{i_win};
                
            case 'PredError'
                effect_title = ['Prediction error (explain p34 activity by absolute p*34-p34 difference [lin reg t-stat]) - ' inputData ' - ' BL_label];
                pval = PredErrorEffect.pval{i_k}{i_win};
                
            case 'SimplePredError'
                effect_title = ['Simple Prediction error (explain p34 activity by absolute p33-p34 difference [lin reg t-stat]) - ' inputData ' - ' BL_label] ;
                pval = SimplePredErrorEffect.pval{i_k}{i_win};
        end
        
        
       
        
        win_title = ['TW = [' w1 ' - ' w2 ' ms]'];
        
        %3.3.2 Set up data struct for plotting
        clear dat
        dat.dimord = 'chan_time';
        dat.time   = 0;
        dat.label  = i_SelectChan;
        dv = eval([effect 'Effect.tval' '{i_k}{i_win}'';']);
        dat.avg = dv; %
        
        if FDR_correct == 1
            pval_FDRcorrect = mafdr(pval,'BHFDR', true);
            dat.sign_elecs = pval_FDRcorrect < pval_plotting;
            sumSignELec = sum(pval_FDRcorrect < pval_plotting);
            FDR_label = 'FDRcorrected';
            dat.sign_elecs = pval_FDRcorrect < pval_plotting;          
        else
            dat.sign_elecs = pval < pval_plotting;
            FDR_label = 'uncorrected';
            sumSignELec = sum(pval < pval_plotting);
            dat.sign_elecs = pval < pval_plotting;           
        end 
                
%         %1)Plot matrix with association measure for meaningfull channel
%         %order (e.g., anterior to posterior)
%         x = (1:length(dat.avg)); %proxy x axis array
%         y = dat.avg';%proxy y axis array
%         c = [-(max([dv_absmax, 1.1])), max([dv_absmax, 1.1])]; %colorbar
%         
%         h = figure;
%         set(gcf,'units','normalized','outerposition',[0 0 1 1])
%         hold on
%         subFigtitle = {[sub ' - ' effect_title] ['ToneDur: ' ToneDur_text 'ms - ERP ' win_title '; p < ' num2str(pval_plotting) ' (uncorrected)']};
%         suptitle(subFigtitle)
%         
%         plot(x,y)
%         xlim([1, length(x)])
%         ylim([-(max([dv_absmax, 1.1])), max([dv_absmax, 1.1])])
%         grid on
%         set(gca,'xtick',1:length(x))
%         ax = gca;
%         ax.XMinorGrid = 'on';
%         ax.XAxis.TickLabel = dat.label;
%         set(gca,'FontSize',10,'XTickLabelRotation',90)
%         hold on
%         plot(x(dat.sign_elecs),y(dat.sign_elecs),'.r','Markersize',30) %highlight sign. channels
%         
%         if save_poststepFigs == 1
%            path_fig1 = [path_fig '2Dline/' effect '/' FDR_label '/pval' num2str(pval_plotting) '/'];
%             mkdir([path_fig1]);
%             filename     = ['2Dline_' effect '_' inputData '_' w1 '-' w2 'msTW_' ToneDur_text 'msTD_' sub '_' BL_label '.png'];
%             figfile      = [path_fig1 filename];
%             
%             saveas(gcf, [figfile], 'png'); %save png version
%             
%             close all;
%         end        
        
        %2) Plot electrodes as spheres on MNI brain with each sphere color-coded
        coords = preprocData_AllTrials.elec.chanpos(i_SelectChan,:); %MNI coordinates for selected electrodes
        vals = dat.avg; %parameter of interest for resp. electrodes
        chanSize = (vals./vals)*3; %electrode size (arbitrary)
        sign_index = dat.sign_elecs;
        clims = [-(max([dv_absmax, 1.1])), max([dv_absmax, 1.1])]; %free scaling
        %                 clims = [-3 +3]; %fixed scaling
        cmap = 'jet';
        view_angle = [270,0]; %(270,0) - LH, (90,0) = RH
        
        NASTD_ECoG_Plot_PlotSignElecsSurf(coords,vals,sign_index,chanSize,clims,cmap,view_angle);
        
        Figtitle =  {[sub ' - ' effect_title] ['ToneDur: ' ToneDur_text 'ms - ERP ' win_title '; p < ' num2str(pval_plotting) '  (' FDR_label ')'] ...
            ['Sign. elecs:' preprocData_AllTrials.label{i_SelectChan(dat.sign_elecs)} '; ' num2str(sumSignELec) ' Elecs sign.']};
        title(Figtitle)
        
        if save_poststepFigs == 1
            path_fig1 = [path_fig '3Dsurf/' effect '/' FDR_label '/pval' num2str(pval_plotting) '/'];
            mkdir([path_fig1]);
            filename     = ['3DsurfLH_' effect '_' inputData '_' w1 '-' w2 'msTW_' ToneDur_text 'msTD_' sub '_' BL_label '.png'];
            figfile      = [path_fig1 filename];
            
            saveas(gcf, [figfile], 'png'); %save png version
            
            close all;
        end
                
        if strcmp(sub,'NY723') %For subject 4 also plot RH
            view_angle = [90,0]; %(270,0) - LH, (90,0) = RH
            NASTD_ECoG_Plot_PlotSignElecsSurf(coords,vals,sign_index,chanSize,clims,cmap,view_angle);
            
        Figtitle =  {[sub ' - ' effect_title] ['ToneDur: ' ToneDur_text 'ms - ERP ' win_title '; p < ' num2str(pval_plotting) '  (' FDR_label ')'] ...
            ['Sign. elecs:' preprocData_AllTrials.label{i_SelectChan(dat.sign_elecs)} '; ' num2str(sumSignELec) ' Elecs sign.']};
        title(Figtitle)
            
            if save_poststepFigs == 1
                path_fig1 = [path_fig '3Dsurf/' effect '/' FDR_label '/pval' num2str(pval_plotting) '/'];
                mkdir([path_fig1]);
                filename     = ['3DsurfRH_' effect '_' inputData '_' w1 '-' w2 'msTW_' ToneDur_text 'msTD_' sub '_' BL_label '.png'];
                figfile      = [path_fig1 filename];
                
                saveas(gcf, [figfile], 'png'); %save png version
                
                close all;
            end
            
        end    
        
    end
end

end