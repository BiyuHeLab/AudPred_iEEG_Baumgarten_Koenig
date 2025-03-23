function NASTD_ECoG_Predict_Plot_AllsubSignElec...
    (subs, ToneDur_text, inputData, ...
    predictive_sequencerange, effect, SamplesTW, Label_TW, BL_label,...
    FDR_correct, pval_plotting, ...
    plot_SubplotperTW, save_poststepFigs, paths_NASTD_ECoG, vars)
%Aim: Plot electrodes showing uncorrected prediciton effect on MNI volume for all subjects

% ToneDur_text = '0.2'
% inputData = 'LF'
% effect = 'Pred'

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/NASTD_ECoG_Matlab/Plotting/')

%Add base dir and own script dir
addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

path_fig = [paths_NASTD_ECoG.Fig_Prediction_SingleSubs 'Allsub/Exp/' Label_TW 'sTW/'];

%% 0.2) Determine analysis-specific parameters (whole-recording)
% set(0, 'DefaultFigureVisible', 'on');

%Tone duration condition for data load-in
if strcmp(ToneDur_text,'0.2')
    tonedur_title = '200';
elseif  strcmp(ToneDur_text,'0.4')
    tonedur_title = '400';
elseif  strcmp(ToneDur_text,'all')
    tonedur_title = 'all';
end

%% 1) Load prediction effect data and preproc data (for electrode info) and copy relevant info into var storing all subs
coords_allsub = [];
label_allsub = [];
sub_index = [];

for i_sub = 1:length(subs)
    
    sub = subs{i_sub};
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    subs_PreProcSettings = NASTD_ECoG_SubPreprocSettings; %load in file with individual preproc infos
    
%     path_inputdata = [paths_NASTD_ECoG.ECoGdata_Prediction sub '/Exp/' Label_TW 'sTW/'];  
    path_inputdata = [paths_NASTD_ECoG.ECoGdata_Prediction sub '/' Label_TW 'sTW/'];  

    load([path_inputdata sub '_Predp33_regstatsEXP_' inputData '_' tonedur_title 'msTD_' BL_label '.mat']);
    
    %Also load ECoG preproc data for channel labels and position
    loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];%path to indiv preprocessed ECoG data
    tic
    disp(['Loading preprocessed data set for sub: ' sub])
    load(loadfile_ECoGpreprocdata);
    disp(['done loading in ' num2str(toc) ' sec'])
    
    %Store electrode coordinates
    i_SelectChan = setdiff(preprocData_AllTrials.cfg.info_elec.selected.index4EDF, subs_PreProcSettings.(sub).rejectedChan_index);
    
    coords_sub{i_sub} = preprocData_AllTrials.hdr.elec.chanpos...
        (i_SelectChan,:);
    coords_allsub = [coords_allsub; coords_sub{i_sub}];
    
    for i_chan = 1:length(i_SelectChan)
        label_allsub{end+1} = [preprocData_AllTrials.label{i_SelectChan(i_chan)} '_' sub];
    end
    %Store t and p-values from analyis
    tval{i_sub} = eval([effect 'Effect.tval']);
    sub_index = [sub_index; ones(length(coords_sub{i_sub}),1)*i_sub]; %used to differentiate subjects
    pval{i_sub} = eval([effect 'Effect.pval']);
    
    %Optional: FDR-correct pvals per sub and win
    for i_win = 1:length(pval{i_sub}{1})
%          [~,~,~,pval_FDRcorrect{i_sub}{1}{i_win}] = fdr_bh(pval{i_sub}{1}{i_win},0.05, 'pdep','no');        
         pval_FDRcorrect{i_sub}{1}{i_win} = mafdr(pval{i_sub}{1}{i_win},'BHFDR', true);
    end
        
    for i_win = 1:length(tval{i_sub}{1})
        if i_sub == 1
            tval_allsubs{i_win} = tval{i_sub}{1}{i_win}';
            pval_allsubs{i_win} = pval{i_sub}{1}{i_win}';   
            pval_FDRcorrect_allsubs{i_win} = pval_FDRcorrect{i_sub}{1}{i_win}';           
        else
            tval_allsubs{i_win} = [tval_allsubs{i_win}; tval{i_sub}{1}{i_win}'];
            pval_allsubs{i_win} = [pval_allsubs{i_win}; pval{i_sub}{1}{i_win}'];
            pval_FDRcorrect_allsubs{i_win} = [pval_FDRcorrect_allsubs{i_win}; pval_FDRcorrect{i_sub}{1}{i_win}'];   
        end
    end
    
end
label_allsub = label_allsub';

%% 2) Define analysis parameters
%1.2 Define analysis parameters
fsample = preprocData_AllTrials.fsample;
if  ~strcmp(ToneDur_text,'all')
    toneDur_inSecs  = str2num(ToneDur_text);
elseif strcmp(ToneDur_text,'all')
    toneDur_inSecs  = 0.2;
end
    
nSamplesPerTone = toneDur_inSecs * fsample;

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

if plot_SubplotperTW == 1 %one subplot per TW
    h = figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    
    DimSubplot = [2,size(windows,1)/2];
    CounterSubplot = 1;
end

for i_k = 1:length(predictive_sequencerange) %loop across tones on which prediction is focused on
    %3.2 find zlimits from t-values for constant plotting across time windows
    dv_absmax = 0;
    for i_win = 1:size(windows,1)
        dv_win_absmax = max(abs(tval_allsubs{i_win}));
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
        dat.label  = label_allsub;
        dv = tval_allsubs{i_win};
        
%         dat.avg = abs(dv); %test stat
%         clims = [0 3]; %fixed scaling

        dat.avg = sub_index; %subject indicator
        clims = [1,7]; %colormap limits

        if FDR_correct == 1
            dat.sign_elecs = pval_FDRcorrect_allsubs{i_win} < pval_plotting;
            FDR_label = 'FDRcorrected';
        else
            dat.sign_elecs = pval_allsubs{i_win} < pval_plotting;
            FDR_label = 'uncorrected';
        end
        
        sign_title = [num2str(sum(dat.sign_elecs)) ' sign. elecs'];

        %1) Plot electrodes as spheres on MNI brain with each sphere color-coded
        coords = coords_allsub; %MNI coordinates for selected electrodes
        vals = dat.avg; %parameter of interest for resp. electrodes
        chanSize = (vals./vals)*5; %electrode size (arbitrary)
        sign_index = dat.sign_elecs;
        cmap = 'hot';
        
%         view_angle = [270,0]; %(270,0) - LH, (90,0) = RH
        if plot_SubplotperTW == 0 %one separate figure per TW
            NASTD_ECoG_Plot_PlotSignElecsSurf_SubplotHemis(coords,vals,sign_index,chanSize,clims,cmap) %Both H

            Figtitle =  {['All sub (n = ' num2str(length(subs)) ') '  effect_title] ...
                ['ToneDur: ' ToneDur_text 'ms - ERP ' win_title '; p < ' num2str(pval_plotting) ' (' FDR_label '); ' sign_title]};
            suptitle(Figtitle)
            
             if save_poststepFigs == 1
                path_fig1 = [path_fig '3Dsurf/' effect '/' FDR_label '/pval' num2str(pval_plotting) '/'];

                mkdir([path_fig1]);
                filename     = ['3DsurfLH_' effect '_' inputData '_' w1 '-' w2 'msTW_' ToneDur_text 'msTD_Allsub_' BL_label '.png'];
                figfile      = [path_fig1 filename];

                saveas(gcf, [figfile], 'png'); %save png version

                close all;
             end
        
        else
            NASTD_ECoG_Plot_SubplotSignElecsSurf_SubplotHemis...
                (coords,vals,sign_index,chanSize,clims,cmap,...
                DimSubplot,CounterSubplot,...
                [win_title '; ' sign_title])
            
            CounterSubplot = CounterSubplot+2;
            
            Figtitle =  {['All sub (n = ' num2str(length(subs)) ') '  effect_title] ...
                ['ToneDur: ' ToneDur_text 'ms; p < ' num2str(pval_plotting) ' (' FDR_label ')']};
            
        end
       
%         %Separate RH figure
%         view_angle = [90,0]; %(270,0) - LH, (90,0) = RH
%         c = NASTD_ECoG_Plot_PlotSignElecsSurf(coords,vals,sign_index,chanSize,clims,cmap,view_angle);
%         
%         Figtitle =  {['All sub (n = ' num2str(length(subs)) ') '  effect_title] ['ToneDur: ' ToneDur_text 'ms - ERP ' win_title '; p < ' num2str(pval_plotting) ' (uncorrected)'] ...
%             ['Sign. elecs:' label_allsub{dat.sign_elecs}]};
%         title(Figtitle)
%         
%         if save_poststepFigs == 1
%            filename     = ['3DsurfRH_' effect '_' inputData '_' w1 '-' w2 'msTW_' ToneDur_text 'msTD_Allsub' FDR_label '_' pval_plotting '.png'];
%             figfile      = [path_fig1 filename];
%             
%             saveas(gcf, [figfile], 'png'); %save png version
%             
%             close all;
%         end
        
    end
    
    if plot_SubplotperTW == 1 %one separate figure per TW

        suptitle(Figtitle)
        
       %Adjust colorbar
        ColorbarPosition = [0.1 0.3 0 3];

        h = colorbar;
        h.Position(1) = ColorbarPosition(1); %sets colorbar to the right
        h.Position(2) = ColorbarPosition(2); %sets colorbar higher
        h.Position(4) = h.Position(4)*ColorbarPosition(4); %makes colorbar shorter
        h.Label.String = ['Subject #'];
        caxis([1 length(subs)])
        h.Ticks = [1:length(subs)];
        h.Label.FontSize = 16;
        
        if save_poststepFigs == 1
            path_fig1 = [path_fig '3Dsurf/' effect '/' FDR_label '/pval' num2str(pval_plotting) '/'];
            mkdir([path_fig1]);
            filename     = [effect '_' inputData '_AllTW_' ToneDur_text 'msTD_Allsub_' BL_label '.png'];
            figfile      = [path_fig1 filename];
            saveas(gcf, [figfile], 'png'); %save png version
            close all;
        end   
    
    end

end

end


