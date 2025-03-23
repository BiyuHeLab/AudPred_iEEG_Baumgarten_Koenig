function NASTD_ECoG_Predict_PlotSignElec_SsubperTW...
    (sub, ToneDur_text, inputData, effect, ...
    predictive_sequencerange, SamplesTW, Label_TW, BL_label,...
    FDR_correct, pval_plotting,...
    plot_SubplotperTW, save_poststepFigs, paths_NASTD_ECoG, vars)

%Aim: Plot electrodes showing sign. prediction and prediction error effects
%per input data & TD condition
%Output: Per Subject and TW, plot 1) 1 surface plot where color-coding indicates p-value

% set(0, 'DefaultFigureVisible', 'on');

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')

%Add base dir and own script dir
addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

path_inputdata = [paths_NASTD_ECoG.ECoGdata_Prediction sub '/' Label_TW 'sTW/'];
path_fig = [paths_NASTD_ECoG.Fig_Prediction_SingleSubs sub '/' Label_TW 'sTW/'];

%% 0.2) Determine subject-specific parameters
NASTD_ECoG_subjectinfo %load subject info file (var: si)
subs_PreProcSettings = NASTD_ECoG_SubPreprocSettings; %load in file with individual preproc infos

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

%% 2) Define analysis parameters (based on shorter TD)
%2.1 Define analysis parameters
fsample = preprocData_AllTrials.fsample;
if strcmp(ToneDur_text,'all')
    toneDur_inSecs  = 0.2;
else
    toneDur_inSecs  = str2num(ToneDur_text);
end
nSamplesPerTone = toneDur_inSecs * fsample;

IndexSelElecs = setdiff(preprocData_AllTrials.cfg.info_elec.selected.index4EDF, ...
    subs_PreProcSettings.(sub).rejectedChan_index); %selected and clean ECoG chan

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

%% 3) find zlimits from t-values for constant plotting across time windows
% disp(['-Plotted effect = ' effect])
for i_k = 1:length(predictive_sequencerange) %loop across tones on which prediction is focused on
    dv_absmax = 0;
    for i_win = 1:size(windows,1)
        dv_win_absmax = eval(['max( abs( ' effect 'Effect.tval{i_k}{i_win} ) )']);
        if dv_win_absmax > dv_absmax
            dv_absmax = dv_win_absmax;
        end
    end
    
    %% 4) Plot figure (surface)
    %4.1 Set up figure for TW-subplotting
    if plot_SubplotperTW == 1
        h = figure;
        set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
        
        if strcmp(sub,'NY723') %both hemispheres
            DimSubplot = [4,size(windows,1)/2];
        else
            DimSubplot = [2,size(windows,1)/2];
        end
        CounterSubplot = 1;
        CounterSurfplot = 1;

        SizeFactor = 4;

    end
    
    %4.2 Set up input data for current TW
    for i_win = 1:size(windows,1)
        if plot_SubplotperTW == 0 %1 Figure per TW
            h = figure;
            set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
            if strcmp(sub,'NY723') %both hemispheres
                DimSubplot = [1,2];
            else
                DimSubplot = [1,1];
            end
            CounterSubplot = 1;
            CounterSurfplot = 1;
            
            SizeFactor = 3;
        end
        
        %4.2.1 Determine sign. electrodes
        switch effect
            case 'Pred'
                pval = PredEffect.pval{i_k}{i_win};
            case 'PredError'
                pval = PredErrorEffect.pval{i_k}{i_win};
            case 'SimplePredError'
                pval = SimplePredErrorEffect.pval{i_k}{i_win};
        end
        
        if FDR_correct == 1
            FDR_label = 'FDRcorrected';
            pval_FDRcorrect = mafdr(pval,'BHFDR', true);
            Index_SuprathresElecs{i_win} = pval_FDRcorrect < pval_plotting;
            sumSignELec = sum(pval_FDRcorrect < pval_plotting);
            pval_inputplot = pval_FDRcorrect;
        else
            FDR_label = 'uncorrected';
            Index_SuprathresElecs{i_win} = pval < pval_plotting;
            sumSignELec = sum(pval < pval_plotting);
            pval_inputplot = pval;
        end
        
        %4.2.2 Set up title information
        w1 = num2str( 1000 * windows(i_win, 1) / fsample , 3 );
        w2 = num2str( 1000 * windows(i_win, 2) / fsample , 3 );
        win_title = ['TW = [' w1 ' - ' w2 ' ms]'];
        
        switch effect
            case 'Pred'
                effect_title = ['Prediction (predict p*34 based on p33 activity) - ' ...
                    tonedur_title 'msTD - '  inputData ' - ' BL_label ' - ' FDR_label];
            case 'PredError'
                effect_title = ['Prediction error (explain p34 activity by absolute p*34-p34 difference) - ' ...
                    tonedur_title 'msTD - '  inputData ' - ' BL_label ' - ' FDR_label];
            case 'SimplePredError'
                effect_title = ['Simple Prediction error (explain p34 activity by absolute p33-p34 difference) - ' ...
                    tonedur_title 'msTD - '  inputData ' - ' BL_label ' - ' FDR_label] ;
        end
        
        %4.2.3 Determine labels of sign. elecs
        index_signelec_allelec = IndexSelElecs(find(Index_SuprathresElecs{i_win}));
        label_signelec = [];
        if ~isempty(index_signelec_allelec)
            for i_signelec = 1:length(index_signelec_allelec)
                label_signelec = [label_signelec ' ' preprocData_AllTrials.label{index_signelec_allelec(i_signelec)}];
            end
             sign_title = [num2str(sumSignELec)...
                ' / ' num2str(length(IndexSelElecs)) ' supra-thresh. elecs']; %1 line
        else
             sign_title = ['0 / ' num2str(length(IndexSelElecs)) ' supra-thresh. elecs']; %1 line 
        end
        
        %4.2.4 Set up data struct for plotting
        clear dat
        dat.dimord = 'chan_time';
        dat.time   = 0;
        dat.label  = IndexSelElecs;
        dv = eval([effect 'Effect.tval' '{i_k}{i_win}'';']);
        dat.avg = dv; %
        dat.sign_elecs = Index_SuprathresElecs{i_win};
        sign_index = Index_SuprathresElecs{i_win};
        
        
        %4.2.5 Plot electrodes as spheres on MNI brain with each sphere color-coded
        %depending on p-value
        coords = preprocData_AllTrials.elec.chanpos(IndexSelElecs,:); %MNI coordinates for selected electrodes
%         vals = dat.avg; %t-value lin reg
%         clims = [-4 4]; %fixed scaling

        vals = pval_inputplot; %p-value
        clims = [0 0.05]; %fixed scaling
       
        chanSize = ones(1,length(dat.avg))*SizeFactor; %electrode size (arbitrary)
%         clims = [max([dv_absmax, 1.1])*-1 max([dv_absmax, 1.1])]; %free scaling
        cmap = 'jet';
        view_angle = [270,0];        
        
        %4.2.6Adjust subplot position according to number of subplots
        if plot_SubplotperTW == 1 %All TW in 1 plot (2 subplot (surface + 2D plot) per TW)
            SubplotPosition = [0 0 0 0];
            ColorbarPosition = [0 0 0 0];
            %         suptitle([sub ' - ' inputData ' - Kprime Across TD'])
        else %1 Plot per TW (2 subplot (surface + 2D plot) per TW)
            SubplotPosition = [0 0 0 DimSubplot(2)*0.5];
            ColorbarPosition = [0.1 0.3 0 0.3];
        end
        
        if strcmp(sub,'NY723') %both hemispheres
            sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_SubplotHemis...
                (coords,preprocData_AllTrials.label(IndexSelElecs),vals,sign_index,...
                chanSize,clims,cmap,DimSubplot,CounterSubplot,SubplotPosition,ColorbarPosition,...
                []);
            sp_handle_surf{CounterSurfplot} = sp_handle_surf_temp.L;
            sp_handle_surf{CounterSurfplot+1} = sp_handle_surf_temp.R;            
            CounterSubplot = CounterSubplot+2;
            CounterSurfplot = CounterSurfplot+2;
            
        else %left hemisphere
            sp_handle_surf{CounterSurfplot} = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label...
                (coords,preprocData_AllTrials.label(IndexSelElecs),vals,sign_index,...
                chanSize,clims,cmap,view_angle,DimSubplot,CounterSubplot,SubplotPosition,ColorbarPosition);
            CounterSubplot = CounterSubplot+1;
            CounterSurfplot = CounterSurfplot+1;
        end
        
        %Adjust header/title
        Figtitle = [sub ' - ' effect_title];
        if plot_SubplotperTW == 1 %All TW in 1 plot (2 subplot (surface + 2D plot) per TW)
        	title({[win_title] [sign_title]},'FontSize',10)
        else %1 Plot per TW (2 subplot (surface + 2D plot) per TW)
            suptitle({[Figtitle] [win_title  ' - ' sign_title]})
        end
        
        %% 5) Save Figure
        if plot_SubplotperTW == 0 %1 plot per TW
            
            %Adjust colorbar 
            ColorbarPosition = [0.1 0.2 0 0.75];
            h = colorbar;
            h.Position(1) = ColorbarPosition(1); %sets colorbar to the right
            h.Position(2) = ColorbarPosition(2); %sets colorbar higher
            h.Position(4) = h.Position(4)*ColorbarPosition(4); %makes colorbar shorter
            h.Label.String = ['t-statistic (lin reg p33 activity on p*34)'];
            h.FontSize = 12;
            caxis(clims)
            %                 h.Ticks = [0, 1 , 2];
            
            if save_poststepFigs == 1
                path_fig1 = [path_fig '3Dsurf/' effect '/' FDR_label '/pval' num2str(pval_plotting) '/'];
                mkdir([path_fig1]);
                filename     = ['3DsurfRH_' effect '_' inputData '_' w1 '-' w2 'msTW_' ToneDur_text 'msTD_' sub '_' BL_label '.png'];
                figfile      = [path_fig1 filename];
                saveas(gcf, [figfile], 'png'); %save png version
                close all;
            end
        end
        
    end %End TW loop
    
    if plot_SubplotperTW == 1 %All TW in 1 plot
        %Adjust subplot sizes
        if strcmp(sub,'NY723') %both hemispheres
            SubplotPosition = [0 0 0 0.135]; %slightly enlarge
            ColorbarPosition = [0.1 0.25 0 8];
        else
            SubplotPosition = [0 0 0 0.1]; %slightly enlarge
            ColorbarPosition = [0.1 0.3 0 1.5];
        end
        
        for i_surfplot = 1:length(sp_handle_surf)
            if strcmp(sub,'NY723') %both hemispheres
                if i_surfplot <= DimSubplot(2)
                    pos = get(sp_handle_surf{i_surfplot},'Position');
                    newpos = pos;
                    newpos(4) = pos(4) + SubplotPosition(4);
                    set(sp_handle_surf{i_surfplot},'Position',newpos)
                else
                    pos = get(sp_handle_surf{i_surfplot},'Position');
                    pos_above = get(sp_handle_surf{i_surfplot-DimSubplot(2)},'Position');
                    newpos = pos;
                    newpos(2) = pos_above(2) - 0.225;
                    newpos(4) = pos(4)  + SubplotPosition(4);
                    set(sp_handle_surf{i_surfplot},'Position',newpos)                    
                end
            else
                pos = get(sp_handle_surf{i_surfplot},'Position');
                newpos = pos + SubplotPosition;
                set(sp_handle_surf{i_surfplot},'Position',newpos)
            end
            
            if i_surfplot == 1
                %Adjust title
                %Adjust colorbar
                h = colorbar;
                h.Position(1) = ColorbarPosition(1); %sets colorbar to the right
                h.Position(2) = ColorbarPosition(2); %sets colorbar higher
                h.Position(4) = h.Position(4)*ColorbarPosition(4); %makes colorbar shorter
%                 h.Label.String = ['t-statistic (lin reg p33 activity on p*34)'];
                h.Label.String = ['p-value (' FDR_label '; t-stat for lin reg p33 activity on p*34)'];
                h.FontSize = 12;
                caxis(clims)
%                 h.Ticks = [0, 1 , 2];
            end            
        end
        suptitle(Figtitle)

        %Save
        if save_poststepFigs == 1
            path_fig1 = [path_fig '3Dsurf/' effect '/' FDR_label '/pval' num2str(pval_plotting) '/'];
            mkdir([path_fig1]);
            filename     = [effect '_' inputData '_allTW_' ToneDur_text 'msTD_' sub '_' BL_label '.png'];
            figfile      = [path_fig1 filename];
            saveas(gcf, [figfile], 'png'); %save png version
            close all;
        end
    end
    
end

end
