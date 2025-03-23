function NASTD_ECoG_TFA_TFAallTrials...
    (sub, FuncInput_ToneDur_text,  ...
    TFA_FOI, ...
    FuncInput_InputData, ...
    plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG)

%Aim: Compute Time-Frequency Representation across all frequencies using
%Wavelet approach (compromise between both low and high freqs).
%Plot output across electrodes and on 2D surface

%% 0.1) Specify vars, paths, and setup fieldtrip
%Add base dir and own script dir
addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));

%% 1) Determine time points of tone presentation
%1.1 Define time window parameters depending on TD
nTrials     = length(FuncInput_InputData.trial);
nSensors    = size(FuncInput_InputData.trial{1},1);
SampleFreq  = FuncInput_InputData.fsample;
ToneDur_Sec  = str2num(FuncInput_ToneDur_text);
nSamples_perTone = ToneDur_Sec * SampleFreq;

%1.2 Determine TP/samples for each tone start+end
TP_Tone_StartStop = NaN(36,2);
TP_Tone_StartStop(1,1) = 0; %p1 set as t = 0 in trial definition
for i_tone = 2:36
    Dist = abs(FuncInput_InputData.time{1} - ((str2num(FuncInput_ToneDur_text)*i_tone) - str2num(FuncInput_ToneDur_text)));
    minDist = min(Dist);
    i_minDist = find(Dist == minDist);
    TP_Tone_StartStop(i_tone,1) = FuncInput_InputData.time{1}(i_minDist);
end
for i_tone = 1:35
    i_LastSampleTone = find(FuncInput_InputData.time{1} == TP_Tone_StartStop(i_tone + 1,1));
    TP_Tone_StartStop(i_tone,2) = FuncInput_InputData.time{1}(i_LastSampleTone);
end
TP_Tone_StartStop = TP_Tone_StartStop(1:34,:);

%1.3 Check if all tones are of equal length
%(if not then choose min length by deleting the last sample of longer trials)
minSeqLength_sec = min(TP_Tone_StartStop(:,2)-TP_Tone_StartStop(:,1));
for i_tone = 1:34
    if (TP_Tone_StartStop(i_tone,2) - TP_Tone_StartStop(i_tone,1)) > minSeqLength_sec
        TP_Tone_StartStop(i_tone,2) = ...
            FuncInput_InputData.time{1}(...
            find(FuncInput_InputData.time{1} == TP_Tone_StartStop(i_tone,2))-1);
    end
end

%1.4 Determine samples corresponding to TP
Sample_Tone_StartStop = NaN(34,2);
for i_tone = 1:34
    Sample_Tone_StartStop(i_tone,1) = ...
        find(FuncInput_InputData.time{1} == TP_Tone_StartStop(i_tone,1));
    Sample_Tone_StartStop(i_tone,2) = ...
        find(TP_Tone_StartStop(i_tone,2) == FuncInput_InputData.time{1});
end

%% 2) Compute TFA across all freqs

%2.1) NaN-interpolation after IED-ramoval
%Determine how many samples are changed to NaN
Found_NaNs = zeros(length(FuncInput_InputData.trial),length(FuncInput_InputData.label));
PercNaNsamples_perChan = zeros(length(FuncInput_InputData.trial),length(FuncInput_InputData.label));
for i_trial = 1:length(FuncInput_InputData.trial)
    PercNaNsamples_perChan(i_trial,:) = (sum(isnan(FuncInput_InputData.trial{i_trial})')./length(FuncInput_InputData.trial{i_trial})) * 100;
    Found_NaNs(i_trial,:) = any(isnan(FuncInput_InputData.trial{i_trial})');
end
AvgPercNaNsamples_perChan = mean(PercNaNsamples_perChan);
disp(['Average/Max percent NaN-samples per trial = ' ...
    num2str(mean(AvgPercNaNsamples_perChan)) '%/' ...
    num2str(max(AvgPercNaNsamples_perChan)) '%']);

%2.2 Interpolate NaN samples
cfg = [];
cfg.method = 'linear'; %'nearest','linear','spline','pchip','cubic','v5cubic'
cfg.prewindow = 0.01; 
cfg.postwindow = 0.01;
FuncInput_InputData = ft_interpolatenan(cfg, FuncInput_InputData);

Found_NaNs = zeros(length(FuncInput_InputData.trial),length(FuncInput_InputData.label));
for i_trial = 1:length(FuncInput_InputData.trial)
    Found_NaNs(i_trial,:) = any(isnan(FuncInput_InputData.trial{i_trial})');
end
if sum(sum(Found_NaNs)) > 0 
    disp(['Remaining NaNs in data (elec: ' num2str(find(any(Found_NaNs)))]);
    pause
end

%Compute TFA
cfg = [];
cfg.output     = 'pow';
cfg.channel    = FuncInput_InputData.label; %All preprocessed elecs
cfg.foi        = TFA_FOI;
cfg.toi        = -0.5 : 0.05 : (TP_Tone_StartStop(end,2) + 0.4); %500ms prestim until 400ms after last tone ended in steps of 50ms

cfg.padtype    = 'zero';
if strcmp(FuncInput_ToneDur_text, '0.2') %0.2s TD
    cfg.pad = 10; %zero pad out to 10s
elseif strcmp(FuncInput_ToneDur_text, '0.4') %0.4s TD
    cfg.pad = 20; %zero pad out to 20s
end

if TFA_FOI(end) <= 40 %For lower freqs use Hanning with freq-dependent window-length
    
    cfg.method     = 'mtmconvol'; 
    cfg.taper      = 'hanning';
    cfg.t_ftimwin  = 5./cfg.foi; %5 cycles per time window
    cfg.toi        = -0.5 : 0.05 : (TP_Tone_StartStop(end,2) + 0.4); %500ms prestim until 400ms after last tone ended

    TFA_outputdata = ft_freqanalysis(cfg, FuncInput_InputData);
    TFA_type = 'HannFlexTW';
    
elseif TFA_FOI(end) > 40 %For higher freqs use wavelet-based approach
    
    cfg.method     = 'wavelet'; %Morlet wavelets
    cfg.width      = 5; %Width of wavelength in number of cycles

    TFA_outputdata = ft_freqanalysis(cfg, FuncInput_InputData);
    TFA_type = 'Wavelet';
    
end

%Compute power values relative to baseline (relative
%increase or decrease with respect to the power in the baseline interval)
cfg = [];
cfg.baseline = [-0.5 0];
cfg.baselinetype = 'relative';
cfg.parameter = 'powspctrm';
TFA_outputdata = ft_freqbaseline(cfg, TFA_outputdata);

%Determine colorlimits
% clim = [0 ...
%     round(median(max(max(TFA_outputdata_BC.powspctrm))),0)];
clim = [0 2.5];

%Determine output folder
TFA_label = [TFA_type '_' num2str(TFA_FOI(1)) '-' num2str(TFA_FOI(end)) 'Hz'];

% path_dataoutput = [paths_NASTD_ECoG.ECoGdata_TFA '/allTrials/' TFA_label '/'...
%     sub '/Data/'];
% if (~exist(path_dataoutput, 'dir')); mkdir(path_dataoutput); end

path_fig = ([paths_NASTD_ECoG.ECoGdata_TFA '/allTrials/' TFA_label '/' ...
    sub '/Figs/']);
if (~exist(path_fig, 'dir')); mkdir(path_fig); end

%% 3) Subplot TFA representations across all electrodes
if plot_poststepFigs == 1
    NumFig = ceil(nSensors/30); %30 elecs per plot
    
    i_elec = 0;
    for i_fig = 1 : NumFig
        h = figure('visible','off'); %ensures figure doesn't pop up during plotting
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        
        for i_subplot = 1:30
            if i_elec < nSensors
                
                subplot(5,6,i_subplot)
                i_elec = i_elec+1;
                
                imagesc(TFA_outputdata.time, TFA_outputdata.freq, ...
                    squeeze(TFA_outputdata.powspctrm(i_elec,:,:)), ...
                    clim);
                axis xy % flip vertically
                colorbar;
                                               
                for i_tone = [1 10 20 30]
                    hold on; plot([TP_Tone_StartStop(i_tone,1) ...
                        TP_Tone_StartStop(i_tone,1)], [5 150], ...
                        'k--', 'LineWidth', 0.5) %Add vertical lines to indicate tone starts
                end
                
                hold on; plot([TP_Tone_StartStop(33,1) ...
                    TP_Tone_StartStop(33,1)], [5 150], ...
                    'w--', 'LineWidth', 0.5) %Add vertical line to indicate p33 start
%                 hold on; plot([TP_Tone_StartStop(33,2) ...
%                     TP_Tone_StartStop(33,2)], [5 150], ...
%                     'w--', 'LineWidth', 0.5) %Add vertical line to indicate p33 stop
                hold on; plot([TP_Tone_StartStop(34,1) ...
                    TP_Tone_StartStop(34,1)], [5 150], ...
                    'r--', 'LineWidth', 0.5) %Add vertical line to indicate p34 start
                
                if i_subplot == 25 %For each lower left subplot
                    xlabel('Time [s]');
                    ylabel('Freq [Hz]');               
                end
                
                title(['Elec: ' FuncInput_InputData.label{i_elec}])
                
            end
        end
        
        sgtitle({[sub ' - ' FuncInput_ToneDur_text...
            'msTD'] ['TFA (' TFA_type '; ' num2str(TFA_FOI(1)) '-'  num2str(TFA_FOI(end)) ...
            'Hz) across tone sequence - Ratio of power change vs. baseline (-0.5 to 0s; ' ...
            num2str(i_fig) '/' num2str(NumFig) ')']},'Interpreter','none')
        
        if save_poststepFigs == 1
            filename     = ['TFA' num2str(TFA_FOI(1)) '-'  num2str(TFA_FOI(end)) ...
                'Hz_AllElecs_' sub '_' ...
                FuncInput_ToneDur_text 'sTD_' num2str(i_fig) '.png'];
            figfile      = [path_fig filename];
            saveas(gcf, [figfile], 'png'); %save png version
            close
        end
    end
end

%% 4) Multiplot for both hemispheres


% To do: Check if this step is necessary, or of FT gets this done by label alone
% %Adjust electrode coordinates & labels (.elec subfield to .label subfield)
% for i_elec = 1:length(FuncInput_InputData.label)%All preprocessed elecs
%     used_elecs_chanposIndex(1,i_elec) = ...
%         find(strcmp(FuncInput_InputData.label{i_elec}, ...
%         FuncInput_InputData.elec.label));    
% end
% FuncInput_InputData.elec.chanpos = FuncInput_InputData.elec.chanpos(used_elecs_chanposIndex,:); 


%Construct surface models for plotting
LH = load([paths_NASTD_ECoG.FieldTrip '/template/anatomy/surface_pial_left.mat']);
RH = load([paths_NASTD_ECoG.FieldTrip '/template/anatomy/surface_pial_right.mat']);

%Prepare elec layout model (works for both sides)
cfg             = [];
if strcmp(sub, 'NY787')
    cfg.headshape   = RH.mesh;
    cfg.viewpoint   = 'right';
else   
    cfg.headshape   = LH.mesh;
    cfg.viewpoint   = 'left';
end

cfg.projection  = 'orthographic';
cfg.channel     = FuncInput_InputData.label; %All preprocessed elecs
cfg.mask        = 'convex';
% cfg.boxchannel = {'G*', 'G**'}; %Channels determining channel box size
cfg.boxchannel  = {TFA_outputdata.label{1}, TFA_outputdata.label{2}}; %Channels determining channel box size
cfg.overlap     = 'shift';
lay = ft_prepare_layout(cfg, TFA_outputdata);

%Plot multiplot on 2D surface
cfg             = [];
cfg.layout      = lay;
cfg.showoutline = 'yes';
cfg.zlim        = clim;
cfg.colorbar    = 'yes';
cfg.showlabels  = 'yes';
cfg.showscale   = 'no';
cfg.showcomment = 'no';
cfg.fontsize    = 8;

ft_multiplotTFR(cfg, TFA_outputdata);
set(gcf,'units','normalized','outerposition',[0 0 1 1])

sgtitle({[sub ' - ' FuncInput_ToneDur_text...
    'msTD'] ['TFA (' TFA_type '; ' num2str(TFA_FOI(1)) '-'  num2str(TFA_FOI(end)) ...
    'Hz) across tone sequence - Ratio of power change vs. baseline (-0.5 to 0s)']},...
    'Interpreter','none')

if save_poststepFigs == 1
    filename     = ['TFA' num2str(TFA_FOI(1)) '-'  num2str(TFA_FOI(end)) ...
        'Hz_MultiSurf_' sub '_' ...
        FuncInput_ToneDur_text 'sTD.png'];
    figfile      = [path_fig filename];
    saveas(gcf, [figfile], 'png'); %save png version
    close
end

%% 5) Topo-plot a specific freqband
%Determine data selection over freqs and time
cfg = [];

if TFA_FOI(end) < 40
    TFA_FOI_topoplot = {[6 8], [8 12], [15 25]};
elseif TFA_FOI(end) > 40
    TFA_FOI_topoplot = {[30 70], [70 150]};
end

for i_TFA_FOI_topoplot = 1:length(TFA_FOI_topoplot)
    
    cfg.frequency               = TFA_FOI_topoplot{i_TFA_FOI_topoplot};
    cfg.avgoverfreq             = 'yes';
    cfg.latency                 = [TP_Tone_StartStop(1,1) TP_Tone_StartStop(34,2)]; %Sequence start to end
    cfg.avgovertime             = 'yes';
    TFA_outputdata_BC_freqband  = ft_selectdata(cfg, TFA_outputdata);

    %plot on topoplot
    cfg                 = [];
    cfg.parameter       = 'powspctrm';
    cfg.marker          = 'on';
    cfg.markersymbol    = '.';
    cfg.markersize      = 8;
    cfg.layout          = lay;
    cfg.showoutline     = 'yes';
    cfg.zlim            = clim;
    cfg.interplimits    = 'electrodes';
    cfg.shading         = 'interp';
    
    ft_topoplotTFR(cfg, TFA_outputdata_BC_freqband);
    set(gcf,'units','normalized','outerposition',[0 0 1 1])

    sgtitle({[sub ' - ' FuncInput_ToneDur_text...
        'msTD'] ['TFA (' TFA_type '; ' ...
        num2str(TFA_FOI_topoplot{i_TFA_FOI_topoplot}(1)) '-'  num2str(TFA_FOI_topoplot{i_TFA_FOI_topoplot}(end)) ...
        'Hz) averaged across tone sequence - Ratio of power change vs. baseline (-0.5 to 0s)']},...
        'Interpreter','none')

    if save_poststepFigs == 1
        filename     = ['TFA' num2str(TFA_FOI_topoplot{i_TFA_FOI_topoplot}(1)) ...
            '-'  num2str(TFA_FOI_topoplot{i_TFA_FOI_topoplot}(end)) ...
            'Hz_TopoAvgFreq_' sub '_' ...
            FuncInput_ToneDur_text 'sTD.png'];
        figfile      = [path_fig filename];
        saveas(gcf, [figfile], 'png'); %save png version
        close
    end

%     %Plot on surface brain
%     cfg                 = [];
%     cfg.funparameter    = 'powspctrm';
%     cfg.funcolorlim     = [0 2];
%     cfg.method          = 'surface';
%     cfg.interpmethod    = 'sphere_weighteddistance';
%     cfg.sphereradius    = 10;
%     cfg.camlight        = 'no';
%     cfg.funcolormap     = 'parula';
%     cfg.renderer        = 'painters';
%     cfg.surfdownsample = 10;  % downsample to speed up processing
%     ft_sourceplot(cfg, TFA_outputdata_BC_freqband, LH.mesh);
%     set(gcf,'units','normalized','outerposition',[0 0 1 1])
% 
%     if strcmp(sub,'NY723')
%         view([0 90]); %Above        
%     elseif strcmp(sub,'NY787')
%         view([-90 20]); %RH lateral
%     else
%         view([-90 20]); %LH lateral
%     end
%     material dull;
%     lighting gouraud;
%     camlight;
    
%     sgtitle({[sub ' - ' FuncInput_ToneDur_text...
%         'msTD'] ['TFA (' TFA_type '; ' ...
%         num2str(TFA_FOI_topoplot{i_TFA_FOI_topoplot}(1)) '-'  num2str(TFA_FOI_topoplot{i_TFA_FOI_topoplot}(end)) ...
%         'Hz) averaged across tone sequence - Ratio of power change vs. baseline (-0.5 to 0s)']},...
%         'Interpreter','none')
% 
%     if save_poststepFigs == 1
%         filename     = ['TFA' num2str(TFA_FOI_topoplot{i_TFA_FOI_topoplot}(1)) ...
%             '-'  num2str(TFA_FOI_topoplot{i_TFA_FOI_topoplot}(end)) ...
%             'Hz_SurfAvgFreq_' sub '_' ...
%             FuncInput_ToneDur_text 'sTD.png'];
%         figfile      = [path_fig filename];
%         saveas(gcf, [figfile], 'png'); %save png version
%         close
%     end

%     %Select processed elecs and add them to surface plot
%     for i_elec = 1:length(TFA_outputdata_BC_freqband.label)
%         elecs_chanposIndex(1,i_elec) = ...
%             find(strcmp(TFA_outputdata_BC_freqband.label{i_elec}, ...
%             TFA_outputdata_BC_freqband.elec.label));
%     end
%    TFA_outputdata_BC_freqband.elec.chanpos = ...
%         TFA_outputdata_BC_freqband.elec.chanpos(elecs_chanposIndex,:);
%    TFA_outputdata_BC_freqband.elec.chantype = ...
%         TFA_outputdata_BC_freqband.elec.chantype(elecs_chanposIndex);
%    TFA_outputdata_BC_freqband.elec.chanunit = ...
%         TFA_outputdata_BC_freqband.elec.chanunit(elecs_chanposIndex);
%    TFA_outputdata_BC_freqband.elec.elecpos = ...
%         TFA_outputdata_BC_freqband.elec.elecpos(elecs_chanposIndex,:);
%    TFA_outputdata_BC_freqband.elec.label = ...
%         TFA_outputdata_BC_freqband.elec.label(elecs_chanposIndex);
%    TFA_outputdata_BC_freqband.elec.tra = ...
%         TFA_outputdata_BC_freqband.elec.tra(elecs_chanposIndex,elecs_chanposIndex);
% 
%     ft_plot_sens(TFA_outputdata_BC_freqband.elec);
    
end

end