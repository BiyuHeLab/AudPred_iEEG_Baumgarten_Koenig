function NASTD_ECoG_TFA_CorrPredPerfPowChanges...
    (sub, FuncInput_ToneDur_text,  ...
    TFA_FOI, ...
    FuncInput_InputData, ...
    plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG)

%Aim: Correlate prediction performance metric with power changes 
%to see if power changes indicate effective prediction

%% 0.1) Specify vars, paths, and setup fieldtrip
%Add base dir and own script dir
addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));

%Create output folder
path_fig_FTPLrating = ([paths_NASTD_ECoG.Analysis_Behavior '/FTPLratings_p34predp34distance/' ...
    sub '/Figs/']);
if (~exist(path_fig_FTPLrating, 'dir')); mkdir(path_fig_FTPLrating); end

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

%% 2) Compute single-trial prediction performance metric and plot performance summary figure
%Idea: 5 FTPLratings, 5 PE values (i.e., p34-p*34 distance). Code Pitch
%difference by inverted rank (1 = largest, 5 smallest). Use Spearman rank 
%correlation to correlate rank with FTPL rating. The smaller the distance
%between rank & FTPL rating, the better the prediction performance.
%Good performance = FTPL 5 linked to rank 5; Bad performance = FTPL 5
%linked to rank 1;

%For each trial, determine absolute and normal difference between p34 and p*34 and FTPL rating
for i_trial = 1:nTrials    
    p34_logf = FuncInput_InputData.behav.stim.logf_final(i_trial);
    predp34_logf = FuncInput_InputData.behav.stim.logf_pred(i_trial);
    
    absdifflogf_p34predp34(i_trial) = abs(p34_logf - predp34_logf);
    difflogf_p34predp34(i_trial) = (p34_logf - predp34_logf);    
    FTPLrating(i_trial) = FuncInput_InputData.behav.resp_prob(i_trial);
end

%Compute mean FTPLrating for each possible difference
%Determine all unique possibilities for absolute and normal p34-p*34 difference
possible_absdifflogf_p34predp34 = unique(round(absdifflogf_p34predp34,4)); %5 possibilities for abs diff
possible_difflogf_p34predp34 = unique(round(difflogf_p34predp34,4)); %10 possibilities for normal diff

%Compute descriptive statistics over trials for absolute and normal p34-p*34 difference
for i_possiblediff = 1:length(possible_absdifflogf_p34predp34)
    filt = round(absdifflogf_p34predp34,4) == possible_absdifflogf_p34predp34(i_possiblediff);
    avgFTPLrating_perabsdiff(i_possiblediff) = mean(FTPLrating(filt));
    SDFTPLrating_perabsdiff(i_possiblediff) = std(FTPLrating(filt));    
end
for i_possiblediff = 1:length(possible_difflogf_p34predp34)
    filt = round(difflogf_p34predp34,4) == possible_difflogf_p34predp34(i_possiblediff);
    avgFTPLrating_perdiff(i_possiblediff) = mean(FTPLrating(filt));
    SDFTPLrating_perdiff(i_possiblediff) = std(FTPLrating(filt));    
end

%Compute quadratic and lin reg between p34-p*34 difference and FTPL rating to assess subject performance
stat_diff = regstats(FTPLrating, difflogf_p34predp34, 'quadratic', 'tstat');
stat_absdiff = regstats(FTPLrating, absdifflogf_p34predp34, 'linear', 'tstat');

%Code inverted rank for abs p34-p*34 difference
for i_absdiff = 1:length(possible_absdifflogf_p34predp34)
    possible_absdifflogf_p34predp34(2,i_absdiff) = ...
    length(possible_absdifflogf_p34predp34) - find(sort(possible_absdifflogf_p34predp34(1,:)) == ...
    possible_absdifflogf_p34predp34(1,i_absdiff)) + 1;
end

%Compute trial-wise abs difference between FTPLrating and inverted rank
for i_trial = 1:nTrials
    
    InvRank_pertrial(i_trial) = ...
        possible_absdifflogf_p34predp34(2, ...
        find(possible_absdifflogf_p34predp34(1,:) ...
        == round(absdifflogf_p34predp34(i_trial),4)));
        
    PredPerformance_pertrial(i_trial) = ...
            abs(FTPLrating(i_trial) - InvRank_pertrial(i_trial));

end

% %plot summary for quadratic (normal) and linear (absolute) analysis
% if plot_poststepFigs == 1    
%     figure;
%     set(gcf,'units','normalized','outerposition',[0 0 1 1])
%     
%     subplot(1,3,1)
%     bar(possible_difflogf_p34predp34, avgFTPLrating_perdiff)
%     hold on;
%     errorbar(possible_difflogf_p34predp34, avgFTPLrating_perdiff, SDFTPLrating_perdiff, 'LineStyle', 'none','Color', [0 0 0])
%     scatter(difflogf_p34predp34, FTPLrating, 16, ...
%         'filled','ko','jitter', 'on', 'jitterAmount', 0.03)
%     [polynom, ~] = fit(difflogf_p34predp34', FTPLrating', 'poly2');
%     p = plot(polynom, difflogf_p34predp34, FTPLrating);
%     p(2).LineWidth = 3;
%     p(1).Visible = 'off';
%     hFig=findall(0,'type','figure');
%     hLeg=findobj(hFig(1,1),'type','legend');
%     set(hLeg,'visible','off')
%     xlim([-1 1]);
%     ylim([0.5 6]);
%     yticks([1:5])
%     xticks(unique(round(difflogf_p34predp34,2)))
%     xlabel('Pitch Difference p34 - p*34 [Log Freq]')
%     ylabel('FTPL rating')
%     text_linreg = ['p = ' num2str(round(stat_diff.tstat.pval(2),2))];
%     title({['QuadReg: Pitch Difference p34 - p*34 / FTPL rating'] text_linreg})
%     
%     subplot(1,3,2)
%     bar(possible_absdifflogf_p34predp34(1,:), avgFTPLrating_perabsdiff)
%     hold on;
%     errorbar(possible_absdifflogf_p34predp34(1,:), avgFTPLrating_perabsdiff, SDFTPLrating_perabsdiff, 'LineStyle', 'none','Color', [0 0 0])
%     scatter(absdifflogf_p34predp34, FTPLrating, 16, ...
%         'filled','ko','jitter', 'on', 'jitterAmount', 0.02)
%     plot(absdifflogf_p34predp34, stat_absdiff.tstat.beta(1) + (absdifflogf_p34predp34 * stat_absdiff.tstat.beta(2)),'r-','LineWidth', 3)
%     xlim([0 1]);
%     ylim([0.5 6]);
%     yticks([1:5])
%     xticks(unique(round(absdifflogf_p34predp34,2)))
%     xlabel('Absolute Pitch Difference p34 - p*34 [Log Freq]')
%     ylabel('FTPL rating')
%     text_linreg = ['beta = ' num2str(round(stat_absdiff.tstat.beta(2),2)) ', p = ' num2str(round(stat_absdiff.tstat.pval(2),2))];
%     title({['LinReg: Pitch Difference p34 - p*34 / FTPL rating'] text_linreg})
% 
%     subplot(1,3,3)
%     histogram(PredPerformance_pertrial)
%     ylim([0 30])
%     xticks(0:4)
%     xticklabels({'lowest', 'low', 'medium', 'high', 'highest'})
%     xlabel({['Prediction Performance'] ['abs(FTPLrating - Inverted ranked p34 - p*34 difference)']})
%     ylabel('Frequency')
%     title('Histogram: Trial-wise Prediction Performance')    
%     
%     sgtitle([sub ' - ' FuncInput_ToneDur_text 's TD - FTPL rating as function of p34-p*34 distance'])
%     
%     if save_poststepFigs == 1
%         filename     = ['BehavPredPerf_' ...
%             sub '_' ...
%             FuncInput_ToneDur_text 'sTD.png'];
%         figfile      = [path_fig_FTPLrating filename];
%         saveas(gcf, figfile, 'png'); %save png version
%         close
%     end
% end

%% 3) Compute time-frequency representation (TFR)
%Compute TFR
cfg = [];
cfg.output     = 'pow';
cfg.channel    = FuncInput_InputData.label; %All preprocessed elecs
cfg.foi        = TFA_FOI;
cfg.toi        = -0.5 : 0.05 : (TP_Tone_StartStop(end,2) + 0.4); %500ms prestim until 400ms after last tone ended
cfg.keeptrials  = 'yes'; %Keep single trials for later statistical analysis

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

    TFA_outputdata = ft_freqanalysis(cfg, FuncInput_InputData);   
    TFA_type = 'HannFlexTW';
    
elseif TFA_FOI(end) > 40 %For higher freqs use wavelet-based approach
    
    cfg.method     = 'wavelet'; %Morlet wavelets
    cfg.width      = 5; %Width of wavelength in number of cycles

    TFA_outputdata = ft_freqanalysis(cfg, FuncInput_InputData);  
    TFA_type = 'Wavelet';
    
end

%Compute power values relative to baseline 
%(relative change with respect to the power in the baseline interval)
% cfg                 = [];
% cfg.baseline        = [-0.5 0];
% cfg.baselinetype    = 'relative';
% cfg.parameter       = 'powspctrm';
% TFA_outputdata   = ft_freqbaseline(cfg, TFA_outputdata);

%Determine output folder
TFA_label = [TFA_type '_' num2str(TFA_FOI(1)) '-' num2str(TFA_FOI(end)) 'Hz'];

path_fig = ([paths_NASTD_ECoG.ECoGdata_TFA '/Stat_Perf/' ...
    TFA_label '/' sub '/']);
if (~exist(path_fig, 'dir')); mkdir(path_fig); end

%% 4) Statistically test relationship between power (per time-freq tile) and prediciton performance metric 
%Statistical comparison via Spearman rank correlation (cluster-corrected)
for i_elec = 1:nSensors
    
    cfg = [];
    cfg.channel          = TFA_outputdata.label{i_elec};
    cfg.latency          = [0 TP_Tone_StartStop(end,2) + 0.4]; %During sequence + 400 ms after
    cfg.frequency        = 'all';
    
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'ft_statfun_correlationT';
    cfg.type             = 'Spearman';
    
    cfg.tail             = 0; %2-tailed
    cfg.alpha            = 0.025; %i.e., alpha = 0.05 for 2-tailed test   
        
    cfg.correctm         = 'cluster';
    cfg.clustertail      = 0;
    cfg.clusteralpha     = 0.1;
    cfg.clusterstatistic = 'maxsum';

    cfg.numrandomization = 1000;

    cfg.design(1,1:nTrials) = PredPerformance_pertrial; %Prediction performance metric per trial as IV
    cfg.ivar = 1;
    
    TFA_stat{i_elec} = ft_freqstatistics(cfg, TFA_outputdata);    
    
    if any(any(squeeze(TFA_stat{i_elec}.mask)))
        disp(['sign. effect in electrode: ' FuncInput_InputData.label{i_elec}])
    end
    
end

    alpha = cfg.alpha;
    clusteralpha = cfg.clusteralpha;
    
%Combine stat output into a joint struct in FT-format
TFA_statdata = TFA_stat{1};
for i_elec = 2:nSensors
    TFA_statdata.label{1,i_elec} = TFA_stat{i_elec}.label{1};
    TFA_statdata.prob(i_elec,:,:) = TFA_stat{i_elec}.prob;
    TFA_statdata.posclusterslabelmat(i_elec,:,:) = TFA_stat{i_elec}.posclusterslabelmat;
    TFA_statdata.posdistribution(i_elec,:) = TFA_stat{i_elec}.posdistribution;
    TFA_statdata.negclusterslabelmat(i_elec,:,:) = TFA_stat{i_elec}.negclusterslabelmat;
    TFA_statdata.negdistribution(i_elec,:) = TFA_stat{i_elec}.negdistribution;    
    TFA_statdata.cirange(i_elec,:,:) = TFA_stat{i_elec}.cirange;
     TFA_statdata.mask(i_elec,:,:) = TFA_stat{i_elec}.mask;
     TFA_statdata.stat(i_elec,:,:) = TFA_stat{i_elec}.stat;
     TFA_statdata.ref(i_elec,:,:) = TFA_stat{i_elec}.ref;
     TFA_statdata.rho(i_elec,:,:) = TFA_stat{i_elec}.rho;
end
TFA_statdata.label = TFA_statdata.label';

% Plot results (test-statistic with highlighted sign-mask)
clim = [-0.5 0.5]; %Color scale limits

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
               
                imagesc(TFA_stat{i_elec}.time, TFA_stat{i_elec}.freq, ...
                    squeeze(TFA_stat{i_elec}.rho(1,:,:)), ...
                    clim);
                axis xy % flip vertically
                c = colorbar;
                c.Label.String = 'Spearmans Rho';
                
                if any(any(squeeze(TFA_stat{i_elec}.mask)))
                    hold on;
                    contour(TFA_stat{i_elec}.time, TFA_stat{i_elec}.freq, ...
                        squeeze(TFA_stat{i_elec}.mask(1,:,:)),'r-','LineWidth', 1);
                    axis xy % flip vertically
                end
                
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
                
                title(TFA_stat{i_elec}.label)
                
            end
        end
        
        sgtitle({[sub ' - ' FuncInput_ToneDur_text 's TD - ' ...
            num2str(TFA_FOI(1)) '-' num2str(TFA_FOI(end)) 'Hz'] ...
            ['Lin. Correlation (Spearmans Rank) - power (rel. to baseline)' ...
            ' & Prediction Performance Metric' ...
            '(p < ' num2str(alpha) '; clusterp < ' num2str(clusteralpha) ...
            '; page: ' num2str(i_fig) '/' num2str(NumFig) ')']}, ...
            'Interpreter','none')
        
        if save_poststepFigs == 1
            filename     = ['TestAbsPow_StatCompTFAPerf_' ...
                num2str(TFA_FOI(1)) '-'  num2str(TFA_FOI(end)) ...
                'Hz_AllElecs_' sub '_' ...
                FuncInput_ToneDur_text 'sTD_' num2str(i_fig) '.png'];
            figfile      = [path_fig filename];
            saveas(gcf, figfile, 'png'); %save png version
            close    
        end
    end
end

%% Topoplot highlighting sign. effects
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
cfg.channel     = TFA_outputdata.label; %All preprocessed elecs
cfg.mask        = 'convex';
% cfg.boxchannel = {'G*', 'G**'}; %Channels determining channel box size
cfg.boxchannel  = {TFA_outputdata.label{1}, TFA_outputdata.label{2}}; %Channels determining channel box size
cfg.overlap     = 'shift';
lay = ft_prepare_layout(cfg, TFA_outputdata);

%Plot multiplot on 2D surface
cfg             = [];
cfg.layout      = lay;
cfg.showoutline = 'yes';
cfg.parameter   = 'rho';
cfg.zlim        = clim;
cfg.colorbar    = 'yes';
cfg.showlabels  = 'yes';
cfg.showscale   = 'no';
cfg.showcomment = 'no';
cfg.fontsize    = 8;
cfg.box         = 'no';

cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';

ft_multiplotTFR(cfg, TFA_statdata);
set(gcf,'units','normalized','outerposition',[0 0 1 1])

sgtitle({[sub ' - ' FuncInput_ToneDur_text 's TD - ' ...
    num2str(TFA_FOI(1)) '-' num2str(TFA_FOI(end)) 'Hz'] ...
    ['Lin. Correlation (Spearmans Rank) - power (rel. to baseline)' ...
    ' & Prediction Performance Metric' ...
    '(p < ' num2str(alpha) '; clusterp < ' num2str(clusteralpha) ')']}, ...
    'Interpreter','none')
                
if save_poststepFigs == 1
    filename     = ['TestAbsPow_StatCompTFAPerf_' ...
        num2str(TFA_FOI(1)) '-'  num2str(TFA_FOI(end)) ...
        'Hz_MultiSurf_' sub '_' ...
        FuncInput_ToneDur_text 'sTD_' num2str(i_fig) '.png'];  
    figfile      = [path_fig filename];
    saveas(gcf, [figfile], 'png'); %save png version
    close   
end

end