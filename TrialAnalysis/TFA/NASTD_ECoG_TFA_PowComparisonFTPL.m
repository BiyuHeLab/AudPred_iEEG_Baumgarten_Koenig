function NASTD_ECoG_TFA_PowComparisonFTPL...
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

%% 2) Differentiate trials based on the FTPL rating
%For each trial, determine absolute and normal difference between p34 and p*34 and FTPL rating
for i_trial = 1:nTrials    
    p34_logf = FuncInput_InputData.behav.stim.logf_final(i_trial);
    predp34_logf = FuncInput_InputData.behav.stim.logf_pred(i_trial);
    
    absdifflogf_p34predp34(i_trial) = abs(p34_logf - predp34_logf);
    difflogf_p34predp34(i_trial) = (p34_logf - predp34_logf);    
    FTPLrating(i_trial) = FuncInput_InputData.behav.resp_prob(i_trial);
end

%Compute mean FTPLrating for each possible difference
possible_absdifflogf_p34predp34 = unique(round(absdifflogf_p34predp34,4));
for i_possiblediff = 1:length(possible_absdifflogf_p34predp34)
    filt = round(absdifflogf_p34predp34,4) == possible_absdifflogf_p34predp34(i_possiblediff);
    avgFTPLrating_perabsdiff(i_possiblediff) = mean(FTPLrating(filt));
    SDFTPLrating_perabsdiff(i_possiblediff) = std(FTPLrating(filt));    
end

possible_difflogf_p34predp34 = unique(round(difflogf_p34predp34,4));
for i_possiblediff = 1:length(possible_difflogf_p34predp34)
    filt = round(difflogf_p34predp34,4) == possible_difflogf_p34predp34(i_possiblediff);
    avgFTPLrating_perdiff(i_possiblediff) = mean(FTPLrating(filt));
    SDFTPLrating_perdiff(i_possiblediff) = std(FTPLrating(filt));    
end

%Compute quadratic and lin reg between p34-p*34 difference and FTPL rating to assess subject performance
stat_diff = regstats(FTPLrating, difflogf_p34predp34, 'quadratic', 'tstat');
stat_absdiff = regstats(FTPLrating, absdifflogf_p34predp34, 'linear', 'tstat');

%plot summary for quadratic (normal) and linear (absolute) analysis
if plot_poststepFigs == 1    
    figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    
    subplot(1,2,1)
    bar(possible_difflogf_p34predp34, avgFTPLrating_perdiff)
    hold on;
    errorbar(possible_difflogf_p34predp34, avgFTPLrating_perdiff, SDFTPLrating_perdiff, 'LineStyle', 'none','Color', [0 0 0])
    scatter(difflogf_p34predp34, FTPLrating, 16, ...
        'filled','ko','jitter', 'on', 'jitterAmount', 0.03)
    [polynom, gof] = fit(difflogf_p34predp34', FTPLrating', 'poly2');
    p = plot(polynom, difflogf_p34predp34, FTPLrating);
    p(2).LineWidth = 3;
    p(1).Visible = 'off';
    hFig=findall(0,'type','figure');
    hLeg=findobj(hFig(1,1),'type','legend');
    set(hLeg,'visible','off')
    xlim([-1 1]);
    ylim([0.5 6]);
    yticks([1:5])
    xticks(unique(round(difflogf_p34predp34,2)))
    xlabel('Pitch Difference p34 - p*34 [Log Freq]')
    ylabel('FTPL rating')
    text_linreg = ['p = ' num2str(round(stat_diff.tstat.pval(2),2))];
    title({['QuadReg: Pitch Difference p34 - p*34 / FTPL rating'] text_linreg})
    
    subplot(1,2,2)
    bar(possible_absdifflogf_p34predp34, avgFTPLrating_perabsdiff)
    hold on;
    errorbar(possible_absdifflogf_p34predp34, avgFTPLrating_perabsdiff, SDFTPLrating_perabsdiff, 'LineStyle', 'none','Color', [0 0 0])
    scatter(absdifflogf_p34predp34, FTPLrating, 16, ...
        'filled','ko','jitter', 'on', 'jitterAmount', 0.02)
    plot(absdifflogf_p34predp34, stat_absdiff.tstat.beta(1) + (absdifflogf_p34predp34 * stat_absdiff.tstat.beta(2)),'r-','LineWidth', 3)
    xlim([0 1]);
    ylim([0.5 6]);
    yticks([1:5])
    xticks(unique(round(absdifflogf_p34predp34,2)))
    xlabel('Absolute Pitch Difference p34 - p*34 [Log Freq]')
    ylabel('FTPL rating')
    text_linreg = ['beta = ' num2str(round(stat_absdiff.tstat.beta(2),2)) ', p = ' num2str(round(stat_absdiff.tstat.pval(2),2))];
    title({['LinReg: Pitch Difference p34 - p*34 / FTPL rating'] text_linreg})
    
    sgtitle([sub ' - ' FuncInput_ToneDur_text 's TD - FTPL rating as function of p34-p*34 distance'])
    
    if save_poststepFigs == 1
        filename     = ['FTPLrating_p34predp34distance_' ...
            sub '_' ...
            FuncInput_ToneDur_text 'sTD.png'];
        figfile      = [path_fig_FTPLrating filename];
        saveas(gcf, figfile, 'png'); %save png version
        close
    end
end

%Differentiate trials with 'correct' vs. 'incorrect' response 
%Correct trials: low diff + high FTPL; high diff + low FTPL
%Incorrect trials: high diff + high FTPL; low diff + low FTPL

filt_lowabsdiff1 = round(absdifflogf_p34predp34,4) == possible_absdifflogf_p34predp34(1);
filt_lowabsdiff2 = round(absdifflogf_p34predp34,4) == possible_absdifflogf_p34predp34(2);
filt_lowabsdiff = or(filt_lowabsdiff1, filt_lowabsdiff2);
ind_lowabsdiff = find(filt_lowabsdiff);
ind_lowabsdiff_lowFTPLrating = ind_lowabsdiff(FTPLrating(filt_lowabsdiff) < median(FTPLrating(filt_lowabsdiff)));
ind_lowabsdiff_highFTPLrating = ind_lowabsdiff(FTPLrating(filt_lowabsdiff) > median(FTPLrating(filt_lowabsdiff)));

filt_highabsdiff1 = round(absdifflogf_p34predp34,4) == possible_absdifflogf_p34predp34(end);
filt_highabsdiff2 = round(absdifflogf_p34predp34,4) == possible_absdifflogf_p34predp34(end-1);
filt_highabsdiff = or(filt_highabsdiff1, filt_highabsdiff2);
ind_highabsdiff = find(filt_highabsdiff);
ind_highabsdiff_lowFTPLrating = ind_highabsdiff(FTPLrating(filt_highabsdiff) < median(FTPLrating(filt_highabsdiff)));
ind_highabsdiff_highFTPLrating = ind_highabsdiff(FTPLrating(filt_highabsdiff) > median(FTPLrating(filt_highabsdiff)));

filt_medabsdiff = round(absdifflogf_p34predp34,4) == possible_absdifflogf_p34predp34(3);
ind_medabsdiff = find(filt_medabsdiff);
ind_medabsdiff_medFTPLrating = ind_medabsdiff(FTPLrating(filt_medabsdiff) == 3);
ind_medabsdiff_extremeFTPLrating = ind_medabsdiff(FTPLrating(filt_medabsdiff) ~= 3);

ind_correctpred = sort([ind_lowabsdiff_highFTPLrating ind_highabsdiff_lowFTPLrating ind_medabsdiff_medFTPLrating]);
ind_incorrectpred = sort([ind_lowabsdiff_lowFTPLrating ind_highabsdiff_highFTPLrating ind_medabsdiff_extremeFTPLrating]);

disp(['Number trials with rather correct FTPL ratings: ' num2str(length(ind_correctpred))])
disp(['Number trials with rather wrong FTPL ratings: ' num2str(length(ind_incorrectpred))])

%% 3) Compute TFA for all conditions 
%Seperate trials and ensure all trials are of equal length
cfg = [];
cfg.toilim = [-1 TP_Tone_StartStop(end,2) + 1];
cfg.trials = ind_correctpred;
InputData_correcttrials = ft_redefinetrial(cfg, FuncInput_InputData);
cfg.trials = ind_incorrectpred;
InputData_incorrecttrials = ft_redefinetrial(cfg, FuncInput_InputData);

%Compute TFA
cfg = [];
cfg.output     = 'pow';
cfg.channel    = FuncInput_InputData.label; %All preprocessed elecs
cfg.foi        = TFA_FOI;
cfg.toi        = -0.5 : 0.05 : (TP_Tone_StartStop(end,2) + 0.4); %500ms prestim until 400ms after last tone ended
cfg.keeptrials  = 'yes';

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

    TFA_outputdata_correcttrials = ft_freqanalysis(cfg, InputData_correcttrials);
    TFA_outputdata_incorrecttrials = ft_freqanalysis(cfg, InputData_incorrecttrials);
    
    TFA_type = 'HannFlexTW';
    
elseif TFA_FOI(end) > 40 %For higher freqs use wavelet-based approach
    
    cfg.method     = 'wavelet'; %Morlet wavelets
    cfg.width      = 5; %Width of wavelength in number of cycles

    TFA_outputdata_correcttrials = ft_freqanalysis(cfg, InputData_correcttrials);
    TFA_outputdata_incorrecttrials = ft_freqanalysis(cfg, InputData_incorrecttrials);
    
    TFA_type = 'Wavelet';
    
end

%Compute power values relative to baseline 
%(relative change with respect to the power in the baseline interval)
cfg                 = [];
cfg.baseline        = [-0.5 0];
cfg.baselinetype    = 'relative';
cfg.parameter       = 'powspctrm';
TFA_outputdata_correcttrials_BC     = ft_freqbaseline(cfg, TFA_outputdata_correcttrials);
TFA_outputdata_incorrecttrials_BC   = ft_freqbaseline(cfg, TFA_outputdata_incorrecttrials);

%Determine output folder
TFA_label = [TFA_type '_' num2str(TFA_FOI(1)) '-' num2str(TFA_FOI(end)) 'Hz'];

path_fig_powcomp = ([paths_NASTD_ECoG.ECoGdata_TFA '/CompFTPLratings/' ...
    TFA_label '/' sub '/Figs/']);
if (~exist(path_fig_powcomp, 'dir')); mkdir(path_fig_powcomp); end

%% 4) Contrast TFAs between correct vs. incorrect trials
%Statistical comparison via t-Test (cluster-corrected)
for i_elec = 1:nSensors
    
    cfg = [];
    cfg.channel          = FuncInput_InputData.label{i_elec};
    cfg.latency          = 'all';
    cfg.frequency        = 'all';
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'ft_statfun_indepsamplesT';
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.tail             = 0;
    cfg.clustertail      = 0;
    cfg.alpha            = 0.025;
    cfg.numrandomization = 100;

    design = zeros(1,size(TFA_outputdata_correcttrials_BC.powspctrm,1) + size(TFA_outputdata_incorrecttrials_BC.powspctrm,1));
    design(1,1:size(TFA_outputdata_correcttrials_BC.powspctrm,1)) = 1;
    design(1,(size(TFA_outputdata_correcttrials_BC.powspctrm,1)+1):(size(TFA_outputdata_correcttrials_BC.powspctrm,1)+...
    size(TFA_outputdata_incorrecttrials_BC.powspctrm,1))) = 2;

    cfg.design           = design;
    cfg.ivar             = 1;

    TFA_stat{i_elec} = ft_freqstatistics(cfg, TFA_outputdata_correcttrials_BC, TFA_outputdata_incorrecttrials_BC);    
    
    if any(any(squeeze(TFA_stat{i_elec}.mask)))
        disp(['sign. effect in electrode: ' FuncInput_InputData.label{i_elec}])
    end
end

% Plot results (test-statistic with highlighted sign-mask)
clim = [-4 4]; %Color scale limits

if plot_poststepFigs == 1
    NumFig = ceil(nSensors/30); %30 elecs per plot
    
    i_elec = 0;
    for i_fig = 1 : NumFig
        h = figure;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        
        for i_subplot = 1:30
            if i_elec < nSensors
                
                subplot(5,6,i_subplot)
                i_elec = i_elec+1;
               
                imagesc(TFA_stat{i_elec}.time, TFA_stat{i_elec}.freq, ...
                    squeeze(TFA_stat{i_elec}.stat(1,:,:)), ...
                    clim);
                axis xy % flip vertically
                colorbar;
                
                if any(any(squeeze(TFA_stat{i_elec}.mask)))
                    hold on;
                    contour(TFA_stat{i_elec}.time, TFA_stat{i_elec}.freq, ...
                        squeeze(TFA_stat{i_elec}.mask(1,:,:)),'r-','LineWidth', 3);
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
                
                title(['Elec: ' FuncInput_InputData.label{i_elec}])
                
            end
        end
        
        sgtitle({[sub ' - ' FuncInput_ToneDur_text 's TD - ' ...
            num2str(TFA_FOI(1)) '-' num2str(TFA_FOI(end)) 'Hz'] ...
            ['Stat. comparison (indep. t-test) of power relative to baseline ' ...
            'between trials with correct (n = ' num2str(size(TFA_outputdata_correcttrials_BC.powspctrm,1)) ') ' ...
            'vs. incorrect (n = ' num2str(size(TFA_outputdata_incorrecttrials_BC.powspctrm,1)) ') FTPLrating (' ...
            num2str(i_fig) '/' num2str(NumFig) ')']},'Interpreter','none')
        
        if save_poststepFigs == 1
            filename     = ['StatCompRelPow_AllElecs_' sub '_' ...
                FuncInput_ToneDur_text 'sTD_' num2str(i_fig) '.png'];
            figfile      = [path_fig_powcomp filename];
            saveas(gcf, figfile, 'png'); %save png version
            close
        end
    end
end


%% To do: add topoplot highlighting sign. effects


end