function InputData = ...
    NASTD_ECoG_Preproc_AmpEnvelope...
    (sub, FilterOrder, highPassFreq, lowPassFreq, Freq_Label, InputData, ...
    paths_NASTD_ECoG, plot_poststepFigs, save_poststepFigs)

%Note: Old filter function without any NaN-interpolation - use only with no
%NaN data!

%Aim: For each trial, filter signal to get frequency-band specific activity,
%then compute Hilbert-transform to get the amplitude envelope.
%Save amplitude value normalized by average amplitude across block and
%log-amplitude normalized by average log-amplitude across block.

%For gamma-band: Smooth activity

%1. Determine parameters and set up struct
SampleFreq     = InputData.fsample;

BPsignalFilt        = [];
BPsignalAmp         = [];
BPsignalAmpSmooth   = [];

FieldLabel = [Freq_Label '_LogAmp'];
% FieldLabel = [Freq_Label '_NormLogAmp'];
InputData.(FieldLabel) = [];

%2. Band-pass filter
for i_trial = 1:length(InputData.trial)
    
    %2.1 Low pass filter
    [lpb,lpa] = butter(FilterOrder,(lowPassFreq*2)/SampleFreq,'low');
    BPsignalFilt{i_trial} = filtfilt(lpb,lpa,InputData.trial{i_trial}')';
    
    %2.2 High pass filter
    [lpb,lpa] = butter(FilterOrder,(highPassFreq*2)/SampleFreq,'high');
    BPsignalFilt{i_trial} = filtfilt(lpb,lpa,BPsignalFilt{i_trial}')';
    
    %2.3 Hilbert transform and extraction of amplitude envelope (absolute part of Hilbert-transform)
    for i_chan = 1:length(InputData.label) %Chan loop, otherwise Hilbert function replicates one channel
        BPsignalAmp{i_trial}(i_chan,:) = abs(hilbert(BPsignalFilt{i_trial}(i_chan,:)));
    end
    %Note: Take absolute hilbert as amplitude envelope (Hilbert transform
    %basically gives you an analytical signal, allowing to determine
    %the amplitude and phase information at a certain point of the signal)
    
    %2.4 Smooth (moving average) amplitude envelope by 100ms for high freqs
    if strcmp(Freq_Label, 'Gamma') || strcmp(Freq_Label, 'HighGamma') 
        samples_100ms = round(InputData.fsample/10);
        for i_chan = 1:length(InputData.label)
            %Chan loop, otherwise Hilbert function replicates one channel
            BPsignalAmp{i_trial}(i_chan,:) = ...
                smooth(BPsignalAmp{i_trial}(i_chan,:), samples_100ms, 'moving')'; 
        end
    end
    
    %2.5 Compute Amplitue Signal normalized by across-trial signal (non-log)
    Temp_NormAmp{i_trial} = ...
        BPsignalAmp{i_trial} - ...
        repmat(mean(BPsignalAmp{i_trial},2),...
        [1 size(BPsignalAmp{i_trial},2)]);
    
    %2.6 Compute log10 amplitude per channel for each TP across the trial.
    InputData.(FieldLabel){i_trial}  = ...
        log10(BPsignalAmp{i_trial});
    
    %2.7 Compute log10 amplitude per channel for each TP across the trial.
    %Normalized by taking the log(base10) of the amplitude envelope and
    %subtracting the mean of the log10 amplitude across the entire trial.
    Temp_LogNormAmp{i_trial}  = ...
        log10(BPsignalAmp{i_trial}) ...
        - repmat(mean(log10(BPsignalAmp{i_trial}),2),...
        [1 size(BPsignalAmp{i_trial},2)]);
    
end

% %3. Visually check output
% figure;
% i_chan = randi(length(InputData.label),1);
% 
% % %3.1 Check if filter worked
% subplot(1,2,1);
% pwelch(BPsignalFilt{i_trial}(i_chan,:),[],[],[],512); %filtered data
% xlim([0 150])
% 
% %3.2 Plot all and compare
% subplot(1,2,2);
% hold on;
% TW = 1:length(BPsignalFilt{i_trial}); %Time window in samples
% TW = 1:1000; %Time window in samples
% 
% plot(1:length(BPsignalFilt{i_trial}(i_chan,TW)),...
%     BPsignalFilt{i_trial}(i_chan,TW),'k-');
% plot(1:length(BPsignalAmp{i_trial}(i_chan,TW)),...
%     BPsignalAmp{i_trial}(i_chan,TW),'b-', 'LineWidth', 2);
% plot(1:length(Temp_NormAmp{i_trial}(i_chan,TW)),...
%     Temp_NormAmp{i_trial}(i_chan,TW),'m-', 'LineWidth', 2);
% plot(1:length(InputData.(FieldLabel){i_trial}(i_chan,TW)),...
%     InputData.(FieldLabel){i_trial}(i_chan,TW),'y-', 'LineWidth', 3);
% plot(1:length(Temp_LogNormAmp{i_trial}(i_chan,TW)),...
%     Temp_LogNormAmp{i_trial}(i_chan,TW),'c-', 'LineWidth', 2);
% legend('FilteredSignal','Amp','NormAmp','LogAmp','LogNormAmp');
% xlabel('Samples')
% 
% sgtitle([sub ' - ' Freq_Label ' - Elec: ' InputData.label{i_chan} ' Trial: ' num2str(i_trial) ...
%   ' - Filter summary']);

% %3.3 Plot filtered data power spectrum for each channel
% if plot_poststepFigs == 1
%     tic
%     for i_chan = 1:length(InputData.label)
%         figure
%         set(gcf,'units','normalized','outerposition',[0 0 1 1]) %fullscreen
%         for i_trial = 1:length(InputData.sampleinfo)
%             
%             subplot(10,6,i_trial);
%             
%             pwelch(BPsignalFilt{i_trial}(i_chan,:),[],[],[],512); %filtered data
%             subFigtitle = strcat({['Trial:' num2str(i_trial)]});
%             title(subFigtitle)
%             %                 xlim([highPassFreq lowPassFreq]);
%             xlim([0 lowPassFreq*1.5]);
%             ylim([-80 40]);
%             ylabel('')
%             yticks([-80 -40 0 40])
%         end
%         
%         Figtitle = [sub ' Power-Spectrum for BP-filtered data (' num2str(highPassFreq) '-' num2str(lowPassFreq) 'Hz) for Elec:' ...
%             InputData.label(i_chan) '(' num2str(i_chan) ')'];
%         sgtitle(strcat(Figtitle{:}))
%         
%         if save_poststepFigs == 1
%             path_fig = [paths_NASTD_ECoG.Preproc_ECoGdata sub...
%                 '/Figs/Preproc/BPfilt_' num2str(highPassFreq)...
%                 '-' num2str(lowPassFreq) 'Hz/' ];
%             mkdir(path_fig);
%             
%             filename     = strcat([sub '_PSD_BP' num2str(highPassFreq) ...
%                 '-' num2str(lowPassFreq) 'Hz_AllTrials_Chan' ...
%                 InputData.label{i_chan} '(' ...
%                 num2str(i_chan) ').png']);
%             figfile      = [path_fig filename];
%             
%             saveas(gcf, [figfile], 'png'); %save png version
%             close
%         end
%     end
%     toc
% end

% %3.4 Plot normalized and log-normalized Gamma-Amplitude per channel
% if plot_poststepFigs == 1
%     tic
%     for i_chan = 1:length(InputData.label)
%         figure
%         set(gcf,'units','normalized','outerposition',[0 0 1 1]) %fullscreen
%         for i_trial = 1:length(InputData.sampleinfo)
%             
%             subplot(10,6,i_trial);
%             hold on;
%             plot(1:length(Temp_NormAmp{i_trial}(i_chan,:)), ...
%                 Temp_NormAmp{i_trial}(i_chan,:));
%             plot(1:length(InputData.(FieldLabel){i_trial}(i_chan,:)), ...
%                 InputData.(FieldLabel){i_trial}(i_chan,:));
%             
%             subFigtitle = strcat({['Trial:' num2str(i_trial)]});
%             title(subFigtitle)
%             %                 xlim([highPassFreq lowPassFreq]);
%             xlim([1 length(InputData.(FieldLabel){i_trial})]);
%             ylim([-20 20]);
%             ylabel('')
%             yticks([-20 0 20])
%         end
%         
%         Figtitle = [sub ' Normalized (non-log & log) Amplitude (' num2str(highPassFreq) '-' num2str(lowPassFreq) 'Hz) for Elec:' ...
%             InputData.label(i_chan) '(' num2str(i_chan) ')'];
%         sgtitle(strcat(Figtitle{:}))
%         
%         if save_poststepFigs == 1
%             path_fig = [paths_NASTD_ECoG.Preproc_ECoGdata sub...
%                 '/Figs/Preproc/Amp_' num2str(highPassFreq)...
%                 '-' num2str(lowPassFreq) 'Hz/' ];
%             mkdir(path_fig);
%             
%             filename     = strcat([sub '_Amp' num2str(highPassFreq) ...
%                 '-' num2str(lowPassFreq) 'Hz_AllTrials_Chan' ...
%                 InputData.label{i_chan} '(' ...
%                 num2str(i_chan) ').png']);
%             figfile      = [path_fig filename];
%             
%             saveas(gcf, [figfile], 'png'); %save png version
%         end
%         close
%     end
%     toc
% end

end