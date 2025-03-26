function data = NASTD_ECoG_Preproc_RemovePulseArtifact(sub, data, subs_PreProcSettings)

%Determine vars and set up struct and figs
EKGchans = subs_PreProcSettings.(sub).EKG_chanindex;
EKGthresh = subs_PreProcSettings.(sub).EKG_thresh;
close all
heartbeat = struct;
numBlocks = length(data.trial);
if numBlocks > 6
    sbX = 3;
    sbY = 4;
else
    sbX = 2;
    sbY = 3;
end

%Determine filter parameters
fs = data.fsample;
filterOrder = 3;
bigTime = fs*2;
hbeatCutoff = 0.5*fs;
numChans = size(data.trial{1},1);

%% create and plot heartbeat signals + threshold
for i_block = 1:numBlocks
    %Get heartbeat by taking bipolar montage of 2 heartbeat channels
    heartbeat(i_block).Signal = zscore(data.trial{i_block}(EKGchans(1),:) - data.trial{i_block}(EKGchans(2),:));
    heartbeat(i_block).lnSignal = length(heartbeat(i_block).Signal);
    
    %Filter heartbeat signal to remove noise and slow drifts
    highPassFreq = 4 %subjectPreProcessing.EKGFilterHP(i_block); %~4 Hz
    [hpb,hpa] = butter(filterOrder,(highPassFreq*2)/fs,'high');
    heartbeat(i_block).Signal = filtfilt(hpb,hpa,double(heartbeat(i_block).Signal)')';
    lowPassFreq = 20 %subjectPreProcessing.EKGFilterLP(i_block); %~20 Hz
    [lpb,lpa] = butter(filterOrder,(lowPassFreq*2)/fs,'low');
    heartbeat(i_block).Signal = filtfilt(lpb,lpa,double(heartbeat(i_block).Signal)')';
    
    %Detect heartbeats as crossing of picked threshold
    %         thresholdIndex = floor((subjectPreProcessing.EKGThresholds(i_block) / 100) * heartbeat(i_block).lnSignal); %~98
    thresholdIndex = floor((EKGthresh(i_block) / 100) * heartbeat(i_block).lnSignal); %~98
    heartbeat(i_block).threshold = sort(heartbeat(i_block).Signal);
    heartbeat(i_block).threshold = heartbeat(i_block).threshold(thresholdIndex);
    heartbeat(i_block).times = find([0 diff(heartbeat(i_block).Signal>heartbeat(i_block).threshold)==1]);
    
    %remove heartbeats less than 2 seconds from beginning or end of block,
    %and heartbeats that are less than 0.5 seconds apart (likely due to
    %noise in heartbeat signal)
    heartbeat(i_block).times = heartbeat(i_block).times(heartbeat(i_block).times > ...
        bigTime & heartbeat(i_block).times < heartbeat(i_block).lnSignal -(bigTime) - 1);
    hTimes = ones(length(heartbeat(i_block).times),1);
    tooFastHeartbeats = find(diff(heartbeat(i_block).times) < hbeatCutoff);
    hTimes(tooFastHeartbeats) = 0;
    heartbeat(i_block).times  = heartbeat(i_block).times(find(hTimes));
    
    %Fig1: Plot Zscored EKG signal on cumulative distribution function (CDF)
    %Shows which Zscored EKG signal values lie above the threshold
    figure(1)
    subplot(sbX,sbY,i_block);
    plot(sort(heartbeat(i_block).Signal),100 * (1/heartbeat(i_block).lnSignal:1/heartbeat(i_block).lnSignal:1));
    %         refline(0,subjectPreProcessing.EKGThresholds(i_block));
    refline(0,98);
    xlabel('zscore(EKG1 - EKG2)');
    ylabel('CDF (%)')
    title(['Block ' num2str(i_block)]);
    axis([-5 5 0 100]);
    
    %Fig2: Plot Zscored bipolar EKG signal across time along with
    %threshold. Shows if threshold is suitable to catch all heartbeats.
    figure(2)
    subplot(numBlocks,1,i_block);
    plot(1/fs:1/fs:heartbeat(i_block).lnSignal/fs,heartbeat(i_block).Signal);
    
    refline(0,heartbeat(i_block).threshold);
    set(gca,'ylim',[-2*heartbeat(i_block).threshold 2*heartbeat(i_block).threshold]);
    title(['Block ' num2str(i_block)]);
    
    %Fig3: Plot histogram of heartbeat intervals.
    %Shows if heartbeat intervals lie within a physiologically sane
    %range (~1s)
    figure(3)
    subplot(sbX,sbY,i_block);
    [y,x] = hist(diff(heartbeat(i_block).times)/fs,[0.4:0.05:1.6]);
    bar(x,log10(0.1 + (100 * (y/sum(y)))));
    xlabel('HeartbeatInterval (seconds)')
    ylabel('log10(heartbeats %)');
    title(['Block ' num2str(i_block)]);
end

%Rescale and save figures
figure(1);
set(gcf,'position',[45         100        1721         716])
sgtitle('threshold percentile')
print('-dpng',['Heartbeat CDF']);

figure(2);
xlabel('Time (seconds)')
set(gcf,'position',[14        0        2220         886])
print('-dpng',['heartbeat_Threshold.png']);

figure(3);
set(gcf,'position',[  138         100        2008         569])
title('heartbeat interval distributions')
print('-dpng',['heartbeat_Intervals.png']);

%%

for i_block = 1:numBlocks
    
    %Get half of average heartbeatInterval in samples
    halfTime = floor(median(diff(heartbeat(i_block).times))/2);
    
    %Get all heartbeat-evoked trials
    trial = zeros(numChans,length(heartbeat(i_block).times),1 + (bigTime*2));
    for i_chan = 1:numChans
        for i_beat = 1:length(heartbeat(i_block).times)
            beatIndices = heartbeat(i_block).times(i_beat) - (bigTime) : heartbeat(i_block).times(i_beat) + (bigTime);
            trial(i_chan,i_beat,:) = data.trial{i_block}(i_chan,beatIndices);
        end
    end
    
    %Average over heartbeats to get a sensor*time matrix
    averageWaveForm = squeeze(mean(trial,2));
    
    %LP filter average waveform
    lowPassFreq = 5;
    [lpb,lpa] = butter(filterOrder,(lowPassFreq*2)/fs,'low');
    filteredAverageWaveForm = filtfilt(lpb,lpa,averageWaveForm')';
    filteredAverageWaveForm = filteredAverageWaveForm(:,bigTime + 1 - halfTime:bigTime + 1 + halfTime);
    
    %Apply tapered window (with 10 % cutoff
    taper = tukeywin(1 + (halfTime*2),0.1);
    filteredAverageWaveForm = filteredAverageWaveForm .* repmat(taper',[numChans 1]);
    
    %% remove average waveform from each detected heartbeat
    data.trial_raw{i_block} = data.trial{i_block};
%     data.trial_heartbeat{i_block} = data.trial{i_block};
    
    for i_beat = 1:length(heartbeat(i_block).times)
        beatIndices = heartbeat(i_block).times(i_beat) - (halfTime):heartbeat(i_block).times(i_beat) + (halfTime);
        data.trial{i_block}(:,beatIndices) =  data.trial_raw{i_block}(:,beatIndices) - filteredAverageWaveForm;
    end
    
    data.trial_raw{i_block} = [];
    
end

data = rmfield(data, 'trial_raw');
close all

end