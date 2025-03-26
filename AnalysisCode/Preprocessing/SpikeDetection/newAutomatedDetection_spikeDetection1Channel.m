function [spikeDetectedHigh, localMaximumHigh, localMaximumHighK, ... 
          spikeDetectedLow,  localMaximumLow, localMaximumLowK ] = ...
          newAutomatedDetection_spikeDetection1Channel...
          (filteredData,fs,windowStartIndices,settings)
    
numWindows = length(windowStartIndices);
lengthDataSamples = size(filteredData,1);

%% Step 1: Get Amplitude Envelope
amplitudeEnvelope=abs(hilbert(filteredData)); % Hilbert's envelope (intense envelope)

%% Step 2: estimation of windows distribution using MLE
windowStatistics = zeros(numWindows,2);
for i_window=1:length(windowStartIndices) % for each segment

windowedEnvelope = amplitudeEnvelope(windowStartIndices(i_window): ...
    windowStartIndices(i_window)+settings.detectionBand_winSizeSamples-1);

windowStatistics(i_window,1) = mean(log(windowedEnvelope)); % median
windowStatistics(i_window,2) = std(log(windowedEnvelope));
%
end

r= lengthDataSamples/length(windowStartIndices);
factor = round(settings.detectionBand_winSizeSamples / r);
filtPoints = ones(factor,1)/factor;
windowStatistics(:,1)=filtfilt(filtPoints,1,windowStatistics(:,1));
windowStatistics(:,2)=filtfilt(filtPoints,1,windowStatistics(:,2));

% interpolation of thresholds value to threshold curve (like background)
interpolatedWindowStatistics=[];
if size(windowStatistics,1)>1
    interpolatedWindowStatistics(:,1) = interp1(windowStartIndices+round(settings.detectionBand_winSizeSamples/2),windowStatistics(:,1), ...
        (windowStartIndices(1):windowStartIndices(end))+round(settings.detectionBand_winSizeSamples/2),'spline');
    
    interpolatedWindowStatistics(:,2) = interp1(windowStartIndices+round(settings.detectionBand_winSizeSamples/2),windowStatistics(:,2), ...
        (windowStartIndices(1):windowStartIndices(end))+round(settings.detectionBand_winSizeSamples/2),'spline');
    
    interpolatedWindowStatistics=[ones(floor(settings.detectionBand_winSizeSamples/2),size(windowStatistics,2)).* ...
        repmat(interpolatedWindowStatistics(1,:),floor(settings.detectionBand_winSizeSamples/2),1); ...
        interpolatedWindowStatistics ; ...
        ones(size(amplitudeEnvelope,1)-(length(interpolatedWindowStatistics)+floor(settings.detectionBand_winSizeSamples/2)), ...
        size(windowStatistics,2)) .* ...
        repmat(interpolatedWindowStatistics(end,:),size(amplitudeEnvelope,1)-(length(interpolatedWindowStatistics)+floor(settings.detectionBand_winSizeSamples/2)),1)];
else
    interpolatedWindowStatistics=windowStatistics.*ones(lengthDataSamples,2);
end

lognormal_mode = exp(interpolatedWindowStatistics(:,1)-interpolatedWindowStatistics(:,2).^2);
lognormal_median = exp(interpolatedWindowStatistics(:,1));

%% Step 3: Compare window Statistics to distribution
interpolatedWindowStatistics(:,1)=settings.detectionBand_k1*(lognormal_mode+lognormal_median);
if ~(settings.detectionBand_k2==settings.detectionBand_k1)
    interpolatedWindowStatistics(:,2)=settings.detectionBand_k2*(lognormal_mode+lognormal_median);
end

% -------------------------------------------------------------------------
% detection of obvious and ambiguous spike
% -------------------------------------------------------------------------
[spikeDetectedHigh, localMaximumHigh] = local_maxima_detection(amplitudeEnvelope,interpolatedWindowStatistics(:,1),fs,settings);
if ~(settings.detectionBand_k2==settings.detectionBand_k1)
    [spikeDetectedLow, localMaximumLow]=local_maxima_detection(amplitudeEnvelope,interpolatedWindowStatistics(:,2),fs,settings);
else
    spikeDetectedLow = spikeDetectedHigh;
    localMaximumLow = localMaximumHigh;
end

localMaximumHighK = amplitudeEnvelope(localMaximumHigh==1) ./ (lognormal_mode(localMaximumHigh==1)+lognormal_median(localMaximumHigh==1));
localMaximumLowK = amplitudeEnvelope(localMaximumLow==1) ./ (lognormal_mode(localMaximumLow==1)+lognormal_median(localMaximumLow==1));


    function [spikeDetected, localMaximum] = local_maxima_detection(amplitudeEnvelope,interpolatedWindowStatistics,fs,settings)
        
        
        spikeDetected=zeros(size(amplitudeEnvelope));
        spikeDetected(amplitudeEnvelope(:)>interpolatedWindowStatistics(:))=1; % crossing of high threshold
        
        thresholdCrossings=[];
        thresholdCrossings(:,1)=find(diff([0;spikeDetected])>0); % strat crossing
        thresholdCrossings(:,2)=find(diff([spikeDetected;0])<0); % end crossing
        
        % union of section, where local maxima are close together <(1/f_low + 0.02 sec.)~ 120 ms
        for i_thresholdCrossing=2:size(thresholdCrossings,1)
            endPreviousSpike = thresholdCrossings(i_thresholdCrossing-1,2);
            startThisSpike = thresholdCrossings(i_thresholdCrossing,1);
            if startThisSpike - endPreviousSpike < settings.combineSpikesSeconds*fs
                spikeDetected(endPreviousSpike:startThisSpike) = 1;
            end
        end
        
        localMaximum = spikeDetected;
        % finding of the highes maxima of the section with local maxima
        thresholdCrossings=[];
        thresholdCrossings(:,1)=find(diff([0;spikeDetected])>0); % start
        thresholdCrossings(:,2)=find(diff([spikeDetected;0])<0); % end
        
        % local maxima
        for i_thresholdCrossing=1:size(thresholdCrossings,1)
            if thresholdCrossings(i_thresholdCrossing,2)-thresholdCrossings(i_thresholdCrossing,1)>0
                [~,max_poz]=max(amplitudeEnvelope(thresholdCrossings(i_thresholdCrossing,1):thresholdCrossings(i_thresholdCrossing,2)));
                localMaximum(thresholdCrossings(i_thresholdCrossing,1):thresholdCrossings(i_thresholdCrossing,2))=false;
                localMaximum(thresholdCrossings(i_thresholdCrossing,1)+max_poz-1)=true;
            end
        end
        
    end

end