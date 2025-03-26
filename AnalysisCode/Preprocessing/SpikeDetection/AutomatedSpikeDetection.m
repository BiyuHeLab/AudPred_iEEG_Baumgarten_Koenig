function [spikeDetectedHigh, localMaximumHigh, spikeDetectedHighK, ...
          spikeDetectedLow,  localMaximumLow, spikeDetectedLowK, ...
          MuscleArtifact] ...
          = AutomatedSpikeDetection(data,fs,settings)
      
%Aim: Detect interictal epileptic discharges (IEDs) in the signal. Provide
%a descriptive summary about these IEDs, mark trials, and replace 
%contaminated sections with NaNs.
%Based on script by Richard Hardstone and on Janca et al. 2015       

% %% Start parallel pool if one not already running
% if license('test','distrib_computing_toolbox')==1
%     myCluster = parcluster('local');
%     numWorkers = myCluster.NumWorkers;
%     
%     poolobj = gcp('nocreate'); % If no pool, do not create new one.
%     if isempty(poolobj)
%         poolobj = parpool('local',numWorkers-1);
%     end
% end


%% Step 0: Get and display size of data
numberChannels = size(data,2);
lengthDataSamples = size(data,1);
lengthSeconds = lengthDataSamples / fs;

disp(['Data has ' num2str(numberChannels) ' Channels and is ' num2str(lengthSeconds,2) ' seconds long']);

%% Step 1: Resample to 200 Hz

if fs>settings.resamplingFrequency
    resampledData=resample(data(:,1),settings.resamplingFrequency,fs,100);
    resampledData = zeros(size(resampledData,1),numberChannels);
    
    if ~license('test','distrib_computing_toolbox')==1
        for i_chan=1:numberChannels
            resampledData(:,i_chan)=resample(data(:,i_chan),settings.resamplingFrequency,fs,100);
        end
    else
        parfor i_chan=1:numberChannels
            resampledData(:,i_chan)=resample(data(:,i_chan),settings.resamplingFrequency,fs,100);
        end
    end
    data = resampledData;
    clear resampledData;
    originalFs = fs;
    fs=settings.resamplingFrequency;
    lengthDataSamples = size(data,1);
end
settings.detectionBand_winSizeSamples = settings.detectionBand_winSizeSeconds*fs;
settings.detectionBand_winOverlapSamples = settings.detectionBand_winOverlapSeconds*fs;

%% Step 2: split data into segments if necessary
indexStart = 1;
indexStop = lengthDataSamples;

%% Step 3: Apply SpikeDetector to each segment

% Get starting sample of each detection window
windowStartIndices = 1:settings.detectionBand_winSizeSamples - settings.detectionBand_winOverlapSamples: ...
    size(data,1)-settings.detectionBand_winSizeSamples+1;

% Filter data (10-60Hz, zero-phase) into detection window and remove line noise

filteredData = newAutomatedDetection_Filtering([settings.detectionBand_highPass settings.detectionBand_lowPass],data,fs,settings.filterType);
filteredData = newAutomatedDetection_LineNoiseFiltering(filteredData,fs,settings.lineNoiseFrequency);

localMaximumHigh = zeros(size(data));
localMaximumLow = zeros(size(data));
spikeDetectedHigh = zeros(size(data));
spikeDetectedLow = zeros(size(data));
localMaximumHighK = struct;
localMaximumLowK = struct;

%if ~license('test','distrib_computing_toolbox')==1
for i_chan = 1:numberChannels
%     i_chan
    [spikeDetectedHigh(:,i_chan), localMaximumHigh(:,i_chan), localMaximumHighK(i_chan).ks, ...
        spikeDetectedLow(:,i_chan),  localMaximumLow(:,i_chan), localMaximumLowK(i_chan).ks,] = ...
        newAutomatedDetection_spikeDetection1Channel(filteredData(:,i_chan),fs,windowStartIndices,settings);
end
%else
%   parfor i_chan = 1:numberChannels
%i_chan
%      [spikeDetectedHigh(:,i_chan), localMaximumHigh(:,i_chan), ...
%         spikeDetectedLow(:,i_chan),  localMaximumLow(:,i_chan)] = newAutomatedDetection_spikeDetection1Channel(filteredData(:,i_chan),fs,windowStartIndices,settings);
% end
%end

%% Step 4: Detect Muscle artifacts and remove Segments

% beta activity detection
if settings.muscleArtifactBand_frequency<fs/2 && settings.muscleArtifactBand_windowSizeSeconds>0
    MuscleArtifact=newAutomatedDetection_DetectBeta(data,fs,settings);
    %spikeDetectedHigh(M_beta)=false;
    %localMaximumHigh(M_beta)=false;
    %spikeDetectedLow(M_beta)=false;
    %localMaximumLow(M_beta)=false;
end


if settings.outputAtOriginalSamplingRate == true
    resampledSpikeDetectedHigh = resample(spikeDetectedHigh(:,1),originalFs,fs,100);
    resampledSpikeDetectedHigh = zeros(size(resampledSpikeDetectedHigh,1),numberChannels);
    resampledMuscleArtifact = zeros(size(resampledSpikeDetectedHigh,1),numberChannels);
    
    if ~license('test','distrib_computing_toolbox')==1
        for i_chan=1:numberChannels
            resampledSpikeDetectedHigh(:,i_chan)=resample(spikeDetectedHigh(:,i_chan),originalFs,fs,100);
            resampledMuscleArtifact(:,i_chan)=resample(MuscleArtifact(:,i_chan),originalFs,fs,100);
        end
    else
        parfor i_chan=1:numberChannels
            resampledSpikeDetectedHigh(:,i_chan)=resample(spikeDetectedHigh(:,i_chan),originalFs,fs,100);
            resampledMuscleArtifact(:,i_chan)=resample(double(MuscleArtifact(:,i_chan)),originalFs,fs,100);
        end
    end
    spikeDetectedHigh = double(resampledSpikeDetectedHigh>0.5);
    MuscleArtifact = resampledMuscleArtifact>0.5;
end

counter = 1;

for i_chan = 1:numberChannels
%     i_chan
    starts = find(diff(double([0 ; spikeDetectedHigh(:,i_chan)]))==1);
    stops = find(diff(double([spikeDetectedHigh(:,i_chan) ; 0]))==-1);
    if ~isempty(starts)
        for i_spike = 1:length(starts)
            spikeDetectedHigh(starts(i_spike):stops(i_spike),i_chan) = counter;
            spikeDetectedHighK(counter) = localMaximumHighK(i_chan).ks(i_spike);
            counter = counter + 1;
        end
    else
        spikeDetectedHighK = NaN;
    end
end

counter = 1;

for i_chan = 1:numberChannels
%     i_chan
    starts = find(diff(double([0 ; spikeDetectedLow(:,i_chan)]))==1);
    stops = find(diff(double([spikeDetectedLow(:,i_chan) ; 0]))==-1);
    if ~isempty(starts)
        for i_spike = 1:length(starts)
            spikeDetectedLow(starts(i_spike):stops(i_spike),i_chan) = counter;
            spikeDetectedLowK(counter) = localMaximumLowK(i_chan).ks(i_spike);
            counter = counter + 1;
        end
    else
       spikeDetectedLowK  = NaN;
    end        
end


% if ~license('test','distrib_computing_toolbox')==1
%     delete(poolobj);
% end