classdef SpikeDetectionSettings
    properties
        detectionBand_highPass = 10;
        detectionBand_lowPass = 60;
        detectionBand_k1 = 3;
        detectionBand_k2 = [];
        detectionBand_winSizeSeconds = 5;
        detectionBand_winSizeSamples double
        detectionBand_winOverlapSeconds = 4;
        detectionBand_winOverlapSamples double
        bufferingSeconds = 300;
        lineNoiseFrequency = 60;
        muscleArtifactBand_frequency = 13;
        muscleArtifactBand_windowSizeSeconds = 5;
        muscleArtifactBand_autoregressiveSamples = 12;
        filterType = 1;
        discharge_tol = 0.005;
        combineSpikesSeconds = 0.12;
        resamplingFrequency = 200;
        outputAtOriginalSamplingRate = true;
    end
    
    
    
end