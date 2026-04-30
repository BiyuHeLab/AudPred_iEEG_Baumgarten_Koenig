Project overview 

File structure:

1) ParadigmCode: Donwload entire folder and run on Matlab 2017 +. Folders for:
    • \CreateWav\ used for creating wav files used in experiment
    • \functions\ contains helper functions used in paradigm
    
2) AnalysisCode:
a) Preprocessing:
    • Main function: NASTD_ECoG_Preproc_Main.m
      ◦ \HelperFunctions\ stores subfunctions for preprocessing
      ◦ \IED_detector_rev3.7.2014\ stores functions for detection of interictal discharges
      ◦ \SpikeDetection\ stores functions for detection of interictal discharges
b) Analysis_and_Plotting
    • \Behavior\contains scripts to analyze behavioral responses
    • \Connectivity\: Granger causality scripts; 
      ◦ Base script: NASTD_ECoG_Connectivity_Main_lk.m
    • \DataPrep\: Scripts for data preparation prior to specific analysis (e.g., filtering, NaN-removal)
    • \HistoryTracking\: Sensory history tracking (SHI) analysis
      ◦ Base script: NASTD_ECoG_HisTrack_Main_lk.m
      ◦ Subfunctions stored in \Helper Functions\
    • \Plotting\: functions used for plotting (function name indicates what is plotted)
    • \Prediction\: Prediction and Prediction error analysis
      ◦ Base script: NASTD_ECoG_Predict_Main_lk.m
      ◦ Subfunctions stored in \Helper Functions\
    • \StimulusCorrelation\: Stimulus tracking analysis
      ◦ Base script: NASTD_ECoG_StimCorr_Main_lk.m   
      ◦ Subfunctions stored in \Helper Functions\ 
c) HelperFunctions
    • Unspecific/non-project-specific functions used for analyses (e.g., statistical tests, color scales, shuffle functions)
                
