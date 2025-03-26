NASTD_iEEG Project overview 

Project aim: 
    • Show how the human brain processes naturalistic auditory tone sequences using iEEG recordings
    • Focus on prediction/prediction error components, sensory history integration comparisons across different tone durations, directed connectivity measures
Project duration: 
    • preparation: 2018
    • Recording: Oct 2018 - Nov 2020
    • Analysis: 2020-2023
Responsible persons:
    • Preparation and planning: Brian Maniscalco, Jennifer L. Lee, Thomas J. Baumarten, Biyu J. He

File structure:
Base directory: 
\\bjhlpdcpvm04.nyumc.org\VMdrive\NaturalisticAuditorySequences_ToneDuration(NAS_TD)\ECoG\
Documents:
    • \Docs\
        ◦ Analysis notes
        ◦ Notes from previous colleagues (Jenn, Brian Maniscalco)
        ◦ Project-related talks
        ◦ Interim overview for respective analysis steps
Paradigm/Experiment files:
    • \Experiment_Versions\
        ◦ Different versions of paradigm
        ◦ Used version: sfa_expt4_Adeen_Current
        ◦ Images of sound sequences
    • \Experiment_Versions\Trigger\
        ◦ Files for used auditory trigger (delineating trial start/end)
Scripts:
    • NASTD_ECoG_paths.m: master script for project paths
    • NASTD_ECoG_setVars.m: master script for analysis variables
    • NASTD_ECoG_subjectinfo.m: master script for subject-specific variables (important for preprocessing of single-subject data)
    • \NASTD_ECoG_Matlab\
        ◦ Separated by analysis topic
            ▪ \Behavior\
                • Scripts to analyze behavioral responses
            ▪ \CreateWav\
                • Scripts to create sound files used in experiment
                • Subfolder \Figs\ stores figures showing tones used in each sound sequence
                • Subfolder \Sequences\ stores sound files for tone sequences (*.wav)
            ▪ \HelperFunctions\
                • Unspecific/non-project-specific functions used for analyses (e.g., statistical tests, color scales, shuffle functions)
            ▪ \Plotting\
                • Functions used for plotting (function name indicates what is plotted)
                • Subfolder \HelperFunctions\ stores unspecific/non-project-specific plotting functions
            ▪ \Preprocessing\
                • Different functions used for data preprocessing
                    ◦ Main function: NASTD_ECoG_Preproc_Main.m
                    ◦ Subfolder \HelperFunctions\ stores subfunctions for preprocessing
                    ◦ Subfolder \IED_detector_rev3.7.2014\ stores functions for detection of interictal discharges
                    ◦ Subfolder \SpikeDetection\ stores functions for detection of interictal discharges
            ▪ \Project_OCNC\
                • Folder with scripts for Thomas’ trip to the OCNC 2019
                • Not project-relevant
            ▪ \TrialAnalysis\
                • Stores all scripts for analyses based on trial-defined data (basically all analyses)
                • \Connectivity\: Granger causality scripts; 
                    ◦ Base script: NASTD_ECoG_Connectivity_Main.m
                • \DataPrep\: Scripts for data preparation prior to specific analysis (e.g., filtering, NaN-removal)
                • \HistoryTracking\: Sensory history tracking (SHI) analysis
                    ◦ Base script: NASTD_ECoG_HisTrack_Main.m
                    ◦ Subfunctions stored in \Helper Functions\
                • \Prediction\: Prediction and Prediction error analysis
                    ◦ Base script: NASTD_ECoG_Predict_Main.m
                    ◦ Subfunctions stored in \Helper Functions\
                • \StimulusCorrelation\: Stimulus tracking analysis
                    ◦ Base script: NASTD_ECoG_StimCorr_Main.m
                    ◦ Subfolder \SequenceTracking\: Scripts from previous project (Maniscalco et al., 2018) used as template
                • \TFA\: Time-frequency analyses (not in final manuscript version)
                    ◦ Base script: NASTD_ECoG_TFA_Main.m
                • \TLA\: Time-locked (ERF) analysis (not in final manuscript version)

Data:
Raw data:
    • \Data\raw\
        ◦ Subject-wise (NY***) data, including iEEG time series data (*.EDF), electrode labels (*.txt), electrode locations (\images\*.png), recording notes (*__NotesMeasurement.docx)
Preprocessed data:
    • \Data\preprocessed\
        ◦ Subject-wise (NY***) data (artifact-cleaned, trial-separated); 2 versions: with (NY***_DataClean_AllTrials.mat) or without pulse-artifact correction (NY***_DataClean_AllTrials_NoPulseIEDCorr.mat)
        ◦ For each subject, subfolder with images for
            ▪ Trigger definition
            ▪ Plots for different steps of preprocessing
Analysis data:
    • \Analysis\
        ◦ \Behavior\: Group-level and single subject analyses of behavioral responses
            ▪ Central group-level results Final tone pitch likelihood rating (Manuscript Fig. 2) 
            ▪ Note: Sub1 (NY622 not added in group due to different paradigm version and different tone sequences)
        ◦ \Connectivity\: Granger Causality (GC) analyzes
            ▪ \ElectSelect\: Files and figures containing target electrodes for (GC) analyzes
            ▪ \GC\: Files and figures containing GC analyzes results
                • Files distinguished by p-value threshold, signal component (broadband vs. high gamma), multiple-comparison correction (uncorr vs. corr), analysis time window (full time window vs. first 100 ms)
                • Subfolder \PEeffect\: GC analysis results for prediction error effect
                • Subfolder \Figs\: result figures (basis for manuscript Fig. 6)
        ◦ \ElecParcel\: Figures showing electrode categorization into different anatomical regions (basis for manuscript Fig. 1C and S1)
        ◦ \Figures\: Outdated figures
        ◦ \HisTrack\: Sensory history tracking results and figures for group-level and single subjects
            ▪ Subfolder \Exp\: Results for experimental/real/recorded data
            ▪ Subfolder \Shuffle\: Results for null distribution based on shuffled data
            ▪ Subfolder \ExpvsShuff\: Final results from comparison experimental/real/recorded vs. null distribution data
            ▪ Subfolder \Figs\: Result figures, separated by surface plots (\Kprime_Surf\), region-specific plots (\Kprime_AnatReg\), surface ratio-plots (\KprimeTDRatio_Surf\) and videos (\vid\, i.e., significant electrodes on surface plot as function of time)
            ▪ Files distinguished by p-value threshold, signal component, tone duration (200ms vs. 400 ms), multiple-comparison correction (uncorr vs. corr), analysis time window (averaged across all time windows vs. per time window); (basis for manuscript Fig. 5)
        ◦ \Prediction\: Prediction and prediction error analyzes files and figures for group-level and single subjects
            ▪ Separated into timelocked/ERF effects (\PredEffects_ERF\; shows recording of sign. Electrodes as function of time) and non-timelocked effects (\PredEffects_Surf\; shows just sign. Electrodes on brain surface)
            ▪ Separated into prediction (\PredEffect\), simple prediction error (\SimplePredErrEffect\) and complex prediction error (\ComplexPredErrEffect\) effects
            ▪ Separated by statistical method applied (cluster corrected over time (\ClusterCorr\), FDR-corrected (\FDR corr\), uncorrected (\uncorr\))
            ▪ Videos (\vid\, i.e., significant electrodes on surface plot as function of time)
            ▪ Files distinguished by signal component, tone duration (200ms vs. 400 ms vs. all), 
            ▪ Basis for manuscript Fig. 3,4
        ◦ \StimCorr\: Sequence tracking analysis result files and figures for group-level and single subjects
            ▪ Note: Sub1 (NY622 not added in group due to different paradigm version and different tone sequences)
            ▪ Basis for manuscript Fig. 2
        ◦ \TFA\: Time-frequency analyses result files and figures for group-level and single subjects
            ▪ Not used in final manuscript
        ◦ \TimelockedResp\: Time-locked/ERF analyses results files and figures for group-level and single subjects
            ▪ Not used in final manuscript
