%% Project: NASTD_ECoG
%Localizing electrodes showing timelocked respose to auditory sitmulation
%(tonal tracking)

%1. Load in pre-processed data
%2. Prepare input data
%3. Optionally compute amplitude envelope for specific freuency band
%4. Compute N100 as index for early auditory processing
%5. Plot output (single-subject and across-subjects)


%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
%Add project base path

NASTD_ECoG_setVars
paths_NASTD_ECoG = NASTD_ECoG_paths;

addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
%Add project base and script dir

%Determine subjects
sub_list = vars.sub_list;

%Load in file with individual preproc infos
subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;

plot_poststepFigs = 1;
save_poststepFigs = 1;

%% 1. Load in pre-processed data
for i_sub = vars.validSubjs
   
    sub = sub_list{i_sub};
    
    timelock_dir = [paths_NASTD_ECoG.ECoGdata_Timelocked sub '/'];
    if (~exist(timelock_dir, 'dir')); mkdir(timelock_dir); end
    
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    
    %Load in preprocessed neural and behavioral data
    loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];%path to indiv preprocessed ECoG data
    
    % 0.3) Load in preproc data
    tic
    disp(['Loading preprocessed data set for sub: ' sub])
    load(loadfile_ECoGpreprocdata);
    disp(['done loading in ' num2str(toc) ' sec'])
    
    % 0.4) Specify analysis info
    Selected_Elecs = ...
        setdiff(DataClean_AllTrials.cfg.info_elec.selected.index4EDF, ...
        subs_PreProcSettings.(sub).rejectedChan_index);
    ToneDur_text = {'0.2' '0.4'};
    SampleFreq = DataClean_AllTrials.fsample;
    
    %% 2) Prepare input data
    %2.1) Select only clean/valid trials
    DataClean_CleanTrials = NASTD_ECoG_Preproc_SelCleanTrials(sub, DataClean_AllTrials);
    
    %2.2) Select only valid electrodes (Grid+Strip, clean, MNI-coordinates present)
    Index_Valid_Elecs = setdiff(DataClean_CleanTrials.cfg.info_elec.selected.index4EDF, ...
        subs_PreProcSettings.(sub).rejectedChan_index); %via index
    Label_Valid_Elecs = setdiff(DataClean_CleanTrials.cfg.info_elec.selected.Label, ...
        subs_PreProcSettings.(sub).rejectedChan_label); %via label
    if ~isempty(find((strcmp(sort(DataClean_CleanTrials.label(Index_Valid_Elecs)), ...
            sort(Label_Valid_Elecs))) ==0))
        disp(['Mismatch in electrode selection']);
    end
    
    %    cd([paths_NASTD_ECoG.FieldTrip 'utilities/']) %In case ft_selectdata not found
    cfg         = [];
    cfg.channel = (Label_Valid_Elecs);
    DataClean_CleanTrialsElecs = ft_selectdata(cfg, DataClean_CleanTrials);
    
    %2.3) Select only trials with specific Tone Duration (TD)
    for i_TD = 1:length(ToneDur_text)
        DataClean_CleanTrialsElecs_perTD{i_TD} = ...
            NASTD_ECoG_Preproc_SelTrialsperTD...
            (ToneDur_text{i_TD}, DataClean_CleanTrialsElecs);
        
        disp(['Total trial number for TD ' ToneDur_text{i_TD} 's: ' ...
            num2str(length(DataClean_CleanTrialsElecs_perTD{i_TD}.trial))])
    end
    
    %Cleanup
    clear DataClean_CleanTrialsElecs DataClean_CleanTrials DataClean_AllTrials
    
    %% 3) Prepare input signals: slow arrhythmic activity & gamma amplitude envelope
    
%     %3.0A Baseline Correction: Normalize whole-trial MEG activity relative to prestimulus baseline 
%     BL_win = [-0.5 0];
%     
%     for i_TD = 1:length(ToneDur_text)
%         for i_trial = 1:length(DataClean_CleanTrialsElecs_perTD{i_TD}.trial)
%             
%             for i_elec = 1:length(DataClean_CleanTrialsElecs_perTD{i_TD}.label)
%                 DataClean_CleanTrialsElecs_perTD{i_TD}.trial{i_trial}(i_elec,:) = ...
%                     DataClean_CleanTrialsElecs_perTD{i_TD}.trial{i_trial}(i_elec,:) - ...
%                     mean(DataClean_CleanTrialsElecs_perTD{i_TD}.trial{i_trial}(i_elec,...
%                     find(DataClean_CleanTrialsElecs_perTD{i_TD}.time{1} == BL_win(1)): ...
%                     find(DataClean_CleanTrialsElecs_perTD{i_TD}.time{1} == BL_win(2))));
%             end
% 
%         end
%     end
    %3.0B Demean: Normalize whole-trial MEG activity for across-trial average     
    for i_TD = 1:length(ToneDur_text)
        for i_trial = 1:length(DataClean_CleanTrialsElecs_perTD{i_TD}.trial)
            
            for i_elec = 1:length(DataClean_CleanTrialsElecs_perTD{i_TD}.label)
                DataClean_CleanTrialsElecs_perTD{i_TD}.trial{i_trial}(i_elec,:) = ...
                    DataClean_CleanTrialsElecs_perTD{i_TD}.trial{i_trial}(i_elec,:) - ...
                    mean(DataClean_CleanTrialsElecs_perTD{i_TD}.trial{i_trial}(i_elec,:));
            end

        end
    end 
    
    %3.1 Slow arrhythmic activity
    FilterOrder = 3;
    LowPassFreq = 35;
    
    for i_TD = 1:length(ToneDur_text)
        for i_trial = 1:length(DataClean_CleanTrialsElecs_perTD{i_TD}.trial)
            
            %LP filter
            [lpb,lpa] = butter(FilterOrder,(LowPassFreq*2)/SampleFreq,'low');
            DataClean_CleanTrialsElecs_perTD{i_TD}.LP35Hz{i_trial} = ...
                filtfilt(lpb,lpa,DataClean_CleanTrialsElecs_perTD{i_TD}.trial{i_trial}')';
            
            %         %Visual check (LP filtered vs. non-LP filtered
            %         figure
            %         plot(DataClean_CleanTrialsElecs_perTD{i_TD}.trial{i_trial}(1,:),'r-')
            %         hold on; plot(DataClean_CleanTrialsElecs_perTD{i_TD}.LP35Hz{i_trial}(1,:),'b-')
            %         legend('Non-LP-filtered data','LP-filtered data')
            
        end
    end
    
    %3.2 Frequency-band specific amplitude envelope
    FrequencyBands = [30, 70; 70, 150]; %Hz [HP,LP]
    FrequencyBand_Labels = {'Gamma', 'HighGamma'};
    
    for i_TD = 1:length(ToneDur_text)
        for i_freqbands = 1:length(FrequencyBands)
            
            DataClean_CleanTrialsElecs_perTD{i_TD} = ...
                NASTD_ECoG_Preproc_AmpEnvelope...
                (sub, FilterOrder, ...
                FrequencyBands(i_freqbands,1), FrequencyBands(i_freqbands,2), ...
                FrequencyBand_Labels{i_freqbands}, ...
                DataClean_CleanTrialsElecs_perTD{i_TD}, ...
                paths_NASTD_ECoG, plot_poststepFigs, save_poststepFigs);
            
        end
    end    
    
    %% 4) Use N100/tonal response as functional localizer
    %Compute N100 windowed activity timelocked to tone presentation
    %to determine electrodes showing increased activity after tone presentation
    
    N1_win = [0.09 0.11]; %should focus on N1-peak at 100 ms
    BL_win = [-0.5 0];
    
    inputData = {'LP35Hz','Gamma_NormLogAmp', 'HighGamma_NormLogAmp'};
    ThresComp = {'BL','AvgTone';};
    
    for i_TD = 1:length(ToneDur_text)
        for i_inputData = 1:length(inputData)
            
            if i_inputData == 1
                STD_thresh = 5; %subs_PreProcSettings.(sub).N1localizer_STDthresh.LF;
            else
                STD_thresh = 5; %subs_PreProcSettings.(sub).N1localizer_STDthresh.GammaAmp;
            end
                        
            %4.1 Use localizer to determine electrodes with tonal tracking
            N100SuperThresElecs{i_sub}.(inputData{i_inputData}){i_TD} = ...
                NASTD_ECoG_TLA_CompN1localizer...
                (sub, inputData{i_inputData}, ToneDur_text{i_TD}, ...
                N1_win, BL_win, ThresComp{1}, STD_thresh, ...
                DataClean_CleanTrialsElecs_perTD{i_TD}, ...
                plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
            
        end          
    end
    
    %% 5. Plot brain volume with N100 superthres elecs for all subjects
    
    %%Stopped here on 12/03/21
    %to do: update plotting for single subject / all electrodes, update
    %plotting across subjects / all electrodes
    
    cmap = 'hot';
    view_angle = [270,0]; %(270,0) - LH, (90,0) = RH
    
    %Read out all coordinates for all elecs for each subject
    for i_sub = 1:length(sub_list)
        
        sub = sub_list{i_sub};
        timelock_dir = [paths_NASTD_ECoG.ECoGdata_Timelocked '/' sub '/'];
        
        NASTD_ECoG_subjectinfo %load subject info file (var: si)
        subs_PreProcSettings = NASTD_ECoG_SubPreprocSettings; %load in file with individual preproc infos
        %Load in preprocessed neural and behavioral data
        loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];%path to indiv preprocessed ECoG data
        tic
        disp(['Loading preprocessed data set for sub: ' sub])
        load(loadfile_ECoGpreprocdata);
        disp(['done loading in ' num2str(toc) ' sec'])
        
        coords_sub{i_sub} = DataClean_AllTrials.hdr.elec.chanpos;
        
        clear DataClean_AllTrials
    end
    
    %Optional: Combine elecs across TD
    temp = [];
    for i_sub = 1:length(sub_list)
        for i_inputData = 1:length(inputData)
            for i_TD = 1:length(ToneDur_text)
                temp = [temp, TonalTrackingElec{i_sub}.(inputData{i_inputData}).Index{i_TD}];
            end
            TonalTrackingElec_AllTD{i_sub}.(inputData{i_inputData}).Index = unique(temp);
            temp = [];
        end
    end
    
    for i_inputData = 1:length(inputData)
        for i_TD = 1:length(ToneDur_text)
            
            coords = [];
            vals = [];
            chanSize = [];
            
            %3.1A Only tonal sensors
            for i_sub = 1:length(sub_list)
                
                clims = [1,length(sub_list)]; %colormap limits
                coords = [coords; TonalTrackingElec{i_sub}.(inputData{i_inputData}).ChanPos{i_TD}];
                vals = [vals; ones(length(TonalTrackingElec{i_sub}.(inputData{i_inputData}).Index{i_TD}),1)*i_sub]; %used to differentiate subjects
                chanSize = [chanSize; ones(length(TonalTrackingElec{i_sub}.(inputData{i_inputData}).ChanPos{i_TD}),1)*3]; %plot tonal elecs
                %             chanSize = [chanSize; ones(length(TonalTrackingElec_AllTD{i_sub}.(inputData{i_inputData}).Index),1)*3]; %plot tonal elecs combined across TD
                
            end
            
            %3.1B All sensors but highlight tonal sensors
            for i_sub = 1:length(sub_list)
                
                coords = [coords; coords_sub{i_sub}];
                vals = [vals; ones(length(coords_sub{i_sub}),1)*i_sub]; %used to differentiate subjects
                chanSize_sub = ones(length(coords_sub{i_sub}),1);
                chanSize_sub(TonalTrackingElec{i_sub}.(inputData{i_inputData}).Index{i_TD}) = 3; %highlight tonal elecs
                %             chanSize_sub(TonalTrackingElec_AllTD{i_sub}.(inputData{i_inputData}).Index) = 3; %highlight tonal elecs combined across TD
                chanSize = [chanSize; chanSize_sub];
                clims = [1,length(sub_list)+1]; %colormap limits
                
            end
            
            NASTD_ECoG_Plot_PlotEleconSurf...
                (coords,vals,chanSize,...
                clims,cmap,view_angle,...
                inputData{i_inputData}, ToneDur_text{i_TD}) %Single hemisphere, defined by view_angle
            
            NASTD_ECoG_Plot_PlotEleconSurf_SubplotHemis...
                (coords,vals,chanSize,...
                clims,cmap,...
                inputData{i_inputData}, ToneDur_text{i_TD}) %Both hemispheres in subplots
            
        end
    end
end