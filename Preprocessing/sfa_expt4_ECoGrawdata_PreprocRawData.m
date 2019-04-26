%TJB: Load in FT-structure raw ECoG data and then preprocess following steps: 
    %1) Create FT-compatible struct with selected (MNI-EDF matched, strip+grid)
        %electrodes (sfa_expt4_ECoGrawdata_CreateFTstruct )
    %2) Read out triggers and define blocks and trials
    %3) Cut continious recording into blocks ()
    %4) Demean and filter ()
    %5) Remove pulse/cardiac artefacts ()
    %6) Compute common reference data ()
    %7) Cut block data into single trials ()
    %8) Detailed preprocessing (i.e., reject noisy trials/channels)

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/data/gogodisk4/thomas/AuditoryPrediction/iEEG/')
% addpath('//gogo.sb.nyumc.org/data/gogodisk4/thomas/AuditoryPrediction/iEEG/')

sfa_expt4_setVars
paths_sfa_expt4 = sfa_expt4_paths(vars.location);

%Add base dir and own script dir
addpath(genpath(paths_sfa_expt4.BaseDir));
addpath(genpath(paths_sfa_expt4.ScriptsDir));

% addpath(genpath('/data/gogodisk4/thomas/AuditoryPrediction/iEEG/Scripts/FromRichard'));

%% 0.2) Determine subject-specific parameters (whole-recording)
i_sub = 1;
sub_list = vars.sub_list; %patients
sub = sub_list{i_sub};

sfa_expt4_subjectinfo %load subject info file (var: si)
subs_PreProcSettings = sfa_expt4_subjectPreProcSettings; %load in file with individual preproc infos

save_rawFTstruct = 0;
%Define directory for preprocessed data
preproc_dir = [paths_sfa_expt4.Preproc_ECoGdata '/' sub '/'];
mkdir(preproc_dir);
cd(preproc_dir);

plot_poststepFigs = 0;
save_poststepFigs = 0;

%% Step 1: Create FT-compatible struct with selected (MNI-EDF matched, strip+grid) electrodes (whole-recording)
data_ECoGraw = sfa_expt4_ECoGrawdata_CreateFTstruct(i_sub, save_rawFTstruct, vars);

% %plot selected channel with time scale in s
% sel_channel_index = 72;
% sel_channel_label = data_ECoGraw.label(sel_channel_index);
% figure 
% set(gcf,'units','normalized','outerposition',[0 0 1 1]) %fullscreen% 
% plot(1/data_ECoGraw.fsample:1/data_ECoGraw.fsample:length(data_ECoGraw.trial{1})/... %Selected channel
%     data_ECoGraw.fsample,data_ECoGraw.trial{1}(sel_channel_index,:)')
% 
% figure 
% set(gcf,'units','normalized','outerposition',[0 0 1 1]) %fullscreen
% plot(1/data_ECoGraw.fsample:1/data_ECoGraw.fsample:length(data_ECoGraw.trial{1})/... %trigger channel
%     data_ECoGraw.fsample,data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,:)')

%% Step 2: read out triggers from auditory trigger channel (whole-recording)
[data_ECoGraw.info_trigger] = sfa_expt4_ECoGrawdata_ReadTriggers(sub, data_ECoGraw, plot_poststepFigs, save_poststepFigs);

%% Step 3: Cut continious recording into blocks (whole-recording to blocks)
%3.1) Create trialfunction in FT-format defining blocks. Block start and 
    %end is based on trigger definition done in previous step and stored in 
    %info_trigger subfield.
    offsetBlock_inSec = 1; %enter if/how long additional time in sec should be cut out before/after block start/end
    offsetBlock_inSamples = offsetBlock_inSec*data_ECoGraw.fsample;
    block = sfa_expt4_trialfun_defineblocks(data_ECoGraw,offsetBlock_inSamples);

%3.2) cut continuous data into trials 
    cfg = [];
    cfg.trl = block;
    data_ECoGraw_blocks = ft_redefinetrial(cfg, data_ECoGraw);
    %Copy info_elec and info_trigger subfields into cfg.subfield (not done by FT-functions)
    data_ECoGraw_blocks.cfg.info_elec = data_ECoGraw.info_elec;
    data_ECoGraw_blocks.cfg.info_trigger = data_ECoGraw.info_trigger;
    
%3.3 Check if data is present for all sensors & timepoints in all blocks
for i_block = 1:length(data_ECoGraw_blocks.trial)
    if any(find(isnan(data_ECoGraw_blocks.trial{i_block}))) == 1
    missing_data = find(isnan(data_ECoGraw_blocks.trial{i_block}(1,:)));
    disp(['Data missing in block ' num2str(i_block) ' from timepoint ' num2str(data_ECoGraw_blocks.time{1}(missing_data(1))) ' sec'])
    end
end

%     clear data_ECoGraw

%     %3.4) Visually inspect data for newly cut blocks   
%     cfg = [];
%      cfg.channel = data_ECoGraw_blocks.label(subs_PreProcSettings.(sub).Trigger_chanindex); %trigger chan
%     cfg.ylim = 'maxabs';
% %       cfg.channel = 68:71; %5 chan
%     cfg.viewmode = 'vertical';
%     cfg = ft_databrowser(cfg,data_ECoGraw_blocks);
    
% %3.4) Plot power spectrum of specific channel for all blocks
%     if plot_poststepFigs == 1
%         for i_chan = data_ECoGraw_blocks.cfg.info_elec.selected.index4EDF' %for all selected chans
%               figure
%               set(gcf,'units','normalized','outerposition',[0 0 1 1]) %fullscreen%   
%               Figtitle = strcat([sub ' Power-Spectrum for Channel: ' data_ECoGraw_blocks.label{i_chan} '(' num2str(i_chan) ')']);
%               suptitle(strcat(Figtitle))              
% 
%               for i_block = 1:length(data_ECoGraw_blocks.sampleinfo)
%                 subplot(length(data_ECoGraw_blocks.sampleinfo)/4,length(data_ECoGraw_blocks.sampleinfo)/3,i_block);  
%                 pwelch(data_ECoGraw_blocks.trial{i_block}(i_chan,:),[],[],[],512);
%                 xlim([0 250]);
%                 ylim([-60 60]);       
%                 subFigtitle = strcat({['Block ' num2str(i_block)]});
%                 title(subFigtitle)
%               end      
%         
%               if save_poststepFigs == 1            
%                 path_fig = [paths_sfa_expt4.Fig_ECoGdataraw_SingleSubs sub '/' 'Preproc/PSD_raw/' ];
%                 mkdir(path_fig);
% 
%                 filename     = strcat([sub '_PSD_raw_AllBlocks_Chan' ...
%                     data_ECoGraw_blocks.label{i_chan} '(' ...
%                     num2str(i_chan) ').png']);                     
%           
%                 figfile      = [path_fig filename];              
%                 saveas(gcf, [figfile], 'png'); %save png version  
%               end
%               close
%         end
%     end
%% Step 4: Basic preprocessing 1: Detrend and filter (block-wise)
%4.1 Define importaned parameters 
    filterOrder = 3; %filter order

% %4.2 Try different filter options and plot them in order to decide on the best
    %test now in sfa_expt4_TEST_LineNoiseFilterVersions
    %Note: FT-filter don't really work - LineNoise remnants remain - presumably
    %due to amplitude chages because of long block duration. Thus, decided to
    %use manual filter)
%4.2.1 Manual version base on filtfilt.m (from Richard - EBST_DetrendFilterData.m)
data_ECoGfilt_blocks = data_ECoGraw_blocks; %copy struct
    for i_block = 1:length(data_ECoGraw_blocks.sampleinfo)
        data_ECoGfilt_blocks.trial{i_block} = ...
            sfa_expt4_ECoGrawdata_DetrendLineNoiseFilter...
            (i_block, sub, subs_PreProcSettings, filterOrder, data_ECoGraw_blocks);
    end

%4.3 Check if data is present for all sensors & timepoints in all blocks
for i_block = 1:length(data_ECoGfilt_blocks.trial)
    if any(find(isnan(data_ECoGfilt_blocks.trial{i_block}))) == 1
    missing_data = find(isnan(data_ECoGfilt_blocks.trial{i_block}(1,:)));
    disp(['Data missing in block ' num2str(i_block) ' from timepoint ' num2str(data_ECoGfilt_blocks.time{1}(missing_data(1))) ' sec'])
    end
end    
    
%4.4 Plot manualy filtered data for selected electrodes
    if plot_poststepFigs == 1
        tic
        for i_chan = data_ECoGfilt_blocks.cfg.info_elec.selected.index4EDF'
            figure
            set(gcf,'units','normalized','outerposition',[0 0 1 1]) %fullscreen%for i_chan = 100 %length(data_ECoGraw_blocks.label)
              for i_block = 1:length(data_ECoGfilt_blocks.sampleinfo)      
                subplot(length(data_ECoGfilt_blocks.sampleinfo)/4,length(data_ECoGfilt_blocks.sampleinfo)/3,i_block);    
                pwelch(data_ECoGfilt_blocks.trial{i_block}(i_chan,:),[],[],[],512); %filtered data
                subFigtitle = strcat({['Block:' num2str(i_block)]});
                title(subFigtitle)
                xlim([0 250]);
        %         xlim([55 65]);
                ylim([-60 60]);       
              end

              Figtitle = [sub ' Power-Spectrum after manual BS-Filter for Channel:' ...
                  data_ECoGfilt_blocks.label(i_chan) '(' num2str(i_chan) ')'];
              suptitle(strcat(Figtitle{:}))

                if save_poststepFigs == 1            
                    path_fig = [paths_sfa_expt4.Fig_ECoGdataraw_SingleSubs sub '/' 'Preproc/PSD_afterFilt/' ];
                    mkdir(path_fig);

                    filename     = strcat([sub '_PSD_LNfilterDetrend_AllBlocks_Chan' ...
                        data_ECoGfilt_blocks.label{i_chan} '(' ...
                        num2str(i_chan) ').png']);            
                    figfile      = [path_fig filename];                

                    saveas(gcf, [figfile], 'png'); %save png version  
                end
                close
        end
        toc
    end

%% 5) Remove pulse/cardiac artefacts
%5.1 Richards Option (remove waveform of estimated heart beat)

%5.2 Mutual information with afferent ECG channels (MI-informed rejection
%of of ICA channels)


%% 6) Compute common reference data (block-wise)
%Idea: compute average signal (per sample) across all ECoG electrodes 
%within a block and then subtract this common average from each individual 
%electrode. However, for that to work, we need to determine if all channels 
%are clean (i.e., low variance of signal), in order not to mess up 
%the common average.

%6.1 Determine Threshold (as multiple of STD) determining which channels 
    %are excluded from CAR 
%     numstd_thresh = 3; %determine threshhold level as multiple of std 
numstd_thresh = subs_PreProcSettings.(sub).numSTDthresh_forCAR; %determine threshhold level as multiple of std

%6.2 Determine suprathresh channels, take them out, compute common average 
    %reference (CAR) pr block based on remaining trials and then subtract that 
    %reference from every channel per block. Plot images for threshold and
    %cross-correlation pre vs post-referencing
data_ECoGfiltref_blocks = sfa_expt4_ECoGrawdata_CommonAvgReference...
    (sub, numstd_thresh, data_ECoGfilt_blocks, subs_PreProcSettings, paths_sfa_expt4, plot_poststepFigs, save_poststepFigs);


%6.3 Plot re-referenced data for selected electrodes
    if plot_poststepFigs == 1
        tic
        for i_chan = data_ECoGfiltref_blocks.cfg.info_elec.selected.index4EDF'
            figure
            set(gcf,'units','normalized','outerposition',[0 0 1 1]) %fullscreen%for i_chan = 100 %length(data_ECoGraw_blocks.label)
              for i_block = 1:length(data_ECoGfiltref_blocks.sampleinfo)      
                subplot(length(data_ECoGfiltref_blocks.sampleinfo)/4,length(data_ECoGfiltref_blocks.sampleinfo)/3,i_block);    
                pwelch(data_ECoGfiltref_blocks.trial{i_block}(i_chan,:),[],[],[],512); %filtered data
                subFigtitle = strcat({['Block:' num2str(i_block)]});
                title(subFigtitle)
                xlim([0 250]);
        %         xlim([55 65]);
                ylim([-60 60]);       
              end

              Figtitle = [sub ' Power-Spectrum after common avg referencing for Channel:' ...
                  data_ECoGfiltref_blocks.label(i_chan) '(' num2str(i_chan) ')'];
              suptitle(strcat(Figtitle{:}))

                if save_poststepFigs == 1            
                    path_fig = [paths_sfa_expt4.Fig_ECoGdataraw_SingleSubs sub '/' 'Preproc/CommonAvgRef/' ];
                    mkdir(path_fig);

                    filename     = strcat([sub '_PSD_CAref_AllBlocks_Chan' ...
                        data_ECoGfiltref_blocks.label{i_chan} '(' ...
                        num2str(i_chan) ').png']);            
                    figfile      = [path_fig filename];                

                    saveas(gcf, [figfile], 'png'); %save png version  
                end
                close
        end
        toc
    end

%6.4 Check post-referencing power spectrum to identify noisy channels -
%note noisy channels in subjectPreProcSettings.m
   
    
%% To do: Channel-wise block-wise plot for after rerencing in order to 
%determine where the additional noise peaks (ratio of power for 99-101 Hz vs small window before/after)
%Aim: Get rid of additional peaks

%% 7) Cut block data into single trials
%Trial definition based on above-defined auditory triggers
%7.1) Create trialfunction in FT-format defining trials. Trial start and 
    %end is based on trigger definition done in previous step and stored in 
    %info_trigger subfield.    
    offsetTrial_inSec = 1; %enter if/how long additional time in sec should be cut out before/after block start/end
    offsetTrial_inSamples = offsetTrial_inSec*data_ECoGfiltref_blocks.fsample;
    trials = sfa_expt4_trialfun_definetrials(data_ECoGfiltref_blocks,offsetTrial_inSamples);

%7.2) Cut block data into trials 
    cfg = [];
    cfg.trl = trials;
    data_ECoGfiltref_trials = ft_redefinetrial(cfg, data_ECoGfiltref_blocks);
  
    %Copy info_elec and info_trigger subfields into cfg.subfield (not done by FT-functions)
    data_ECoGfiltref_trials.cfg.info_elec = data_ECoGfiltref_blocks.cfg.info_elec;
    data_ECoGfiltref_trials.cfg.info_trigger = data_ECoGfiltref_blocks.cfg.info_trigger;
    data_ECoGfiltref_trials.cfg.info_ref = data_ECoGfiltref_blocks.cfg.info_ref;
    %Clear old fields from cfg.previous
    data_ECoGfiltref_trials.cfg.previous = rmfield(data_ECoGfiltref_trials.cfg.previous,'info_elec');
    data_ECoGfiltref_trials.cfg.previous = rmfield(data_ECoGfiltref_trials.cfg.previous,'info_trigger');
    data_ECoGfiltref_trials.cfg.previous = rmfield(data_ECoGfiltref_trials.cfg.previous,'info_ref');

    
%7.3 Check if data is present for all sensors & timepoints in all trials
for i_trial = 1:length(data_ECoGfiltref_trials.trial)
    if any(find(isnan(data_ECoGfiltref_trials.trial{i_trial}))) == 1
    missing_data = find(isnan(data_ECoGfiltref_trials.trial{i_trial}(1,:)));
    disp(['Data missing in trial ' num2str(i_trial) ' from timepoint ' num2str(data_ECoGfiltref_trials.time{1}(missing_data(1))) ' sec'])
    end
end
    
    
% %7.4) Visually inspect data for newly cut trials   
    cfg = [];
    cfg.channel = data_ECoGfiltref_trials.label(subs_PreProcSettings.(sub).Trigger_chanindex); %trigger chan
    cfg.channel = data_ECoGfiltref_trials.label(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF); %selected channels
    
    cfg.ylim = 'maxabs';
%      cfg.channel = 71; %5 chan
    cfg.viewmode = 'vertical';
    cfg = ft_databrowser(cfg,data_ECoGfiltref_trials);
    
    
%7.4) Save summary file to indiv. preproc path
    savefile = [preproc_dir sub '_data_preproc.mat'];
    save(savefile, 'data_ECoGfiltref_trials','-v7.3');
    
% %% 8) Detailed preprocessing (i.e., reject noisy trials/channels)
% 
% %8. Preprocessing Step 1 - remove jump and muscle artefacts (important to do before ICA)
% %Aim; Remove big/prominent jump and muscle artifacts - be liberal, since we will do the thorough rejection afterwards
%     %8.1Jumps
%     cfg = [];
%     cfg.artfctdef.reject = 'partial';
%     cfg.trl = trials;
%     cfg.continuous = 'no'; %because we have trial-wise data now
%     
%     cfg.artfctdef.jump.channel = data_ECoGfiltref_trials.cfg.info_elec.selected.Label; %only selected ECoG channels
%     cfg.artfctdef.jump.cutoff = 20;
%     cfg.artfctdef.jump.trlpadding = 0;
%     cfg.artfctdef.jump.artpadding = 0;
%     cfg.artfctdef.jump.fltpadding = 0;
%     cfg.artfctdef.jump.cumulative = 'yes';
%     cfg.artfctdef.jump.medianfilter = 'yes';
%     cfg.artfctdef.jump.medianfilterord = 9;
%     cfg.artfctdef.jump.absdiff = 'yes';
%     cfg.artfctdef.jump.interactive = 'yes';
% 
%     [cfg, artifcat_jump] = ft_artifact_jump(cfg, data_ECoGfiltref_trials);
% 
%     %8.2Muscle
%     cfg.artfctdef.muscle.channel = data_ECoGfiltref_trials.cfg.info_elec.selected.Label;  %only selected ECoG channels
%     cfg.artfctdef.muscle.cutoff = 4;
%     cfg.artfctdef.muscle.trlpadding = 0;
%     cfg.artfctdef.muscle.artpadding = 0;
%     cfg.artfctdef.muscle.fltpadding = 0.1;
%     cfg.artfctdef.muscle.bpfilter = 'yes';
%     cfg.artfctdef.muscle.bpfreq = [110 140];
%     cfg.artfctdef.muscle.bpfiltord = 9;
%     cfg.artfctdef.muscle.bpfilttype = 'but';
%     cfg.artfctdef.muscle.hilbert = 'yes';
%     cfg.artfctdef.muscle.boxcar = 0.2;
%     cfg.artfctdef.muscle.interactive = 'yes';
% 
%     [cfg, artifcat_muscle] = ft_artifact_muscle(cfg, data_ECoGfiltref_trials);
% 
%     %EOG
%     cfg = [];
%     cfg.trl = trl_struct;
%     cfg.continuous = 'yes';
%     
%     cfg.artfctdef.zvalue.channel = data_ECoGfiltref_trials.cfg.info_elec.selected.Label;  %only selected ECoG channels
%     cfg.artfctdef.zvalue.cutoff = 4;
%     cfg.artfctdef.zvalue.trlpadding = 0;
%     cfg.artfctdef.zvalue.artpadding = 0.1;
%     cfg.artfctdef.zvalue.fltpadding = 0;
%     
%     cfg.artfctdef.zvalue.bpfilter = 'yes';
%     cfg.artfctdef.zvalue.bpfilttype = 'but';
%     cfg.artfctdef.zvalue.bpfreq = [1 15];
%     cfg.artfctdef.zvalue.bpfiltord = 4;
%     cfg.artfctdef.zvalue.hilbert = 'yes';
%     
%     cfg.artfctdef.zvalue.interactive = 'yes';
%     
%     [cfg, artifcat_EOG] = ft_artifact_zvalue(cfg, data_ECoGfiltref_trials);
% 
%     %6.3Remove semi-automatically determined artifacts
%     cfg.artfctdef.reject = 'partial';
%     cfg.artfctdef.minaccepttim = 1; %Minimum length in s of remaining trial
%     % cfg.artfctdef.feedback        = 'yes';
%     dataAll_EEGcleanedSemiauto = ft_rejectartifact(cfg,dataAll_preproc1);



    