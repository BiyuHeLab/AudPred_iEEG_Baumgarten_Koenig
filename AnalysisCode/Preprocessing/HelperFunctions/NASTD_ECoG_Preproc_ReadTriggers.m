function [info_trigger] = NASTD_ECoG_Preproc_ReadTriggers...
    (sub, data_ECoGraw, plot_poststepFigs, save_poststepFigs)

%Aim: Read out auditory triggers from raw data file, choose max/first
%marker per deflection as trigger-defining data point and then
%distinguish trigger-types (3 pos. deflections = block start, 2 = trial
%start, 1 = trial end). Optional plotting and saving of process.

NASTD_ECoG_setVars
paths_NASTD_ECoG = NASTD_ECoG_paths;
NASTD_ECoG_subjectinfo

%1.1) Determine threshold for trigger selection
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));

subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings; %load in file with individual preproc infos
relthresh_level = subs_PreProcSettings.(sub).Triggerthresh; %Set Triggerthreshold (good default = )

%1.2) Determine absolute max/min thresh_level and timepoints of above/below thresh signal
triggerval_max = max(data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,:));
triggerval_maxthresh = triggerval_max*relthresh_level;
index_DPabovemaxthresh = ... %index of data points (relative to .trial) where data points are above thresh
    find(data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,:) ...
    >= triggerval_max*relthresh_level);
TP_DPabovemaxthresh = data_ECoGraw.time{1}... %time point of above thresh data points
    (find(data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,:) ...
    >= triggerval_max*relthresh_level));

triggerval_min = ...
    min(data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,:));
triggerval_minthresh = ...
    triggerval_min*relthresh_level;
index_DPbelowminthresh = ...
    find(data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,:) ...
    <= triggerval_min*relthresh_level);
TP_DPbelowminthresh = ...
    data_ECoGraw.time{1}...
    (find(data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,:) ...
    >= triggerval_min*relthresh_level));

% disp([num2str(size(TP_DPabovemaxthresh,2)) ' data points above threshold; ' ...
%     num2str(size(TP_DPabovemaxthresh,2)) ' data points below threshold']);

%plot trigger chan, max+min thres line, and markers for abovethresh data points
maxthres_line = ...
    [(data_ECoGraw.trial{1}(1,:)./data_ECoGraw.trial{1}(1,:)) ...
    .*triggerval_maxthresh]; %create row vector for plotting thres line
minthres_line = ...
    [(data_ECoGraw.trial{1}(1,:)./data_ECoGraw.trial{1}(1,:)) ...
    .*triggerval_minthresh]; %create row vector for plotting thres line

if plot_poststepFigs == 1
    figure
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %fullscreen%
    plot(1/data_ECoGraw.fsample:1/data_ECoGraw.fsample:length(data_ECoGraw.trial{1})/data_ECoGraw.fsample,... %trigger channel
        data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,:)',['k'])
    xlabel('Recording Time [sec]')
    ylabel('Trigger Channel Amplitude [au]')
    hold on;
    plot(1/data_ECoGraw.fsample:1/data_ECoGraw.fsample:length(data_ECoGraw.trial{1})/data_ECoGraw.fsample,... %max line
        maxthres_line',['r'])
    plot(1/data_ECoGraw.fsample:1/data_ECoGraw.fsample:length(data_ECoGraw.trial{1})/data_ECoGraw.fsample,... %min line
        minthres_line',['r'])
    plot(data_ECoGraw.time{1}(index_DPabovemaxthresh), ...
        data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,index_DPabovemaxthresh),['b','.'])
    plot(data_ECoGraw.time{1}(index_DPbelowminthresh), ...
        data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,index_DPbelowminthresh),['b','.'])
    title([sub ' TriggerChan WholeRecording TriggerThresh = ' num2str(relthresh_level*100) '% of max/min TriggerVal'])
    
    if save_poststepFigs == 1
        path_fig = [paths_NASTD_ECoG.Preproc_ECoGdata sub '/Figs/' 'TriggerDef/' ];        
        mkdir(path_fig);
        
        filename     = [sub '_TriggerThresh' num2str(relthresh_level*100) '%maxmin_AllTrigger.png'];
        figfile      = [path_fig filename];
        
        saveas(gcf, [figfile], 'png'); %save png version
    end
end

%1.3) Discriminate different trigger types
%6 square waves/deflections = block start - ~25 samples/deflection for fsample 512 Hz, ~102 samples/deflection for fsample 2048 Hz, duration ~ 0.3s
%2 deflections = trial start - ~25 samples/deflection for fsample 512 Hz, ~102 samples/deflection for fsample 2048 Hz, duration ~ 0.2s
%1 deflection = trial end - ~25 samples/deflection for fsample 512 Hz, ~102 samples/deflection for fsample 2048 Hz, duration ~ 0.1s

%1.3.1) Select all timepoints per deflection where signal is over threshold
%and group them per deflection
clear counter_deflection index_DPabovemaxthresh_perdeflection index_DPabovemaxthresh_perdeflection
counter_deflection = 1; %for using new cell per deflection
index_DPabovemaxthresh_perdeflection = [];
for i_abovethresh_val = 1 :(length(index_DPabovemaxthresh)-1) %loop across all found above deflection samples
    %Group them per deflection by means of temporal adjacency - select
    %all above thresh samples where the next (+1) sample is
    %really the next neighbouring sample
    if index_DPabovemaxthresh(i_abovethresh_val+1) == index_DPabovemaxthresh(i_abovethresh_val)+1 ...
            || index_DPabovemaxthresh(i_abovethresh_val+1) == index_DPabovemaxthresh(i_abovethresh_val)+2 %+2 needed for counter_deflection+1, other wise empty
        index_DPabovemaxthresh_perdeflection = [index_DPabovemaxthresh_perdeflection,index_DPabovemaxthresh(i_abovethresh_val)];
    else
        index_DPabovemaxthresh_alldeflections{counter_deflection} =  index_DPabovemaxthresh_perdeflection; %holds indices for all data points within one deflection (1 deflection per cell)
        index_DPabovemaxthresh_perdeflection = [];
        counter_deflection = counter_deflection+1;
    end
end
%and also take the final deflection into account
index_DPabovemaxthresh_alldeflections{counter_deflection} =  index_DPabovemaxthresh_perdeflection; %holds indices for all data points within one deflection (1 deflection per cell)
index_DPabovemaxthresh_perdeflection = [];

%1.3.2 Select first and max data point per deflection (because we either want the earliest or max val per trigger deflection)
if ~isempty(index_DPabovemaxthresh_alldeflections{end}) && ~isempty(index_DPabovemaxthresh_alldeflections{1})
    for i_neighbours = 1:length(index_DPabovemaxthresh_alldeflections) %loop over all thresh deflections
        index_firstDPabovemaxthresh(1,i_neighbours) = min(index_DPabovemaxthresh_alldeflections{i_neighbours}); %restricts indices to first val in deflection
        [M, I] = ...
            max(data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,index_DPabovemaxthresh_alldeflections{i_neighbours})); %restricts indices to first val in deflection
        index_maxDPabovemaxthresh(1,i_neighbours) = index_DPabovemaxthresh_alldeflections{i_neighbours}(I);
    end
elseif isempty(index_DPabovemaxthresh_alldeflections{1}) %if first deflection is missing
    for i_neighbours = 2:length(index_DPabovemaxthresh_alldeflections) %loop over above thresh deflections from 2 on
        index_firstDPabovemaxthresh(1,i_neighbours) = min(index_DPabovemaxthresh_alldeflections{i_neighbours}); %restricts indices to first val in deflection
        [M, I] = ...
            max(data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,index_DPabovemaxthresh_alldeflections{i_neighbours})); %restricts indices to first val in deflection
        index_maxDPabovemaxthresh(1,i_neighbours) = index_DPabovemaxthresh_alldeflections{i_neighbours}(I);
    end
    index_firstDPabovemaxthresh(1) = [];
    index_maxDPabovemaxthresh(1) = [];
elseif isempty(index_DPabovemaxthresh_alldeflections{end}) %if last deflection is missing
    for i_neighbours = 1:length(index_DPabovemaxthresh_alldeflections)-1 %loop over above thresh deflections eccept last one
        index_firstDPabovemaxthresh(1,i_neighbours) = min(index_DPabovemaxthresh_alldeflections{i_neighbours}); %restricts indices to first val in deflection
        [M, I] = ...
            max(data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,index_DPabovemaxthresh_alldeflections{i_neighbours})); %restricts indices to first val in deflection
        index_maxDPabovemaxthresh(1,i_neighbours) = index_DPabovemaxthresh_alldeflections{i_neighbours}(I);
    end
end

%plot trigger chan, max+min thresh line, and markers for first valid marker
%per deflection (red) and max value marker per deflection (cyan)
if plot_poststepFigs == 1
    figure
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %fullscreen%
    plot(1/data_ECoGraw.fsample:1/data_ECoGraw.fsample:length(data_ECoGraw.trial{1})/data_ECoGraw.fsample,... %trigger channel
        data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,:)',['k'])
    xlabel('Recording Time [sec]')
    ylabel('Trigger Channel Amplitude [au]')
    hold on;
    plot(1/data_ECoGraw.fsample:1/data_ECoGraw.fsample:length(data_ECoGraw.trial{1})/data_ECoGraw.fsample,... %max line
        maxthres_line',['r'])
    plot(1/data_ECoGraw.fsample:1/data_ECoGraw.fsample:length(data_ECoGraw.trial{1})/data_ECoGraw.fsample,... %min line
        minthres_line',['r'])
    plot(data_ECoGraw.time{1}(index_firstDPabovemaxthresh), ...
        data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,index_firstDPabovemaxthresh),['r','.'])    %everything larger than min distance results results in 1 hit per
    %         plot(data_ECoGraw.time{1}(index_maxDPabovemaxthresh), ...
    %             data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,index_maxDPabovemaxthresh),['c','.'])    %everything larger than min distance results results in 1 hit per
    title([sub ' TriggerChan WholeRecording TriggerThresh = ' num2str(relthresh_level*100) '% First Marker per Deflection'])
    
    if save_poststepFigs == 1
        path_fig = [paths_NASTD_ECoG.Preproc_ECoGdata sub '/Figs/' 'TriggerDef/' ];        
        
        filename     = [sub '_TriggerThresh' num2str(relthresh_level*100) '%maxmin_FirstMarker.png'];
        figfile      = [path_fig filename];
        
        saveas(gcf, [figfile], 'png'); %save png version
    end
        
    figure
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %fullscreen%
    plot(1/data_ECoGraw.fsample:1/data_ECoGraw.fsample:length(data_ECoGraw.trial{1})/data_ECoGraw.fsample,... %trigger channel
        data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,:)',['k'])
    xlabel('Recording Time [sec]')
    ylabel('Trigger Channel Amplitude [au]')
    hold on;
    plot(1/data_ECoGraw.fsample:1/data_ECoGraw.fsample:length(data_ECoGraw.trial{1})/data_ECoGraw.fsample,... %max line
        maxthres_line',['r'])
    plot(1/data_ECoGraw.fsample:1/data_ECoGraw.fsample:length(data_ECoGraw.trial{1})/data_ECoGraw.fsample,... %min line
        minthres_line',['r'])
    plot(data_ECoGraw.time{1}(index_maxDPabovemaxthresh), ...
        data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,index_maxDPabovemaxthresh),['c','.'])    %everything larger than min distance results results in 1 hit per
    title([sub ' TriggerChan WholeRecording TriggerThresh = ' num2str(relthresh_level*100) '% Max Val Marker per Deflection'])
    
    if save_poststepFigs == 1
        path_fig = [paths_NASTD_ECoG.Preproc_ECoGdata sub '/Figs/' 'TriggerDef/' ];        
        
        filename     = [sub '_TriggerThresh' num2str(relthresh_level*100) '%maxmin_MaxValMarker.png'];
        figfile      = [path_fig filename];
        
        saveas(gcf, [figfile], 'png'); %save png version
    end
end

%1.3.3 Now check how close the selected indices per deflection are,
%in order to determine amount of deflections and thus the trigger type
if data_ECoGraw.fsample == 512
    samples_between_deflections = 53; %Max number of samples between subsequent deflections
    %Idea: take (number of possible samples per deflection *2)+1
elseif data_ECoGraw.fsample == 2048
    samples_between_deflections = 207;
end

index_DP_3deflections = [];
index_DP_2deflections_uncleaned = [];
index_DP_1deflection_uncleaned = [];
for i_firstTPovermaxthresh = 1:length(index_firstDPabovemaxthresh) %loop over selected data points
    %if difference between adjacent vals is small than selected distance in
    %samples
    if i_firstTPovermaxthresh < length(index_firstDPabovemaxthresh)-2
        if index_firstDPabovemaxthresh(i_firstTPovermaxthresh+2)-index_firstDPabovemaxthresh(i_firstTPovermaxthresh) < samples_between_deflections*3 %3 deflections
            index_DP_3deflections = [index_DP_3deflections, index_firstDPabovemaxthresh(i_firstTPovermaxthresh)];
        elseif index_firstDPabovemaxthresh(i_firstTPovermaxthresh+1)-index_firstDPabovemaxthresh(i_firstTPovermaxthresh) < samples_between_deflections*2 %2 deflections
            index_DP_2deflections_uncleaned = [index_DP_2deflections_uncleaned, index_firstDPabovemaxthresh(i_firstTPovermaxthresh)];
        else
            index_DP_1deflection_uncleaned = [index_DP_1deflection_uncleaned, index_firstDPabovemaxthresh(i_firstTPovermaxthresh)];
        end
    elseif i_firstTPovermaxthresh < length(index_firstDPabovemaxthresh)-1
        if index_firstDPabovemaxthresh(i_firstTPovermaxthresh+1)-index_firstDPabovemaxthresh(i_firstTPovermaxthresh) < samples_between_deflections*2 %2 deflections
            index_DP_2deflections_uncleaned = [index_DP_2deflections_uncleaned, index_firstDPabovemaxthresh(i_firstTPovermaxthresh)];
        else
            index_DP_1deflection_uncleaned = [index_DP_1deflection_uncleaned, index_firstDPabovemaxthresh(i_firstTPovermaxthresh)];
        end
    elseif i_firstTPovermaxthresh == length(index_firstDPabovemaxthresh)-1 %Added 09.02.21 for sub NY794 to get final trial end trigger
        index_DP_1deflection_uncleaned = [index_DP_1deflection_uncleaned, index_firstDPabovemaxthresh(i_firstTPovermaxthresh)];
    elseif i_firstTPovermaxthresh == length(index_firstDPabovemaxthresh)
        index_DP_1deflection_uncleaned = [index_DP_1deflection_uncleaned, index_firstDPabovemaxthresh(i_firstTPovermaxthresh)];
    end
end
% disp([num2str(length(index_DP_3deflections)) ' Block start triggers, ' num2str(length(index_DP_2deflections_uncleaned)) ' Trial start triggers (double entries), ' num2str(length(index_DP_1deflection_uncleaned)) ' Trial end triggers (double entries)'])
%Remaining problem: double markers (e.g., within the 3 deflections (block
%start) there is also 1 2 deflections marker (trial start))

%2.3.4 Clean double indeces
%2.3.4.1 Clean double Trial start and trial end markers within block start
%deflections (idea: set minimum temporal distance between markers to
%3*(between deflection sample number)
index_DP_2deflections = [];
index_DP_1deflection_uncleaned2 = [];
for i_blocks = 1:length(index_DP_3deflections)
    if i_blocks < length(index_DP_3deflections)
        index_DP_2deflections = [index_DP_2deflections, ...
            index_DP_2deflections_uncleaned([index_DP_2deflections_uncleaned > index_DP_3deflections(i_blocks)+(samples_between_deflections*3) ...
            & index_DP_2deflections_uncleaned < index_DP_3deflections(i_blocks+1)])];
        index_DP_1deflection_uncleaned2 = [index_DP_1deflection_uncleaned2, ...
            index_DP_1deflection_uncleaned([index_DP_1deflection_uncleaned > index_DP_3deflections(i_blocks)+(samples_between_deflections*3) ...
            & index_DP_1deflection_uncleaned < index_DP_3deflections(i_blocks+1)])];
    else
        index_DP_2deflections = [index_DP_2deflections, ...
            index_DP_2deflections_uncleaned([index_DP_2deflections_uncleaned > index_DP_3deflections(i_blocks)+(samples_between_deflections*3)])];
        index_DP_1deflection_uncleaned2 = [index_DP_1deflection_uncleaned2, ...
            index_DP_1deflection_uncleaned([index_DP_1deflection_uncleaned > index_DP_3deflections(i_blocks)+(samples_between_deflections*3)])];
    end
end

%2.3.4.2 Clean double Trial end triggers (clear all trial end triggers that
%are within a second before/after the next/last trial start trigger

%NOTE: Here, we use the first sample in a single trigger deflection
%as indicator of trial end. However, the last sample in said trigger
%deflection would correlate better with the real tone sequence end. We
%correct this now be adding an additional ~100ms/6 samples) to each
%trial end in the trial definition script.
index_DP_1deflection = [];
for i_trialstarts = 1:length(index_DP_2deflections)
    if i_trialstarts < length(index_DP_2deflections)
        index_DP_1deflection = [index_DP_1deflection, ...
            index_DP_1deflection_uncleaned2([index_DP_1deflection_uncleaned2 > index_DP_2deflections(i_trialstarts)+data_ECoGraw.fsample ...
            & index_DP_1deflection_uncleaned2 < index_DP_2deflections(i_trialstarts+1)])];
    else
        index_DP_1deflection = [index_DP_1deflection, ...
            index_DP_1deflection_uncleaned2([index_DP_1deflection_uncleaned2 > index_DP_2deflections(i_trialstarts)+data_ECoGraw.fsample])];
    end
end

disp(['Total triggers found: ' num2str(length(index_DP_3deflections)) ' Block start triggers, ' num2str(length(index_DP_2deflections)) ' Trial start triggers, ' num2str(length(index_DP_1deflection)) ' Trial end triggers'])

%Plot triggerchannel overview of entire recording with all distinct markers
if plot_poststepFigs == 1
    figure
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %fullscreen%
    plot(1/data_ECoGraw.fsample:1/data_ECoGraw.fsample:length(data_ECoGraw.trial{1})/data_ECoGraw.fsample,... %trigger channel
        data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,:)',['k'])
    xlabel('Recording Time [sec]')
    ylabel('Trigger Channel Amplitude [au]')
    hold on;
    plot(1/data_ECoGraw.fsample:1/data_ECoGraw.fsample:length(data_ECoGraw.trial{1})/data_ECoGraw.fsample,... %max line
        maxthres_line',['r'])
    plot(1/data_ECoGraw.fsample:1/data_ECoGraw.fsample:length(data_ECoGraw.trial{1})/data_ECoGraw.fsample,... %min line
        minthres_line',['r'])
    
    plot(data_ECoGraw.time{1}(index_DP_1deflection), ...
        data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,index_DP_1deflection),...
        ['k','<'],'MarkerFaceColor','y','Markersize',5)%1 deflection = trial end
    plot(data_ECoGraw.time{1}(index_DP_2deflections), ...
        data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,index_DP_2deflections),...
        ['b','>'],'MarkerFaceColor','b','Markersize',5)%2 deflections = trial start
    plot(data_ECoGraw.time{1}(index_DP_3deflections), ...
        data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,index_DP_3deflections),...
        ['r','^'],'MarkerFaceColor','r','Markersize',10)%3deflections = block start
    title([sub ' FinalTriggerDefinition: ' ...
        num2str(length(index_DP_3deflections)) ' Block start triggers (red), '...
        num2str(length(index_DP_2deflections)) ' Trial start triggers (blue), '...
        num2str(length(index_DP_1deflection)) ' Trial end triggers (yellow)'])
    
    if save_poststepFigs == 1
        path_fig = [paths_NASTD_ECoG.Preproc_ECoGdata sub '/Figs/' 'TriggerDef/' ];
        
        filename     = [sub '_FinalTriggerDef_WholeRecording.png'];
        figfile      = [path_fig filename];
        
        saveas(gcf, [figfile], 'png'); %save png version
    end    
end

%Plot blockwise triggerchannel overview of entire recording with all markers
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings; %load in file with individual preproc infos
index_FirstExpBlock = subs_PreProcSettings.(sub).index_FirstExpBlock;

if plot_poststepFigs == 1
    for i_block = index_FirstExpBlock:length(index_DP_3deflections)
        
        index_blockstart = index_DP_3deflections(i_block);
        if i_block < length(index_DP_3deflections)
            index_trialstart_inblock = index_DP_2deflections([index_DP_2deflections >...
                index_blockstart & index_DP_2deflections < index_DP_3deflections(i_block+1)]);
            index_trialend_inblock = index_DP_1deflection([index_DP_1deflection >...
                index_blockstart & index_DP_1deflection < index_DP_3deflections(i_block+1)]);
        else
            index_trialstart_inblock = index_DP_2deflections([index_DP_2deflections > index_blockstart]);
            index_trialend_inblock = index_DP_1deflection([index_DP_1deflection > index_blockstart]);
        end
        
        figure
        set(gcf,'units','normalized','outerposition',[0 0 1 1]) %fullscreen%
        plot(data_ECoGraw.time{1}([index_DP_3deflections(i_block):index_trialend_inblock(end)+512]),...
            data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,...
            [index_DP_3deflections(i_block):index_trialend_inblock(end)+512]),...
            ['k']); %plot triggerchan data from block start to final trial of block end + 1sec
        xlabel('Recording Time [sec]')
        ylabel('Trigger Channel Amplitude [au]')
        hold on;
        plot(data_ECoGraw.time{1}(index_DP_3deflections(i_block)), ...
            data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,...
            index_DP_3deflections(i_block)),...
            ['r','^'],'MarkerFaceColor','r','Markersize',20)%3deflections = block start
        plot(data_ECoGraw.time{1}(index_trialstart_inblock), ...
            data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,...
            index_trialstart_inblock),...
            ['b','>'],'MarkerFaceColor','b','Markersize',10)%2 deflections = trial start
        plot(data_ECoGraw.time{1}(index_trialend_inblock), ...
            data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,...
            index_trialend_inblock),...
            ['k','<'],'MarkerFaceColor','y','Markersize',5)%1 deflection = trial end
        title([sub ' FinalTriggerDefinition for Block: ' ...
            num2str(i_block-(index_FirstExpBlock-1))...
            ' | Total TrialStart-Trigger: ' num2str(length(index_trialstart_inblock)) ...
            ' | Total TrialEnd-Trigger: ' num2str(length(index_trialend_inblock))])
        
        if save_poststepFigs == 1
            path_fig = [paths_NASTD_ECoG.Preproc_ECoGdata sub '/Figs/' 'TriggerDef/' ];        
            
            filename     = [sub '_FinalTriggerDef_Block' num2str(i_block-(index_FirstExpBlock-1)) '.png'];
            figfile      = [path_fig filename];
            
            saveas(gcf, [figfile], 'png'); %save png version
        end        
    end
end

%1.4 Copy trigger indices and time points in output file
info_trigger = struct;
num_nonexpblocks = index_FirstExpBlock-1;
for i_block = 1:length(index_DP_3deflections)
    if i_block < index_FirstExpBlock %for training blocks
        info_trigger.BlockLabel{i_block} = 'Training Block';
        info_trigger.Training.BlockStart_index(i_block) = ...
            index_DP_3deflections(i_block);
        info_trigger.Training.BlockStart_TP(i_block) = ...
            data_ECoGraw.time{1}(info_trigger.Training.BlockStart_index(i_block));
        
        info_trigger.Training.TrialStart_index{i_block} = ...
            index_DP_2deflections([index_DP_2deflections > info_trigger.Training.BlockStart_index(i_block) & ...
            index_DP_2deflections < index_DP_3deflections(i_block+1)]);
        info_trigger.Training.TrialStart_TP{i_block} = ...
            data_ECoGraw.time{1}(info_trigger.Training.TrialStart_index{i_block});
        
        info_trigger.Training.TrialEnd_index{i_block} = ...
            index_DP_1deflection([index_DP_1deflection > info_trigger.Training.BlockStart_index(i_block) & ...
            index_DP_1deflection < index_DP_3deflections(i_block+1)]);
        info_trigger.Training.TrialEnd_TP{i_block} = ...
            data_ECoGraw.time{1}(info_trigger.Training.TrialEnd_index{i_block});
        
    else %for experimental block
        info_trigger.BlockLabel{i_block} = ['Exp. Block ' num2str(i_block-num_nonexpblocks)];
        info_trigger.Exp.BlockStart_index(i_block-num_nonexpblocks) = ...
            index_DP_3deflections(i_block);
        info_trigger.Exp.BlockStart_TP(i_block-num_nonexpblocks) = ...
            data_ECoGraw.time{1}(info_trigger.Exp.BlockStart_index(i_block-num_nonexpblocks));
        
        if i_block < length(index_DP_3deflections) %for non-final block
            info_trigger.Exp.TrialStart_index{i_block-num_nonexpblocks} = ...
                index_DP_2deflections([index_DP_2deflections > info_trigger.Exp.BlockStart_index(i_block-num_nonexpblocks) & ...
                index_DP_2deflections < index_DP_3deflections(i_block+1)]);
            info_trigger.Exp.TrialStart_TP{i_block-num_nonexpblocks} = ...
                data_ECoGraw.time{1}(info_trigger.Exp.TrialStart_index{i_block-num_nonexpblocks});
            
            info_trigger.Exp.TrialEnd_index{i_block-num_nonexpblocks} = ...
                index_DP_1deflection([index_DP_1deflection > info_trigger.Exp.BlockStart_index(i_block-num_nonexpblocks) & ...
                index_DP_1deflection < index_DP_3deflections(i_block+1)]);
            info_trigger.Exp.TrialEnd_TP{i_block-num_nonexpblocks} = ...
                data_ECoGraw.time{1}(info_trigger.Exp.TrialEnd_index{i_block-num_nonexpblocks});
            
        else %for final block
            info_trigger.Exp.TrialStart_index{i_block-num_nonexpblocks} = ...
                index_DP_2deflections([index_DP_2deflections > info_trigger.Exp.BlockStart_index(i_block-num_nonexpblocks)]);
            info_trigger.Exp.TrialStart_TP{i_block-num_nonexpblocks} = ...
                data_ECoGraw.time{1}(info_trigger.Exp.TrialStart_index{i_block-num_nonexpblocks});
            info_trigger.Exp.TrialEnd_index{i_block-num_nonexpblocks} = ...
                index_DP_1deflection([index_DP_1deflection > info_trigger.Exp.BlockStart_index(i_block-num_nonexpblocks)]);
            info_trigger.Exp.TrialEnd_TP{i_block-num_nonexpblocks} = ...
                data_ECoGraw.time{1}(info_trigger.Exp.TrialEnd_index{i_block-num_nonexpblocks});
        end        
    end
end

%1.5 %Check if respective trigger-defined trial-lengths agree with stimulus
%file determined file lengths: long tone dur = 0.4s a 34 tones = 13.6 s,
%short tone dur = 0.2s a 34 tones = 6.8 s
load ([si.path_behavioraldata_sub]);%Load raw behavioral data

for i_block = 1:length(info_trigger.Exp.BlockStart_index)
    %copy continuous trials from behav matrix into trial per block struct
    data_behav_trialtonedur_perblock{i_block} = ...
        data.stim.toneDur((10*(i_block-1))+1:(10*(i_block)));
    %add potential repeated trials to behavioral data matrix
    if ~isempty(data.redoneTrialNums{i_block})
        for i_repeatedtrial = 1:length(data.redoneTrialNums{i_block})
            if length(data.missedTrialNums{i_block}) == 1 %When only one trial was misses (even multiple times)
                data_behav_trialtonedur_perblock{i_block}(end+1) ...
                    = data.stim.toneDur(data.missedTrialNums{i_block}(1));
            else %when multiple trials were missed                
                data_behav_trialtonedur_perblock{i_block}(end+1) ...
                    = data.stim.toneDur(...
                    data.missedTrialNums{i_block}(data.missedTrialNums{i_block}...
                    == data.redoneTrialNums{i_block}(i_repeatedtrial)));
            end            
        end
    end
    
    for i_trial = 1:length(info_trigger.Exp.TrialStart_index{i_block})
        if info_trigger.Exp.TrialEnd_TP{i_block}(i_trial) - ... %if trigger-def trialdur is short (smaller than 7s)
                info_trigger.Exp.TrialStart_TP{i_block}(i_trial) <= 7 && ...
                data_behav_trialtonedur_perblock{i_block}(i_trial) ~= 0.2 %and stim-def trialdur is not short
            disp(['TrialDur mismatch in Block: ' num2str(i_block) ', Trial: ' num2str(i_trial)])
        elseif info_trigger.Exp.TrialEnd_TP{i_block}(i_trial) - ... %if trigger-def trialdur is long (larger than 10s)
                info_trigger.Exp.TrialStart_TP{i_block}(i_trial) >= 13 && ...
                data_behav_trialtonedur_perblock{i_block}(i_trial) ~= 0.4 %and stim-def trialdur is not long
            disp(['TrialDur mismatch in Block: ' num2str(i_block) ', Trial: ' num2str(i_trial)])
        end
    end
end


