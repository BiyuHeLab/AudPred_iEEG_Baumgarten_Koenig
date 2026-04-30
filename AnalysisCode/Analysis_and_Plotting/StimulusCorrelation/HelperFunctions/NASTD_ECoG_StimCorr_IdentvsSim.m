function NASTD_ECoG_StimCorr_IdentvsSim...
    (sub, FuncInput_DataType, FuncInput_ToneDur_text,  ...
    param,...
    FuncInput_InputData, ...
    save_poststepFigs, paths_NASTD_ECoG)

%Aim: Compute tone sequence tracking in neural data for identical vs similar
%(i.e., same trend/p*34) tone sequences.

%% 0.1) Specify vars, paths, and setup fieldtrip
%Add base dir and own script dir
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));

path_dataoutput = [paths_NASTD_ECoG.ECoGdata_StimCorr '/' sub '/Data/'];
if (~exist(path_dataoutput, 'dir')); mkdir(path_dataoutput); end

path_fig = ([paths_NASTD_ECoG.ECoGdata_StimCorr ...
    '/' sub '/Figs/IdentvsSim/']);
if (~exist(path_fig, 'dir')); mkdir(path_fig); end

%% 1. Select current input data in FT-struct
FuncInput_InputData.trial = [];
FuncInput_InputData.trial = FuncInput_InputData.(FuncInput_DataType);

%% 2. Determine time points and samples for each tone
nTrials             = length(FuncInput_InputData.trial);
nSensors            = size(FuncInput_InputData.trial{1},1);
SampleFreq          = FuncInput_InputData.fsample;
ToneDur_Sec         = str2num(FuncInput_ToneDur_text);
nSamples_perTone    = ToneDur_Sec * SampleFreq;
SensorLabels        = FuncInput_InputData.label;

%2.2 Determine TP/samples for each tone start+end
TP_Tone_StartStop       = NaN(36,2);
TP_Tone_StartStop(1,1)  = 0; %p1 set as t = 0 in trial definition
for i_tone = 2:36
    Dist = abs(FuncInput_InputData.time{1} - ...
        ((str2num(FuncInput_ToneDur_text)*i_tone) - str2num(FuncInput_ToneDur_text)));
    minDist = min(Dist);
    i_minDist = find(Dist == minDist);
    TP_Tone_StartStop(i_tone,1) = FuncInput_InputData.time{1}(i_minDist);
end
for i_tone = 1:35
    i_LastSampleTone = find(FuncInput_InputData.time{1} == TP_Tone_StartStop(i_tone + 1,1));
    TP_Tone_StartStop(i_tone,2) = FuncInput_InputData.time{1}(i_LastSampleTone);
end
TP_Tone_StartStop = TP_Tone_StartStop(1:34,:);

%Check if all tones are of equal length
%(if not then choose min length by deleting the last sample of longer trials)
minSeqLength_sec = min(TP_Tone_StartStop(:,2)-TP_Tone_StartStop(:,1));
for i_tone = 1:34
    if (TP_Tone_StartStop(i_tone,2) - TP_Tone_StartStop(i_tone,1)) > minSeqLength_sec
        TP_Tone_StartStop(i_tone,2) = ...
            FuncInput_InputData.time{1}(...
            find(FuncInput_InputData.time{1} == TP_Tone_StartStop(i_tone,2))-1);
    end
end

%Determine samples corresponding to TP
Sample_Tone_StartStop = NaN(34,2);
for i_tone = 1:34
    Sample_Tone_StartStop(i_tone,1) = ...
        find(FuncInput_InputData.time{1} == TP_Tone_StartStop(i_tone,1));
    Sample_Tone_StartStop(i_tone,2) = ...
        find(TP_Tone_StartStop(i_tone,2) == FuncInput_InputData.time{1});
end

%% 3. Calculate correlations between neural time series data for identical tone sequences for each sensor
%I.e., each identical sequence is presented 4 times in each TD, resulting
%in 6 pairs.
tic
disp(['-- Computing correlations between identical sequences for sub: ' sub ' --'])

%Determine unique trends/p*34 IDs
predIDs = unique(FuncInput_InputData.behav.stim.predID);
%Determine unique sequences
seqIDs = unique(FuncInput_InputData.behav.stim.uSeqID);

for i_predID = 1:length(predIDs)
    
    %Find all trials with sequences matching current predID
    ind_allpredIDtrials = find(...
        FuncInput_InputData.behav.stim.predID == predIDs(i_predID));
    %Find all useqIDs used in trials with current predID
    ind_allpredID_useqIDs = unique(...
        FuncInput_InputData.behav.stim.uSeqID(ind_allpredIDtrials));
    
    for i_seqID = 1:length(ind_allpredID_useqIDs)
        
        % extract data for trials at a specific p*34 and sequence ID
        filter = ...
            FuncInput_InputData.behav.stim.predID == predIDs(i_predID) & ...
            FuncInput_InputData.behav.stim.uSeqID == ind_allpredID_useqIDs(i_seqID);
        field_length = length(FuncInput_InputData.trial);
        
        d               = myft_trialfilter(FuncInput_InputData, filter, field_length);
        d.behav         = myft_trialfilter(d.behav, filter, field_length);
        d.behav.stim    = myft_trialfilter(d.behav.stim, filter, field_length);
        
        % now that the data for this tone sequence is selected,
        % conduct all possible correlations at each sensor...
        
        nTrials = length(d.trial);
        if nTrials > 1 %Across-trial correlation doesnt work for 1 trial
            nPerm   = nchoosek(nTrials, 2);
            
            corr_within{i_predID,i_seqID}.r                 = zeros(nSensors, nPerm);
            corr_within{i_predID,i_seqID}.p                 = zeros(nSensors, nPerm);
            corr_within{i_predID,i_seqID}.trial_pair        = zeros(nPerm, 2);
            corr_within{i_predID, i_seqID}.DataType_label   = FuncInput_DataType;
            corr_within{i_predID, i_seqID}.ToneDur_label    = FuncInput_ToneDur_text;
            
            ind = 0;
            %Correlate each trial with all other trials
            for i_trial1 = 1 : nTrials-1 %1st trial for correlation
                for i_trial2 = i_trial1+1 : nTrials %2nd trial for correlation
                    
                    ind = ind + 1;
                    
                    corr_within{i_predID,i_seqID}.trial_pair(ind, :) = [i_trial1,  i_trial2];
                    
                    for i_sensor = 1 : nSensors
                        
                        %compute correlation after removing the final tone
%                         %Option A: Use all samples (no window-based average)
%                         s1 = d.trial{i_trial1}...
%                             (i_sensor, Sample_Tone_StartStop(1,1) : ...
%                             Sample_Tone_StartStop(param.StimCorr_Seqrange,2));
%                         s2 = d.trial{i_trial2}...
%                             (i_sensor, Sample_Tone_StartStop(1,1) : ...
%                             Sample_Tone_StartStop(param.StimCorr_Seqrange,2));
                        
                        %Option B: Average neural signal across time windows (50 ms) 
                        win_size    = param.TW_Samples; 
                        s_per_win = (1/SampleFreq)*win_size;
                        win_overlap = 0;

                        windows = [Sample_Tone_StartStop(1,1) Sample_Tone_StartStop(1,1) + win_size];
                        while windows(end,end) < Sample_Tone_StartStop(param.StimCorr_Seqrange,2)                       
                            windows = ...
                                [windows; windows(end,:) + ...
                                (1 + win_size - win_overlap)];
                        end
                        if windows(end,end) > Sample_Tone_StartStop(param.StimCorr_Seqrange,2)
                            windows(end,:) = [];
                        end

                        s1 = [];
                        s2 = [];
                        for i_win = 1:size(windows,1)
                            ind_start = windows(i_win, 1);
                            ind_end   = windows(i_win, 2);
                            s1 = [s1 squeeze( mean( d.trial{i_trial1}(i_sensor, ind_start:ind_end, :), 2) )];
                            s2 = [s2 squeeze( mean( d.trial{i_trial2}(i_sensor, ind_start:ind_end, :), 2) )];
                        end    
                        
                        [r p] = corr(s1', s2', 'type', 'Pearson');
                        %Note: non-phase correlation computed here
                        
                        corr_within{i_predID,i_seqID}.r(i_sensor, ind) = r;
                        corr_within{i_predID,i_seqID}.p(i_sensor, ind) = p;
                        
                    end
                end
            end
        end
    end
end
disp(['-- Finished Computing correlations between identical sequences for sub: ' ...
    sub ' after ' num2str(round(toc/60,2)) ' min --'])

%% 4. Calculate correlations between neural time series data for similar tone sequences for each sensor
%For similar trials, use trial swith same p*34/trend, but remove those
%trials that have identical sequences. %I.e., each identical sequence is
%presented 4 times in each TD and has 16 similar sequences with which it
%can be correlated, resulting in 64 pairs.

tic
disp(['-- Computing correlations between similar sequences for sub: ' sub ' --'])

for i_predID = 1:length(predIDs)
    
    %Find all trials with sequences matching current predID
    ind_allpredIDtrials = find(...
        FuncInput_InputData.behav.stim.predID == predIDs(i_predID));
    %Find all useqIDs used in trials with current predID
    ind_allpredID_useqIDs = unique(...
        FuncInput_InputData.behav.stim.uSeqID(ind_allpredIDtrials));
    
    for i_seqID = 1:length(ind_allpredID_useqIDs)
        
        %For current seqID, determine all trials with identical p*34 and
        %identical seqID (i.e., identical trials). These will be
        %our 'within' variable, containing identical trials which we will
        %correlate with all possible similar trials
        filter = ...
            FuncInput_InputData.behav.stim.predID == predIDs(i_predID)...
            & FuncInput_InputData.behav.stim.uSeqID == ind_allpredID_useqIDs(i_seqID);
        
        field_length = length(FuncInput_InputData.trial);
        
        d_wtn               = myft_trialfilter(FuncInput_InputData, filter, field_length);
        d_wtn.behav         = myft_trialfilter(d_wtn.behav, filter, field_length);
        d_wtn.behav.stim    = myft_trialfilter(d_wtn.behav.stim, filter, field_length);
        
        %For current seqID, determine all trials with identical p*34, but
        %different seqID (i.e., only similar but non-identical sequences).
        %This will be our 'across' variable, containing all similar
        %(non-identical) trials for the current within variable
        filter = FuncInput_InputData.behav.stim.predID == predIDs(i_predID)...
            & FuncInput_InputData.behav.stim.uSeqID ~= ind_allpredID_useqIDs(i_seqID);
        
        field_length = length(FuncInput_InputData.trial);
        
        d_acr               = myft_trialfilter(FuncInput_InputData, filter, field_length);
        d_acr.behav         = myft_trialfilter(d_acr.behav, filter, field_length);
        d_acr.behav.stim    = myft_trialfilter(d_acr.behav.stim, filter, field_length);
        
        nTrials_wtn = length(d_wtn.trial);
        nTrials_acr = length(d_acr.trial);
        
        nPerm   = nTrials_wtn * nTrials_acr;
        
        corr_across{i_predID, i_seqID}.r                = zeros(nSensors, nPerm);
        corr_across{i_predID, i_seqID}.p                = zeros(nSensors, nPerm);
        corr_across{i_predID, i_seqID}.trial_pair       = zeros(nPerm, 2);
        corr_across{i_predID, i_seqID}.DataType_label   = FuncInput_DataType;
        corr_across{i_predID, i_seqID}.ToneDur_label    = FuncInput_ToneDur_text;
        
        ind = 0;
        for i_trial_wtn = 1 : nTrials_wtn
            for i_trial_acr = 1 : nTrials_acr
                
                ind = ind + 1;
                
                corr_across{i_predID, i_seqID}.trial_pair(ind, :) = ...
                    [i_trial_wtn,  i_trial_acr];
                
                for i_sensor = 1 : nSensors
                    
                    % compute correlation after removing the final tone
%                   %Option A: Use all samples (no window-based average)
%                     s1 = d_wtn.trial{i_trial_wtn}...
%                         (i_sensor, Sample_Tone_StartStop(1,1) : ...
%                         Sample_Tone_StartStop(param.StimCorr_Seqrange,2));
%                     s2 = d_acr.trial{i_trial_acr}...
%                         (i_sensor, Sample_Tone_StartStop(1,1) : ...
%                         Sample_Tone_StartStop(param.StimCorr_Seqrange,2));
                    
                    %Option B: Average neural signal across time windows (50 ms)
                    win_size    = param.TW_Samples;
                    s_per_win = (1/SampleFreq)*win_size;
                    win_overlap = 0;
                    
                    windows = [Sample_Tone_StartStop(1,1) Sample_Tone_StartStop(1,1) + win_size];
                    while windows(end,end) < Sample_Tone_StartStop(param.StimCorr_Seqrange,2)
                        windows = ...
                            [windows; windows(end,:) + ...
                            (1 + win_size - win_overlap)];
                    end
                    if windows(end,end) > Sample_Tone_StartStop(param.StimCorr_Seqrange,2)
                        windows(end,:) = [];
                    end
                    
                    s1 = [];
                    s2 = [];
                    for i_win = 1:size(windows,1)
                        ind_start = windows(i_win, 1);
                        ind_end   = windows(i_win, 2);
                        s1 = [s1 squeeze( mean( d_wtn.trial{i_trial_wtn}(i_sensor, ind_start:ind_end, :), 2) )];
                        s2 = [s2 squeeze( mean( d_acr.trial{i_trial_acr}(i_sensor, ind_start:ind_end, :), 2) )];
                    end
                    
                    [r p] = corr(s1', s2', 'type', 'Pearson');
                    
                    corr_across{i_predID, i_seqID}.r(i_sensor, ind) = r;
                    corr_across{i_predID, i_seqID}.p(i_sensor, ind) = p;
                    
                end
                
            end
        end        
    end
end

disp(['-- Finished Computing correlations between similar sequences for sub: ' ...
    sub ' after ' num2str(round(toc/60,2)) ' min --'])

%% 5) Statistically compare identical vs. similar correlations
%5.1 Aggregate Fisher-z-transformed correlation coefficients for identical vs similar
%for each sensor and then compare distributions via one-sample ttest (one-sided)
corr_ttest          = [];
corrcoeff_identical = [];
corrcoeff_similar   = [];

for i_sensor = 1:nSensors
    
    rs_within = [];
    rs_across = [];
    
    %For each predID-useqID combination, z-score identical and similar
    %correlation coefficients, average coefficients across trial pairs,
    %and aggregate them across predID-useqID for each sensor.
    for i_predID = 1:length(predIDs)
        for i_seqID = 1:length(corr_within)
            if ~isempty(corr_within{i_predID,i_seqID})
                rs_within = ...
                    [rs_within, mean(r2z(corr_within{i_predID,i_seqID}.r(i_sensor,:)))]; %Fisher-z-transform
                rs_across = ...
                    [rs_across, mean(r2z(corr_across{i_predID,i_seqID}.r(i_sensor,:)))];
            end
        end
    end
    
    %Compare identical vs. similar correlation coefficients for each electrode
%     [h, p, ci, stats] = ttest2(rs_within', rs_across');
    [h, p, ci, stats] = ttest(rs_within', rs_across', 'tail', 'right');
    
    %Store results for saving & plotting
    corr_ttest.t(i_sensor,1)    = stats.tstat;
    corr_ttest.df(i_sensor,1)   = stats.df;
    corr_ttest.sd(i_sensor,1)   = stats.sd;
    corr_ttest.p(i_sensor,1)    = p;
    corr_ttest.DataType_label   = FuncInput_DataType;
    corr_ttest.ToneDur_label    = FuncInput_ToneDur_text;
    
    corrcoeff_identical = [corrcoeff_identical; mean(rs_within)];
    corrcoeff_similar   = [corrcoeff_similar; mean(rs_across)];
    
end

%5.2 Aggregate Fisher-z-transformed correlation coefficients for identical vs similar
%for each sensor and then compare distributions via ANOVA, including stimulus condition factors

%Create empty proxies
zvals  = [];
predID = [];
seqID  = [];
Within = [];

for i_sensor = 1:nSensors
    
    zvals_i     = [];
    predID_i    = [];
    seqID_i     = [];
    Within_i    = [];
    
    for i_predID = 1:length(predIDs)
        
        %Find all trials with sequences matching current predID
        ind_allpredIDtrials = find(...
            FuncInput_InputData.behav.stim.predID == predIDs(i_predID));
        %Find all useqIDs used in trials with current predID
        ind_allpredID_useqIDs = unique(...
            FuncInput_InputData.behav.stim.uSeqID(ind_allpredIDtrials));
        
        for i_seqID = 1:length(ind_allpredID_useqIDs)
            
            if ~isempty(corr_within{i_predID,i_seqID})
                % within-stimulus correlations            
                rs_within = [corr_within{i_predID,i_seqID}.r(i_sensor,:)'] ;
                zvals_i = [zvals_i; r2z(rs_within)];

                %Markers
                predID_i    = [predID_i; predIDs(i_predID) * ones(size(rs_within))];
                seqID_i     = [seqID_i; ind_allpredID_useqIDs(i_seqID) * ones(size(rs_within))];
                Within_i    = [Within_i; ones(size(rs_within))];

                % across-stimulus correlations
                rs_across = [corr_across{i_predID,i_seqID}.r(i_sensor,:)'] ;
                zvals_i = [zvals_i; r2z( rs_across )];

                %Markers
                predID_i    = [predID_i; predIDs(i_predID) * ones(size(rs_across))];
                seqID_i     = [seqID_i; ind_allpredID_useqIDs(i_seqID) * ones(size(rs_across))];
                Within_i    = [Within_i; zeros(size(rs_across))];
            end
            
        end
    end
    
    zvals(i_sensor, :)  = zvals_i';
    Pred(i_sensor, :)   = predID_i';
    seqID(i_sensor, :)  = seqID_i';
    Within(i_sensor, :) = Within_i';
    
    [p, anova_tab{i_sensor}] = anovan(zvals_i, ...
        {predID_i seqID_i Within_i}, ...
        'model', 3, 'sstype', 2, 'display', 'off', ...
        'varnames', char('predID', 'seqID', 'Within'));
    
    anova_F(i_sensor,1) = anova_tab{i_sensor}{4,6};
    anova_p(i_sensor,1) = anova_tab{i_sensor}{4,7};
    
end

nFactors = 7;
for i_factor = 1 : nFactors
    factor_labels{i_factor} = anova_tab{1}{i_factor+1, 1};
end

%% 6. Save data
savefile = [path_dataoutput sub '_StimCorr_' FuncInput_DataType '_' ...
    FuncInput_ToneDur_text 'sTD.mat'];

save(savefile, ...
    'corr_within', 'corr_across', 'corrcoeff_identical', 'corrcoeff_similar', ...
    'corr_ttest', 'anova_F', 'anova_p', 'anova_tab', 'factor_labels', ...
    'zvals', 'Pred', 'seqID', 'Within', ...
    'SensorLabels', ...
    '-v7.3');

%% 7. Plot stats of comparison on surface brain and highlight sign. electrodes
%Plot summary figure containing:
%1) Correlation coefficient to identical sequences
%2) Correlation coefficient to similar sequences
%3) Statistic of indentical vs. similar comparison
%4) Table with anatomical labels for sign. electrodes

clear plot_struct

%7.0 Set up subplot structure
h = figure('visible','off'); %ensures figure doesn't pup during plotting
set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen

DimSubplot          = [2 2];
SubplotPosition     = [0 -0.05 0 0];
ColorbarPosition    = [0 0 0 0];
SizeFactor          = 4;
CounterSubplot      = 1;

sgtitle({...
    ['Across-trial correlation of neural data for identical vs. similar sequences'], ...
    [sub ' - ' FuncInput_DataType ' - ' FuncInput_ToneDur_text ' - (RH elecs projected on LH)']}, ...
    'Interpreter','none')

%Find index of each electrode label in cleaned data .elec struct
%and read out elec-specific MNI coordinates for plotting
Index_Elecs2ChanPos = NaN(1,nSensors);
for i_elec = 1:length(FuncInput_InputData.label)
    Index_Elecs2ChanPos(1,i_elec) = ...
        find(strcmp(FuncInput_InputData.label{i_elec}, ...
        FuncInput_InputData.elec.label));
end

plot_struct.coords = FuncInput_InputData.elec.chanpos(Index_Elecs2ChanPos,:);
%Project all electrodes on left hemisphere,
plot_struct.coords(:,1) = abs(plot_struct.coords(:,1)) * -1;

%7.1 Plot surface plot with correlation coefficient for identical sequences
plot_struct.dimord          = 'chan_time';
plot_struct.time            = 0;
plot_struct.sign_elecs      = logical(ones(length(corrcoeff_identical),1));
plot_struct.clims           = [-0.25 0.25]; %free symmetric scaling
plot_struct.chanSize        = ...
    ones(1,length(plot_struct.sign_elecs))*SizeFactor; %electrode size (arbitrary)
plot_struct.cmap            = 'jet';
plot_struct.textcolor_rgb   = [0 0 0];

for i_elec = 1:length(FuncInput_InputData.label)
    plot_struct.label{i_elec} = '';
end

plot_struct.avg         = corrcoeff_identical;

%Plot surface with highlighted sign electrodes
sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH...
    (plot_struct.coords, plot_struct.label, plot_struct.avg, plot_struct.sign_elecs,...
    plot_struct.chanSize, plot_struct.clims, plot_struct.cmap, plot_struct.textcolor_rgb, ...
    DimSubplot, CounterSubplot, SubplotPosition, ColorbarPosition, []);

sp_handle_surf{CounterSubplot} = sp_handle_surf_temp.L;
CounterSubplot = CounterSubplot + 1;
title('Identical sequences','FontSize',14)

%7.2 Plot surface plot with correlation coefficient for similar sequences
plot_struct.avg         = corrcoeff_similar;

%Plot surface with highlighted sign electrodes
sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH...
    (plot_struct.coords, plot_struct.label, plot_struct.avg, plot_struct.sign_elecs,...
    plot_struct.chanSize, plot_struct.clims, plot_struct.cmap, plot_struct.textcolor_rgb, ...
    DimSubplot, CounterSubplot, SubplotPosition, ColorbarPosition, []);

sp_handle_surf{CounterSubplot} = sp_handle_surf_temp.L;
CounterSubplot = CounterSubplot + 1;
title('Similar sequences','FontSize',14)

%Add colorbar
ColorbarPosition = [sp_handle_surf{CounterSubplot-1}.Position(1)-0.05 ...
    sp_handle_surf{CounterSubplot-1}.Position(2)+0.05 0 0];
h = colorbar;
h.Position(1) = ColorbarPosition(1); %sets colorbar to the right
h.Position(2) = ColorbarPosition(2); %sets colorbar higher
h.Position(4) = h.Position(4)*0.8; %makes colorbar shorter
h.Label.String = ['Pearsons r'];
h.FontSize = 14;
caxis(plot_struct.clims)

%7.3 Plot surface plot showing stats for identical vs. similar contrast
%Determine labels of sign. electrodes
SignElecs.index = find(corr_ttest.p < 0.05);
for i_signelecs = 1:length(SignElecs.index)
    temp_elecindex = find(strcmp(...
        FuncInput_InputData.label{SignElecs.index(i_signelecs)}, ...
        FuncInput_InputData.elec.label));
    
    SignElecs.anatlabel{i_signelecs,1} = ...
        FuncInput_InputData.elec.T1AnatLabel{temp_elecindex};
end

plot_struct.label           = FuncInput_InputData.label;
plot_struct.avg             = corr_ttest.t;
plot_struct.sign_elecs      = logical(corr_ttest.p < 0.05);
plot_struct.clims           = [0 5]; %free symmetric scaling
plot_struct.textcolor_rgb   = [0 0 0];

%Plot surface with highlighted sign electrodes
sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH...
    (plot_struct.coords, plot_struct.label, plot_struct.avg, plot_struct.sign_elecs,...
    plot_struct.chanSize, plot_struct.clims, plot_struct.cmap, plot_struct.textcolor_rgb, ...
    DimSubplot, CounterSubplot, SubplotPosition, ColorbarPosition, []);
title({['Identical vs. similar sequences (p < 0.05 uncorr)'], ...
    ['# sign. elecs / all elecs: ' ...
    num2str(length(find(corr_ttest.p < 0.05))) ' / ' num2str(nSensors)]}, ...
    'FontSize',14)

sp_handle_surf{CounterSubplot} = sp_handle_surf_temp.L;
CounterSubplot = CounterSubplot + 1;

%Add colorbar
ColorbarPosition = [sp_handle_surf{CounterSubplot-1}.Position(1)-0.05 ...
    sp_handle_surf{CounterSubplot-1}.Position(2)+0.02 0 0];
h = colorbar;
h.Position(1) = ColorbarPosition(1); %sets colorbar to the right
h.Position(2) = ColorbarPosition(2); %sets colorbar higher
h.Position(4) = h.Position(4)*0.8; %makes colorbar shorter
h.Label.String = ['t-value (identical vs. different)'];
h.FontSize = 14;
caxis(plot_struct.clims)

%7.4 Table with analysis & electrode information
subplot(DimSubplot(1), DimSubplot(2), 4)
textbox_info1 = {...
    [sub] ...
    ['# sign. elecs / all elecs: ' ...
    num2str(length(SignElecs.index)) ' / ' num2str(nSensors)] ...
    ''};
textbox_info2 = {};
if length(SignElecs.index) < 30 %one row
    for i_elec = 1:length(SignElecs.index)
        sel_elec = SignElecs.index(i_elec);
        textbox_info2{i_elec} = ...
            [FuncInput_InputData.label{sel_elec} ' = ' SignElecs.anatlabel{i_elec}];
    end
    textbox = [textbox_info1 textbox_info2];
    t = text(0, 0.5, 0, textbox, 'FontSize',12,'Interpreter','none');
else
    for i_elec = 1:25 %two rows
        sel_elec = SignElecs.index(i_elec);
        textbox_info2{i_elec} = ...
            [FuncInput_InputData.label{sel_elec} ' = ' SignElecs.anatlabel{i_elec}];
    end
    for i_elec = 26:length(SignElecs.index)
        sel_elec = SignElecs.index(i_elec);
        textbox_info3{i_elec} = ...
            [FuncInput_InputData.label{sel_elec} ' = ' SignElecs.anatlabel{i_elec}];
    end
    textbox = [textbox_info1 textbox_info2];
    t = text(-0.3, 0.5, 0, textbox, 'FontSize',12,'Interpreter','none');
    t2 = text(0.6, 1, 0, textbox_info3, 'FontSize',12,'Interpreter','none');
end
set(gca,'visible','off')

%7.5 Save Figure
if save_poststepFigs == 1
    filename     = ['StimCorr_Surf1HSignElec_' sub '_' ...
        FuncInput_DataType '_' FuncInput_ToneDur_text 'sTD.png'];
    figfile      = [path_fig filename];
    saveas(gcf, figfile, 'png'); %save png version
    close all;
end

end