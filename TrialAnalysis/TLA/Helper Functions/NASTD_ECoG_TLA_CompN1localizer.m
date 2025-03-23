function Index_N100SuperThresElecs = NASTD_ECoG_TLA_CompN1localizer...
    (sub, InputDataType, tonedur_label, ...
    N1_win, BL_win, ThresComp, STD_thresh, ...
    InputData, ...
    plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG)

%Aim: Compute timelocked response (N100) to each tone to find electrodes that
%react to tone presentation. Can be done based on ERFs computed for only
%first tone or all tones.(i.e., functional localizer)

%Method: Compute timelocked response in N100 time window corresponding to 
%early auditory processing for each electrode. Compare N100 ERF against 
%A) Prestimulus baseline activity + Threshold STD or 
%B) average ERF across whole tone duration. Select channels that show an 
%N100 response that is larger than comaprison metric + a priori defined 
%threshold 

%% 1. Determine start+end time points (s) for each tone
TP_Tone_StartStop = NaN(36,2);
TP_Tone_StartStop(1,1) = 0; %p1 set as t = 0 in trial definition
for i_tone = 2:36
    Dist = abs(InputData.time{1} - ((str2num(tonedur_label)*i_tone) - str2num(tonedur_label)));
    minDist = min(Dist);
    i_minDist = find(Dist == minDist);
    TP_Tone_StartStop(i_tone,1) = InputData.time{1}(i_minDist);
end
for i_tone = 1:35
    i_LastSampleTone = find(InputData.time{1} == TP_Tone_StartStop(i_tone + 1,1));
    TP_Tone_StartStop(i_tone,2) = InputData.time{1}(i_LastSampleTone);
end
TP_Tone_StartStop = TP_Tone_StartStop(1:34,:);

%Check if all tones are of equal length
%(if not then choose min length by deleting the last sample of longer trials)
minSeqLength_sec = min(TP_Tone_StartStop(:,2)-TP_Tone_StartStop(:,1));
for i_tone = 1:34
    if (TP_Tone_StartStop(i_tone,2) - TP_Tone_StartStop(i_tone,1)) > minSeqLength_sec
        TP_Tone_StartStop(i_tone,2) = ...
            InputData.time{1}(...
            find(InputData.time{1} == TP_Tone_StartStop(i_tone,2))-1);
    end
end

%Determine samples corresponding to TP
Sample_Tone_StartStop = NaN(34,2);
for i_tone = 1:34
    Sample_Tone_StartStop(i_tone,1) = ...
        find(InputData.time{1} == TP_Tone_StartStop(i_tone,1));
    Sample_Tone_StartStop(i_tone,2) = ...
        find(TP_Tone_StartStop(i_tone,2) == InputData.time{1});
end

%Determine N100 start/stop samples
N1_win_samples = ...
    [round(N1_win(1)*InputData.fsample) round(N1_win(2)*InputData.fsample)];

%% 2. Select current input data in FT-struct
InputData.trial = [];
InputData.trial = InputData.(InputDataType);

%% 3. Average timelocked activity across trials
ERF_avgTrials = ft_timelockanalysis([],InputData);

%% 4. Separate GAvg trial into tone-epochs and demean epochs
tone_range = 1:34; %which tones to include (all or first only)

for i_tone = tone_range
    ERF_perTone{i_tone} = ERF_avgTrials.avg...
        (:,Sample_Tone_StartStop(i_tone,1):Sample_Tone_StartStop(i_tone,2));
end

% %Visualize ERFs for each tone per electrode
% if plot_poststepFigs == 1
%     for i_elec = 1:length(InputData.label)
%         
%         figure; hold on;
%         set(gcf,'units','normalized','outerposition',[0 0 1 1])
%         
%         for i_tone = tone_range
%             subplot(6,6,i_tone)
%             plot(1:size(ERF_perTone{i_tone},2), ...
%                 ERF_perTone{i_tone}(i_elec,:), ...
%                 'LineWidth',2);
%             
%             hold on;
%             proxy = ones(1,length(N1_win_samples(1):N1_win_samples(2)));
%             area(N1_win_samples(1):N1_win_samples(2), ...
%                 proxy * min(ERF_perTone{i_tone}(i_elec,:)),...
%                 'basevalue',0,'FaceColor',[0.1, 0.1, 0.1],...
%                 'FaceAlpha', 0.5,'LineStyle','none');
%             area(N1_win_samples(1):N1_win_samples(2), ...
%                 proxy * max(ERF_perTone{i_tone}(i_elec,:)),...
%                 'basevalue',0,'FaceColor',[0.1, 0.1, 0.1],...
%                 'FaceAlpha', 0.5,'LineStyle','none');
%             
%             title(['Tone ' num2str(i_tone)])
%         end
%         
%         sgtitle(['ERF (elec: ' ERF_avgTrials.label{i_elec} ...
%             ') per tone ' num2str(tone_range(1)) '-' num2str(tone_range(end))])
%         
%         if save_poststepFigs == 1
%             path_fig = ([paths_NASTD_ECoG.ECoGdata_Timelocked ...
%                 '/N1localizer/' sub '/Figs/ERF_perTone/' ...
%                 InputDataType '/' tonedur_label 'msTD/']);
%             if (~exist(path_fig, 'dir')); mkdir(path_fig); end
%             
%             filename     = ['ERFperTone_' InputData.label{i_elec} ...
%                 '_'  sub '_TD' tonedur_label  's_' InputDataType '.png'];
%             figfile      = [path_fig filename];
%             saveas(gcf, [figfile], 'png'); %save png version
%             close
%         end
%     end
% end

%% 5. Average across all tones
tone_range = 1:34; %All tones

for i_tone = tone_range
    proxy_tonematrix(i_tone,:,:) = ERF_perTone{i_tone};
end
ERF_avgTone = squeeze(mean(proxy_tonematrix,1));
clear proxy*

% %Visualize ERF for first and averaged across all tones
% for i_elec = 1:length(InputData.label)
%     figure;
%     set(gcf,'units','normalized','outerposition',[0 0 1 1])
%     
%     subplot(1,2,1) %Add average across tones
%     plot(1:length(ERF_avgTone), ...
%         ERF_perTone{1}(i_elec,:), ...
%         'LineWidth',2);
%     title('First tone')
%     
%     hold on;
%     proxy = ones(1,length(N1_win_samples(1):N1_win_samples(2)));
%     area(N1_win_samples(1):N1_win_samples(2), ...
%         proxy * min(ERF_perTone{1}(i_elec,:)),...
%         'basevalue',0,'FaceColor',[0.1, 0.1, 0.1],...
%         'FaceAlpha', 0.5,'LineStyle','none');
%     area(N1_win_samples(1):N1_win_samples(2), ...
%         proxy * max(ERF_perTone{1}(i_elec,:)),...
%         'basevalue',0,'FaceColor',[0.1, 0.1, 0.1],...
%         'FaceAlpha', 0.5,'LineStyle','none');
%     
%     array_xtick = [0:10:length(ERF_perTone{1})];
%     for i_tick = 1:length(array_xtick)
%         labels_x{i_tick} = ...
%             num2str(round(array_xtick(i_tick)./InputData.fsample,2));
%     end
%     xticks(0:10:length(ERF_avgTone))
%     xticklabels(labels_x)
%     xlabel('Time [sec]')
%     xlim([0 length(ERF_avgTone)])
%     ylabel('Amplitude')
%     
%     subplot(1,2,2) %Add average across tones
%     plot(1:length(ERF_avgTone), ...
%         ERF_avgTone(i_elec,:), ...
%         'LineWidth',2);
%     title(['Averaged across tones (' ...
%         num2str(tone_range(1)) '-' num2str(tone_range(end)) ')'])
%     hold on
%     
%     hold on;
%     proxy = ones(1,length(N1_win_samples(1):N1_win_samples(2)));
%     area(N1_win_samples(1):N1_win_samples(2), ...
%         proxy * min(ERF_avgTone(i_elec,:)),...
%         'basevalue',0,'FaceColor',[0.1, 0.1, 0.1],...
%         'FaceAlpha', 0.5,'LineStyle','none');
%     area(N1_win_samples(1):N1_win_samples(2), ...
%         proxy * max(ERF_avgTone(i_elec,:)),...
%         'basevalue',0,'FaceColor',[0.1, 0.1, 0.1],...
%         'FaceAlpha', 0.5,'LineStyle','none');
%     
%     array_xtick = [0:10:length(ERF_avgTone)];
%     for i_tick = 1:length(array_xtick)
%         labels_x{i_tick} = ...
%             num2str(round(array_xtick(i_tick)./InputData.fsample,2));
%     end
%     xticks(0:10:length(ERF_avgTone))
%     xticklabels(labels_x)
%     xlabel('Time [sec]')
%     xlim([0 length(ERF_avgTone)])
%     ylabel('Amplitude')
%     
%     sgtitle(['ERF - elec: ' InputData.label{i_elec}]);
%     
%     if save_poststepFigs == 1
%         path_fig = ([paths_NASTD_ECoG.ECoGdata_Timelocked ...
%             '/N1localizer/' sub '/Figs/ERF_FirstvsAvgTones/' ...
%             InputDataType '/' tonedur_label 'msTD/']);
%         if (~exist(path_fig, 'dir')); mkdir(path_fig); end
%         
%         filename     = ['ERF_FirstvsAvgTones_' InputData.label{i_elec} ...
%             '_'  sub '_TD' tonedur_label  's_' InputDataType '.png'];
%         figfile      = [path_fig filename];
%         saveas(gcf, [figfile], 'png'); %save png version
%         close
%     end
% end

%% 6. For LF-ERFs, Square raw amplitude
%Don't do this for gamma envelope though, since those are only positive values
if strcmp(InputDataType, 'LP35Hz')
    for i_elec = 1:size(ERF_avgTone,1)
        for i_sample = 1:size(ERF_avgTone,2)
            ERF_avgTone(i_elec,i_sample) = ...
                ERF_avgTone(i_elec,i_sample)^2;
            ERF_firstTone(i_elec,i_sample) = ...
                ERF_perTone{1}(i_elec,i_sample)^2;
        end
    end
else
    for i_elec = 1:size(ERF_avgTone,1)
        for i_sample = 1:size(ERF_avgTone,2)
            ERF_firstTone(i_elec,i_sample) = ...
                ERF_perTone{1}(i_elec,i_sample);
        end
    end    
end

%% 7. Compare avg N100 ERF against threshold

%7.0 Decide which N1 marker to use (first tone only or avg across all tones)
Input_ERF = ERF_firstTone; 
Input_Label_N100ERF = 'T1';
disp('Tonal tracking based on first tone only')
% Input_ERF = ERF_avgTone;
% Input_Label_N100ERF = 'AllT';
% disp('Tonal tracking based on avg across all tones')

%7.1 Compute avg ERF activity during N100 window
avgERF_acrossSamples_N100win = ... %Across samples
    mean(Input_ERF(:,N1_win_samples(1):N1_win_samples(2)),2);
avgERF_acrossElecs_N100win = ... %Across electrodes
    mean(Input_ERF(:,N1_win_samples(1):N1_win_samples(2)),1);

%7.2A Option 1: Threshold determined per electrode for baseline
%Compute AVG and STD of ERF activity during prestim baseline (Per electrode, across baseline samples):
cfg = [];
cfg.latency = [BL_win(1) BL_win(end)];
temp_ERF_BLwin = ft_timelockanalysis(cfg,InputData);

if strcmp(InputDataType, 'LP35Hz')
    temp_ERF_BLwin.avg = temp_ERF_BLwin.avg.^2; %Square ERF for LFs
end

avgERF_acrossSamples_BLwin = mean(temp_ERF_BLwin.avg,2); %Mean across samples for each electrode
for i_elec = 1:length(temp_ERF_BLwin.label)
    stdERF_acrossSamples_BLwin(i_elec,1) = std(temp_ERF_BLwin.avg(i_elec,:)); %STD across samples for each channel
end

%7.2B Option 2: Threshold determined across electrodes for tone presentation window
%Compute AVG across tone presentation samples for each electrode and STD across electrodes
avgERF_acrossElec_Tonewin = mean(Input_ERF,2);%Mean across samples
stdERF_acrossElec_Tonewin = std(avgERF_acrossElec_Tonewin);%STD across electrodes

%7.3 Read out elecs whose avg N1 ERF is higher than its avg BL ERF - (STD_thresh * STD)
Index_SuperThresElec_BLwin = [];
Index_SuperThresElec_Tonewin = [];

for i_elec = 1:length(InputData.label)
    %Option 1: Threshold determined per electrode for baseline
    if avgERF_acrossSamples_N100win(i_elec) > ...
            avgERF_acrossSamples_BLwin(i_elec) +(STD_thresh * stdERF_acrossSamples_BLwin(i_elec))
        Index_SuperThresElec_BLwin = ...
            [Index_SuperThresElec_BLwin i_elec];
    end
    %Option 2: Threshold determined across electrodes for tone presentation window
    if avgERF_acrossSamples_N100win(i_elec) > ...
            avgERF_acrossElec_Tonewin(i_elec) + (STD_thresh * stdERF_acrossElec_Tonewin)
        Index_SuperThresElec_Tonewin = ...
            [Index_SuperThresElec_Tonewin i_elec];
    end
end
InputData.label(Index_SuperThresElec_BLwin)
InputData.label(Index_SuperThresElec_Tonewin)

%Visually compare threshold options
figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1])
plot(avgERF_acrossSamples_N100win,'r-', 'LineWidth', 2)
hold on;
shadedErrorBar([], avgERF_acrossSamples_BLwin, ...
    (STD_thresh * stdERF_acrossSamples_BLwin), ...
    'lineprops',{'-b','linewidth',2},'transparent',1)
hold on;
shadedErrorBar([], avgERF_acrossElec_Tonewin, ...
    ones(1,length(avgERF_acrossElec_Tonewin)) * ...
    (STD_thresh * stdERF_acrossElec_Tonewin), ...
    'lineprops',{'-k','linewidth',2},'transparent',1)
legend('N100 ERF', 'Avg ERF Baseline', 'Avg ERF tone')
xlabel('Electrode number');
ylabel('Squared ERF amplitude')
temp_ymax = max([max(avgERF_acrossSamples_N100win) ...
    max(avgERF_acrossSamples_BLwin) ...
    max(avgERF_acrossElec_Tonewin)]);
% ylim([0 temp_ymax])

%Add text inset listing superthres elecs
if ~isempty(Index_SuperThresElec_BLwin)
    text_SuperThresElec_BLwin = ['Superthres elecs (BLwin):'];
    for i_elec = Index_SuperThresElec_BLwin
        text_SuperThresElec_BLwin = ...
            strcat(text_SuperThresElec_BLwin , {' '}, InputData.label{i_elec});
    end
text(1,max(ylim)*0.9, ...
    text_SuperThresElec_BLwin, 'FontSize',12,'Color','b')    
end
if ~isempty(Index_SuperThresElec_Tonewin)
    text_SuperThresElec_Tonewin = ['Superthres elecs (Tonewin):'];
    for i_elec = Index_SuperThresElec_Tonewin
        text_SuperThresElec_Tonewin  = ...
            strcat(text_SuperThresElec_Tonewin , {' '}, InputData.label{i_elec});
    end
text(1,max(ylim)*0.8, ...
    text_SuperThresElec_Tonewin, 'FontSize',12,'Color','k')
end

title([sub ' - ' InputDataType ' - ' tonedur_label 's - ' ...
    'N100ERF ' Input_Label_N100ERF ...
    ' - Superthres Elecs (' ...
    ThresComp ' - ' num2str(STD_thresh) 'STD'], ...
    'Interpreter', 'none')

if save_poststepFigs == 1
    path_fig = ([paths_NASTD_ECoG.ECoGdata_Timelocked ...
        '/N1localizer/' sub '/Figs/SuperThres_Elecs/' ...
        InputDataType '/' tonedur_label 'msTD/']);
    if (~exist(path_fig, 'dir')); mkdir(path_fig); end
    
    filename     = ['ThresComp_' ThresComp num2str(STD_thresh) 'STD_N100ERF_'  ...
        sub '_' InputDataType '_' tonedur_label '.png'];
    figfile      = [path_fig filename];
    saveas(gcf, [figfile], 'png'); %save png version
    close
end

%% 8 Plot summary figure for each supra-thresh electrode
if strcmp(ThresComp,'BL')
    Index_SuperThresElec_Selected = Index_SuperThresElec_BLwin;
    Threshold_Selected = avgERF_acrossSamples_BLwin(Index_SuperThresElec_Selected) ...
                + (STD_thresh * stdERF_acrossSamples_BLwin(Index_SuperThresElec_Selected));
elseif strcmp(ThresComp,'AvgTone')
    Index_SuperThresElec_Selected = Index_SuperThresElec_Tonewin;
    Threshold_Selected = avgERF_acrossElec_Tonewin(Index_SuperThresElec_Selected) ...
                + (STD_thresh * stdERF_acrossElec_Tonewin);
end        
        
if plot_poststepFigs == 1
    for i_elec = 1:length(Index_SuperThresElec_Selected)
        
        %set up marker to shade N100 window
        marker_N1_win = zeros(1,size(Input_ERF,2));
        marker_N1_win(N1_win_samples(1):N1_win_samples(2)) ...
            = max(Input_ERF(Index_SuperThresElec_Selected(i_elec),:));
        
        h = figure;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        
        %8.1 ERF (squared) for first and averaged across all tones
        s1 = subplot(2,2,1);
        plot(1:size(Input_ERF,2),... %ERF (squared) for first tone
            Input_ERF(Index_SuperThresElec_Selected(i_elec),:),'b','LineWidth',2);
        hold on
        plot(1:size(ERF_firstTone,2),...
            ones(1,size(Input_ERF,2)) * ...
            (Threshold_Selected(i_elec)), ...
            'k--');        
        hold on
        area(1:size(Input_ERF,2), ...%Mark N100 window
            marker_N1_win,... 
            'basevalue',0,'FaceColor',[0.1, 0.1, 0.1], ...
            'FaceAlpha', 0.5,'LineStyle','none');
        
        legend(...
            'ERF',...
            ['Thresh (' num2str(STD_thresh) '* STD)'],...
            ['N1 window (' num2str(N1_win(1)) '-' num2str(N1_win(2)) 's)'], ...
            'Interpreter','none')
        
        xlim([0 size(ERF_avgTone,2)])
        xticks([0 (size(ERF_avgTone,2)/4)*1 ...
            (size(ERF_avgTone,2)/4)*2 ...
            (size(ERF_avgTone,2)/4)*3 ...
            (size(ERF_avgTone,2))])
        xlabel('Time')
        ylabel('Amplitude (squared)')
        
        if length(InputData.trial{1})/InputData.fsample - 2 < 10
            xticklabels({'0 ms', '50 ms', '100 ms', '150 ms', '200 ms'})
        else
            xticklabels({'0 ms', '100 ms', '200 ms', '300 ms', '400 ms'})
        end
        title(['ERF used for threshold determination (' ...
            Input_Label_N100ERF ')'], 'Interpreter', 'none');
        
        %8.2 Plot Suprathresh electrode highlighted on MNI template brain
        addpath(paths_NASTD_ECoG.Freesurfer);
        
        [lhvtx,lhfaces]=read_surf('lh.pial'); %pial surface - 3D surface file used to estimate cortical gray matter volume
        [rhvtx,rhfaces]=read_surf('rh.pial');
        t.faces = [lhfaces+1 ; rhfaces+1+length(lhvtx)];
        t.vertices = [lhvtx ; rhvtx];
        t.EdgeColor = 'none';
        t.faceColor = 'interp';
        t.facevertexCData = 0.9*ones(length(t.vertices),3);
        t.facealpha = 1;
        
        s2 = subplot(2,2,2)
        
        patch(t);
        hold on
        h = rotate3d;
        h.ActionPostCallback = @RotationCallback;
        h.Enable = 'on';
        
        axis square;
        axis equal;
        
        view([-90 20])
        c = camlight('headlight','infinite');
        lighting gouraud
        material dull;
        axis off
        
        coords = InputData.elec.chanpos(...
            find(strcmp(InputData.label{Index_SuperThresElec_Selected(i_elec)}, ...
            InputData.elec.label)),...
            :);        
        nearestVertex = nearestVertices(t.vertices,coords); %uses subfunction nearestVertices to find nearest vertex for each chan position
        [xs,ys,zs] = sphere;
        temp = surf2patch(xs,ys,zs,ones(size(zs)));
        counter = 0;
        electrodePatch = struct;
        
        numChans = 1;
        chanSize = 5;
        clims = [-3 +3]; %fixed scaling
        cmap = 'jet';
        colormap(cmap);
        set(gca,'clim',clims);
        vals = 2;
        
        electrodePatch.faces = zeros(numChans * size(temp.faces,1),4);
        electrodePatch.vertices = zeros(numChans * size(temp.vertices,1),3);
        electrodePatch.facevertexcdata = zeros(numChans * size(temp.facevertexcdata,1),1);
        electrodePatch.edgecolor = 'none';
        electrodePatch.facecolor = 'interp';
        
        for cn = 1:numChans
            x = (chanSize(cn)*xs) + t.vertices(nearestVertex(cn),1);
            y = (chanSize(cn)*ys) + t.vertices(nearestVertex(cn),2);
            z = (chanSize(cn)*zs) + t.vertices(nearestVertex(cn),3);
            temp = surf2patch(x,y,z,vals * ones(size(z)));
            faceInds = 1 + (cn-1) * size(temp.faces,1):cn * size(temp.faces,1);
            verticesInds =  1 + (cn-1) * size(temp.vertices,1):cn * size(temp.vertices,1);
            electrodePatch.faces(faceInds,:) = temp.faces + verticesInds(1) - 1;
            electrodePatch.vertices(verticesInds,:) = temp.vertices;
            electrodePatch.facevertexcdata(verticesInds,:) = temp.facevertexcdata;
        end
        patch(electrodePatch);
        material dull;
        lighting gouraud
        
        axis([-100 100 -115 85 -100 100]);
        title('Suprathres electrode location on MNI brain')
        
        %8.3 ERF avg across trials, but for each tone
        s3 = subplot(2,2,3);
        plot(1:length(ERF_avgTrials.avg), ... %ERF trace
            ERF_avgTrials.avg(Index_SuperThresElec_Selected(i_elec),:),...
            'b','LineWidth',2)
        hold on;
        plot(1:length(ERF_avgTrials.avg),zeros(1,length(ERF_avgTrials.avg)), 'k-');
        
        for i_tone = 1:length(Sample_Tone_StartStop)
            hold on; ... %Tone time window marks
                plot([Sample_Tone_StartStop(i_tone,1) Sample_Tone_StartStop(i_tone,1)],...
                [min(ERF_avgTrials.avg(Index_SuperThresElec_Selected(i_elec),:)) ...
                max(ERF_avgTrials.avg(Index_SuperThresElec_Selected(i_elec),:))], 'k-');
        end
        xticks(Sample_Tone_StartStop(:,1))
        xticklabels(num2str(round(TP_Tone_StartStop(:,1),2)))
        xlim([find(ERF_avgTrials.time == -0.5) Sample_Tone_StartStop(end,2)+0.5*512])
        s3.XTickLabelRotation = 270;
        ylim([min(ERF_avgTrials.avg(Index_SuperThresElec_Selected(i_elec),...
            find(ERF_avgTrials.time == -0.5):Sample_Tone_StartStop(end,2)+0.5*512))*0.9 ...
            max(ERF_avgTrials.avg(Index_SuperThresElec_Selected(i_elec),...
            find(ERF_avgTrials.time == -0.5):Sample_Tone_StartStop(end,2)+0.5*512))*1.1])
        xlabel('Time [s]')
        ylabel('Amplitude')
        
        title('ERF across tone sequence (averaged across trials)')
        
        pos3 = get(s3, 'Position');
        new_pos3 = pos3 + [0 0 0.4 0];
        set(s3, 'Position', new_pos3);
        
        %8.4 Save figure
        sgtitle([InputData.label{Index_SuperThresElec_Selected(i_elec)} ' - ' ...
            sub ' (' InputDataType '; ' tonedur_label 'TD) - '...
            Input_Label_N100ERF ' ERF - SuperThresh (' ThresComp '-' num2str(STD_thresh) 'STD)'], ...
            'Interpreter', 'none')
        
        if save_poststepFigs == 1
            path_fig = ([paths_NASTD_ECoG.ECoGdata_Timelocked ...
            '/N1localizer/' sub '/Figs/SuperThres_Elecs/' ...
            InputDataType '/' tonedur_label 'msTD/']);    
            if (~exist(path_fig, 'dir')); mkdir(path_fig); end
            
            filename = [InputData.label{Index_SuperThresElec_Selected(i_elec)} '_'...
                Input_Label_N100ERF 'ERF_SuperThres' ThresComp  num2str(STD_thresh) 'STD_'  ...
                sub '_' InputDataType '_TD' tonedur_label  's.png'];
            figfile      = [path_fig filename];
            saveas(gcf, [figfile], 'png'); %save png version
            close
        end
        
    end
end

%% 9 Provide output
disp(['N100 Superthres elecs (' InputDataType ', ' Input_Label_N100ERF ', ' ...
    num2str(STD_thresh) 'STD, ' tonedur_label 'TD):'])

    Index_N100SuperThresElecs.Info.N100ERF = Input_Label_N100ERF;
    Index_N100SuperThresElecs.Info.ThresComp = [ThresComp num2str(STD_thresh)];   
 
if ~isempty(Index_SuperThresElec_Selected)
    for i_elec = 1:length(Index_SuperThresElec_Selected)

        Index_N100SuperThresElecs.Label{i_elec} = ...
            InputData.label{Index_SuperThresElec_Selected(i_elec)};
        Index_N100SuperThresElecs.Index{i_elec} = ...
            Index_SuperThresElec_Selected(i_elec);    

        disp(InputData.label{Index_SuperThresElec_Selected(i_elec)})
    end
else
    Index_N100SuperThresElecs.Label = [];
    Index_N100SuperThresElecs.Index = [];
    disp('No N100 SuperThres electrodes found.')
end