function NASTD_ECoG_TLA_PlotSignPredp34...
    (sub, i_TD, fun_inputData, i_inputData, fun_inputStat, ...
    fun_TimelockData_perTDPredp34, fun_TimelockStat_Predp34_perTD, fun_Predp34DiffElec, ...
    save_fig, paths_NASTD_ECoG)

%% 1. Setup
%1.1 Set up figure and plotting param
figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1])
hold on
subplot_position = 0;

num_plots = length(fun_Predp34DiffElec.Index)+2;

ToneDur_text = {'0.2' '0.4'};
label_Predp34 = [-1 0 1];% low, medium, high;
  
%1.2 relabel input data
if strcmp(fun_inputData,'LF')
    inputData_Avg = 'Avg_LF';
    inputData_STD = 'STD_LF';
elseif strcmp(fun_inputData,'GammaAmp')
    inputData_Avg = 'Avg_GammaAmp';
    inputData_STD = 'STD_GammaAmp';    
elseif strcmp(fun_inputData,'LogGammaAmp')
    inputData_Avg = 'Avg_LogGammaAmp';
    inputData_STD = 'STD_LogGammaAmp';
end

%1.3. Read out TP and samples corresponding to tone presentation and response window 
if i_TD == 1
    ToneDur_text = '0.2';
    ToneDur_title = '200msTD';
    
elseif i_TD == 2
    ToneDur_text = '0.4';
    ToneDur_title = '400msTD';
end

TP_Tone_StartStop = NaN(35,2);
TP_Tone_StartStop(1,1) = 0;
for i_tone = 2:35
    Dist = abs(fun_TimelockData_perTDPredp34{1}.time - ((str2num(ToneDur_text)*i_tone)-str2num(ToneDur_text)));
    minDist = min(Dist);
    i_minDist = find(Dist == minDist);
    TP_Tone_StartStop(i_tone,1) = fun_TimelockData_perTDPredp34{1}.time(i_minDist);
end
for i_tone = 1:34
    i_LastSampleTone = find(fun_TimelockData_perTDPredp34{1}.time == TP_Tone_StartStop(i_tone+1,1));
    TP_Tone_StartStop(i_tone,2) = fun_TimelockData_perTDPredp34{1}.time(i_LastSampleTone);
end
TP_Tone_StartStop(35,1) = TP_Tone_StartStop(34,2)+0.4;
TP_Tone_StartStop(35,2) = NaN;

%check if all tones are of equal length, if not then choose min length
%by deleting the last sample
minSeqLength_sec = min(TP_Tone_StartStop(1:34,2)-TP_Tone_StartStop(1:34,1));
for i_tone = 1:34
    if TP_Tone_StartStop(i_tone,2)-TP_Tone_StartStop(i_tone,1) > minSeqLength_sec
        TP_Tone_StartStop(i_tone,2) = ...
            fun_TimelockData_perTDPredp34{1}.time(find(fun_TimelockData_perTDPredp34{1}.time == ...
            TP_Tone_StartStop(i_tone,2))-1);
    end
end

Sample_Tone_StartStop = NaN(35,2);
for i_tone = 1:34
    Sample_Tone_StartStop(i_tone,1) = ...
        find(fun_TimelockData_perTDPredp34{1}.time == TP_Tone_StartStop(i_tone,1));
    Sample_Tone_StartStop(i_tone,2) = ...
        find(TP_Tone_StartStop(i_tone,2) == fun_TimelockData_perTDPredp34{1}.time);
end
    Sample_Tone_StartStop(35,1) = ...
        nearest(fun_TimelockData_perTDPredp34{1}.time,TP_Tone_StartStop(35,1));

    
%1.4 Compute min/max vals across conditions and TD for common y-axis    
for i_chan = 1:length(fun_Predp34DiffElec.Index)
    for i_Predp34 = 1:length(label_Predp34)
        minAmp_perPredp34(i_chan,i_Predp34) = min(fun_TimelockData_perTDPredp34{i_Predp34}.(inputData_Avg)...
            (fun_Predp34DiffElec.Index(i_chan),Sample_Tone_StartStop(1):Sample_Tone_StartStop(34,2)));
        maxAmp_perPredp34(i_chan,i_Predp34) = max(fun_TimelockData_perTDPredp34{i_Predp34}.(inputData_Avg)...
            (fun_Predp34DiffElec.Index(i_chan),Sample_Tone_StartStop(1):Sample_Tone_StartStop(34,2)));
    end
    
    if strcmp(fun_inputData,'LF') %For LF, round up to next decimal
        minAmp(i_chan) = floor(floor(min(minAmp_perPredp34(i_chan,:)))*0.1)*10;
        maxAmp(i_chan) = ceil(ceil(max(maxAmp_perPredp34(i_chan,:)))*0.1)*10;
    else %For Gamma, don't round up
        minAmp(i_chan) = min(minAmp_perPredp34(i_chan,:))+(0.1*min(minAmp_perPredp34(i_chan,:)));
        maxAmp(i_chan) = max(maxAmp_perPredp34(i_chan,:))+(0.1*max(maxAmp_perPredp34(i_chan,:)));
    end
end
%% 2 Plot
%2.1 Plot surface with sign. electrodes
coords = fun_Predp34DiffElec.coords;
for i_Chan = 1:length(fun_Predp34DiffElec.Index)
    vals(1,i_Chan) = ...
        sum(fun_TimelockStat_Predp34_perTD{fun_Predp34DiffElec.Index(i_Chan)}{i_TD}{i_inputData}.mask)./...
        length(fun_TimelockStat_Predp34_perTD{fun_Predp34DiffElec.Index(i_Chan)}{i_TD}{i_inputData}.mask);
    chanSize(i_Chan) = 3;
end
clims = [min(vals),max(vals)]; %colormap limits            
cmap = 'parula';

%2.1.1 Left hemisphere
subplot_position = subplot_position + 1;
sp_handle1 = subplot(ceil(num_plots/2),2,subplot_position);
pos1 = get(sp_handle1,'Position');
newpos1 = pos1 + [0 0 0 0.15];
set(sp_handle1,'Position',newpos1)

addpath('/isilon/LFMI/VMdrive/Thomas/toolboxes/fieldtrip-20190314/external/freesurfer/');

[lhvtx,lhfaces]=read_surf('lh.pial'); %pial surface - 3D surface file used to estimate cortical gray matter volume
[rhvtx,rhfaces]=read_surf('rh.pial');
t.faces = [lhfaces+1 ; rhfaces+1+length(lhvtx)];
t.vertices = [lhvtx ; rhvtx];
t.EdgeColor = 'none';
t.faceColor = 'interp';
t.facevertexCData = 0.9*ones(length(t.vertices),3);
t.facealpha = 1;

patch(t);
hold on
h = rotate3d;
h.ActionPostCallback = @RotationCallback;
h.Enable = 'on';

axis square;
axis equal;

view([270,0])
c = camlight('headlight','infinite');
lighting gouraud
material dull;
axis off

colormap(cmap); 
h = colorbar;
ylabel(h, 'Ratio sign. samples / all tested samples')
set(gca,'clim',clims);

nearestVertex = nearestVertices(t.vertices,coords); %uses subfunction nearestVertices to find nearest vertex for each chan position
[xs,ys,zs] = sphere;
temp = surf2patch(xs,ys,zs,ones(size(zs)));
electrodePatch = struct;
numChans = length(vals);

electrodePatch.faces = zeros(numChans * size(temp.faces,1),4);
electrodePatch.vertices = zeros(numChans * size(temp.vertices,1),3);
electrodePatch.facevertexcdata = zeros(numChans * size(temp.facevertexcdata,1),1);
electrodePatch.edgecolor = 'none';
electrodePatch.facecolor = 'interp';

for cn = 1:numChans
    x = (chanSize(cn)*xs) + t.vertices(nearestVertex(cn),1);
    y = (chanSize(cn)*ys) + t.vertices(nearestVertex(cn),2);
    z = (chanSize(cn)*zs) + t.vertices(nearestVertex(cn),3);
    temp = surf2patch(x,y,z,vals(cn) * ones(size(z)));
    faceInds = 1 + (cn-1) * size(temp.faces,1):cn * size(temp.faces,1);
    verticesInds =  1 + (cn-1) * size(temp.vertices,1):cn * size(temp.vertices,1);
    electrodePatch.faces(faceInds,:) = temp.faces + verticesInds(1) - 1;
    electrodePatch.vertices(verticesInds,:) = temp.vertices;
    electrodePatch.facevertexcdata(verticesInds,:) = temp.facevertexcdata;
    
    text(t.vertices(nearestVertex(cn),1)-25,t.vertices(nearestVertex(cn),2),t.vertices(nearestVertex(cn),3)+5,...
        fun_Predp34DiffElec.Label{cn},...
        'FontSize',10)
end

patch(electrodePatch);
material dull;
lighting gouraud

% axis([-100 100 -115 85 -100 100]);
axis([-80 80 -120 70 -50 80]);

title('Left Hemisphere');
  
%2.1.2 Right hemisphere
subplot_position = subplot_position + 1;
sp_handle2 = subplot(ceil(num_plots/2),2,subplot_position);
pos2 = get(sp_handle2,'Position');
newpos2 = pos2 + [0 0 0 0.15];
set(sp_handle2,'Position',newpos2)

addpath('/isilon/LFMI/VMdrive/Thomas/toolboxes/fieldtrip-20190314/external/freesurfer/');

[lhvtx,lhfaces]=read_surf('lh.pial'); %pial surface - 3D surface file used to estimate cortical gray matter volume
[rhvtx,rhfaces]=read_surf('rh.pial');
t.faces = [lhfaces+1 ; rhfaces+1+length(lhvtx)];
t.vertices = [lhvtx ; rhvtx];
t.EdgeColor = 'none';
t.faceColor = 'interp';
t.facevertexCData = 0.9*ones(length(t.vertices),3);
t.facealpha = 1;

patch(t);
hold on
h = rotate3d;
h.ActionPostCallback = @RotationCallback;
h.Enable = 'on';

axis square;
axis equal;

view([90,0])
c = camlight('headlight','infinite');
lighting gouraud
material dull;
axis off

colormap(cmap); colorbar off
set(gca,'clim',clims);

nearestVertex = nearestVertices(t.vertices,coords); %uses subfunction nearestVertices to find nearest vertex for each chan position
[xs,ys,zs] = sphere;
temp = surf2patch(xs,ys,zs,ones(size(zs)));
electrodePatch = struct;
numChans = length(vals);

electrodePatch.faces = zeros(numChans * size(temp.faces,1),4);
electrodePatch.vertices = zeros(numChans * size(temp.vertices,1),3);
electrodePatch.facevertexcdata = zeros(numChans * size(temp.facevertexcdata,1),1);
electrodePatch.edgecolor = 'none';
electrodePatch.facecolor = 'interp';

for cn = 1:numChans
    x = (chanSize(cn)*xs) + t.vertices(nearestVertex(cn),1);
    y = (chanSize(cn)*ys) + t.vertices(nearestVertex(cn),2);
    z = (chanSize(cn)*zs) + t.vertices(nearestVertex(cn),3);
    temp = surf2patch(x,y,z,vals(cn) * ones(size(z)));
    faceInds = 1 + (cn-1) * size(temp.faces,1):cn * size(temp.faces,1);
    verticesInds =  1 + (cn-1) * size(temp.vertices,1):cn * size(temp.vertices,1);
    electrodePatch.faces(faceInds,:) = temp.faces + verticesInds(1) - 1;
    electrodePatch.vertices(verticesInds,:) = temp.vertices;
    electrodePatch.facevertexcdata(verticesInds,:) = temp.facevertexcdata;
    
    text(t.vertices(nearestVertex(cn),1)+25,t.vertices(nearestVertex(cn),2),t.vertices(nearestVertex(cn),3)+5,...
        fun_Predp34DiffElec.Label{cn},...
        'FontSize',10)
end

patch(electrodePatch);
material dull;
lighting gouraud

% axis([-100 100 -115 85 -100 100]);
axis([-80 80 -120 70 -50 80]);
title('Right Hemisphere');
        
    
%3.2 Plot timelocked p*34 traces with sign. different samples
for i_chan = fun_Predp34DiffElec.Index
    
    subplot_position = subplot_position + 1;
    subplot(ceil(num_plots/2),2,subplot_position);
    hold on;
    
    %Vertical line for each tone/event
    Array_Tone_Pos = zeros(1,length(fun_TimelockData_perTDPredp34{1}.(inputData_Avg)));
    Array_Tone_Neg = zeros(1,length(fun_TimelockData_perTDPredp34{1}.(inputData_Avg)));
    
    for i_SampleToneStart = [Sample_Tone_StartStop(1:34,1)' Sample_Tone_StartStop(34,2) Sample_Tone_StartStop(35,1)]
        Array_Tone_Pos(i_SampleToneStart) = maxAmp(subplot_position-2);
        Array_Tone_Neg(i_SampleToneStart) = minAmp(subplot_position-2);
    end
    
    line(1:length(Array_Tone_Pos),...
        Array_Tone_Pos,...
        'color',[0.8 0.8 0.8],'LineWidth',0.01) %plot mark for each tone-start-sample
    line(1:length(Array_Tone_Neg),...
        Array_Tone_Neg,...
        'color',[0.8 0.8 0.8],'LineWidth',0.01) %plot mark for each tone-start-sample
    
    Labels_ToneStart = {'T1','','','','','','','','','T10',...
        '','','','','','','','','','T20',...
        '','','','','','','','','','',...
        '','','T33(STABLE)','','', 'RespDisp'};
    
    
    plot_param.legend_label = {'predp34 = low','predp34 =  med','predp34 =  high'};
    plot_param.color = {[0, 0.28, 0.73],[1, 0.75, 0],[0.69, 0, 0.16]};
    plot_param.linestyle = {'-','-','-'};
    plot_param.lineThicknessToneMark = 2;
        
    %Low p*34
    plot(1:length(fun_TimelockData_perTDPredp34{1}.(inputData_Avg)), ...
        fun_TimelockData_perTDPredp34{1}.(inputData_Avg)(i_chan,:),...
        'color',plot_param.color{1},'LineStyle',plot_param.linestyle{1},'LineWidth',plot_param.lineThicknessToneMark)
%     plot(1:length(fun_TimelockData_perTDPredp34{1}.(inputData_Avg)), ...
%         zscore(fun_TimelockData_perTDPredp34{1}.(inputData_Avg)(i_chan,:)),...
%         'color',plot_param.color{1},'LineStyle',plot_param.linestyle{1},'LineWidth',plot_param.lineThicknessToneMark)
    hold on;
    
    %Med p*34
    plot(1:length(fun_TimelockData_perTDPredp34{2}.(inputData_Avg)), ...
        fun_TimelockData_perTDPredp34{2}.(inputData_Avg)(i_chan,:),...
        'color',plot_param.color{2},'LineStyle',plot_param.linestyle{2},'LineWidth',plot_param.lineThicknessToneMark)
%     plot(1:length(fun_TimelockData_perTDPredp34{2}.(inputData_Avg)), ...
%         zscore(fun_TimelockData_perTDPredp34{2}.(inputData_Avg)(i_chan,:)),...
%         'color',plot_param.color{2},'LineStyle',plot_param.linestyle{2},'LineWidth',plot_param.lineThicknessToneMark)
    hold on;    
    
    %High p*34
    plot(1:length(fun_TimelockData_perTDPredp34{3}.(inputData_Avg)), ...
        fun_TimelockData_perTDPredp34{3}.(inputData_Avg)(i_chan,:),...
        'color',plot_param.color{3},'LineStyle',plot_param.linestyle{3},'LineWidth',plot_param.lineThicknessToneMark)
%     plot(1:length(fun_TimelockData_perTDPredp34{3}.(inputData_Avg)), ...
%         zscore(fun_TimelockData_perTDPredp34{3}.(inputData_Avg)(i_chan,:)),...
%         'color',plot_param.color{3},'LineStyle',plot_param.linestyle{3},'LineWidth',plot_param.lineThicknessToneMark)
    hold on;
    
    %Sign. samples
    area(1:length(fun_TimelockData_perTDPredp34{1}.(inputData_Avg)),...
        fun_TimelockStat_Predp34_perTD{i_chan}{i_TD}{i_inputData}.plot_statp*maxAmp(subplot_position-2),...
        'basevalue',0,'FaceColor',[0.3, 0.3, 0.3],'FaceAlpha', 0.5,'LineStyle','none');
    hold on;
    area(1:length(fun_TimelockData_perTDPredp34{1}.(inputData_Avg)),...
        fun_TimelockStat_Predp34_perTD{i_chan}{i_TD}{i_inputData}.plot_statp*minAmp(subplot_position-2),...
        'basevalue',0,'FaceColor',[0.3, 0.3, 0.3],'FaceAlpha', 0.8,'LineStyle','none');
    hold on;    
    
    %x-axis shows tones [tone number]
    xlim([1  length(fun_TimelockData_perTDPredp34{1}.(inputData_Avg))])
    ylim([minAmp(subplot_position-2) maxAmp(subplot_position-2)])
%     minAmpElec = min([min(fun_TimelockData_perTDPredp34{1}.(inputData_Avg)(i_chan,Sample_Tone_StartStop(1):Sample_Tone_StartStop(34,2))),...
%         min(fun_TimelockData_perTDPredp34{2}.(inputData_Avg)(i_chan,Sample_Tone_StartStop(1):Sample_Tone_StartStop(34,2))),...
%         min(fun_TimelockData_perTDPredp34{3}.(inputData_Avg)(i_chan,Sample_Tone_StartStop(1):Sample_Tone_StartStop(34,2)))]); 
%     maxAmpElec = max([max(fun_TimelockData_perTDPredp34{1}.(inputData_Avg)(i_chan,Sample_Tone_StartStop(1):Sample_Tone_StartStop(34,2))),...
%         max(fun_TimelockData_perTDPredp34{2}.(inputData_Avg)(i_chan,Sample_Tone_StartStop(1):Sample_Tone_StartStop(34,2))),...
%         min(fun_TimelockData_perTDPredp34{3}.(inputData_Avg)(i_chan,Sample_Tone_StartStop(1):Sample_Tone_StartStop(34,2)))]);
%     ylim([minAmpElec maxAmpElec ])

    set(gca,'FontSize',2,'XTickLabelRotation',0)
    set(gca,'xtick',Sample_Tone_StartStop(:,1)')
    x_label = Labels_ToneStart;
    set(gca,'FontSize',10)
    xticklabels(x_label);
    
%     if subplot_position == 3
%         legend('Low Predp34', 'Med Predp34', 'High Predp34','Location','NorthWest')
%     end

        title({[fun_TimelockData_perTDPredp34{1}.label{i_chan} '; ' ...
            num2str(sum(fun_TimelockStat_Predp34_perTD{i_chan}{i_TD}{i_inputData}.mask)) ' sign. samples (~' ...
            num2str(round(sum(fun_TimelockStat_Predp34_perTD{i_chan}{i_TD}{i_inputData}.mask)./...
        length(fun_TimelockStat_Predp34_perTD{i_chan}{i_TD}{i_inputData}.mask)*100,2)) '%)']});        

end
       
suptitle({[sub '; ' fun_inputData '; ' ToneDur_title ],...
    ['Electrodes with sign. difference (' fun_inputStat ') between timelocked p*34 (GAvg across trials)']});


if save_fig == 1
    path_fig = ([paths_NASTD_ECoG.Fig_Timelocked_SingleSubs sub '/Stat/' fun_inputStat '/']);
    mkdir(path_fig);
    filename     = [sub '_' fun_inputData '_' ToneDur_title '_TL' fun_inputStat '_Summary.png'];
    figfile      = [path_fig filename];
    saveas(gcf, [figfile], 'png'); %save png version
    close
end

end
