function sp_handle_surf = NASTD_ECoG_Plot_SubplotSignElecsSurf_ColorAnatAllLobes_LH_Pairs...
    (coords_1, coords_2, coords_3, ...
    vals_1, vals_2, vals_3, ....
    sign_index_1, sign_index_2, sign_index_3, ...
    pairings_1, pairings_2, pairings_3, ...
    Pairing2plot, ...
    chanSize_1, chanSize_2, chanSize_3, ...
    clims, cmap, ...
    DimSubplot, CounterSubplot, SubplotPosition, ...
    win_title)


%% Plot 3D Cortex surface
addpath('/isilon/LFMI/VMdrive/Thomas/toolboxes/fieldtrip-20190314/external/freesurfer/');

[lhvtx,lhfaces]=read_surf('lh.pial'); %pial surface - 3D surface file used to estimate cortical gray matter volume
[rhvtx,rhfaces]=read_surf('rh.pial');
t.faces = [lhfaces+1 ; rhfaces+1+length(lhvtx)];
t.vertices = [lhvtx ; rhvtx];
t.EdgeColor = 'none';
t.faceColor = 'interp';
t.facevertexCData = 0.9*ones(length(t.vertices),3);
t.facealpha = 1;

%% LH
sp_handle_surf = subplot(DimSubplot(1),DimSubplot(2),CounterSubplot);
pos1 = get(sp_handle_surf,'Position');
newpos1 = pos1 + SubplotPosition;
set(sp_handle_surf,'Position',newpos1)

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

%% colormap
colormap(cmap);
temp_colormap = colormap;
set(gca,'clim',clims);
colorbar off

%Set 1
coords_signElecs_1 = coords_1((sign_index_1),:);
nearestVertex_1 = nearestVertices(t.vertices,coords_1(find(sign_index_1),:));
[xs,ys,zs] = sphere;
temp = [];
temp = surf2patch(xs,ys,zs,ones(size(zs)));
counter = 0;
electrodePatch = struct;

sign_vals_1 = vals_1(find(sign_index_1)); %select only vals with sign. p-val
numChans_1 = length(sign_vals_1);

electrodePatch.faces = zeros(numChans_1 * size(temp.faces,1),4);
electrodePatch.vertices = zeros(numChans_1 * size(temp.vertices,1),3);
electrodePatch.facevertexcdata = zeros(numChans_1 * size(temp.facevertexcdata,1),1);
electrodePatch.edgecolor = 'none';
electrodePatch.facecolor = 'interp';

for cn = 1:numChans_1
    x = (chanSize_1(cn)*xs) + t.vertices(nearestVertex_1(cn),1);
    
    %CAVE: Changing x-axis to highlight extreme values
    if numChans_1 > (length(coords_1)/2)
        x = x - (abs(sign_vals_1(cn)*22)); %Places extreme values in front of plot
        disp('-- Warning: changing x-axis to place high values in front --')
    end
    
    y = (chanSize_1(cn)*ys) + t.vertices(nearestVertex_1(cn),2);
    z = (chanSize_1(cn)*zs) + t.vertices(nearestVertex_1(cn),3);
    temp = surf2patch(x,y,z,sign_vals_1(cn) * ones(size(z)));
    faceInds = 1 + (cn-1) * size(temp.faces,1):cn * size(temp.faces,1);
    verticesInds =  1 + (cn-1) * size(temp.vertices,1):cn * size(temp.vertices,1);
    electrodePatch.faces(faceInds,:) = temp.faces + verticesInds(1) - 1;
    electrodePatch.vertices(verticesInds,:) = temp.vertices;
    electrodePatch.facevertexcdata(verticesInds,:) = temp.facevertexcdata;
    
end
patch(electrodePatch);

%Set 2
coords_signElecs_2 = coords_2((sign_index_2),:);
nearestVertex_2 = nearestVertices(t.vertices,coords_2(find(sign_index_2),:));
[xs,ys,zs] = sphere;
temp = [];
temp = surf2patch(xs,ys,zs,ones(size(zs)));
counter = 0;
electrodePatch = struct;

sign_vals_2 = vals_2(find(sign_index_2)); %select only vals with sign. p-val
numChans_2 = length(sign_vals_2);

electrodePatch.faces = zeros(numChans_2 * size(temp.faces,1),4);
electrodePatch.vertices = zeros(numChans_2 * size(temp.vertices,1),3);
electrodePatch.facevertexcdata = zeros(numChans_2 * size(temp.facevertexcdata,1),1);
electrodePatch.edgecolor = 'none';
electrodePatch.facecolor = 'interp';

for cn = 1:numChans_2
    x = (chanSize_2(cn)*xs) + t.vertices(nearestVertex_2(cn),1);
    
    %CAVE: Changing x-axis to highlight extreme values
    if numChans_2 > (length(coords_2)/2)
        x = x - (abs(sign_vals_2(cn)*20)); %Places extreme values in front of plot
        disp('-- Warning: changing x-axis to place high values in front --')
    end
    
    y = (chanSize_2(cn)*ys) + t.vertices(nearestVertex_2(cn),2);
    z = (chanSize_2(cn)*zs) + t.vertices(nearestVertex_2(cn),3);
    temp = surf2patch(x,y,z,sign_vals_2(cn) * ones(size(z)));
    faceInds = 1 + (cn-1) * size(temp.faces,1):cn * size(temp.faces,1);
    verticesInds =  1 + (cn-1) * size(temp.vertices,1):cn * size(temp.vertices,1);
    electrodePatch.faces(faceInds,:) = temp.faces + verticesInds(1) - 1;
    electrodePatch.vertices(verticesInds,:) = temp.vertices;
    electrodePatch.facevertexcdata(verticesInds,:) = temp.facevertexcdata;
    
end
patch(electrodePatch);

%Set 3
coords_signElecs_3 = coords_3((sign_index_3),:);
nearestVertex_3 = nearestVertices(t.vertices,coords_3(find(sign_index_3),:));
[xs,ys,zs] = sphere;
temp = [];
temp = surf2patch(xs,ys,zs,ones(size(zs)));
counter = 0;
electrodePatch = struct;

sign_vals_3 = vals_3(find(sign_index_3)); %select only vals with sign. p-val
numChans_3 = length(sign_vals_3);

electrodePatch.faces = zeros(numChans_3 * size(temp.faces,1),4);
electrodePatch.vertices = zeros(numChans_3 * size(temp.vertices,1),3);
electrodePatch.facevertexcdata = zeros(numChans_3 * size(temp.facevertexcdata,1),1);
electrodePatch.edgecolor = 'none';
electrodePatch.facecolor = 'interp';

for cn = 1:numChans_3
    x = (chanSize_3(cn)*xs) + t.vertices(nearestVertex_3(cn),1);
    
    %CAVE: Changing x-axis to highlight extreme values
    if numChans_3 > (length(coords_3)/2)
        x = x - (abs(sign_vals_3(cn)*20)); %Places extreme values in front of plot
        disp('-- Warning: changing x-axis to place high values in front --')
    end
    
    y = (chanSize_3(cn)*ys) + t.vertices(nearestVertex_3(cn),2);
    z = (chanSize_3(cn)*zs) + t.vertices(nearestVertex_3(cn),3);
    temp = surf2patch(x,y,z,sign_vals_3(cn) * ones(size(z)));
    faceInds = 1 + (cn-1) * size(temp.faces,1):cn * size(temp.faces,1);
    verticesInds =  1 + (cn-1) * size(temp.vertices,1):cn * size(temp.vertices,1);
    electrodePatch.faces(faceInds,:) = temp.faces + verticesInds(1) - 1;
    electrodePatch.vertices(verticesInds,:) = temp.vertices;
    electrodePatch.facevertexcdata(verticesInds,:) = temp.facevertexcdata;
    
end
patch(electrodePatch);

%Add lines between electrode pairs specified by var: pairings, visualizing connectivity computations
if strcmp(Pairing2plot, 'Frontal_Temporal')
    hold on;
    for i_sourceelec = unique(pairings_1(:,1),'stable')'
        if ~isnan(i_sourceelec)
            %Coordinates of current elec
            x = t.vertices(nearestVertex_1(i_sourceelec),1);
            y = t.vertices(nearestVertex_1(i_sourceelec),2);
            z = t.vertices(nearestVertex_1(i_sourceelec),3);
            %Find electrodes linked to current electrodes
            connected_elecs = unique(pairings_1(pairings_1(:,1) == i_sourceelec,2),'stable');
            if length(connected_elecs) > 0
                for i_connected_elecs = 1:length(connected_elecs)
                    x_other = t.vertices(nearestVertex_1(connected_elecs(i_connected_elecs)),1);
                    y_other = t.vertices(nearestVertex_1(connected_elecs(i_connected_elecs)),2);
                    z_other = t.vertices(nearestVertex_1(connected_elecs(i_connected_elecs)),3);
                    hold on;
                    l = line([(x-20) (x_other-20)],[y y_other],[z z_other], 'LineWidth', 0.5);
                    stepsize_colormap = 1;
                    if any(sum(~isnan(pairings_1)) > 5000)
                        l.Color = [0.3, 0.3, 0.3, 0.1]; %Color & Transparency
                    elseif any(sum(~isnan(pairings_1)) > 1000)
                        l.Color = [0.3, 0.3, 0.3, 0.5]; %Color & Transparency
                    else
                        l.Color = [0.2, 0.2, 0.2, 0.5];%[0, 0, 0]; %Color
                    end
                    %                 l.Color = temp_colormap(1 + ((unique(sub_index) - 1) * stepsize_colormap),:);
                end
            end
        end
    end
elseif strcmp(Pairing2plot, 'Frontal_Parietal')
    hold on;
    for i_sourceelec = unique(pairings_2(:,1),'stable')'
        if ~isnan(i_sourceelec)
            %Coordinates of current elec
            x = t.vertices(nearestVertex_2(i_sourceelec),1);
            y = t.vertices(nearestVertex_2(i_sourceelec),2);
            z = t.vertices(nearestVertex_2(i_sourceelec),3);
            %Find electrodes linked to current electrodes
            connected_elecs = unique(pairings_2(pairings_2(:,1) == i_sourceelec,2),'stable');
            if length(connected_elecs) > 0
                for i_connected_elecs = 1:length(connected_elecs)
                    x_other = t.vertices(nearestVertex_2(connected_elecs(i_connected_elecs)),1);
                    y_other = t.vertices(nearestVertex_2(connected_elecs(i_connected_elecs)),2);
                    z_other = t.vertices(nearestVertex_2(connected_elecs(i_connected_elecs)),3);
                    hold on;
                    l = line([(x-20) (x_other-20)],[y y_other],[z z_other], 'LineWidth', 0.5);
                    stepsize_colormap = 1;
                    if any(sum(~isnan(pairings_2)) > 5000)
                        l.Color = [0.3, 0.3, 0.3, 0.1]; %Color & Transparency
                    elseif any(sum(~isnan(pairings_2)) > 1000)
                        l.Color = [0.3, 0.3, 0.3, 0.5]; %Color & Transparency
                    else
                        l.Color = [0.2, 0.2, 0.2, 0.5];%[0, 0, 0]; %Color
                    end
                    %                 l.Color = temp_colormap(1 + ((unique(sub_index) - 1) * stepsize_colormap),:);
                end
            end
        end
    end
elseif strcmp(Pairing2plot, 'Parietal_Temporal')
    hold on;
    for i_sourceelec = unique(pairings_3(:,1),'stable')'
        if ~isnan(i_sourceelec)
            %Coordinates of current elec
            x = t.vertices(nearestVertex_3(i_sourceelec),1);
            y = t.vertices(nearestVertex_3(i_sourceelec),2);
            z = t.vertices(nearestVertex_3(i_sourceelec),3);
            %Find electrodes linked to current electrodes
            connected_elecs = unique(pairings_3(pairings_3(:,1) == i_sourceelec,2),'stable');
            if length(connected_elecs) > 0
                for i_connected_elecs = 1:length(connected_elecs)
                    x_other = t.vertices(nearestVertex_3(connected_elecs(i_connected_elecs)),1);
                    y_other = t.vertices(nearestVertex_3(connected_elecs(i_connected_elecs)),2);
                    z_other = t.vertices(nearestVertex_3(connected_elecs(i_connected_elecs)),3);
                    hold on;
                    l = line([(x-20) (x_other-20)],[y y_other],[z z_other], 'LineWidth', 0.5);
                    stepsize_colormap = 1;
                    if any(sum(~isnan(pairings_3)) > 5000)
                        l.Color = [0.3, 0.3, 0.3, 0.1]; %Color & Transparency
                    elseif any(sum(~isnan(pairings_3)) > 1000)
                        l.Color = [0.3, 0.3, 0.3, 0.5]; %Color & Transparency
                    else
                        l.Color = [0.2, 0.2, 0.2, 0.5];%[0, 0, 0]; %Color
                    end
                    %                 l.Color = temp_colormap(1 + ((unique(sub_index) - 1) * stepsize_colormap),:);
                end
            end
        end
    end
end


patch(electrodePatch);
material dull;
lighting gouraud

axis([-150 150 -110 70 -50 80]); %Focus
title([win_title],'FontSize',8);

