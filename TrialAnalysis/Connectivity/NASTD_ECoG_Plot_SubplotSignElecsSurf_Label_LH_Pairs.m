function sp_handle_surf = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH_Pairs...
    (coords, SelElecLabes, vals, sign_index, typeindex, ...
    sub_index, pairings, effect_comp, ...
    chanSize, clims, cmap, textcolor_rgb, ...
    DimSubplot, CounterSubplot, SubplotPosition, ColorbarPosition,...
    win_title)

% DimSubplot = [1 1];
% CounterSubplot = 1;
% coords = elec_coords
% SelElecLabes = elec_labels
% vals = elec_data
% sign_index = elec_signindex
% sub_index = ones(length(elec_coords),1)*i_sub;
% pairings = elec_pairings
% effect_comp = effect_comp;
% typeindex = elec_typeindex
% SubplotPosition = [0 0 0 0];

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

coords_signElecs = coords((sign_index),:);
label_signElecs = SelElecLabes(sign_index);

nearestVertex = nearestVertices(t.vertices,coords(find(sign_index),:));
%uses subfunction nearestVertices to find nearest vertex for selected chan position
[xs,ys,zs] = sphere;
temp = surf2patch(xs,ys,zs,ones(size(zs)));
counter = 0;
electrodePatch = struct;

sign_vals = vals(find(sign_index)); %select only vals with sign. p-val
numChans = length(sign_vals);

electrodePatch.faces = zeros(numChans * size(temp.faces,1),4);
electrodePatch.vertices = zeros(numChans * size(temp.vertices,1),3);
electrodePatch.facevertexcdata = zeros(numChans * size(temp.facevertexcdata,1),1);
electrodePatch.edgecolor = 'none';
electrodePatch.facecolor = 'interp';

for cn = 1:numChans
    x = (chanSize(cn)*xs) + t.vertices(nearestVertex(cn),1);
    
    %CAVE: Changing x-axis to highlight extreme values
    if numChans > (length(coords)/2)
        x = x - (abs(sign_vals(cn)*10)); %Places extreme values in front of plot
        disp('-- Warning: changing x-axis to place high values in front --')
    end
    
    y = (chanSize(cn)*ys) + t.vertices(nearestVertex(cn),2);
    z = (chanSize(cn)*zs) + t.vertices(nearestVertex(cn),3);
    temp = surf2patch(x,y,z,sign_vals(cn) * ones(size(z)));
    faceInds = 1 + (cn-1) * size(temp.faces,1):cn * size(temp.faces,1);
    verticesInds =  1 + (cn-1) * size(temp.vertices,1):cn * size(temp.vertices,1);
    electrodePatch.faces(faceInds,:) = temp.faces + verticesInds(1) - 1;
    electrodePatch.vertices(verticesInds,:) = temp.vertices;
    electrodePatch.facevertexcdata(verticesInds,:) = temp.facevertexcdata;
    
    %Add electrode labels add respective position
    u = text(coords_signElecs(cn,1) - 150,... %ensures that labels are not behind surface
        coords_signElecs(cn,2),... %little to the left so that it overlaps with elec symbol
        coords_signElecs(cn,3),...
        label_signElecs(cn),'Interpreter', 'none');
    if typeindex(cn) == 1
        set(u,'Color', [0 0 0.3],'FontSize',9,'FontWeight', 'bold'); %Source - blue
    elseif typeindex(cn) == 2
        set(u,'Color', [0.3 0 0],'FontSize',9,'FontWeight', 'bold');  %Target - red
    elseif typeindex(cn) == 3
        set(u,'Color', [0 0 0],'FontSize',9,'FontWeight', 'bold');  %Both - black
    elseif typeindex(cn) == 4
        set(u,'Color', [0 0.5 0],'FontSize',9,'FontWeight', 'bold');  %Sign. Pred & PE Effect - Green        
    end
end
patch(electrodePatch);

%Add lines between electrode pairs specified by var: pairings, visualizing connectivity computations
%Re-number pairings to increasing indices
if ~isnan(pairings)
    if strcmp(effect_comp, 'Pred_Pred') || strcmp(effect_comp, 'PE_PE') || strcmp(effect_comp, 'All_All')  %If same effect, sorting elec indices can run across both source and target elecs
        orig_num = unique(pairings);
        [~,sort_num] = sort(orig_num);
        [rows, columns] = size(pairings);
        for i_row = 1:rows
            for i_columns = 1:columns
                for i_entry = 1:length(orig_num)
                    if pairings(i_row,i_columns) == orig_num(i_entry)
                        pairings(i_row,i_columns) = sort_num(i_entry);
                    end
                end
            end
        end
    elseif strcmp(effect_comp, 'Pred_PE') || strcmp(effect_comp, 'PE_Pred') %If different effect, sorting elec indices must run separately for source and target elecs
        orig_num1 = unique(pairings(:,1));
        [~,sort_num1] = sort(orig_num1);
        orig_num2 = unique(pairings(:,2));
        [~,sort_num2] = sort(orig_num2);
        [rows, ~] = size(pairings);
        for i_row = 1:rows
            for i_entry = 1:length(orig_num1)
                if pairings(i_row,1) == orig_num1(i_entry)
                    pairings(i_row,1) = sort_num1(i_entry);
                end
            end
            for i_entry = 1:length(orig_num2)
                if pairings(i_row,2) == orig_num2(i_entry)
                    pairings(i_row,2) = sort_num2(i_entry) + max(sort_num1);
                end
            end
        end
    end
    
    hold on;
    for i_sourceelec = unique(pairings(:,1))'
        %Coordinates of current elec
        x = t.vertices(nearestVertex(i_sourceelec),1);
        y = t.vertices(nearestVertex(i_sourceelec),2);
        z = t.vertices(nearestVertex(i_sourceelec),3);
        %Find electrodes linked to current electrodes
        connected_elecs = unique(pairings(pairings(:,1) == i_sourceelec,2));
        if length(connected_elecs) > 0
            for i_connected_elecs = 1:length(connected_elecs)
                x_other = t.vertices(nearestVertex(connected_elecs(i_connected_elecs)),1);
                y_other = t.vertices(nearestVertex(connected_elecs(i_connected_elecs)),2);
                z_other = t.vertices(nearestVertex(connected_elecs(i_connected_elecs)),3);
                hold on;
                l = line([(x-50) (x_other-50)],[y y_other],[z z_other], 'LineWidth', 1);
                stepsize_colormap = 1;
                if length(pairings) > 100
                    l.Color = [0.3, 0.3, 0.3, 0.1]; %Color & Transparency
                else
                    l.Color = [0, 0, 0]; %Color                   
                end                
%                 l.Color = [0 1 1]; %Cyan
%                 l.Color = temp_colormap(1 + ((unique(sub_index) - 1) * stepsize_colormap),:);
            end
        end
    end
end

patch(electrodePatch);
material dull;
lighting gouraud

axis([-150 150 -110 70 -50 80]); %Focus
title([win_title],'FontSize',8);

