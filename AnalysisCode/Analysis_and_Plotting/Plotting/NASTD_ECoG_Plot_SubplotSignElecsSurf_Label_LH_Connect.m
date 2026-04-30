function sp_handle_surf = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH_Connect...
    (coords, SelElecLabes, vals, sign_index, ...
    sub_index, anatcat_index, ...
    chanSize, clims, cmap, textcolor_rgb, ...
    DimSubplot, CounterSubplot, SubplotPosition, ColorbarPosition,...
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

%all channels (not only sign ones)
nearestVertex = nearestVertices(t.vertices,coords); %uses subfunction nearestVertices to find nearest vertex for each chan position
[xs,ys,zs] = sphere;
temp = surf2patch(xs,ys,zs,ones(size(zs)));
counter = 0;
electrodePatch = struct;
numChans = length(vals);

electrodePatch.faces = zeros(numChans * size(temp.faces,1),4);
electrodePatch.vertices = zeros(numChans * size(temp.vertices,1),3);
electrodePatch.facevertexcdata = zeros(numChans * size(temp.facevertexcdata,1),1);
electrodePatch.edgecolor = 'none';
electrodePatch.facecolor = 'k'; %black

chanSize_nonsignElec = chanSize/3;

for cn = 1:numChans
    x = (chanSize_nonsignElec(cn)*xs) + t.vertices(nearestVertex(cn),1);
    y = (chanSize_nonsignElec(cn)*ys) + t.vertices(nearestVertex(cn),2);
    z = (chanSize_nonsignElec(cn)*zs) + t.vertices(nearestVertex(cn),3);
    temp = surf2patch(x,y,z,vals(cn) * ones(size(z)));
    faceInds = 1 + (cn-1) * size(temp.faces,1):cn * size(temp.faces,1);
    verticesInds =  1 + (cn-1) * size(temp.vertices,1):cn * size(temp.vertices,1);
    electrodePatch.faces(faceInds,:) = temp.faces + verticesInds(1) - 1;
    electrodePatch.vertices(verticesInds,:) = temp.vertices;
    electrodePatch.facevertexcdata(verticesInds,:) = temp.facevertexcdata;
end
patch(electrodePatch);

%sign channels
if any(sign_index) %if there are any sign. channels
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
        
%         %Add electrode labels add respective position
%         u = text(coords_signElecs(cn,1) - 50,... %ensures that labels are not behind surface
%             coords_signElecs(cn,2),... %little to the left so that it overlaps with elec symbol
%             coords_signElecs(cn,3),...
%             label_signElecs(cn),'Interpreter', 'none');
%         set(u,'Color',textcolor_rgb,'FontSize',6,'FontWeight', 'bold'); %black
    end
    patch(electrodePatch);
    
    %Add lines between electrodes in differing anatomical parcellations for the same subject, visualizing potential connectivity computations
    subject_currsignelec = sub_index(sign_index); %Determine subject of sign elecs
    anatparcel_currsignelec = min(anatcat_index(sign_index,:),[],2);
    
    hold on;
    for i_sub = unique(subject_currsignelec)'
        for cn = find(subject_currsignelec == i_sub)'
            
            %Coordinates of current elec
            x = t.vertices(nearestVertex(cn),1);
            y = t.vertices(nearestVertex(cn),2);
            z = t.vertices(nearestVertex(cn),3);
            
            if length(find(subject_currsignelec == i_sub)) > 1
                selected_elecs = setdiff(find(subject_currsignelec == i_sub)', cn);               
                for i_otherelecs = selected_elecs
                    
                    if anatparcel_currsignelec(i_otherelecs) ...
                            ~= anatparcel_currsignelec(cn)  
                        
                        x_other = t.vertices(nearestVertex(i_otherelecs),1);
                        y_other = t.vertices(nearestVertex(i_otherelecs),2);
                        z_other = t.vertices(nearestVertex(i_otherelecs),3);
                        hold on;
                        l = line([(x-50) (x_other-50)],[y y_other],[z z_other], 'LineWidth', 1);
%                         stepsize_colormap = floor((length(temp_colormap)/size(unique(subject_currsignelec),1)));
                        stepsize_colormap = 1;
                        l.Color = temp_colormap(1 + ((i_sub - 1) * stepsize_colormap),:);
                        
                    end
                end
            end
        end
    end
end

patch(electrodePatch);
material dull;
lighting gouraud
 
axis([-150 150 -110 70 -50 80]); %Focus
title([win_title ' Both H projected on LH'],'FontSize',8);

end