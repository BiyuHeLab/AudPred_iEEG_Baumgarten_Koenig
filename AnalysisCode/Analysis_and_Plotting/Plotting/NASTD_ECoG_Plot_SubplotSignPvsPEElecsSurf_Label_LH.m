function sp_handle_surf = NASTD_ECoG_Plot_SubplotSignPvsPEElecsSurf_Label_LH...
    (coords, SelElecLabes, vals, sign_index,...
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
sp_handle_surf.L = subplot(DimSubplot(1),DimSubplot(2),CounterSubplot);
pos1 = get(sp_handle_surf.L,'Position');
newpos1 = pos1 + SubplotPosition;
set(sp_handle_surf.L,'Position',newpos1)

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
set(gca,'clim',clims);
% if CounterSubplot == 1
%     h = colorbar;
%     h.Position(1) = ColorbarPosition(1); %sets colorbar to the right
%     h.Position(2) = ColorbarPosition(2); %sets colorbar higher
%     h.Position(4) = h.Position(4)*ColorbarPosition(4); %makes colorbar shorter
%     ylabel(h, 'Kprime Ratio (# integrated tones in short TD/ # integrated tones in long TD)')
%     h.Ticks = [0, 1 , 2];
%     h.TickLabels = {'0','1','2'};
% else
colorbar off
% end

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

%sign prediction effect channels
chanSize = chanSize/3;

if any(sign_index(:,1)) %if there are any sign. channels
    coords_signElecs = coords((sign_index(:,1)),:);
    label_signElecs = SelElecLabes(sign_index(:,1));
    
    nearestVertex = nearestVertices(t.vertices,coords(find(sign_index(:,1)),:));
    %uses subfunction nearestVertices to find nearest vertex for selected chan position
    [xs,ys,zs] = sphere;
    temp = surf2patch(xs,ys,zs,ones(size(zs)));
    electrodePatch = struct;
    
    sign_vals = vals(find(sign_index(:,1))); %select only vals with sign. p-val
    numChans = length(sign_vals);
    
    electrodePatch.faces = zeros(numChans * size(temp.faces,1),4);
    electrodePatch.vertices = zeros(numChans * size(temp.vertices,1),3);
    electrodePatch.facevertexcdata = zeros(numChans * size(temp.facevertexcdata,1),1);
    electrodePatch.edgecolor = 'interp';
    electrodePatch.facecolor = 'none';

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
        u = text(coords_signElecs(cn,1) - 50,... %ensures that labels are not behind surface
            coords_signElecs(cn,2),... %little to the left so that it overlaps with elec symbol
            coords_signElecs(cn,3),...
            label_signElecs(cn),'Interpreter', 'none');
        set(u,'Color',textcolor_rgb,'FontSize',6,'FontWeight', 'bold'); %black
        
    end
    patch(electrodePatch, 'Marker', '+', 'LineWidth', 2);
end

%sign PE effect channels
if any(sign_index(:,2)) %if there are any sign. channels
    coords_signElecs = coords((sign_index(:,2)),:);
    label_signElecs = SelElecLabes(sign_index(:,2));
    
    nearestVertex = nearestVertices(t.vertices,coords(find(sign_index(:,2)),:));
    %uses subfunction nearestVertices to find nearest vertex for selected chan position
    [xs,ys,zs] = sphere;
    temp = surf2patch(xs,ys,zs,ones(size(zs)));
    electrodePatch = struct;
    
    sign_vals = vals(find(sign_index(:,2))); %select only vals with sign. p-val
    numChans = length(sign_vals);
    
    electrodePatch.faces = zeros(numChans * size(temp.faces,1),4);
    electrodePatch.vertices = zeros(numChans * size(temp.vertices,1),3);
    electrodePatch.facevertexcdata = zeros(numChans * size(temp.facevertexcdata,1),1);
    electrodePatch.edgecolor = 'interp';
    electrodePatch.facecolor = 'none';

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
        u = text(coords_signElecs(cn,1) - 50,... %ensures that labels are not behind surface
            coords_signElecs(cn,2),... %little to the left so that it overlaps with elec symbol
            coords_signElecs(cn,3),...
            label_signElecs(cn),'Interpreter', 'none');
        set(u,'Color',textcolor_rgb,'FontSize',6,'FontWeight', 'bold'); %black
        
    end
    patch(electrodePatch, 'Marker', '^', 'LineWidth', 1);
end

material dull;

axis([-150 150 -110 70 -50 80]); %Focus
title([win_title ' Both H projected on LH'],'FontSize',8);

end