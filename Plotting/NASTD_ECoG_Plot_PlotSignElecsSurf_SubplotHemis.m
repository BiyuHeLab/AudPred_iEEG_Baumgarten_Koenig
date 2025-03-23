function c = NASTD_ECoG_Plot_PlotSignElecsSurf_SubplotHemis...
    (coords, vals, sign_index, ...
    chanSize, clims, cmap)

% %% Test input data
% coords = data_ECoGfiltref_trials.elec.chanpos(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF,:); %electrode MNI coordinates
% vals = dat.avg; %parameter of interest for resp. electrodes
% chanSize = (vals./vals)*2.5; %electrode size (arbitrary)
% clims = [1,max(vals)]; %colormap limits
% cmap = 'jet'; 

%% Plot 3D Cortex surface LH
addpath('/isilon/LFMI/VMdrive/Thomas/toolboxes/fieldtrip-20190314/external/freesurfer/');

[lhvtx,lhfaces]=read_surf('lh.pial'); %pial surface - 3D surface file used to estimate cortical gray matter volume
[rhvtx,rhfaces]=read_surf('rh.pial');
t.faces = [lhfaces+1 ; rhfaces+1+length(lhvtx)];
t.vertices = [lhvtx ; rhvtx];
t.EdgeColor = 'none';
t.faceColor = 'interp';
t.facevertexCData = 0.9*ones(length(t.vertices),3);
t.facealpha = 1;

figure('visible','on'); %ensures figure doesn't pop up during plotting
set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
subplot(1,2,1)

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
h = colorbar;
h.Position(1) = 0.1;
h.Position(2) = 0.4; %sets colorbar higher
% h.Position(4) = h.Position(4)*1.5; %makes colorbar shorter
ylabel(h, 'Exp Kprime (threshold p < 0.05 - exp vs. shuff Kprime)')
yticklabels('auto');
set(gca,'FontSize',16);


%% plot channels
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
    nearestVertex = nearestVertices(t.vertices,coords(find(sign_index),:)); %uses subfunction nearestVertices to find nearest vertex for selected chan position
    [xs,ys,zs] = sphere;
    temp = surf2patch(xs,ys,zs,ones(size(zs)));
    counter = 0;
    electrodePatch = struct;

    i_signelecs = find(sign_index);
    sign_vals   = vals(i_signelecs); %select only vals with sign. p-val
    numChans    = length(sign_vals);

    electrodePatch.faces = zeros(numChans * size(temp.faces,1),4);
    electrodePatch.vertices = zeros(numChans * size(temp.vertices,1),3);
    electrodePatch.facevertexcdata = zeros(numChans * size(temp.facevertexcdata,1),1);
    electrodePatch.edgecolor = 'none';
    electrodePatch.facecolor = 'interp';

    for cn = 1:numChans
        x = (chanSize(i_signelecs(cn))*xs) + t.vertices(nearestVertex(cn),1);
        y = (chanSize(i_signelecs(cn))*ys) + t.vertices(nearestVertex(cn),2);
        z = (chanSize(i_signelecs(cn))*zs) + t.vertices(nearestVertex(cn),3);
        temp = surf2patch(x,y,z,sign_vals(cn) * ones(size(z)));
        faceInds = 1 + (cn-1) * size(temp.faces,1):cn * size(temp.faces,1);
        verticesInds =  1 + (cn-1) * size(temp.vertices,1):cn * size(temp.vertices,1);
        electrodePatch.faces(faceInds,:) = temp.faces + verticesInds(1) - 1;
        electrodePatch.vertices(verticesInds,:) = temp.vertices;
        electrodePatch.facevertexcdata(verticesInds,:) = temp.facevertexcdata;
    end
    patch(electrodePatch);
end

material dull;
lighting gouraud

axis([-100 100 -115 85 -100 100]);
title('Left','FontSize',16);

%% RH
subplot(1,2,2)

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

%% colormap
colormap(cmap); colorbar off
set(gca,'clim',clims);

%% plot channels
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
    nearestVertex = nearestVertices(t.vertices,coords(find(sign_index),:)); %uses subfunction nearestVertices to find nearest vertex for selected chan position
    [xs,ys,zs] = sphere;
    temp = surf2patch(xs,ys,zs,ones(size(zs)));
    counter = 0;
    electrodePatch = struct;

    i_signelecs = find(sign_index);
    sign_vals   = vals(i_signelecs); %select only vals with sign. p-val
    numChans    = length(sign_vals);

    electrodePatch.faces = zeros(numChans * size(temp.faces,1),4);
    electrodePatch.vertices = zeros(numChans * size(temp.vertices,1),3);
    electrodePatch.facevertexcdata = zeros(numChans * size(temp.facevertexcdata,1),1);
    electrodePatch.edgecolor = 'none';
    electrodePatch.facecolor = 'interp';


    for cn = 1:numChans
        x = (chanSize(i_signelecs(cn))*xs) + t.vertices(nearestVertex(cn),1);
        y = (chanSize(i_signelecs(cn))*ys) + t.vertices(nearestVertex(cn),2);
        z = (chanSize(i_signelecs(cn))*zs) + t.vertices(nearestVertex(cn),3);
        temp = surf2patch(x,y,z,sign_vals(cn) * ones(size(z)));
        faceInds = 1 + (cn-1) * size(temp.faces,1):cn * size(temp.faces,1);
        verticesInds =  1 + (cn-1) * size(temp.vertices,1):cn * size(temp.vertices,1);
        electrodePatch.faces(faceInds,:) = temp.faces + verticesInds(1) - 1;
        electrodePatch.vertices(verticesInds,:) = temp.vertices;
        electrodePatch.facevertexcdata(verticesInds,:) = temp.facevertexcdata;
    end
    patch(electrodePatch);
end

material dull;
lighting gouraud

axis([-100 100 -115 85 -100 100]);
title('Right','FontSize',16);

end