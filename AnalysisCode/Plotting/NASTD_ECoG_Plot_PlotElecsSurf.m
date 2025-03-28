function c = NASTD_ECoG_Plot_PlotElecsSurf...
    (coords,vals,chanSize,...
    clims,cmap,view_angle,...
    inputData,ToneDur_text)

% %% Test input data
% coords = data_ECoGfiltref_trials.elec.chanpos(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF,:); %electrode MNI coordinates
% vals = dat.avg; %parameter of interest for resp. electrodes
% chanSize = (vals./vals)*2.5; %electrode size (arbitrary)
% clims = [1,max(vals)]; %colormap limits
% cmap = 'jet'; 

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

h = figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
        
patch(t);
hold on
h = rotate3d;
h.ActionPostCallback = @RotationCallback;
h.Enable = 'on';

axis square;
axis equal;

view(view_angle)
c = camlight('headlight','infinite');
lighting gouraud
material dull;
axis off

%% colormap
colormap(cmap); 
set(gca,'clim',clims);

h = colorbar;
h.Position(2) = 0.31; %sets colorbar higher
h.Position(4) = h.Position(4)*0.6; %makes colorbar shorter
ylabel(h, 'Exp Kprime (no threshold - all elecs)','FontSize',16)
yticklabels('auto');
set(gca,'FontSize',16);

%% plot channels
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
end
patch(electrodePatch);
material dull;
lighting gouraud

axis([-100 100 -115 85 -100 100]);

% title(['N1-localized electrodes - ' inputData '; ' ToneDur_text 'ms TD']);

    function RotationCallback(~,~)
        c = camlight(c,'headlight');
    end
end