function nii_nii2atlas_rois(source, lutname, roiList)
% Convert NIfTI atlas to MZ3 atlas
% source : indexed NIfTI atlas
% lutname : (optional) ImageJ-format LUT file
% roiList : (optional) vector of label indices to convert (e.g., [1 3 5])

if ~exist('source','var') || isempty(source) %atlas not specified
    [file,path] = uigetfile('*.nii;*.nii.gz', 'Select a NIfTI atlas');
    if isnumeric(file), return; end
    source = fullfile(path,file);
end
if ~exist(source, 'file')
    error('Source file not found: %s', source);
end
if ~exist('lutname','var') || isempty(lutname)
    [file,path] = uigetfile('*.lut', '(optional) select an ImageJ-format color table');
    if isequal(file,0)
        lutname = '';
        fprintf('Using default color table\n');
    else
        lutname = fullfile(path,file);
    end
end
if isempty(which('spm'))
    error('Please get SPM and add to your MATLAB path'); 
end
if isempty(which('nii_reslice_target'))
    error('Please get spmScripts from GitHub and add to your path'); 
end

prefix = '';
reduce = 0.5;
interp = 4; % interpolation level
isReverseFace = true;
target = ''; % optional target NIfTI for reslicing
mask = '';   % optional mask file
smoothVox = 2;
isResize = false;
thresh = 0.555;

% Load source image
[hdrS, imgS] = loadNiiSub(source);
if ~isempty(mask)
    [~, imgM] = loadNiiSub(mask);
    imgS(imgM ~= 0) = 0; 
end
if ~isempty(target)
    isResize = true;
    [hdrT, ~] = loadNiiSub(target);
end
imgS = uint16(imgS);
nROI = max(imgS(:)); 

fprintf('Converting up to %d regions of interest\n', nROI);

for i = 1 : nROI
    if ~isempty(roiList) && ~ismember(i, roiList)
        continue; % Skip labels not in list
    end

    imgI = (imgS == i);
    if max(imgI(:)) == 0
        continue; % Empty region
    end
    imgI = double(imgI);
    if smoothVox > 0
        presmooth = imgI + 0; % force copy
        spm_smooth(presmooth, imgI, smoothVox, 0);
    end
    outnm = sprintf('%so%04d.mz3', prefix, i);
    if ~isResize
        img2meshSub(hdrS, imgI, outnm, reduce, thresh, isReverseFace);
    else
        [~, outimg] = nii_reslice_target(hdrS, imgI, hdrT, interp);
        if smoothVox > 0
            presmooth = outimg + 0;
            spm_smooth(presmooth, outimg, smoothVox, 0);
        end
        hdr = hdrS;
        hdr.pinfo = [1;0;0]; 
        img2meshSub(hdrT, outimg, outnm, reduce, thresh, isReverseFace);
    end
end

obj2mz3(prefix, nROI, lutname, true);

end % nii_nii2atlas()


function img2meshSub(hdr, img, outnm, reduce, thresh, isReverse)
thresh = min(img(:)) + (thresh * (max(img(:)) - min(img(:))));
if (max(img(:)) < thresh) || (min(img(:)) > thresh)
    if min(img(:)) == max(img(:))
        warning('Range %g..%g - No voxels survive for %s\n', min(img(:)), max(img(:)), outnm);
        return; 
    end
    thresh = min(img(:)) + (0.5 * (max(img(:)) - min(img(:))));
    warning('Resetting threshold to %g for %s\n', thresh, outnm); 
end
FV = isosurface(img, thresh);
r = reduce;
isManifold = true;
if size(FV.faces,1) < 160
    r = 1;
end
while (~isManifold) && (r < 1)
    nR = round(r * size(FV.faces,1));
    if nR < 160
       nR = 160; 
    end
    if mod(nR,2), nR = nR + 1; end
    FVr = reducepatch(FV, nR);
    isManifold = manifoldSub(FVr);
    if ~isManifold
        r = r + 0.025;
    end
end
if (r < 1) && exist('FVr', 'var')
    FV = FVr;
end
FV.vertices = FV.vertices(:, [2,1,3]); % swap X/Y
vx = [FV.vertices ones(size(FV.vertices,1),1)];
vx = mtimes(hdr.mat, vx')';
FV.vertices = vx(:,1:3);
if exist('isReverse', 'var') && isReverse
    FV.faces = fliplr(FV.faces); % reverse winding
end
writeMz3(outnm, FV.faces, FV.vertices, [], [], false);
end


function isManifold = manifoldSub(fv)
edges = sort(cat(1, fv.faces(:,1:2), fv.faces(:,2:3), fv.faces(:,[3 1])), 2);
[unqEdges, ~, ~] = unique(edges, 'rows');
isManifold = (size(edges,1) == size(unqEdges,1)*2);
end


function [hdr, img] = loadNiiSub(fnm)
hdr = spm_vol(fnm);
img = spm_read_vols(hdr);
end