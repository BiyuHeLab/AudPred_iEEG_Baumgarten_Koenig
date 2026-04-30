atlas_path = 'brodmann.nii'; 

% Define ROIs: Brodmann label → Region name
roi_map = containers.Map( ...
    {41, 42, 44, 45, 8, 40, 6, 22}, ...
    {'Auditory Cortex (BA41)', ...
     'Auditory Cortex (BA42)', ...
     'IFC Dorsal (BA44)', ...
     'IFC Ventral (BA45)', ...
     'dlPFC (BA8)', ...
     'IPL (BA40)', ...
     'PMC (BA6)', ...
     'STS (BA22)'});

% Load volume
V = spm_vol(atlas_path);
Y = spm_read_vols(V);

% Get all unique non-zero labels
labels = unique(Y(:));
labels(labels == 0) = [];  % remove background

% Preallocate results
label_coords = struct();

% Loop through labels
for i = 1:length(labels)
    label_val = labels(i);
    
    % Find voxels with this label
    [x, y, z] = ind2sub(size(Y), find(Y == label_val));
    voxels = [x y z];
    
    % Convert to MNI space
    mni_coords = V.mat * [voxels'; ones(1, size(voxels, 1))];
    mni_coords = mni_coords(1:3, :)';
    
    % Store
    fieldname = ['BA' num2str(label_val)];
    label_coords.(fieldname) = mni_coords;
    
    % Optional: save to CSV
    csv_filename = sprintf('MNI_coords_%s.csv', fieldname);
    writematrix(mni_coords, csv_filename);
    
    fprintf('Saved %d coordinates for %s\n', size(mni_coords, 1), fieldname);
end

% Save all data to .mat file
save('all_BA_MNI_coords.mat', 'label_coords');