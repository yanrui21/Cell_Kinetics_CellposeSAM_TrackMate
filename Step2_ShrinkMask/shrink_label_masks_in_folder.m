%Mask shrinking for trackmate
function shrink_label_masks_in_folder(inDir, outDir, scaleFactor)
% shrink_label_masks_in_folder
% - Reads .tif/.tiff label images from inDir
% - For each unique value v (excluding background=0), treats (I==v) as one mask
% - Shrinks each mask by scaleFactor about its centroid
% - Writes new label images into outDir
%
% Example:
%   shrink_label_masks_in_folder('C:\data\labels', 'C:\data\labels_shrunk', 0.5);

    if nargin < 3
        scaleFactor = 0.5;
    end
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    files = [dir(fullfile(inDir, '*.tif')); dir(fullfile(inDir, '*.tiff'))];
    if isempty(files)
        error('No .tif/.tiff files found in: %s', inDir);
    end

    for k = 1:numel(files)
        inPath = fullfile(inDir, files(k).name);
        I = imread(inPath);

        % Output label matrix (same size/type)
        out = zeros(size(I), 'like', I);

        % Decide what is background. Commonly 0.
        vals = unique(I(:));
        vals(vals == 0) = [];   % remove background label (edit if your background differs)

        ref = imref2d(size(I));

        for i = 1:numel(vals)
            v = vals(i);
            mask = (I == v);

            if ~any(mask(:))
                continue;
            end

            % Centroid of ALL pixels of this value (union), in pixel coordinates:
            % rows = y, cols = x
            [r, c] = find(mask);
            cx = mean(c);
            cy = mean(r);

            % Build affine transform that scales about (cx, cy)
            % x' = s*x + (1-s)*cx
            % y' = s*y + (1-s)*cy
            s = scaleFactor;
            T = [ s  0  0
                  0  s  0
                 (1-s)*cx  (1-s)*cy  1 ];

            tform = affine2d(T);

            % Warp the binary mask back onto same canvas
            mask2 = imwarp(mask, tform, ...
                'OutputView', ref, ...
                'InterpolationMethod', 'nearest');

            % Assign label back
            % NOTE: if overlaps ever happen, later labels overwrite earlier ones.
            out(mask2) = v;
        end

        outPath = fullfile(outDir, files(k).name);
        imwrite(out, outPath);
        fprintf('Wrote: %s\n', outPath);
    end
end
