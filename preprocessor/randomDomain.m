function pore = randomDomain(NY, NX, windowSize, scalingFactor)

    %% Creating random porous media from noise
    close all
    % create noise
    noisyImage = poissrnd(50, NY, NX);
    % filter with Gaussian
    windowSize = 8;
    blurredImage = imgaussfilt(noisyImage, windowSize);
    blurredImage = mat2gray(blurredImage);
    % binarize image
    scalingFactor = 1.99;
    threshold = scalingFactor * mean2(blurredImage);
    pore = blurredImage < threshold;
    % invert the image
       
    pore(:, 1:50) = 0;
    pore(:, end-51:end) = 0;
    pore(1, :) = 1;
    pore(end, :) = 1;
    
    % cleaning and smoothing
    pore_regions = regionprops(pore, 'all');
    for reg = 1:length(pore_regions)
        region = pore_regions(reg);
        if region.Area < 100
            pixels = region.PixelList;
            for pix = 1:size(pixels, 1)
                pore(pixels(pix, 2), pixels(pix, 1)) = 0;
            end
        end
    end
    
    pore_cp = pore;
    pore(pore_cp == 0) = 1;
    pore(pore_cp == 1) = 0;

    % invert image for premature cleaning
    pore_cp = pore;
    pore(pore_cp == 1) = 0;
    pore(pore_cp == 0) = 1;
    pore = bwmorph(pore, 'bridge'); % bridge unconnected pixels
    pore = bwmorph(pore, 'fill');   % fill isolated pixels
    pore = bwmorph(pore, 'majority');   % set a pixel to 1 if five or more pixels in its 3-by-3 neighborhood are 1; otherwise, set the pixel to 0.
    % invert image back to original
    pore_cp = pore;
    pore(pore_cp == 1) = 0;
    pore(pore_cp == 0) = 1;

    % Remove disconnected areas
    % We can either get rid of all the pores that is not connected to the main connected chunk
    % Or we can put a threshold and remove everything below that threshold
    area_thr = 1000; % threshold area
    pore_regions = regionprops(pore, 'all');
    region_areas = [pore_regions.Area];
    max_area = area_thr; % region_areas(region_areas == max(region_areas));
    
    for reg = 1:length(pore_regions)
        if region_areas(reg) < max_area
            % Get the bounding box of the current region
            bbox = pore_regions(reg).BoundingBox;
            
            % Extract the row and column indices of the bounding box
            row_indices = round(bbox(2)) : round(bbox(2) + bbox(4) - 1);
            col_indices = round(bbox(1)) : round(bbox(1) + bbox(3) - 1);
            
            % Set the corresponding pixels to zero
            pore(row_indices, col_indices) = 0;
        end
    end
    pore(1, :) = 0;
    pore(end, :) = 0;
end
