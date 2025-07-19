function [intensity, reason] = calculateSingleCellIntensity(image, x, y, r, bgInnerFactor, bgOuterFactor, bgMethod)
    %CALCULATESINGLECELLINTENSITY calculates the background-subtracted intensity.
    %   This is a standalone helper function, not a method of the app class.

    % This is a function for the RelativeBaselineChangeApp
    
    % Initialize default output
    intensity = -Inf;
    reason = "Calculation failed";
    
    [rows, cols] = size(image);
    
    % Define local window around the cell
    xmin = max(1, floor(x - r*bgOuterFactor));
    xmax = min(cols, ceil(x + r*bgOuterFactor));
    ymin = max(1, floor(y - r*bgOuterFactor));
    ymax = min(rows, ceil(y + r*bgOuterFactor));
    
    if xmin >= xmax || ymin >= ymax
        reason = "Invalid window";
        return;
    end
    
    [X_local, Y_local] = meshgrid(xmin:xmax, ymin:ymax);
    
    % Create masks for the soma and background annulus
    distSq = (X_local - x).^2 + (Y_local - y).^2;
    maskSoma = distSq <= r^2;
    maskBG = distSq > (r*bgInnerFactor)^2 & distSq <= (r*bgOuterFactor)^2;
    
    if ~any(maskSoma(:))
        reason = "Empty soma mask";
        return;
    end
    
    if ~any(maskBG(:))
        reason = "Empty BG mask";
        return;
    end
    
    localPatch = image(ymin:ymax, xmin:xmax);
    
    % Get soma and background pixel values
    somaPixels = localPatch(maskSoma);
    bgPixels = localPatch(maskBG);
    bgPixels = bgPixels(~isnan(bgPixels));
    
    if isempty(bgPixels)
        reason = "No BG pixels";
        return;
    end
    
    % Calculate soma value
    if strcmpi(bgMethod, 'Median')
        somaVal = median(somaPixels, 'all', 'omitnan');
    else
        somaVal = mean(somaPixels, 'all', 'omitnan');
    end
    
    % Calculate robust background value (ignoring top 10% of pixels)
    prctile_90 = prctile(bgPixels, 90);
    robustBgPixels = bgPixels(bgPixels < prctile_90);
    if isempty(robustBgPixels)
        robustBgPixels = bgPixels; % Fallback if all pixels are above 90th percentile
    end
    
    if strcmpi(bgMethod, 'Median')
        bgVal = median(robustBgPixels, 'omitnan');
    else
        bgVal = mean(robustBgPixels, 'omitnan');
    end
    
    % Final intensity calculation
    if isnan(somaVal) || isnan(bgVal)
        reason = "NaN result";
        return;
    end
    
    intensity = somaVal - bgVal;
    
    if intensity <= 0
        reason = "Intensity <= 0";
    else
        reason = "Valid"; % Mark as valid if intensity is positive
    end

end