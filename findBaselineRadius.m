function newRadius = findBaselineRadius(image_patch, expectedRadius)
    % FINDBASELINERADIUS attempts to find a better radius for a dim cell.
    % This version is updated to suppress the 'smallRadius' warning.

    % This function is for the RelativeBaselineChangeApp
    
    % Default value in case of failure
    newRadius = expectedRadius;
    
    % Define a plausible search range
    minR = max(3, floor(expectedRadius * 0.5));
    maxR = min(20, ceil(expectedRadius * 1.5));
    
    % If the range is invalid, exit early
    if minR >= maxR
        return;
    end

    % --- WRAP THE CALL TO imfindcircles WITH WARNING SUPPRESSION ---
    
    % Capture the original warning state
    original_warning_state = warning('query', 'images:imfindcircles:smallRadius');
    
    % Turn off the specific warning
    warning('off', 'images:imfindcircles:smallRadius');
    
    try
        % Use imfindcircles, a powerful but sensitive function.
        [centers, radii] = imfindcircles(image_patch, [minR maxR], ...
            'ObjectPolarity', 'bright', 'Sensitivity', 0.92, 'EdgeThreshold', 0.05);
        
        if ~isempty(radii)
            % If circles are found, find the one closest to the center of the patch.
            patch_center = size(image_patch)/2;
            distances = sqrt((centers(:,1) - patch_center(2)).^2 + (centers(:,2) - patch_center(1)).^2);
            [~, closest_idx] = min(distances);
            newRadius = radii(closest_idx);
        end
    catch ME
        % If imfindcircles fails, don't crash. Just use the original radius.
        fprintf(2, 'Error in findBaselineRadius: %s. Using original radius.\n', ME.message);
    end
    
    % --- Restore the original warning state ---
    warning(original_warning_state.state, 'images:imfindcircles:smallRadius');
    % --- END OF WRAPPER ---
    
end