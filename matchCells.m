function [MatchedCells, cellMatchesFig] = matchCells(folderpath, image1_folder, image2_folder, tform, label1, label2)

% Function to match the cells identified by NeuroCa between multiple
% timepoints for one sample (currently written for 4 timepoints)

% This function is for the RelativeBaselineChangeApp

% If labels are not provided, use default names for backward compatibility.
    if nargin < 5 || isempty(label1)
        % Create a default label from the folder name if needed
        label1 = image1_folder;
        if startsWith(label1, '\') || startsWith(label1, '/')
            label1 = label1(2:end);
        end
    end
    if nargin < 6 || isempty(label2)
        % Create a default label from the folder name if needed
        label2 = image2_folder;
        if startsWith(label2, '\') || startsWith(label2, '/')
            label2 = label2(2:end);
        end
    end
%% Load data from NeuroCa calculations

% Load data from first timepoint
centroids1 = importdata([folderpath,image1_folder,'\center.mat']);
radii1 = importdata([folderpath,image1_folder,'\radii.mat']);
fdata1 = importdata([folderpath,image1_folder,'\fdata.mat']);

% Load data from the second timepoint being compared
centroids2 = importdata([folderpath,image2_folder,'\center.mat']);
radii2 = importdata([folderpath,image2_folder,'\radii.mat']);
fdata2 = importdata([folderpath,image2_folder,'\fdata.mat']);
image1 = im2double(imread([folderpath,image1_folder,'\beforeraw.tif']));
image2 = im2double(imread([folderpath,image2_folder,'\afterraw.tif']));

%% Create transform matrices by manually matching point between the 0minus image and all other timepoints
if isempty(tform)
    %% --- Step 1: Manual Point Match ---
    % Get matching points from time point 1 and time point 2 images to account for any shifts between
    % before and after images
    numPoints = 6;  % Minimum 3 points required
    fprintf('Select %d points in Image 1 (left) first, then %d in Image 2 (right).\n', numPoints, numPoints);
    [points1, points2] = cpselect(image1, image2, 'Wait', true);

   
    % --- Step 2: Compute Transform & Validate ---
    tform = fitgeotrans(points2, points1, 'affine');

    % Calculate reprojection error on MANUAL POINTS (not centroids)
    points2_aligned = transformPointsForward(tform, points2);  % Transform manual points
    error_stretch = mean(sqrt(sum((points1 - points2_aligned).^2, 2)));
    fprintf('Mean alignment error (time1-time2): %.2f pixels\n', error_stretch);
    if error_stretch > 5
        warning('High alignment error between time point 1 and 2 - recheck manual points.');
    end
    
end

%% --- Step 3: Transform images and centroids ---
image2_aligned = imwarp(image2, tform, 'OutputView', imref2d(size(image1)));
centroids2_aligned = transformPointsForward(tform, centroids2);  % Applies to all centroids

%% --- Step 4: Match Somas ---
% maxDist = 10;
maxDist = 5;

% Comparing timepoints 1 and 2
% idx = sample 2 cell indices identified as the closest neighbor to each
% cell in sample 1
% dist = distance between the identified neighbors
[idx, dist] = knnsearch(centroids2_aligned, centroids1, 'K', 1);
validMatches = dist < maxDist; % only neighbors within the set max distance will be counted

% Extract matched pairs
matchedCentroids1 = centroids1(validMatches, :);
findIndices1 = 1:size(centroids1, 1);
matchedIndices(:,1) = findIndices1(validMatches)'; % cell numbers from 0minus

matchedIndices(:,2) = idx(validMatches); % Indices in ORIGINAL centroids2
matchedCentroids2 = centroids2(matchedIndices(:,2), :); % Original Image2 coordinates
matchedRadii1 = radii1(validMatches);
matchedRadii2 = radii2(matchedIndices(:,2)); % From original Image2

%% --- Step 5: Visualize Alignment ---
% Visualize alignment and matched somas
cellMatchesFig = figure;
imshowpair(image1, image2_aligned, 'falsecolor');
hold on;

% Plot Image1 centroids (before stretch)
plot(matchedCentroids1(:,1), matchedCentroids1(:,2), 'ro', ...
    'MarkerSize', 10, 'LineWidth', 1.5);

% Plot Image2 centroids (after stretch, aligned)
% Use idx(validMatches) to get correct indices into centroids2_aligned
matchedIndices2 = idx(validMatches);
plot(centroids2_aligned(matchedIndices2,1), centroids2_aligned(matchedIndices2,2), 'g+', ...
    'MarkerSize', 10, 'LineWidth', 1.5);

% Draw lines between matched pairs
for i = 1:size(matchedCentroids1,1)
    plot([matchedCentroids1(i,1), centroids2_aligned(matchedIndices2(i),1)], ...
         [matchedCentroids1(i,2), centroids2_aligned(matchedIndices2(i),2)], ...
         'b-', 'LineWidth', 1);
end

hold off;
title(['Aligned Images and Soma Centroids (',num2str(size(matchedCentroids1,1)),' Cells Matched)']);
legend(label1, [label2,' (Aligned)'], 'Matches');

%% Save variables from matched pairs

MatchedCells.matchedCellNumbers = matchedIndices;
MatchedCells.matchedCentroidCoordinates_0minus = matchedCentroids1;
MatchedCells.matchedRadii_0minus = matchedRadii1;
MatchedCells.matchedCentroidCoordinates_timepoint2 = matchedCentroids2; % Original Image2 coordinates
MatchedCells.matchedRadii_timepoint2 = matchedRadii2; % From original Image2
MatchedCells.transformMatrix = tform;
MatchedCells.validMatches_image1 = validMatches; % logical whether each cell in sample 1 found a valid match
MatchedCells.GCaMP_Image1 = image1;
MatchedCells.GCaMP_Image2 = image2;

end