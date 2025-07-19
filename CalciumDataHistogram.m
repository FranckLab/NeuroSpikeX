function [CalciumCalculationHist] = CalciumDataHistogram(datasets_cell, dataset_labels, xAxisName, xUnit, binWidth, distribution)
% CalciumDataHistogram: Generates histograms and fitted distributions for calcium imaging data.
%
% Inputs:
%   datasets_cell  - Cell array, where each cell contains a vector of data for one timepoint/condition.
%   dataset_labels - Cell array of strings, labels for each dataset corresponding to datasets_cell.
%   xAxisName      - String, label for the x-axis (e.g., 'Spike Rate').
%   xUnit          - String, unit for the x-axis data (e.g., 'spikes/sec').
%   binWidth       - Scalar, width of histogram bins.
%   distribution   - String, type of distribution to fit (e.g., 'Normal', 'Lognormal', 'Kernel').

CalciumCalculationHist = []; % Initialize to empty

% Number of datasets to plot
num_datasets = numel(datasets_cell);

% Check if any datasets were provided
if num_datasets == 0
    warning('CalciumDataHistogram: No datasets provided for plotting.');
    CalciumCalculationHist = figure;
    text(0.5, 0.5, 'No data provided for histogram.', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    title('No Data');
    return;
end

% Check if each dataset has enough data points (more than 1 non-NaN value)
sufficient_data_flags = cellfun(@(d) sum(~isnan(d(:))) > 1, datasets_cell);
if ~all(sufficient_data_flags)
    warning('CalciumDataHistogram: One or more datasets have insufficient data (<=1 non-NaN point). Histogram cannot be generated reliably.');
    CalciumCalculationHist = figure;
    text(0.5, 0.5, 'Insufficient data in one or more sets.', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    title('Insufficient Data');
    return;
end

datasets = datasets_cell; % Use the validated cell array

%% Discount any zero values in the data if the distribution is set to Lognormal
if strcmpi(distribution, 'Lognormal')
    for i = 1:num_datasets
        if ~isempty(datasets{i})
            original_label = 'Unknown';
            if i <= numel(dataset_labels), original_label = dataset_labels{i}; end

            datasets{i}(datasets{i} == 0) = NaN;
            if sum(~isnan(datasets{i}(:))) <= 1
                 warning('CalciumDataHistogram: Dataset "%s" (index %d) became insufficient after removing zeros for Lognormal fit.', original_label, i);
                 CalciumCalculationHist = figure;
                 text(0.5, 0.5, sprintf('Dataset "%s" insufficient after Lognormal zero removal.', original_label), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                 title('Insufficient Data Post-Processing');
                 return;
            end
        end
    end
end

%% SET HISTOGRAM COLORS and USE PROVIDED LABELS

% Use provided labels, with a fallback for mismatched counts
labels = dataset_labels;
if numel(labels) ~= num_datasets
    warning('CalciumDataHistogram: Mismatch between number of datasets (%d) and provided labels (%d). Using generic labels for some/all.', num_datasets, numel(labels));
    temp_labels = cell(1, num_datasets);
    for i = 1:num_datasets
        if i <= numel(labels) && ~isempty(labels{i})
            temp_labels{i} = labels{i}; % Keep original if available
        else
            temp_labels{i} = sprintf('Dataset %d', i); % Fallback
        end
    end
    labels = temp_labels;
end

% Define colors
color_palette = { % A list of distinguishable colors
    [0 0 1],         % Blue
    [1 0 0],         % Red
    [0 0.7 0],       % Green (darker)
    [0.75 0.75 0],   % Yellow (darker)
    [0 0.75 0.75],   % Cyan (darker)
    [0.75 0 0.75],   % Magenta (darker)
    [0.5 0.5 0.5],   % Gray
    [1 0.5 0],       % Orange
    [0.5 0 0],       % Maroon
    [0 0.5 0],       % Dark Green
    [0 0 0.5],       % Navy
    [0.6 0.2 0.8]    % Purple
};
colors = cell(1, num_datasets);
for i = 1:num_datasets
    colors{i} = color_palette{mod(i-1, numel(color_palette)) + 1};
end

% If you want to keep your specific colors for N=1 to 4, you can add that logic here:
if num_datasets == 1
    colors = {[0 0 1]}; % Blue
elseif num_datasets == 2
    colors = {[0 0 1], [1 0 0]}; % Blue, Red
elseif num_datasets == 3
    colors = {[0 0 1], [1 0 0], [0 0.7 0]}; % Blue, Red, Green
elseif num_datasets == 4
    colors = {[0.5 0.5 0.5], [0 0 1], [1 0 0], [0 0.7 0]}; % Gray, Blue, Red, Green
end


%% CORE PLOTTING LOGIC (largely unchanged from previous version)

% Step 1: Convert each dataset to a column vector
cellColumns = cellfun(@(x) x(:), datasets, 'UniformOutput', false);
allData = vertcat(cellColumns{:});
allData_finite = allData(~isnan(allData) & ~isinf(allData));

if isempty(allData_finite)
    warning('CalciumDataHistogram: All data is NaN or Inf after processing. Cannot determine range.');
    CalciumCalculationHist = figure;
    text(0.5, 0.5, 'All data is NaN/Inf.', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    title('Invalid Data Range');
    return;
end

combinedMin = min(allData_finite);
combinedMax = max(allData_finite);

if combinedMin == combinedMax
    if binWidth > 0
        padding = binWidth * 2; % Ensure some width
        if padding == 0; padding = max(eps(combinedMin), eps('double')); end
        binEdges = [combinedMin - padding/2, combinedMin + padding/2];
    else
        warning('CalciumDataHistogram: binWidth is not positive, and all data points are identical. Cannot create bins.');
        CalciumCalculationHist = figure; text(0.5,0.5,'Cannot bin identical data with given binWidth.','HorizontalAlignment','center');
        return;
    end
else
    binEdges = combinedMin:binWidth:combinedMax;
    if numel(binEdges) < 2
        binEdges = [combinedMin, combinedMax];
        if binEdges(1) == binEdges(2)
             delta = max(eps(combinedMin), eps('double'));
             binEdges = [combinedMin - delta, combinedMax + delta];
        end
    end
end

% Step 2: Compute histogram counts and max counts
counts = cell(num_datasets, 1);
maxCounts = zeros(num_datasets, 1);
for i = 1:num_datasets
    [counts{i}, ~] = histcounts(datasets{i}, binEdges);
    maxCounts(i) = max(counts{i});
end

% Step 3: Sort datasets by max counts (plot largest first for visibility)
[~, sortedIndices] = sort(maxCounts, 'descend');
plottingOrder = sortedIndices;

% Step 4: Fit distributions
pds = cell(num_datasets, 1);
for i = 1:num_datasets
    data_for_fit = datasets{i};
    data_for_fit = data_for_fit(~isnan(data_for_fit) & ~isinf(data_for_fit));
    current_label = labels{i}; % Use the actual label for warning

    if isempty(data_for_fit) || numel(data_for_fit) <= 1 || (numel(unique(data_for_fit)) == 1 && numel(data_for_fit) > 1 && ~strcmpi(distribution, 'Kernel'))
        warning('CalciumDataHistogram: Dataset "%s" has insufficient or constant data for reliable %s distribution fitting. Skipping fit.', current_label, distribution);
        pds{i} = []; continue;
    end
    
    try
        pd = fitdist(data_for_fit, distribution);
        pds{i} = pd;
    catch ME
        warning('CalciumDataHistogram: Could not fit %s distribution to dataset "%s". Error: %s. Skipping fit.', distribution, current_label, ME.message);
        pds{i} = [];
    end
end

% X-values for PDF curves
x_padding = (combinedMax - combinedMin) * 0.05;
if x_padding == 0, x_padding = max(binWidth,1); end
x_min_plot = combinedMin - x_padding;
x_max_plot = combinedMax + x_padding;

if strcmpi(distribution, 'Lognormal') && x_min_plot <= 0
    positive_data = allData_finite(allData_finite > 0);
    if ~isempty(positive_data)
        x_min_plot = min(positive_data)/2;
    else
        x_min_plot = eps('double'); % Smallest positive if no positive data
    end
    if x_min_plot <=0 ; x_min_plot = eps('double'); end
end
if x_min_plot >= x_max_plot
    delta_adjust = max(abs(combinedMin*0.1) + eps, 1);
    x_min_plot = combinedMin - delta_adjust;
    x_max_plot = combinedMax + delta_adjust;
    if x_min_plot == x_max_plot
        x_min_plot = x_min_plot -1; x_max_plot = x_max_plot+1;
    end
end
x_curve = linspace(x_min_plot, x_max_plot, 1000);

%% Step 5: Plot histograms and curves
CalciumCalculationHist = figure;
hold on;
plottedHandles = gobjects(num_datasets * 2, 1);
legendLabelsList = {};
handle_idx = 0;

for i_plot_order = 1:num_datasets
    original_idx = plottingOrder(i_plot_order); % Actual index in 'datasets', 'labels', 'colors'
    
    h_hist = histogram('BinEdges', binEdges, 'BinCounts', counts{original_idx}, ...
        'FaceColor', colors{original_idx}, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    
    is_hist_plotted = any(counts{original_idx} > 0);
    if is_hist_plotted
        handle_idx = handle_idx + 1;
        plottedHandles(handle_idx) = h_hist;
        legendLabelsList{handle_idx} = labels{original_idx}; % Use the passed/determined label
    end
    
    if ~isempty(pds{original_idx})
        pd_current = pds{original_idx};
        data_for_scaling = datasets{original_idx};
        data_for_scaling = data_for_scaling(~isnan(data_for_scaling) & ~isinf(data_for_scaling));
        if ~isempty(data_for_scaling)
            y_curve = pdf(pd_current, x_curve) * numel(data_for_scaling) * binWidth;
            h_fit = plot(x_curve, y_curve, 'Color', colors{original_idx}, 'LineWidth', 3);
             % Add fit to legend only if its corresponding hist was plotted or if fit itself is substantial
            if is_hist_plotted || any(y_curve > max(y_curve)*0.01) % Only add if fit is meaningful
                handle_idx = handle_idx + 1;
                plottedHandles(handle_idx) = h_fit;
                legendLabelsList{handle_idx} = [labels{original_idx} ' (' distribution ' Fit)']; % Use label
            end
        end
    end
end
plottedHandles = plottedHandles(1:handle_idx); % Trim unused handles

%% Configure legend and axes
if ~isempty(plottedHandles) && ~isempty(legendLabelsList)
    legend(plottedHandles, legendLabelsList, 'Location', 'best', 'Interpreter', 'none'); % Added Interpreter none for safety with tags
    legend boxoff;
end

xlim([x_min_plot, x_max_plot]);
box off;
xlabel([xAxisName,' (',xUnit,')']);
ylabel('Number of Active Cells');
title([xAxisName ' Distribution'], 'Interpreter', 'none');
hold off;

end