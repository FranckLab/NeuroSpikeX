function [CalciumCalculationHist] = NormalizedCalciumDataHistogram(data1, data2, data3, data4, num_timepoints, xAxisName, xUnit, binWidth, distribution)

%% Handle NaN values for Lognormal
if strcmp(distribution, 'Lognormal')
    data1(data1 == 0) = NaN;
    if num_timepoints >= 2, data2(data2 == 0) = NaN; end
    if num_timepoints >= 3, data4(data4 == 0) = NaN; end
    if num_timepoints == 4, data3(data3 == 0) = NaN; end
end

%% Organize datasets
datasets = {}; colors = {}; labels = {};
gray = [0.5 0.5 0.5];
switch num_timepoints
    case 1
        datasets = {data1}; colors = {'b'}; labels = {'0-'};
    case 2
        datasets = {data1, data2}; colors = {'b', 'r'}; labels = {'0-', '0+'};
    case 3
        datasets = {data1, data2, data4}; colors = {'b', 'r', 'g'}; labels = {'0-', '0+', '24 hr'};
    case 4
        datasets = {data1, data2, data3, data4}; colors = {gray, 'b', 'r', 'g'}; labels = {'0-', '0+', '1 hr', '24 hr'};
end

%% Compute global bin edges
allData = cellfun(@(x) x(~isnan(x(:))), datasets, 'UniformOutput', false);
allData = vertcat(allData{:});

% --- Critical Fix 1: Check for empty datasets ---
if isempty(allData)
    error('All input data is NaN/empty. Check data preprocessing.');
end

% --- Critical Fix 2: Validate binWidth ---
if ~isscalar(binWidth) || ~isnumeric(binWidth) || binWidth <= 0 || isnan(binWidth)
    error('binWidth must be a positive, finite, numeric scalar.');
end

% Extend range to cover ±4σ of all datasets
sigmaVals = cellfun(@(x) std(x(~isnan(x))), datasets);
sigmaVals(isnan(sigmaVals)) = 0;  % Handle NaN std from empty datasets
maxSigma = max(sigmaVals);

% --- Critical Fix 3: Ensure finite range ---
combinedMin = min(allData) - 4*maxSigma;
combinedMax = max(allData) + 4*maxSigma;

% --- Critical Fix 4: Handle edge cases for bin alignment ---
adjustedMin = floor(combinedMin / binWidth) * binWidth;
adjustedMax = ceil(combinedMax / binWidth) * binWidth;

% Force valid bin range if min/max collapse (e.g., identical data)
if adjustedMin >= adjustedMax
    adjustedMin = adjustedMin - binWidth;
    adjustedMax = adjustedMax + binWidth;
end

% --- Critical Fix 5: Final validation of bin edges ---
if isnan(adjustedMin) || isnan(adjustedMax)
    error('Adjusted bin edges are NaN. Check combinedMin/Max calculations.');
end

binEdges = adjustedMin:binWidth:adjustedMax;
if isempty(binEdges) || numel(binEdges) < 2
    binEdges = [adjustedMin, adjustedMin + binWidth]; % Force valid edges
end

%% Compute histograms and fits
CalciumCalculationHist = figure;
hold on;
legendEntries = {};

for i = 1:numel(datasets)
    data = datasets{i}(:);
    data = data(~isnan(data));
    if isempty(data), error('Dataset %d is empty.', i); end
    
    % Histogram with global bins
    counts = histcounts(data, binEdges, 'Normalization', 'pdf');
    binCenters = binEdges(1:end-1) + binWidth/2; 
    
    % Plot bars (adjacent, no gaps)
    bar(binCenters, counts, 1, 'FaceColor', colors{i}, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'BarWidth', 1);
    
    % Fit distribution and plot PDF
    pd = fitdist(data, distribution);
    x = linspace(min(binEdges), max(binEdges), 1000);
    y = pdf(pd, x);
    plot(x, y, 'Color', colors{i}, 'LineWidth', 2);
    
    legendEntries = [legendEntries, labels{i}, [distribution ' Fit']];
end

%% Finalize plot
xline([-1.0474, 1.4633], '--r', {'2.5%', '97.5%'},'LineWidth',1.5) % thresholds to classify increases and decreases (based on control data)
xlabel([xAxisName ' (' xUnit ')']);
ylabel('Probability Density');
legend(legendEntries, 'Location', 'best', 'Interpreter', 'none');
xlim([min(binEdges), max(binEdges)]);
box off;
hold off;
end

% function [CalciumCalculationHist] = NormalizedCalciumDataHistogram(data1, data2, data3, data4, num_timepoints, xAxisName, xUnit, binWidth, distribution)
% 
% %% discount any zero values in the data if the distribution is set to Lognormal
% if strcmp(distribution, 'Lognormal')
%     data1(data1==0) = NaN;
% 
%     if num_timepoints ~= 1
%         data2(data2==0) = NaN;
%     end
% 
%     if num_timepoints == 4
%         data3(data3==0) = NaN;
%     end
% 
%     if (num_timepoints == 3) || (num_timepoints == 4)
%         data4(data4==0) = NaN;
%     end
% end
% 
% %% Combine datasets into a cell array and define colors/labels
% gray = [0.5 0.5 0.5];
% 
% if num_timepoints == 1
%     datasets = {data1};
%     labels = {'0-'};
%     colors = {'b'};
% elseif num_timepoints == 2
%     datasets = {data1, data2};
%     labels = {'0-', '0+'};
%     colors = {'b', 'r'};
% elseif num_timepoints == 3
%     datasets = {data1, data2, data4};
%     labels = {'0-', '0+', '24 hr'};
%     colors = {'b', 'r', 'g'};
% elseif num_timepoints == 4
%     datasets = {data1, data2, data3, data4};
%     labels = {'0-', '0+', '1 hr', '24 hr'};
%     colors = {gray, 'b', 'r', 'g'};
% end
% 
% %% Remove NaNs from combined data to compute bin edges
% cellColumns = cellfun(@(x) x(:), datasets, 'UniformOutput', false);
% allData = vertcat(cellColumns{:});
% allData = allData(~isnan(allData)); % Exclude NaNs
% 
% if isempty(allData)
%     error('All data points are NaN. Check input data.');
% end
% 
% combinedMin = min(allData);
% combinedMax = max(allData);
% binEdges = combinedMin:binWidth:combinedMax;
% 
% %% Compute histogram counts and fit distributions
% numDatasets = numel(datasets);
% counts = cell(numDatasets, 1);
% maxCounts = zeros(numDatasets, 1);
% pds = cell(numDatasets, 1);
% 
% for i = 1:numDatasets
%     % Remove NaNs from each dataset before processing
%     data = datasets{i}(:);
%     data = data(~isnan(data));
% 
%     % Compute histogram counts
%     counts{i} = histcounts(data, binEdges);
%     maxCounts(i) = max(counts{i});
% 
%     % Fit distribution if data is not empty
%     if ~isempty(data)
%         pds{i} = fitdist(data, distribution);
%     else
%         error(['Dataset ', num2str(i), ' is empty after removing NaNs.']);
%     end
% end
% 
% %% Determine plotting order based on max counts
% [~, sortedIndices] = sort(maxCounts, 'ascend');
% plottingOrder = sortedIndices(end:-1:1);
% 
% %% Extend x-axis for distribution curves
% x_min_extended = combinedMin;
% x_max_extended = combinedMax;
% 
% for i = 1:numDatasets
%     pd = pds{i};
% 
%     if isa(pd, 'prob.NormalDistribution')
%         x_min = pd.mu - 4 * pd.sigma;
%         x_max = pd.mu + 4 * pd.sigma;
%     elseif isa(pd, 'prob.LognormalDistribution')
%         x_min = exp(pd.mu - 4 * pd.sigma);
%         x_max = exp(pd.mu + 4 * pd.sigma);
%     elseif isa(pd, 'prob.KernelDistribution')
%         x_min = combinedMin - 4 * pd.BandWidth;
%         x_max = combinedMax + 4 * pd.BandWidth;
%     else
%         padding = 0.1 * (combinedMax - combinedMin);
%         x_min = combinedMin - padding;
%         x_max = combinedMax + padding;
%     end
% 
%     x_min_extended = min(x_min_extended, x_min);
%     x_max_extended = max(x_max_extended, x_max);
% end
% 
% x = linspace(x_min_extended, x_max_extended, 1000);
% 
% %% Plot histograms and distribution fits
% CalciumCalculationHist = figure;
% hold on;
% 
% legendHandles = [];
% legendEntries = {};
% 
% for i = 1:numDatasets
%     idx = plottingOrder(i);
%     data = datasets{idx}(:);
%     data = data(~isnan(data));
%     color = colors{idx};
% 
%     % Plot histogram
%     h = histogram('BinEdges', binEdges, 'BinCounts', counts{idx}, ...
%         'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'Normalization', 'pdf');
% 
%     % Plot distribution fit
%     pd = pds{idx};
%     y = pdf(pd, x);
%     p = plot(x, y, 'Color', color, 'LineWidth', 2);
% 
%     % Store legend entries
%     legendHandles = [legendHandles, h, p];
%     legendEntries = [legendEntries, labels{idx}, [distribution ' Fit']];
% end
% 
% xlim([x_min_extended, x_max_extended]);
% box off;
% xlabel([xAxisName, ' (', xUnit, ')']);
% ylabel('Probability Density');
% legend(legendHandles, legendEntries, 'Location', 'best', 'Interpreter', 'none');
% legend boxoff;
% hold off;
% 
% end

% function [CalciumCalculationHist] = NormalizedCalciumDataHistogram(data1, data2, data3, data4, num_timepoints, xAxisName, xUnit, binWidth, distribution)
% 
% %% discount any zero values in the data if the distribution is set to Lognormal
% if strcmp(distribution, 'Lognormal')
%     data1(data1==0) = NaN;
% 
%     if num_timepoints ~= 1
%         data2(data2==0) = NaN;
%     end
% 
%     if num_timepoints == 4
%         data3(data3==0) = NaN;
%     end
% 
%     if (num_timepoints == 3) || (num_timepoints == 4)
%         data4(data4==0) = NaN;
%     end
% end
% 
% %% The following code was written with help from DeepSeek
% 
% %% SET POSSIBLE HISTOGRAM COLORS and combine datasets into a cell array
% 
% gray = [0.5 0.5 0.5]; %creating gray color
% 
% if num_timepoints == 1
%     color1 = 'b'; % 0-minus
%     datasets = {data1};
%     % labels = {'0-',[distribution,' Dist.']};
%     labels = {'0-'};
%     colors = {color1};
% elseif num_timepoints == 2
%     color1 = 'b';
%     color2 = 'r'; % 0-plus
%     datasets = {data1, data2};
%     labels = {'0-', '0+'};
%     colors = {color1,color2};
% elseif num_timepoints == 3
%     color1 = 'b';
%     color2 = 'r'; % 0-plus
%     color4 = 'g'; % 24hr
%     datasets = {data1, data2, data4};
%     labels = {'0-', '0+', '24 hr'};
%     colors = {color1,color2,color4};
% elseif num_timepoints == 4
%     color1 = gray;
%     color2 = 'b';
%     color3 = 'r'; % 1-hr
%     color4 = 'g';
%     datasets = {data1, data2, data3, data4};
%     labels = {'0-', '0+', '1 hr', '24 hr'};
%     colors = {color1,color2,color3,color4};
% end
% 
% %%
% 
% % Step 1: Convert each dataset to a column vector (even if already a column)
% cellColumns = cellfun(@(x) x(:), datasets, 'UniformOutput', false);
% 
% % Vertically concatenate all column vectors into a single numeric array
% allData = vertcat(cellColumns{:});
% 
% combinedMin = min(allData);
% combinedMax = max(allData);
% binEdges = combinedMin:binWidth:combinedMax;
% 
% % Step 2: Compute histogram counts and max counts for each dataset
% numDatasets = numel(datasets);
% counts = cell(numDatasets, 1);
% maxCounts = zeros(numDatasets, 1);
% 
% for i = 1:numDatasets
%     [counts{i}, ~] = histcounts(datasets{i}, binEdges);
%     maxCounts(i) = max(counts{i});
% end
% 
% % Step 3: Sort datasets by max counts (smallest last)
% [~, sortedIndices] = sort(maxCounts, 'ascend');
% plottingOrder = sortedIndices(end:-1:1); % Reverse to plot largest first
% 
% % Step 4: Fit distributions and extend x-axis for full curve visibility
% pds = cell(numDatasets, 1);
% x_min_extended = combinedMin;
% x_max_extended = combinedMax;
% 
% % Fit distributions and compute extended x-axis range
% for i = 1:numDatasets
%     data = datasets{i};
% 
%     % Fit a distribution (e.g., Normal, Lognormal, or Kernel)
%     pd = fitdist(data, distribution); % Replace with 'Lognormal' or 'Kernel' as needed
%     pds{i} = pd;
% 
%     % Extend x-axis based on distribution type
%     if isa(pd, 'prob.NormalDistribution')
%         % Normal: extend by ±4σ
%         x_min = pd.mu - 4 * pd.sigma;
%         x_max = pd.mu + 4 * pd.sigma;
%     elseif isa(pd, 'prob.LognormalDistribution')
%         % Lognormal: extend using log-space parameters
%         x_min = exp(pd.mu - 4 * pd.sigma);
%         x_max = exp(pd.mu + 4 * pd.sigma);
%     elseif isa(pd, 'prob.KernelDistribution')
%         % Kernel: extend by 4x bandwidth
%         x_min = combinedMin - 4 * pd.BandWidth;
%         x_max = combinedMax + 4 * pd.BandWidth;
%     else
%         % Default extension (10% padding)
%         padding = 0.1 * (combinedMax - combinedMin);
%         x_min = combinedMin - padding;
%         x_max = combinedMax + padding;
%     end
% 
%     % Update global extended range
%     x_min_extended = min(x_min_extended, x_min);
%     x_max_extended = max(x_max_extended, x_max);
% end
% 
% % Generate x values over the extended range
% x = linspace(x_min_extended, x_max_extended, 1000);
% 
% %% Step 5: Plot histograms and curves in sorted order
% CalciumCalculationHist = figure;
% hold on
% 
% % Loop through datasets in plotting order (largest max count first)
% for i = 1:numDatasets
%     idx = plottingOrder(i);
%     color = colors{idx};
% 
%     % Plot histogram (smallest max count will be on top)
%     histogram('BinEdges', binEdges, 'BinCounts', counts{idx}, ...
%         'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', 0.5,'Normalization', 'pdf');
% 
%     % Plot distribution curve scaled to counts
%     % data = datasets{idx};
%     pd = pds{idx};
%     y = pdf(pd, x);
%     plot(x, y, 'Color', color, 'LineWidth', 3);
% 
%     legend_order{1,i} = labels{idx};
%     % Store handles and labels for legend
%     % legendHandles = [legendHandles, h, p];
%     % legendEntries = [legendEntries, labels{idx}, ['Fit ' num2str(idx)]];
% end
% 
% if num_timepoints == 1
%     final_labels = {'Data',[distribution,' Dist.']};
% elseif num_timepoints == 2
%     final_labels = {legend_order{1,1},[distribution,' Dist.'],legend_order{1,2},[distribution,' Dist.']};
% elseif num_timepoints == 4
%     final_labels = {legend_order{1,1},[distribution,' Dist.'],legend_order{1,2},[distribution,' Dist.'],legend_order{1,3},[distribution,' Dist.'],legend_order{1,4},[distribution,'Dist.']};
% end
% 
% xlim([x_min_extended, x_max_extended]);
% box off
% xlabel([xAxisName,' (',xUnit,')'])
% ylabel('Probability Density (1/a.u.)')
% legend(final_labels, 'Location', 'best');
% legend boxoff
% hold off
% 
% end
% 
% 
% %% First DeepSeek code
% % % Step 1: Define consistent bin edges for both datasets
% % if num_timepoints == 1
% %     allData = data1(:);
% % elseif num_timepoints == 2
% %     allData = [data1(:); data2(:)];
% % elseif num_timepoints == 3
% %     allData = [data1(:); data2(:); data4(:)];
% % elseif num_timepoints == 4
% %     allData = [data1(:); data2(:); data3; data4(:)];
% % end
% % 
% % combinedMin = min(allData);
% % combinedMax = max(allData);
% % binEdges = combinedMin:binWidth:combinedMax;
% % 
% % % Step 2: Compute histogram counts (raw counts)
% % [counts1, ~] = histcounts(data1, binEdges);
% % 
% % % Step 3: Fit distributions (e.g., Normal)
% % pd1 = fitdist(data1, distribution);
% % 
% % % Step 4: Generate distribution curves SCALED TO COUNTS
% % % Extend x-axis by 10% of the data range
% % padding = 0.1 * (combinedMax - combinedMin);
% % x_min_extended = combinedMin - padding;
% % 
% % x = linspace(x_min_extended, combinedMax, 1000);
% % y1 = pdf(pd1, x) * numel(data1) * binWidth; % Scale PDF to match counts
% % 
% % if num_timepoints ~= 1
% %     [counts2, ~] = histcounts(data2, binEdges);
% %     pd2 = fitdist(data2, distribution);
% %     y2 = pdf(pd2, x) * numel(data2) * binWidth;
% % end
% % 
% % if num_timepoints == 3 || num_timepoints == 4
% %     [counts4, ~] = histcounts(data4, binEdges);
% %     pd4 = fitdist(data4, distribution);
% %     y4 = pdf(pd4, x) * numel(data4) * binWidth;
% % end
% % 
% % if num_timepoints == 4
% %     [counts3, ~] = histcounts(data3, binEdges);
% %     pd3 = fitdist(data3, distribution);
% %     y3 = pdf(pd3, x) * numel(data3) * binWidth;
% % end
% % 
% % % Step 5: Plot histograms and curves (no normalization)
% % figure;
% % hold on;
% % 
% % % Plot histograms (raw counts)
% % histogram('BinEdges', binEdges, 'BinCounts', counts1, ...
% %     'FaceColor', color1, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
% % % Overlay scaled distribution curves
% % plot(x, y1, color1, 'LineWidth', 3);
% % 
% % if num_timepoints ~= 1
% %     histogram('BinEdges', binEdges, 'BinCounts', counts2, ...
% %         'FaceColor', color2, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
% %     plot(x, y2, color2, 'LineWidth', 3);
% % end
% % 
% % if num_timepoints == 3 || num_timepoints == 4
% %     histogram('BinEdges', binEdges, 'BinCounts', counts4, ...
% %         'FaceColor', color4, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
% %     plot(x, y4, color4, 'LineWidth', 3);
% % end
% % 
% % if num_timepoints == 4
% %     histogram('BinEdges', binEdges, 'BinCounts', counts3, ...
% %         'FaceColor', color3, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
% %     plot(x, y3, color3, 'LineWidth', 3);
% % end
% % 
% % hold off;
% % legend('0-',[distribution,' Dist.'],'0+',[distribution,' Dist.'],'1 hr',[distribution,' Dist.'],'24 hr',[distribution,' Dist.'],'location','northeast')
% % legend boxoff
% % 
% % xlim([x_min_extended, combinedMax]);
% % xlabel([xAxisName,' (',xUnit,')'])
% % ylabel('Number of Active Cells')
% % 
% % end
% 
% %% My own original code
% % %% SET POSSIBLE HISTOGRAM COLORS
% % 
% % gray = [0.5 0.5 0.5]; %creating gray color
% % 
% % if num_timepoints == 1
% %     color0 = 'b'; % 0-minus
% % elseif num_timepoints == 2
% %     color0 = 'b';
% %     color1 = 'r'; % 0-plus
% % elseif num_timepoints == 3
% %     color0 = 'b';
% %     color1 = 'r';
% %     color3 = 'g'; % 24hr
% % elseif num_timepoints == 4
% %     color0 = gray;
% %     color1 = 'b';
% %     color2 = 'r'; % 1-hr
% %     color3 = 'g';
% % end
% % 
% % %% discount any zero values in the data if the distribution is set to Lognormal
% % if strcmp(distribution, 'Lognormal')
% %     data_0minus(data_0minus==0) = NaN;
% % 
% %     if num_timepoints ~= 1
% %         data_0plus(data_0plus==0) = NaN;
% %     end
% % 
% %     if num_timepoints == 4
% %         data_1hr(data_1hr==0) = NaN;
% %     end
% % 
% %     if (num_timepoints == 3) || (num_timepoints == 4)
% %         data_24hr(data_24hr==0) = NaN;
% %     end
% % end
% % 
% % % Calculating number of bins needed to achieve the user-set bin width
% % b0 = histogram(data_0minus, 'binWidth', bin_width);
% % binCount0 = b0.NumBins;
% % 
% % if num_timepoints ~= 1
% %     b1 = histogram(data_0plus, 'binWidth',bin_width);
% %     binCount1 = b1.NumBins;
% % end
% % 
% % if num_timepoints == 4
% %     b2 = histogram(data_1hr, 'binWidth',bin_width);
% %     binCount2 = b2.NumBins;
% % end
% % 
% % if (num_timepoints == 3) || (num_timepoints == 4)
% %     b3 = histogram(data_24hr, 'binWidth',bin_width);
% %     binCount3 = b3.NumBins;
% % end
% % 
% % %% Plot comparison of the histogram of the data, and the fit
% % CalciumCalculationHist = figure;
% % 
% %     h0 = histfit(data_0minus,binCount0,distribution);
% %     h0(1).FaceColor = color0;
% %     h0(1).FaceAlpha = .5;
% %     h0(1).EdgeColor = 'none';
% %     h0(2).Color = color0;
% %     h0(2).LineWidth = 3;
% % 
% % if num_timepoints ~= 1
% %     hold on
% %     h1 = histfit(data_0plus,binCount1,distribution);
% %     h1(1).FaceColor = color1;
% %     h1(2).Color = color1;
% %     h1(1).FaceAlpha = .5;
% %     h1(1).EdgeColor = 'none';
% %     h1(2).LineWidth = 3;
% % end
% % 
% % if num_timepoints == 4    
% %     hold on
% %     h2 = histfit(data_1hr, binCount2,distribution);
% %     h2(1).FaceColor = color2;
% %     h2(1).FaceAlpha = .5;
% %     h2(1).EdgeColor = 'none';
% %     h2(2).Color = color2;
% %     h2(2).LineWidth = 3;
% % end
% % 
% % if (num_timepoints == 3) || (num_timepoints == 4)
% %     hold on
% %     h3 = histfit(data_24hr, binCount3,distribution);
% %     h3(1).FaceColor = color3;
% %     h3(1).FaceAlpha = .5;
% %     h3(1).EdgeColor = 'none';
% %     h3(2).Color = color3;
% %     h3(2).LineWidth = 3;
% % end
% % 
% % box off
% % legend('0-',[distribution,' Dist.'],'0+',[distribution,' Dist.'],'1 hr',[distribution,' Dist.'],'24 hr',[distribution,' Dist.'],'location','northeast')
% % legend boxoff
% % 
% % xlim([0 inf])
% % xlabel([xAxisName,' (',xUnit,')'])
% % ylabel('Number of Active Cells')
% % hold off