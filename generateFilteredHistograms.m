% generateFilteredHistograms.m

function generateFilteredHistograms(all_calcium_spike_calculations, all_decay_rates, timepointTags_loaded, rootFolder)
% generateFilteredHistograms: Creates filtered histograms for spike rate and decay rate.
% This function takes the results of a primary calcium imaging analysis,
% filters out statistical outliers from the data using the IQR method, and
% then generates new histograms. This is useful for improving the fit of
% distributions to the central tendency of the data.
%
% SYNTAX:
%   generateFilteredHistograms(all_calcium_spike_calculations, all_decay_rates, timepointTags_loaded, rootFolder)
%
% INPUTS:
%   all_calcium_spike_calculations - Cell array from main analysis containing spike rate info.
%   all_decay_rates                - Cell array from main analysis containing decay rate info.
%   timepointTags_loaded           - Cell array of strings with timepoint names/labels.
%   rootFolder                     - String path to the sample's root folder for saving plots.

    %% --- Initial Checks and User Prompt ---
    fprintf('\n---------------------------------------------------------\n');
    fprintf('OPTIONAL: FILTERED HISTOGRAM GENERATION\n');
    fprintf('---------------------------------------------------------\n');

    % Check for necessary variables (internal validation)
    if nargin < 4
        error('This function requires 4 input arguments. Please call it as shown in the documentation.');
    end

    % Ask the user what they want to do
    choice = menu('Which filtered histogram do you want to generate?', ...
                  'Spike Rate Histogram', ...
                  'Decay Rate Histogram', ...
                  'Both', ...
                  'Cancel');

    if choice == 0 || choice == 4
        fprintf('Operation cancelled.\n');
        return;
    end

    %% --- Main Logic ---

    % Spike Rate Processing
    if choice == 1 || choice == 3
        fprintf('\n--- Processing Filtered Spike Rate ---\n');
        config.data_field = 'cell_spike_rate';
        config.plot_title = 'Spike Rate (Filtered)';
        config.x_label = 'Spike Rate';
        config.x_unit = 'spikes/sec';
        config.bin_width = 0.006;
        config.distribution = 'Lognormal';
        config.save_filename_base = 'SpikeRateHistogram_Filtered';
        processAndPlot(all_calcium_spike_calculations, timepointTags_loaded, rootFolder, config);
    end

    % Decay Rate Processing
    if choice == 2 || choice == 3
        fprintf('\n--- Processing Filtered Decay Rate ---\n');
        config.data_field = 'cell_avg_rate_cst';
        config.plot_title = 'Decay Rate Constant (Filtered)';
        config.x_label = 'Decay Rate Constant';
        config.x_unit = 's^{-1}';
        config.bin_width = 0.02;
        config.distribution = 'Normal';
        config.save_filename_base = 'DecayRateConstantHistogram_Filtered';
        processAndPlot(all_decay_rates, timepointTags_loaded, rootFolder, config);
    end

    fprintf('\n---------------------------------------------------------\n');
    fprintf('FILTERED HISTOGRAM GENERATION COMPLETE.\n');
    fprintf('---------------------------------------------------------\n');

end


%% --- LOCAL FUNCTION FOR THE CORE LOGIC ---
function processAndPlot(all_results_cell, timepoint_labels, sample_root_folder, config)
    % This local function extracts data, filters outliers, and plots the histogram.

    num_timepoints = numel(all_results_cell);

    % Step 1: Extract valid data
    plot_inputs = cell(1, num_timepoints);
    labels_for_hist = cell(1, num_timepoints);

    for k = 1:num_timepoints
        if ~isempty(all_results_cell{k}) && isfield(all_results_cell{k}, config.data_field)
            current_data = all_results_cell{k}.(config.data_field);
            if ~isempty(current_data) && sum(~isnan(current_data(:))) > 1
                plot_inputs{k} = current_data;
                labels_for_hist{k} = timepoint_labels{k};
            end
        end
    end
    
    valid_indices = ~cellfun('isempty', plot_inputs);
    if ~any(valid_indices)
        warning('No valid data found for "%s". Cannot generate plot.', config.plot_title);
        return;
    end
    valid_plot_inputs = plot_inputs(valid_indices);
    valid_labels = labels_for_hist(valid_indices);
    
    % Step 2: Filter outliers using IQR
    fprintf('Filtering outliers...\n');
    plot_inputs_filtered = cell(size(valid_plot_inputs));
    for i = 1:numel(valid_plot_inputs)
        data = valid_plot_inputs{i};
        data_nan_removed = data(~isnan(data));

        if numel(data_nan_removed) < 4
            plot_inputs_filtered{i} = data;
            fprintf('  - Timepoint "%s": Too few data points to filter, using original.\n', valid_labels{i});
            continue;
        end

        Q1 = prctile(data_nan_removed, 25);
        Q3 = prctile(data_nan_removed, 75);
        IQR = Q3 - Q1;
        lowerBound = Q1 - 1.5 * IQR;
        upperBound = Q3 + 1.5 * IQR;

        original_count = numel(data_nan_removed);
        filtered_data = data_nan_removed(data_nan_removed >= lowerBound & data_nan_removed <= upperBound);
        filtered_count = numel(filtered_data);

        plot_inputs_filtered{i} = filtered_data;
        fprintf('  - Timepoint "%s": Removed %d outlier(s). Kept %d of %d points.\n', ...
                valid_labels{i}, original_count - filtered_count, filtered_count, original_count);
    end

    % Step 3: Plot histogram using the original CalciumDataHistogram function
    [FilteredHist] = CalciumDataHistogram(plot_inputs_filtered, valid_labels, ...
        config.plot_title, config.x_unit, config.bin_width, config.distribution);
                                       
    if ~isempty(FilteredHist) && isgraphics(FilteredHist, 'figure')
        fig_path = fullfile(sample_root_folder, [config.save_filename_base, '.fig']);
        png_path = fullfile(sample_root_folder, [config.save_filename_base, '.png']);
        savefig(FilteredHist, fig_path);
        saveas(FilteredHist, png_path);
        fprintf('Filtered histogram saved to:\n  %s\n', fig_path);
    else
        warning('Filtered histogram for "%s" was not generated.', config.plot_title);
    end
end