function [calcium_spike_calculations, decay_rates, CWT, calcium_spikes_direct_output] = processTimepoint(rootFolder, tag, fdata, settings,  timepoint1_rate, timepoint1_intensity, timepoint1_decayrate, timepoint1_CWT_max)

% ----- FUNCTION THAT RUNS ALL ANALYSES ON DATA FROM ONE TIMEPOINT ----- %
% INPUTS:
% tag = name of timepoint (e.g. '0minus','0plus', etc.)
%
% fdata = dF/F exported by NeuroCa
%
% settings = contains all settings determined in the
% Preprocessing_Settings_App - signal strength threshold (i.e. noise
% floor), smoothing factor (n), signal boundaries, minimum peak prominence
%
% timepoint1_rate = average network spike rate of baseline (timepoint 1)
%
% timepoint1_intensity = average network spike intensity of baseline (timepoint 1)
%
% timepoint1_decayrate = average network decay rate constant of baseline (timepoint 1)
%
% timepoint1_CWT_max = max CWT magnitude of baseline data (timepoint 1)

    %% Step 1: Extract settings set previously in the Preprocessing_Settings_App
    desired_duration_seconds = settings.analysisDurationSeconds;
    n = settings.n;
    prominence = settings.ProminenceValue;
    signal_strength_threshold = settings.signal_strength_threshold;
    useBounds = settings.useBounds;
    minBound = settings.minBound;
    maxBound = settings.maxBound;
    time_axis = fdata(1,:);
    
    %% Step 2: Pre-process signals and identify calcium spikes (i.e. transients)
    signal_info = ProcessNeuroCa(fdata, desired_duration_seconds);
    smooth_signal = smoothCalciumSignals(signal_info,n);

    % Determine active vs. inactive cells, and detect calcium transient peaks in the active cells
    calcium_spikes = findCalciumSpikes(signal_info,smooth_signal,signal_strength_threshold,useBounds,minBound,maxBound,prominence);
    calcium_spikes_direct_output = calcium_spikes;
    
    % Save main data into the Current Folder
    analysis_struct = struct();
    analysis_struct.(['signal_info_' tag]) = signal_info;
    analysis_struct.(['smooth_signal_' tag]) = smooth_signal;
    analysis_struct.(['calcium_spikes_' tag]) = calcium_spikes;
    analysis_struct.(['signal_strength_threshold_' tag]) = signal_strength_threshold; % Corrected

    save(fullfile(rootFolder,['analyzed_fdata_' tag '.mat']), '-struct', 'analysis_struct');

    %% Step 3: Quantify all the parameters from the calcium transients

    calcium_spike_calculations = quantifyNetworkCalciumSpikes(signal_info,calcium_spikes,timepoint1_rate,timepoint1_intensity);
    [decay_rates, decay_fig] = decayRateConstant(signal_info,calcium_spikes,smooth_signal,timepoint1_decayrate,signal_strength_threshold);
    savefig(decay_fig, fullfile(rootFolder,['RateConstant_Hist_' tag '.fig']));

    %% Step 4: Calculate and plot raster and detect network bursts
    [raster_data, raster_fig] = RasterPlot(smooth_signal,calcium_spikes,signal_info,tag);
    savefig(raster_fig, fullfile(rootFolder,['Raster_' tag '.fig']))
    saveas(raster_fig, fullfile(rootFolder,['Raster_' tag '.png']))

    [burst_data, active_fig, burst_fig] = countNetworkBursts(raster_data,time_axis,tag);
    savefig(active_fig, fullfile(rootFolder,['ActiveCellsPerFrame_' tag '.fig']))
    saveas(active_fig, fullfile(rootFolder,['ActiveCellsPerFrame_' tag '.png']))
    savefig(burst_fig, fullfile(rootFolder,['Burst_raster_' tag '.fig']))
    saveas(burst_fig, fullfile(rootFolder,['Burst_raster_' tag '.png']))

    %% Step 5: Calculate continuous wavelet transform for all active cells in the network & normalize to 0minus
    [CWT, scalogram_fig, norm_scalogram_fig] = wavelet(smooth_signal,signal_info, calcium_spikes, tag, timepoint1_CWT_max);
    savefig(scalogram_fig, fullfile(rootFolder,['MeanNetworkScalogram_' tag '.fig']))
    saveas(scalogram_fig, fullfile(rootFolder,['MeanNetworkScalogram_' tag '.png']))
    savefig(norm_scalogram_fig, fullfile(rootFolder,['Normalized_MeanNetworkScalogram_' tag '.fig']))
    saveas(norm_scalogram_fig, fullfile(rootFolder,['Normalized_MeanNetworkScalogram_' tag '.png']))

    %% Step 6: Save remaining outputs with tagged timepoint name in the Current Folder
    tag_struct = struct();
    tag_struct.(['calcium_spike_calculations_' tag]) = calcium_spike_calculations;
    tag_struct.(['rate_constants_' tag]) = decay_rates;
    tag_struct.(['raster_' tag]) = raster_data;
    tag_struct.(['network_bursts_' tag]) = burst_data;
    tag_struct.(['CWT_' tag]) = CWT;

    save(fullfile(rootFolder,['all_quantitative_analyses_' tag '.mat']), '-struct', 'tag_struct', '-v7.3');

    %% Step 7: Summary Report Box
    ResultSummaryBox(calcium_spike_calculations,calcium_spikes,decay_rates,burst_data,upper(tag))

    close all % close all figures after saving
end