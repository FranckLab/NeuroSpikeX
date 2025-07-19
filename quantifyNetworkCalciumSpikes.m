function [calcium_spike_calculations] = quantifyNetworkCalciumSpikes(signal_info, calcium_transients, avg_network_spike_rate_0minus, avg_network_spike_intensity_0minus)

%% Function to calculate average cell and network spike intensity and firing rate

% INPUTS:
% true_pks = amplitude of each calcium transient peak in each active cell
%
% time_axis = time (seconds)
%
% f_end = final frame # at 120 seconds
%
% avg_network_spike_rate_0minus =
%
% avg_network_spike_intensity = 
%
% OUTPUTS - STRUCTURE "CALCIUM_SPIKE_CALCULATIONS:
% cell_spike_freq = # of spikes from each active cell
%
% cell_spike_rate = spiking rate of each active cell
%
% avg_network_spike_rate = average spiking rate of all active cells
%
% avg_cell_spike_intensity = average spike peak intensity of each active cell
% (from smoothed dF/F signal)
%
% avg_network_spike_intensity = average spike peak intensity of all active cells

%% EXTRACT RELEVANT VARIABLES FROM STRUCT

time_axis = signal_info.time;
f_end = signal_info.analysis_window_frames;
true_pks = calcium_transients.true_peaks;

%% CALCULATE SPIKE RATE

% Spike rate calculations
active_cells = size(true_pks,1);

% Create matrices
cell_spike_freq = zeros(active_cells,1); % total number of spikes of each active cell
cell_spike_rate = zeros(active_cells,1); % spiking rate of each active cell

for i = 1:active_cells
    cell_spike_freq(i,1) = sum( ~isnan(true_pks(i,:)) );
    cell_spike_rate(i,1) = cell_spike_freq(i,1)/time_axis(f_end);
end 

avg_network_spike_rate = mean(cell_spike_rate);

if isnan(avg_network_spike_rate)
    avg_network_spike_rate = 0;
end

%% CACLULATE AVERAGE SPIKE INTENSITY FOR EACH CELL AND FOR THE ENTIRE NETWORK

avg_cell_spike_intensity = mean(true_pks,2, 'omitnan');
avg_network_spike_intensity = mean(avg_cell_spike_intensity,'omitnan');

if isnan(avg_network_spike_intensity)
    avg_network_spike_intensity = 0;
end

%% CALCULATE THE PERCENT CHANGE IN NETWORK OUTPUTS COMPARED TO VALUES AT 0-MINUS (IF APPLICABLE)

% If the function is being run for the initial timepoints data, 
% or if there is only one timepoint, 
% this calculation won't be run
if avg_network_spike_intensity_0minus == 0 && avg_network_spike_rate_0minus == 0
    percent_change_rate = NaN;
    percent_change_intensity = NaN;
else
    percent_change_rate = ((avg_network_spike_rate - avg_network_spike_rate_0minus) / avg_network_spike_rate_0minus) * 100;
    percent_change_intensity = ((avg_network_spike_intensity - avg_network_spike_intensity_0minus) / avg_network_spike_intensity_0minus) * 100;
end

%% CREATE STRUCTURE CONTAINING ALL OUTPUTS

calcium_spike_calculations.cell_spike_count = cell_spike_freq;
calcium_spike_calculations.cell_spike_rate = cell_spike_rate;
calcium_spike_calculations.network_avg_spike_rate = avg_network_spike_rate;
calcium_spike_calculations.cell_avg_spike_intensity = avg_cell_spike_intensity;
calcium_spike_calculations.network_avg_spike_intensity = avg_network_spike_intensity;
calcium_spike_calculations.percent_change_rate = percent_change_rate;
calcium_spike_calculations.percent_change_intensity = percent_change_intensity;

end