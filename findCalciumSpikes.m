function [Calcium_Transients] = findCalciumSpikes(signal_info, smooth_signal, signal_strength_threshold, useBounds, minBound, maxBound, prm)

%% Function to identify the peak of true calcium spikes and find the end of the signal decay
%
% Create a decay windown for each spike
%
% INPUTS:
% frame_rate = image acquisition frame rate (differs between samples)
%
% time_axis = time vector (seconds) (created in "ProcessNeuroCa.m")
%
% smooth_signal = smoothed fdata signals (created using "smoothCalciumSignals.m")
%
% signal_strength_threshold = defining the absolute magnitude of a generalized noise floor (determined empirically)
%
% OUTPUTS - STRUCTURE "CALCIUM_TRANSIENTS" WITH THE FOLLOWING VALUES:
% true_peaks = value of calcium transient peaks for each active cell
%
% peak_locs = location of calcium transient peaks for each active cell
%
% smooth_diff = smoothed derivative of active cells only
%
% decay_min_locs = location of the end of the transient decay for each active cell
% 
% active_cells = list of all active cells
% 
% inactive_cells = list of all inactive cells (noise)

% Define all field names of the output structure
fieldNames = {'true_peaks','peak_locs','active_cells','inactive_cells','smooth_diff','decay_min_locs'};
    
% Initialize all fields to 0 (scalar)
for i = 1:numel(fieldNames)
    Calcium_Transients.(fieldNames{i}) = []; 
end
%% EXTRACT RELEVANT VARIABLES FROM STRUCT

frame_rate = signal_info.fps;
time_axis = signal_info.time;

%% USE CERTAIN CRITERIA TO IDENTIFY THE CELL AS ACTIVE (HAVING TRUE CALCIUM SPIKES) OR INACTIVE (NOISE)

active_count = 0;
inactive_count = 0;
initial_active_cells = [];

% Initialize variables that might not be created if no active cells are found
true_pks = [];
pk_locs = [];
min_locs = []; % For decay_min_locs
smooth_diff = []; % For smooth_diff
active_cells = []; % For active_cells

% Handle case where smooth_signal might be empty or have no rows
if isempty(smooth_signal) || size(smooth_signal,1) == 0
    % If no signal data at all, all cells are inactive
    Calcium_Transients.inactive_cells = 1:size(smooth_signal,1); % Or just 0 if empty
    return;
end

max_signal_all = max(smooth_signal,[],2);
min_signal_all = min(smooth_signal,[],2);

% Run filter for each cell (row in smooth_signal)
for cell = 1:size(smooth_signal,1)

    % Using the global maxiumum and minimum of the signal, identify artifacts 
    % introduced usually in NeuroCa that would indicate the signal as noise
    signal_strength = max_signal_all(cell,1) - min_signal_all(cell,1);

    % The signal_strength_threshold may be manually altered depending on the user's calcium indicator and camera/microscope noise
    % The extreme boundaries of -50 and 50 can also be adjusted by the user
    if signal_strength > signal_strength_threshold && ...
       (~useBounds || (max_signal_all(cell,1) < maxBound && min_signal_all(cell,1) > minBound))

        active_count = active_count + 1;
        initial_active_cells(active_count,1) = cell; % list of active cells        
    else
        inactive_count = inactive_count + 1;
        inactive_cells(inactive_count,1) = cell; % list of inactive cells (noise)
    end
end

% Store the identified inactive cells in the output structure
Calcium_Transients.inactive_cells = inactive_cells;

% In the case where no calcium transients are identified in the entire
% network, exit the function
if isempty(initial_active_cells)
    return
end

%% FIND PEAKS (CALCIUM TRANSIENTS) FROM SMOOTHED SIGNAL OF EACH ACTIVE CELL

i_active = 1; % index used in for loop
inactive_cells = [];
active_cells = [];

% Identify the maximum signal intensity of all cell signals, determine 
% the outlier range to not include spike peak values that lie in that range
% find_outliers = isoutlier(max_signal_all,'gesd'); % use generalized ESD method
find_outliers = isoutlier(max_signal_all,'percentile',[0 95]); % use percentile threshold
max_signal_outliers = max_signal_all(find_outliers);
peak_cutoff = min(max_signal_outliers);

% Run through each possible active cell, based on initial list from above section
for a = 1:size(initial_active_cells,1)

    cell = initial_active_cells(a,1); % the cell # (a.k.a. the row in smooth_signal)
   
    all_pks = []; all_locs = []; % reset vector for each cell
    [all_pks,all_locs] = findpeaks(smooth_signal(cell,:),'MinPeakProminence',prm);
   
    if isempty(all_pks)
        % Add a new cell to the inactive list created in above section
        inactive_count = inactive_count + 1;
        inactive_cells(inactive_count,1) = cell; 
        continue
    end

    % Run through each peak identified -> check if there are any peaks
    % left after excluding above the cutoff
    check = length(all_pks);
    for i = 1:length(all_pks)
        if all_pks(i) >= peak_cutoff % to be excluded
            check = check - 1;
        end
    end
                
    if check == 0
        % Add a new cell to the inactive list created in above section
        inactive_count = inactive_count + 1;
        inactive_cells(inactive_count,1) = cell;
    else
        active_cells(i_active,1) = cell;
        k = 1; % index used in for-loop below

        % Run through each peak identified
        for i = 1:length(all_pks)
            if all_pks(i) < peak_cutoff
                true_pks(i_active,k) = all_pks(i);
                pk_locs(i_active,k) = all_locs(i);
                k = k + 1;
            else
                % if the peak is greater than the cutoff value, set it to zero 
                % to be deleted later
                true_pks(i_active,k) = 0;
                pk_locs(i_active,k) = 0;
                k = k + 1;
            end
        end
        i_active = i_active+1;
    end
end

% To account for each cell having a different number of peaks
% true_pk => amplitude of peak (from filtered signal)
% true_locs => location of peak
true_pks(true_pks==0) = NaN;
pk_locs(pk_locs==0) = NaN;

%% (Using Savitsky-Galoy Filer) FIND THE BOTTOM OF EACH CALCIUM TRANSIENT DECAY IN THE ACTIVE CELLS
% Use derivitative of smoothed signal
% Detect when the signal crosses 0 (from negative to positive)
    
frame_totals = size(smooth_signal,2)-1;
diff_smooth_signal = zeros(size(active_cells,1), frame_totals);
smooth_diff = zeros(size(active_cells,1), frame_totals);


for a = 1:size(active_cells,1)

    cell = active_cells(a,1); % the cell # (a.k.a. the row in smooth_signal matrix)
    % if cell == 105
    %     dummy = 1;
    % end
    
    % Find the derivative of the smoothed calcium dF/F signal
    diff_smooth_signal(a,:) = diff(smooth_signal(cell,:))./diff(time_axis(1,1:frame_totals+1)); % signal / time
    
    % Smooth the derivative signal using Savitsky-Galoy Filter, which
    % preserves the minimums and maximums extra well
    polynomial_order = 3; % Typically 2-4
    window_size = 21; % Must be odd
    smooth_diff(a,:) = sgolayfilt(diff_smooth_signal(a,:), polynomial_order, window_size);

    % Run through each transient peak in the cell being analyzed
    for peak = 1:size(true_pks,2)%sum(~isnan(true_pks(a,:)),2)

        % if cell == 408 && peak == 12
        %     dummy = 1;
        % end

        % some rows may look like [1 2 NaN 3] or [1 NaN NaN Nan], etc. - skip the NaN
        if isnan(true_pks(a,peak))
            min_locs(a,peak) = NaN;
            continue
        end
        % Starting at the frame 0.5 seconds from the peak frame, run
        % through every frame after until the last frame
        search_start = round(0.5*frame_rate); % 0.5 seconds

        % If the search starting frame doesn't have a negative
        % "smooth_diff" value, find the next frame that does have a
        % negative value to be the start
        % start_frame = pk_locs(a,peak) + search_start;
        % if smooth_diff(a, start_frame) > 0


        if (pk_locs(a,peak) + search_start) >= size(smooth_diff,2)
            % If the location of the peak is within 119.5 and 120 seconds
            min_locs(a,peak) = size(smooth_signal,2);
        else
            for frame_num = (pk_locs(a,peak) + search_start) : (size(smooth_diff,2)-3)

                if (smooth_diff(a,frame_num-2)<0 && smooth_diff(a,frame_num-1)<0 && smooth_diff(a,frame_num)<0) && (smooth_diff(a,frame_num+1)>=0 && smooth_diff(a,frame_num+2)>=0 && smooth_diff(a,frame_num+3)>=0)
                    min_locs(a,peak) = frame_num + 1; % location of decay end
                    break;

                % if the decay doesn't reach the asymptote b/c it's too close to the end of the signal,
                % min_loc for the peak will default to the end of the entire signal (120 s)
                elseif frame_num == (size(smooth_diff,2)-3)
                    min_locs(a,peak) = size(smooth_signal,2); 
                end
            end
        end
    end
end

% To account for each cell having a different number of peaks
% min_locs = location of decay minimums after calcium transient peak
min_locs(min_locs==0) = NaN;

%% CREATE STRUCTURE CONTAINING ALL OUTPUTS

Calcium_Transients.true_peaks = true_pks;
Calcium_Transients.peak_locs = pk_locs;
Calcium_Transients.active_cells = active_cells;
Calcium_Transients.inactive_cells = inactive_cells;
Calcium_Transients.smooth_diff = smooth_diff;
Calcium_Transients.decay_min_locs = min_locs;

end