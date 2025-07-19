function [network_bursts, activationFig, burstFig] = countNetworkBursts(raster, time, timePointName)

% Parameters

% USER CAN CHANGE:
percentage = .2;

f_end = size(raster,2); % number of frames in two minutes
num_cells = size(raster,1);
num_active_cells = sum(any(raster,2)); % sum total number of active cells
min_burst_duration = 5; % minimum number of frames in a burst duration
threshold = percentage*num_active_cells; % minimum number of unique active cells in a burst (based on all signaling cells)
% threshold = percentage*num_cells; % minimum number of unique active cells in a burst (based on all detected cells)

if num_active_cells == 0
    network_bursts.burst_rate = 0;
    activationFig = figure;
    burstFig = figure;
    return;
end

% Identify possible burst windows (where at least one cell is active)
active_timepoints = any(raster,1);

% Determine start and end of active period (the frames with 0 active cells
% immediately before and after the active period)
burst_starts = find(diff([0, active_timepoints]) == 1);
burst_ends = find(diff([active_timepoints, 0]) == -1);

% Ensure minimum burst duration threshold set above is met
burst_durations = burst_ends - burst_starts + 1;
valid_bursts = burst_durations >= min_burst_duration;
burst_starts = burst_starts(valid_bursts);
burst_ends = burst_ends(valid_bursts);

% Check each potential burst period to see if it meets the active cell
% threshold
final_burst_starts = [];
final_burst_ends = [];
numBurstCells = [];

for i = 1:length(burst_starts)
    % isolate raster data within the burst window being checked
    burst_window = raster(:,burst_starts(i):burst_ends(i));

    % count number of unique active cells within the burst window
    numActiveCells = sum(any(burst_window,2));

    % validate threshold criteria
    if numActiveCells >= threshold
        final_burst_starts(end+1) = burst_starts(i);
        final_burst_ends(end+1) = burst_ends(i);
        numBurstCells(end+1) = numActiveCells;
    end
end

% % Display results
% disp([timePointName,' Detected Network Bursts:']);
% for i = 1:length(final_burst_starts)
%     fprintf('Burst %d: Start = %d, End = %d, Duration = %d\n',i,final_burst_starts(i),final_burst_ends(i), final_burst_ends(i)-final_burst_starts(i)+1);
% end

%% Calculate burst rate
total_bursts = length(final_burst_starts);
burst_rate = total_bursts/time(1,f_end); % burst/second
%% Plot results
activationFig = figure;
plot(sum(raster,1),'k')
hold on
for i = 1:length(final_burst_starts)
    xline(final_burst_starts(i),'g','Start','LineWidth',2)
    xline(final_burst_ends(i),'m','End','LineWidth',2)
end
xlabel('Frame #')
ylabel('# of Active Cells')
title(['Network Burst Detection: ',timePointName],'(Highlighted Network Bursts - 20% of Active Cells)')
hold off


burstFig = figure; % raster plot

for i = 1:num_cells
    plot(time(1,1:f_end), i*raster(i,1:f_end),'k.','MarkerSize',4)
    hold on
end
for i = 1:length(final_burst_starts)
    xline(time(1,final_burst_starts(i)),'g','Start','LineWidth',2)
    xline(time(1,final_burst_ends(i)),'m','End','LineWidth',2)
end

xlabel('Time (s)')
ylabel('Cell #')
xlim([0 120])
ylim([0.5 size(raster,1)])
title(['Raster - ',timePointName],'(Highlighted Network Bursts - 20% of Active Cells)')
hold off

%% Create structure with network burst data

network_bursts.num_cells_per_burst = numBurstCells; % total number of cells participating in each burst window
network_bursts.burst_starts = final_burst_starts;
network_bursts.burst_ends = final_burst_ends;
network_bursts.total_bursts = total_bursts;
network_bursts.burst_rate = burst_rate;

end