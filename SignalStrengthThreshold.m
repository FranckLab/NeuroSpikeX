function [signal_strength_threshold] = SignalStrengthThreshold(smooth_signal,fdata)

% Function to empirically determine the average noise floor of a samples.
% The noise floor is a function of cells per given day and imaging
% parameters.
% Threshold should remain the same for the timepoints of one sample but can
% be different between samples.
% This code looks at the noise floor on the smoothed signal run through
% "smoothCalciumSignals.m" and not the raw signal
%
% INPUTS:
%
%
% OUTPUTS:
%

%% Have 10 samples pop onto the screen, one by one, and have the user click on the top and bottom points of the noise floor for each

% Choose random cells each time
num_cells = size(smooth_signal, 1);
noise_floor = zeros(10,1); 
valid_count = 0;
attempts = 0;

while valid_count < 10 && attempts < 40  % prevent infinite loops
    cell_num = randi(num_cells);
    attempts = attempts + 1;

    y_axis1 = smooth_signal(cell_num,:);
    y_axis2 = fdata(cell_num,1:size(smooth_signal,2));

    if max(y_axis1) < 50 && min(y_axis1) > -50 && sum(y_axis1) ~= 0
        figure;
        x_axis = 1:size(smooth_signal,2);

        raw = plot(x_axis, y_axis2, '-b');
        hold on
        plot(x_axis, y_axis1, 'LineWidth', 6);
        xlim([0 inf])
        xlabel('Frame #')
        ylabel('dF/F')
        title('Select noise floor',['Cell #', num2str(cell_num)])
        legend('Raw Data','Smoothed Signal')
        hold off

        [~,~,ydata] = selectdata('selectionmode','rect','Ignore',raw);

        noise_floor(valid_count+1) = max(ydata) - min(ydata);
        valid_count = valid_count + 1;

        close all
        clear data
    end
end


if all(isnan(noise_floor))
    signal_strength_threshold = 5; % fallback default
else
    signal_strength_threshold = mean(noise_floor, 'omitnan');
end

end
