function [CWT, scalogram_fig, norm_scalogram_fig] =  wavelet(smooth_signal,signal_info, calcium_transients, timePointName, Max_0minus)

% FUNCTION TO CREATE WAVELET TRANSFORM OF EACH ACTIVE CELL IN NETWORK AND THE MEAN CWT FOR THE NETWORK
% 
% INPUTS:
% fdata = dF/F signals after it was run through "ProcessNeuroCa"
% time_axis = time axis created in "ProcessNeuroCa"
% timePointName = string containing the name of the timepoint ('0minus', '0plus', '1hr', 24hr')
% Max_0minus = maxiumum scalogram magnitude of 0minus
% frame_rate = image acquistion rate output by "ProcessNeuroCa"
% f_end = frame # at 120 sec
% active_cells = list of cell #s of active cells
% 
% OUTPUTS - STRUCTURE "CWT" CONTAINING THE VARIABLES:
% mean_network_scalogram
% tms
% frq
% all_scalograms
% coi
% M
%
% INDEPENDENT OUTPUTS:
% scalogram_fig
% norm_scalogram_fig

%% EXTRACT RELEVANT VARIABLES FROM STRUCT

time_axis = signal_info.time;
frame_rate = signal_info.fps;
f_end = signal_info.analysis_window_frames;
active_cells = calcium_transients.active_cells;

%%

if active_cells == 0
    warning(['No active cells in ',timePointName,'. CWT not run.'])
    CWT = [];
    scalogram_fig = figure;
    norm_scalogram_fig = figure;
    return;
end

% Taking into account that some datasets are less than 2 minutes
if size(smooth_signal,2) < f_end
    f_end = size(smooth_signal,2);
end
fb = cwtfilterbank('SignalLength',f_end,'SamplingFrequency',frame_rate);

% Create 3-D arrays that contains the CWT of all active cells in the network
% and the time power spectrum of each over 120 seconds (2760 frames)
[cfs,frq,coi] = cwt(smooth_signal(1,1:f_end),frame_rate); % cfs = continuous wavelet transform; frq = frequencies of the cwt (Hz); coi = cone of influence (Hz)

all_scalograms = zeros(size(cfs,1),size(cfs,2),size(active_cells,1));
% all_timePower = zeros(size(cfs,1),1,length(fdata(2:end,1)));
% p_f = zeros(size(cfs,1),size(cfs,2),length(fdata(2:end,1)));
% p_int = zeros(size(cfs,1),1,length(fdata(2:end,1)));

% Calculate continuous wavelet transform (CWT) for each cell and their power spectrum [rows = ?,
% columns = frames, 3rd dimension = cell]
for i = 1:size(active_cells,1)
    cell_num = active_cells(i,1);
    
    signal = smooth_signal(cell_num,1:f_end);

    all_scalograms(:,:,i) = wt(fb,signal);
%     all_timePower(:,1,cell_num) = timeSpectrum(fb,fdata(i+1,1:f_end));
%     p_f(:,:,cell_num) = abs(all_scalograms(:,:,i)).^2; % |W(v,t)|^2 function inside integral
%     p_int(:,:,cell_num) = trapz(p_f(:,:,i),2); % Integral of p_f with respect to the second dimension (time)
end

% Calculate average, population scalogram
mean_network_scalogram = mean(all_scalograms,3);

% Calculate Power Spectrum of Network Scalogram
% p_f = abs(mean_network_scalogram).^2; % |W(v,t)|^2 function inside integral
% p_int = trapz(p_f,2); % Integral of p_f with respect to the second dimension (time)
% mean_timePower = mean(all_timePower,3);

% Normalize scalogram magnitude for the plot
if strcmp(timePointName,'0minus')
    M = max(abs(mean_network_scalogram),[],'all');
else
    M = Max_0minus;
end


% Plot average, population scalogram and time power spectrum

tms = (0:numel(time_axis(1,1:f_end))-1)/frame_rate; % vector representing the sample times


scalogram_fig = figure;
PlotScalogram(tms, frq, mean_network_scalogram, coi, 1);
title(['Mean Network Scalogram of Active Cells in ',timePointName])

norm_scalogram_fig = figure;
PlotScalogram(tms, frq, mean_network_scalogram, coi, M);
clim([0 1])
title(['Normalized Magnitude Scalogram - ',timePointName])

% p_fig = figure;
% plot(frq,mean_timePower)
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% xlim([min(frq) max(frq)])
% xlabel('Frequency (Hz)')
% ylabel('Power Spectrum P')

%% CREATE STRUCTURE CONTAINING ALL OUTPUT VARIABLES (EXCEPT THE FIGURES)

CWT.all_scalograms = all_scalograms;
CWT.mean_network_scalogram = mean_network_scalogram;
CWT.time_vector = tms;
CWT.frequencies = frq;
CWT.coi = coi;
CWT.max_0minus = M;

end

%% Function to plot scalogram and normalize the colorbar for the 0-/0+ samples (FROM MATLAB SOURCE)
function hf = PlotScalogram(t,period,wt,coi,M)
% Input dataset time (t = tms), frequency (period = frq), continuous wave transform (wt = scalogram), cone of influence (coi), maximum magnitude of wt of 0minus sample or 24 hour sample (M) 
% Code adapted from parts of Matlab's built in cwt function

    
    % Normalize the magnitude of the wt
    norm_scalogram = wt/M;
    
    % Use magnitude limits in both the analytic and antianalytic parts
    % for the colormap
    cmin = min(abs(norm_scalogram(:)));
    cmax = max(abs(norm_scalogram(:)));
    if cmax <= cmin
        cmax = cmin+eps('single');
    end

    ylbl = [getString(message('Wavelet:getfrequnitstrs:Frequency')) ' (Hz) '];
    xlbl = [getString(message('Wavelet:wcoherence:Time'))  ' (s)'];
    titleString = getString(message('Wavelet:cwt:ScalogramTitle'));
    hf = gcf;
    clf;
    AX = axes('parent',hf);
    
    % The following axes hold must occur before any plotting or the default
    % interactivity is restored
    hold(AX,'on');
    hs = image('Parent',AX,...
        'XData',t,'YData',period,...
        'CData',abs(norm_scalogram), ...
        'CDataMapping','scaled');
    
    
    AX.YLim = [min(period),max(period)];
    AX.XLim = [min(t) max(t)];
    AX.Layer = 'top';
    AX.YDir = 'normal';
    AX.YScale = 'log';
    
    title(AX, titleString);
    ylabel(AX, ylbl)
    xlabel(AX, xlbl)
    
    hcol = colorbar('peer', AX);
    hcol.Label.String = getString(message('Wavelet:cwt:Magnitude'));
    
    
  
    plot(AX,t,coi,'w--','linewidth',2);
    baselevel = min([min(AX.YLim) min(coi)]);
    A1 = area(AX,t,coi,baselevel);
    A1.EdgeColor = 'none';
    A1.FaceColor = [0.8 0.8 0.8];
    alpha(A1,0.4);
    A1.PickableParts = 'none';
    
    hold(AX,'off');
    hf.NextPlot = 'replace';
end