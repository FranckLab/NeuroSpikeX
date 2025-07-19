function [rate_constants, RateConstantHist] = decayRateConstant(signal_info, calcium_transients, smooth_signal, network_rate_cst_0minus, signal_strength_threshold)

%% Function to fit an exponential decay to the signal window and calculate the rate constant
% Using the decay windows detected in "findCalciumSpikes"
%
% INPUTS:
% pk_locs and min_locs and active_cells = outputs from "findCalciumSpikes.m" which give
% the location of the peaks in every cell, decay minimums, and list of active cell #s
%
% t_axis and fdata = outputs from "ProcessNeuroCa.m" (time axis and calcium signals)
%
% true_pks = 
%
% smooth_signal = 
%
% OUTPUTS:
% exp_fit = the cell #, exponential fit, and goodness of fit for each
% transient in each active cell
%
% all_rate_cst = rate constant for each calcium transient in each cell
% calculated from the exponential fit
%
% cell_avg_rate_cst = the average rate constant of every cell
%
% network_avg_rate_cst = the average rate constant of the entire active network
%
% excluded_rate_cst = matrix of cell # (column 1) and transient # (column 2)
% of the transients excluded from the average rate constant calculation
%
% RateConstantHist = histogram (distribution) figure of the average rate constants of every active cell

%% EXTRACT RELEVANT VARIABLES FROM STRUCT

time_axis = signal_info.time;
fdata = signal_info.fdata;
pk_locs = calcium_transients.peak_locs;
min_locs = calcium_transients.decay_min_locs;
active_cells = calcium_transients.active_cells;

%%
% fit the same type of exponential function used in NeuroCa. Get the decay constant for every transient in every cell

% Create cell arrays that contain every calcium transient from every active cell
% Indices for the 3rd dimension:
%   1 - cell #
%   2 - exponential fit
%   3 - goodness of fit (gof)

all_fits = cell(size(pk_locs,1),1);
all_gof = cell(size(pk_locs,1),1);
% rate_cst_not_included = zeros(size(pk_locs));
rate_cst = zeros(size(pk_locs));
time = time_axis(1,:);

%% Attempt to use parallels
if active_cells == 0
    warning('No active cells found in current timepoint. Therefore, no decay rate constants can be calculated.')
    rate_constants.exp_fit = [];
    rate_constants.gof = [];
    rate_constants.all_rate_cst = [];
    rate_constants.cell_avg_rate_cst = [];
    rate_constants.network_avg_rate_cst = 0;
    rate_constants.percent_change_decay_cst = 0;
    RateConstantHist = figure;
    return;
end

parfor c = 1:size(active_cells,1)

    cell_num = active_cells(c,1);
    k = 1;

    temp_r = zeros(1,size(pk_locs(c,:),2));
    cell_fit = cell(1,size(pk_locs(c,:),2));
    cell_gof = cell(1,size(pk_locs(c,:),2));
    cell_pk_locs = pk_locs(c,:);
    cell_min_locs = min_locs(c,:);
    cell_fdata = fdata(cell_num,:);  % decay signal (starting at peak and shifted to end at y=0)
   

    % for loop runs through each transient
    for j = 1 : size(cell_pk_locs,2)%sum(~isnan(pk_locs(c,:)),2)

        % Since not every cell has the same number of pks, skip any cells
        % containting NaN
        if isnan(cell_pk_locs(1,j))
            temp_r(j) = NaN;
            k = k+1;
            continue
        end

        peak_frame = cell_pk_locs(1,j);
        end_frame = cell_min_locs(1,j);

        % % Exclude any transients where the peak and minimum values are 
        % % less than or equal to the set noise floor (default: 2).
        % if isempty(signal_strength_threshold)
        %     thresh = 2;
        % else
        %     thresh = signal_strength_threshold;
        % end
        fit_info = [];
        gof = [];
        
        if isnan(end_frame) % the decay end was not found for this transient
            temp_r(j) = NaN;
            k = k+1;
        elseif j <= (size(pk_locs(c,:),2) - 1) % all peaks except the last one
            % find the column where the next peak is recorded in pk_locs
            find_next_pk = find(~isnan(cell_pk_locs(1,j+1:end))==1, 1, 'first');
            if ~isempty(find_next_pk)
                next_pk = j+find_next_pk;
            else
                next_pk = j+1;
            end

            % Sometimes the code identifies the same minimum for two
            % consecutive peaks that are located close together, or the minimum
            % occurs after the following peak. In these cases, exclude the 
            % first transient of the two from the rate decay average
            if cell_min_locs(1,j) == cell_min_locs(1,next_pk) || cell_min_locs(1,j) >= cell_pk_locs(1,next_pk) % if the decay end location of two consecutive spikes are the same, or if the decay end location appears after the next peak
                % rate_cst_not_included_temp(j) = j;
                % rate_cst_not_included(row_index,1) = cell_num;
                % rate_cst_not_included(row_index,1+k) = j;
                % row_index = row_index + 1;
                temp_r(j) = NaN;
                k = k+1;
            else
                % generate data window of the decay, starting at the peak 
                x = time(1,1:(end_frame-peak_frame+1)); % shifting all exponentials to begin at time 0 
                y = cell_fdata(1,peak_frame:end_frame) - cell_fdata(1,end_frame);  % decay signal (starting at peak and shifted to end at y=0)
   
                % exponential fit method 
                [fit_info, gof] = fit(x', y', 'exp1');
                coeffs = coeffvalues(fit_info);
                temp_r(j) = -coeffs(2); % The decay rate constants (fitted) for every signal in each cell
                % Note: larger rate constant means a faster decay
                k = k+1;
            end
        
        else
            % generate data window of the decay, starting at the peak 
            x = time(1,1:(end_frame-peak_frame+1)); % shifting all exponentials to begin at time 0 
            y = cell_fdata(1,peak_frame:end_frame) - cell_fdata(1,end_frame);  % decay signal (starting at peak and shifted to end at y=0)
   
            % exponential fit method 
            [fit_info, gof] = fit(x', y', 'exp1');
            coeffs = coeffvalues(fit_info);
            temp_r(j) = -coeffs(2); % The decay rate constants (fitted) for every signal in each cell
            % Note: larger rate constant means a faster decay
            k = k+1;
        end
        cell_fit{1,j} = fit_info;
        cell_gof{1,j} = gof;
    
    end
    rate_cst(c,:) = temp_r;
    % rate_cst_not_included(c,:) = rate_cst_not_included_temp;
    [all_fits{c,1}] = cell_fit;
    [all_gof{c,1}] = cell_gof;

end  

%% CALCULATE AVERAGE RATE CONSTANT FOR EACH CELL AND FOR THE ENTIRE NETWORK

% Exclude negative rate constants from calculation (because that would
% indicate an exponential rise, which indicates there was an issue with
% that calculation
neg = rate_cst<0;
rate_cst(neg) = 0;

% To account for each cell having a different number of peaks
rate_cst(rate_cst==0) = NaN;
% rate_cst_not_included(rate_cst_not_included==0) = NaN;

% Average the decay rate constants in every cell
avg_cell_rate_cst = mean(rate_cst, 2, 'omitnan');

% Average decay rate of the entire network
avg_network_rate_cst = mean(avg_cell_rate_cst,'omitnan');

if isnan(avg_network_rate_cst)
    avg_network_rate_cst = 0;
end

%% Create histogram of all decay rate constants

if sum(~isnan(avg_cell_rate_cst)) > 1
    RateConstantHist = figure;

    hist_rate_cst = histfit(avg_cell_rate_cst);
    hist_rate_cst(1).FaceColor = 'r';
    hist_rate_cst(1).FaceAlpha = .5;
    hist_rate_cst(1).EdgeColor = 'none';
    hist_rate_cst(2).Color = 'r';
    hist_rate_cst(2).LineWidth = 3;
    xlabel('Average Decay Rate Constant (s^{-1})')
    ylabel('# of Cells')
else
    RateConstantHist = [];
end

%% CALCULATE THE PERCENT CHANGE IN NETWORK OUTPUTS COMPARED TO VALUES AT 0-MINUS (IF APPLICABLE)

% If the function is being run for the initial timepoints data, 
% or if there is only one timepoint, 
% this calculation won't be run
if network_rate_cst_0minus == 0
    percent_change_decay_cst = NaN;
else
    percent_change_decay_cst = ((avg_network_rate_cst - network_rate_cst_0minus) / network_rate_cst_0minus) * 100;
end

%% STRUCTURE CONTAINING ALL OUTPUTS

rate_constants.exp_fit = all_fits;
rate_constants.gof = all_gof;
rate_constants.all_rate_cst = rate_cst;
rate_constants.cell_avg_rate_cst = avg_cell_rate_cst;
rate_constants.network_avg_rate_cst = avg_network_rate_cst;
rate_constants.percent_change_decay_cst = percent_change_decay_cst;

end