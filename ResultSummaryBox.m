function [] = ResultSummaryBox(calcium_spike_calculations, calcium_transients, rate_constants, network_bursts, timepoint)

meanRate = calcium_spike_calculations.network_avg_spike_rate;
meanIntensity = calcium_spike_calculations.network_avg_spike_intensity;
meanRateConstant = rate_constants.network_avg_rate_cst;
changeRate = calcium_spike_calculations.percent_change_rate;
burstRate = network_bursts.burst_rate;
changeIntensity = calcium_spike_calculations.percent_change_intensity;
changeConstant = rate_constants.percent_change_decay_cst;
activeCells = size(calcium_transients.active_cells,1);
totalCells = activeCells + size(calcium_transients.inactive_cells,1);

summary = msgbox(sprintf('Mean Spike Rate = %2.3g \nSpike Rate Percent Change = %2.3g \n\nNetwork Burst Rate = %2.3g \n\nMean Spike Intensity = %2.2f \nSpike Intensity Percent Change = %2.3g \n\nMean Decay Constant = %2.3g \nDecay Constant Percent Change = %2.3g \n\n# of Active Cells = %2.0f \n# of Total Cells = %2.0f',meanRate,changeRate,burstRate,meanIntensity,changeIntensity,meanRateConstant,changeConstant,activeCells,totalCells),[timepoint,': RESULT SUMMARY'], 'non-modal');

set(summary, 'position', [300 300 350 280]); %makes box bigger
ah = get( summary, 'CurrentAxes' );
ch = get( ah, 'Children' );
set( ch, 'FontSize', 14 ); %makes text bigger

end