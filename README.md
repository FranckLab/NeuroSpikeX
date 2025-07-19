# NeuroSpikeX
Comprehensive Detection and Characterization of Neural Calcium Transients
> A MATLAB toolbox for detecting and quantifying calcium transient spikes

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)
  - [Preparing Data](#preparing-data)
  - [Script Mode](#script-mode)
  - [Graphical User Interfaces](#graphical-user-interfaces)
    - [PreProcessSettingsApp](#preprocesssettingsapp)
    - [PostProcessingApp](#postprocessingapp)
    - [RelativeBaselineChangeApp](#relativebaselinechangeapp)
- [Citation](#citation)
- [License](#license)
- [Acknowledgements](#acknowledgements)
- [Contact](#contact)

## Overview
NeuroSpikeX is a user-friendly MATLAB tool for quantitative analysis of neural calcium dynamics. It provides robust calcium spike detection, comprehensive network metrics, and intuitive graphical interfaces, significantly reducing false positives and simplifying parameter tuning. NeuroSpikeX seamlessly integrates into existing workflows by utilizing outputs from the established algorithm NeuroCa, thereby enhancing accuracy and reproducibility. Validated with hundreds of in vitro datasets, it effectively analyzes calcium dynamics across multiple experimental time points before and after mechanical or chemical stimuli. NeuroSpikeX facilitates detailed cell- and network-level analyses in large datasets, making rigorous calcium transient characterization accessible to researchers with minimal coding expertise.

## Installation
Download MATLAB 2024b or later (mathworks.com/downloads). When installing, ensure the following toolboxes are included: Signal Processing, Image Processing, Statistics and Machine Learning, Curve Fitting, and Parallel Computing. 

## Usage
### Preparing Data
It is recommended that microscope recording timelapses are first processed in ImageJ to adjust brightness and contrast to give better spatial resolution. Timelapses must then be saved as 'Image Sequence' to be compatible with NeuroCa.

This program is a transient detection and calculation software for calcium analysis. An outside software must first be used to detect cell somas, extract fluorescent traces, and normalize the traces. NeuroSpikeX was designed for use with NeuroCa (doi: 10.1117/1.NPh.2.3.035003). If a different software is used, ensure you have 3 MATLAB matrices saved in time point subfolders within a root sample folder (e.g. '0minus', '0plus', '1hr', '24hr' subfolders in '7.9.24.sample1'):
- 'fdata': First row is the time axis (frames/frame_rate). All subsequent rows are the normalized fluorescent intensity traces for each cell in the network.
- 'center': nx2 matrix containing x and y coordinates of all detected cell somas
- 'radii': nx1 matrix containing radii of all detected cell somas
- **Note:** 'center' and 'radii' are only needed if planning to run the Relative Baseline Change App

---

### Script Mode
1. **Run the Pre-Processing App (GUI) first** so that 'preprocessing_settings.mat' exists. The app only needs to be run once, and then the script can be run for many samples. At any time, if a new dataset is to be analyzed, the Pre-Processing App can be run again to change the settings.
2. **Run the main analysis script**
   ```matlab
   run Main_Calcium_Code.m
   ```
3. Follow the prompts:
   - **Enter the number of experimental time points for the sample** e.g. 4
   - **Select the root sample folder:** the folder containing the time point subfolders (e.g. '7.9.24.sample1')

The script then automatically:
- Loads NeuroCa outputs for all experimental time points ('fdata', 'center', 'radii')
- Applies saved noise-floor and signal smoothing settings
- Smooths calcium fluorescent signals, detects peaks, and computes calcium transient metrics
- Exports cell- and network-level results into the root sample folder:
    - 'all_quantitative_analyses_time1.mat' (calcium spike rate and intensity calculations, raster, continuous wavelet transform, network burst, decay rate constants)
    - 'analyzed_fdata_time1.mat' (smoothed signals, time axis, detected peak amplitudes and locations, list of active and inactive cells, derivative of all signals, locations of decay ends)
    - 'ActiveCellsPerFrame_time1.fig' (plot of number of active cells per frame, used to detect network bursts)
    - 'Burst_raster_time1.fig' (raster plot with network bursts highlighted)
    - 'DecayRateConstantHistogram.fig' (histogram plotting the distribution of decay rate constants from all time points)
    - 'MeanNetworkScalogram_time1.fig' (scalogram of the average continuous wavelet transform of all active cells in a network)
    - 'Normalized_MeanNetworkScalogram_time1.fig' (scalogram with the magnitudes normalized to the maximum magnitude of the baseline time point network)
    - 'SpikeRateHistogram.fig' (histogram plotting the distribution of average cell spike rates from all time points)
    - All data and calculations are also saved into the time point subfolders, which is used by the Post-Processing App and Relative Baseline Intensity Change App
    - All figures are also saved as .png files

While most customizations are set in the Pre-Processing App, some optional changes can be made in the following location:
- countNetworkBursts.m function: Network bursts are calculated by finding a period of time where a certain percentage of all cells are firing nonstop. This percentage can be changed on line 6 of the function (default: 20%).
- Main_Calcium_Code.m script: Histograms can be customized with the fit probability distribution functions and histogram bin widths. By default, spike rates are fit with a 'lognormal' distribution and have a bin width of 0.006 (line 400), and decay rate constants are fit with a 'normal' distribution and have a bin width of 0.02 (line 443).

---

### Graphical User Interfaces

#### PreProcessSettingsApp

To launch:
```matlab
PreProcessSettingsApp
```
> **Note:** Run this once **before** 'Main_Calcium_Code.m'

**Tab 1: Loading, Smoothing, and Noise Floor**
1. Click **Select Timepoint Folder(s)** and choose each subfolder for one sample.
2. Define timelapse analysis duration (e.g. '120' seconds).
3. Check or uncheck **Use Fixed Noise Floor** box: If checked, the noise floor you set in the app will be applied to all samples until changed. If unchecked, you will set the noise floor for each sample every time Main_Calcium_Code.m is run (suggested for calcium indicators with noisier signals)
4. Click **Calculate & Set Noise Floor**: About 10 random signals will pop up from each sample. Drag a box around a part of the signal where there are **no** transients. This tells the code the average signal range of the noise. This is then averaged from all time points and rounded up to give the noise threshold. This number can also be manually overridden.
5. Set the smoothing factor (n) signals will be smoothed by, and click **Preview Smoothing** to see what it looks like on your signals. If the value is too low, the noise won't be smoothed enough and may be detected as false calcium transient peaks. If the value is too high, signals will be oversmoothed and true calcium peaks may be flattened.
6. Click **Save All Settings** -> generates 'preprocessing_settings.mat' in the NeuroSpikeX folder

**Tab 2: Spike Detection Parameters**
1. Optionally, check **Use Max/Min Signal Bounds** and set the bounds: some noisier calcium indicators (e.g. fluo-4 AM) may have extreme signal artifacts that need to be excluded
2. Adjust the **Minimum Peak Prominence** used to find calcium peaks. Learn more about peak prominence at https://www.mathworks.com/help/signal/ug/prominence.html
3. Click **Preview Peak Detection** to view the detection of peaks on any smoothed signal in the network.
4. Click **Save All Settings** to update the saved parameters

---

#### PostProcessingApp

To launch:
```matlab
PostProcessingApp
```
> **Note:** Run this **after** 'Main_Calcium_Code.m'

**Tab 1: Visualize Cell Analysis**
Visualize the calcium transient peak detection in any cell and the exponential decay fits to the transients.

**Tab 2: Manual Accuracy Check**
Calculate the accuracy of NeuroSpikeX on your analyses. Set the percentage of all cells to check.

Two accuracy calculatin modes:
- **Mode 1: Cell Active/Inactive** When a signal shows on the screen, select whether it contains transients (signal) or no transients (noise). Accuracy metrics will appear when complete and be saved in the root sample folder as 'PostProcessApp_Cell AccuracyMetrics.mat'.
- **Mode 2: Detailed Transient Count** When a signal appears, manually count how many transients (peaks) you see. Type 0 if no peaks. Error metrics and a histogram of the error distribution will appear when complete and be saved in the root sample folder as 'PostProcessApp_TransientCountMetrics.mat'.

---

#### RelativeBaselineChangeApp

To launch:
```matlab
RelativeBaselineChangeApp
```
> **Note:** Run this **after** 'Main_Calcium_Code.m'

**Tab 1: Frame Selection**
1. Click **Select Timepoint Folders & Load Signals** and choose subfolders of two time points to compare the baseline fluorescent intensity change.
2. Select a frame in each graph where the signals are almost all at baseline. This ensures we won't be comparing the fluorescence during a calcium transient. I suggest choosing a frame close to the end of time point 1 and close to the beginning of time point 2.
3. Instructions will appear on which frames of the timelapse to be saved in the time point subfolders. I suggest using ImageJ to save the frame as .tif and adjust the brightness/contrast to improve spatial resolution if needed. **Note:** If brightness/contrast is adjusted, it must be adjusted the same for both images.

**Tab 2: Analysis Settings (optional)**
1. Check the box to **Dynamically recalculate cell radii at baseline fluorescence** to recalculate all radii instead of using the NeuroCa provided metric. This may improve accuracy as NeuroCa often detects cells when they are at their brightest during a transient. The radii may differ slightly at baseline fluorescence
2. Click **Find Optimal Background Settings** to run optimization to find the best method to calculate the fluorescence of every cell.
3. If satisfied with highlighted settings, click **Apply & Save Optimal Settings**. Or, manually set the background subtraction method and inner and outer factors (related to annulus radii surrounding each cell), then click **Apply & Save Manual Settings**.

**Tab 3: Analysis Pipeline**
Calculate the baseline fluorescence change of every cell as log2(Intensity2/Intensity1) (I2/I1).

Two analysis modes:
- **Mode 1: Control Sample Mode** This must be run first if control samples exist. This mode analyzes all control samples together to establish a threshold for what constitutes a substantial baseline change. The middle 95th percentile of log2(I2/I1) values is characterized as no change. This definition will be applied to all the following experimental samples.
  1. Click **Select Control Folders** and select each control sample folder. Click cancel on the pop-up window once all samples have been selected
  2. Type the correct time point subfolder names (e.g., '0minus' and '0plus')
  3. Click **RUN ANALYSIS**. For each sample, the user must manually select six matching cells in the two images to align them if any shifts have occurred.
  4. When the analysis is complete, navigate to the folder where you want to save the threshold ('Control_BaselineChange_Thresholds.mat'). The distribution of log2(I2/I1) for all control samples will be saved in the same location ('GLOBAL_control_baselineChange_histogram.fig').
  5. In the bottom text window, a summary will show the average log2(I2/I1) and the number of cells characterized as increased/didn't change/decreased baseline fluorescence for all samples.
- **Mode 2: Experimental Sample Mode** Calculate baseline change for an experimental sample.
  1. Click **Select Experimental Folder** to load the root sample folder containing all time point subfolders (e.g., 7.9.24.sample1)
  2. Click **Load Control Thresholds** to load the 'Control_BaselineChange_Thresholds.mat' file calculated in Mode 1.
  3. Type the correct time point subfolder names (e.g., '0minus' and '0plus').
  4. Click **RUN ANALYSIS**. The user will manually select six matching cells in each image to determine any spatial shift.
  5. When the analysis is complete, the summary of the results will appear at the bottom of the text box.
 
A new subfolder ('BaselineChange') will be created in the root sample folder, containing all results:
- 'IntensityChange.mat' (fold change (I2/I1) and log2(fold change) for all cells and the averages)
- 'MatchedCells.mat' (cell numbers, coordinates, and radii matched between the two images, transformation matrix applied to image 2, both images, and a list of the valid matches)
- 'diagnostic_cell_matching.fig' (sample image showing matched cells)
- 'diagnostic_recalculated_radii.fig (sample image showing the original and recalculated radii (if checked in tab 2)
- 'excluded_cells.fig' (sample image showing cells excluded from calculations due to their calculated intensity not being valid or negative)
- 'baselineChange_histogram.fig' (histogram of log2(I2/I1) for the sample with the control thresholds labeled)
- 'baselineChange_spatialDistribution.fig' (sample image labeling each cell as increased, no change, or decreased fluorescence)

---

## Citation
If you use NeuroSpikeX, please cite:

## License
This project is released under the MIT License.

## Acknowledgements
Thank you to Noah Meltzer for signal-processing guidance.
Supported by ONR PANTHER award N00014-22-1-2828

## Contact
Jamie Sergay: jamiea416@gmail.com
Lab Webpage: https://francklab.me.wisc.edu/
