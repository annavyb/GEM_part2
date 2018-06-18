The new data that Margherita and Serge gave to us contain

1) Resting state data along with the artifacts (fs = 1000)
2) 1-s spikes epochs (fs = 250)

—— 

This data was retrieved from the EGI format. However I use eeglab functions to visualize the topographies of the GFT functions and the EEG maps. 

The order of the electrodes differ in the EGI format (xyz file) and when this .xyz file is read by the reallocs function of the EEG lab

Therefore I need to find the mapping from the .xyz to the .locs format, to make sure that the EEG signal corresponds to the correct location in space 

The function convertEGItoChanlocs.m provides the mapping from EGI .xyz to .locs 

—

In our analysis we only care for the electrodes that are positioned at the head and not at the cheeks, therefore we should also go from 257 electrode configuration to the 204 electrode configuration 

_

In the pre-processing phase of the analysis the script preprocessing_spikes.m performs the preprocessing for the spikes, and the preprocessing_resting.m performs the preprocessing for the resting state data ( both clean and with artifacts )



