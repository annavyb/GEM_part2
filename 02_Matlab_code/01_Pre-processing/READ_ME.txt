The new data that Margherita and Serge gave to us contain

1) Resting state data along with the artifacts (fs = 1000)
2) 1-s spikes epochs (fs = 250)

—— 

This data was retrieved from the EPI format. However I use eeglab functions to visualize the topographies of the GFT functions and the EEG maps. 

The order of the electrodes differ in the EPI format (xyz file) and when this .xyz file is read by the reallocs function of the EEG lab

Therefore I need to find the mapping from the .xyz to the .locs format, to make sure that the EEG signal corresponds to the correct location in space 

The function convertEGItoChanlocs.m provides the mapping from EGI .xyz to .locs 

-

I keep all the 257 electrodes, as we need them for the further analysis. The problem is the central electrode in the Slepian analysis for most of the patient will disappear if we only use the 204 configuration
_

In the pre-processing phase of the analysis the script preprocessing_spikes.m performs the preprocessing for the spikes, and the preprocessing_resting.m performs the preprocessing for the resting state data ( both clean and with artifacts )

1. filtering
2. bad channels handling 

-

2 main functions in this folder are: 

1. preprocessing_spikes.m 
2. preprocessing_resting.m 







