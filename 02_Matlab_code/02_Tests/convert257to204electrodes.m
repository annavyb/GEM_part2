% Test the conversion of 257 to 204 electrodes 

clc 
clear all
close all

addpath("/Users/annavybornova/EPFL/Master_4/GEM/GEM_part2/02_Matlab_code/01_Pre-processing")

FileNameChan257 = '/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_Seizure_Epochs_V04_05/EGI 257.Geneva Average 13.10-10.xyz';
FileNameChan204 = '/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_Seizure_Epochs_V04_05/EGI 204.Geneva Average 13.10-10.xyz';

% Step 1: find the correspondance between 257 and 204 channels of the EGI file 
[IndSameChannels, Labels] = findSameChannels(FileNameChan257, FileNameChan204, 257, 204);

% Step 2: find the mapping between the EGI file and the chanlocs 

chanlocs = readlocs(FileNameChan204);
indSame = 1;
for iChL = 1:length(chanlocs)
    for iChO = 1:length(Labels)
        if strcmp(chanlocs(iChL).labels, Labels(iChO))
            IndReorderToChanlocs(iChL) = iChO;
        end
    end
end

% example read file EGI 

[position, channame] = readFileChanloc(FileNameChan257, 257); 


ind_reordered = convertEGItoChanlocs({chanlocs.labels}, Labels); 
isequal(ind_reordered, IndReorderToChanlocs) % sanity check done --> 
% the function convertEGItoChanlocs works

