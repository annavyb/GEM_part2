% Preprocessing of the resting state EEG data 
% 1) Going from 257 electrodes to 204 electrodes 
% 2) Converting EGI channel order into chanlocs (eeglab) channels order
% 3) Filtering the data between 1 Hz and 35 Hz (according to Margherita)
% 4) Bad channels handling should be done once I have the information about
% the bad channels 
% 5) Dividing data into epochs 
% 6) Data saving 
% -------------------------------------------------------------------------
%
%  Version: V1 2018 05 11
%  Author: VYBORNOVA Anna 
% 
% ------------------------- BEGIN CODE ------------------------------------

clc 
clear all
close all

% data filtering params 
ORDER = 8; 
LOW_CF = 1; 
HIGH_CF = 35;
FS = 1000; 
% STEP 1: Loading the data

Path = "/Users/annavybornova/EPFL/Master_4/GEM/GEM_part2/01_Data/20180424EEGresting/Clean/"; 
PathArtefacts = "/Users/annavybornova/EPFL/Master_4/GEM/GEM_part2/01_Data/20180424EEGresting/Artifacts/"

subjects = ["03", "07", "11", "14", "16", "19"]; 

FileNameChan257 = '/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_Seizure_Epochs_V04_05/EGI 257.Geneva Average 13.10-10.xyz';
FileNameChan204 = '/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_Seizure_Epochs_V04_05/EGI 204.Geneva Average 13.10-10.xyz';

[~ , chan_egi_257] = readFileChanloc(FileNameChan257, 257); 
[~ , chan_egi_204] = readFileChanloc(FileNameChan204, 204); 

chanlocs_257 = readlocs(FileNameChan257);
chanlocs_204 = readlocs(FileNameChan204);

for iSub = 1:length(subjects)
    chanlocs_257 = readlocs(FileNameChan257);
    chanlocs_204 = readlocs(FileNameChan204);
    [~ , chan_egi_257] = readFileChanloc(FileNameChan257, 257); 
    [~ , chan_egi_204] = readFileChanloc(FileNameChan204, 204); 
    strcat(" .. Now processing subject ", subjects(iSub))
    filename = strcat(Path, "sub-", subjects(iSub), "_rest_clean.mat"); 
    filename_artefacts = strcat(PathArtefacts, "sub-", subjects(iSub), "_rest_art.mat"); 
    load(filename);
    load(filename_artefacts); 
    
    % going from 257 to 204 channels 
    
    ch_to_delete =  find(~ismember(chan_egi_257, chan_egi_204));
    data_clean(ch_to_delete, :) = [];
    data_art(ch_to_delete, :) = []; 
    chan_egi_257(ch_to_delete) = []; 
    
    % STEP 2: convert egi into chanlocs 
    
    ind_reorder = convertEGItoChanlocs({chanlocs_204.labels},chan_egi_257); 
    data_clean = data_clean(ind_reorder, :);
    data_art = data_art(ind_reorder, :);
    
    % STEP 3: filtering the data 
    
    h1=fdesign.highpass('N,F3dB',ORDER, LOW_CF,...
    FS);
    d1=design(h1,'Butter');
    h2=fdesign.lowpass('N,F3dB',ORDER, HIGH_CF,...
    FS);
    d2=design(h2,'Butter');
    
    data_clean = filterEEG(data_clean, d1, d2);
    data_art = filterEEG(data_art, d1, d2); 
    
    % STEP 4: bad channels handling
    % I'm still missing the infrmation about the bad channel 
    % Complete the bad channels handling once I have the information about
    % the bad channels 
    
    % STEP 5: dividing the data into the 1s epochs
    fs = 1000; 
    data_clean_epochs = epoching(data_clean, fs); 
    data_art_epochs = epoching(data_art, fs);
    
    % STEP 6: save data
    save(strcat(Path, "sub-", subjects(iSub), "_rest_clean_preprocessed.mat"), "data_clean")
    save(strcat(Path, "sub-", subjects(iSub), "_rest_clean_preprocessed_epochs.mat"), "data_clean_epochs")
    save(strcat(PathArtefacts, "sub-", subjects(iSub), "_rest_art_preprocessed.mat"), "data_art")
    save(strcat(PathArtefacts, "sub-", subjects(iSub), "_rest_art_preprocessed_epochs.mat"), "data_art_epochs")
    
end
