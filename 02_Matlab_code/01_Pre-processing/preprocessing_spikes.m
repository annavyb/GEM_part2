
% Preprocessing of the spike epochs data. 
% The following steps are performed on the data 
% 1) Going from 257 electrodes to 204 electrodes 
% 2) Converting EGI channel order into chanlocs (eeglab) channels order
% 3) Finding the bad channels in the chanlocs (eeglab) configuration
% 4) Interpolating bad channels using the "spherical" option of the eeglab 
% built-in function pop-interp
% 5) Saving the pre-processed data

% -------------------------------------------------------------------------
%
%  Version: V1 2018 05 11
%  Author: VYBORNOVA Anna 
% 
% ------------------------- BEGIN CODE ------------------------------------

clc 
clear all
close all

% STEP 1: loading the data

% SPIKE DATA VARIABLES  
% - channels_to_keep - those are the good channels 
% - data - the actual data in the format: time x electrodes x epoch
% - fs - sampling frequency

addpath("../01_Pre-processing")

Path = "/Users/annavybornova/EPFL/Master_4/GEM/GEM_part2/01_Data/20180424EEGspikes/"; 
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
    filename = strcat(Path, "sub-", subjects(iSub), ".mat"); 
    load(filename); 
    
    % Bad_channels structure contains the following fields 
    % indEGI - index of the bad channels in the EGI configuration
    % label - label of the channel (string)
    % indChanloc - in the chanlocs configuration 
    % labelsChanlocs - labels in the chanloc configuration
    
    Bad_channels.indEGI = setdiff(1:size(data, 2), channels_to_keep); 
    Bad_channels.label = chan_egi_257(Bad_channels.indEGI); 
    
    % going from 257 to 204 channels 
    
    ch_to_delete =  find(~ismember(chan_egi_257, chan_egi_204));
    data(:, ch_to_delete, :) = [];
    chan_egi_257(ch_to_delete) = []; 
    
    % STEP 2: convert egi into chanlocs 
    
    ind_reorder = convertEGItoChanlocs({chanlocs_204.labels},chan_egi_257); 
    data = data(:,ind_reorder, :);

    
    % STEP 3: find the indeces of the bad chanels in the chanloc configuration 
    Bad_channels.indChanlocs = zeros(1,3); 
    Bad_channels.labelsChanlocs = {}
    ind = 1; 
    for iChan = 1:length(Bad_channels.indEGI)
        if ~isempty(find(ismember({chanlocs_204.labels},Bad_channels.label{iChan})))
            Bad_channels.indChanlocs(ind)  = find(ismember({chanlocs_204.labels},Bad_channels.label{iChan}));
            Bad_channels.labelsChanlocs{ind} = Bad_channels.label{iChan}; 
            ind = ind + 1; 
        end
    end
    
    Bad_channels.indChanlocs(Bad_channels.indChanlocs == 0) = []; 
    
    % STEP 4: Bad channels handling - interpolation
    % load a datastructure 
    data_prep = []
    for iEpoch = 1:size(data, 3)
        load("/Users/annavybornova/EPFL/Master_4/GEM/GEM_part2/01_Data/DataStructure.mat")
        ALLEEG.data = data(:, :, iEpoch)';
        ALLEEG.chanlocs = chanlocs_204;
        EEGOUT = pop_interp(ALLEEG, Bad_channels.indChanlocs, 'spherical');
        data_prep(:, :, iEpoch) = EEGOUT.data;
    end
    
    clear Bad_channels
    
    % STEP 5: save the data 
    name_save = strcat(Path, "sub-", subjects(iSub), "pre-processed.mat"); 
    save(name_save, "data_prep")
end 
 