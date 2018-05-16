% Statistical tests Spikes vs Clean resting data vs Artifact Resting data
%
% Questions we assess here:
%
%   1) Are there the EEG channels for which the spike data is statistically
%       different from the resting state data? (and same for the spike data
%       vs the artefact data )
%
%   2) Are there the GFT frequencies for which the spike data is
%       statistically different from the resting state data? (and same for
%       the spike data vs the artefact data)
%
% -------------------------------------------------------------------------
%
%  Version: V1 2018 05 15
%  Author: VYBORNOVA Anna 
% 
% ------------------------- BEGIN CODE ------------------------------------

clc 
clear all
close all

addpath("../GFTpackage"); 
R = 80; % adjacency matrix parameter (radius)

RestingClean.Path = "../../01_Data/20180424EEGresting/Clean/"; 
RestingClean.Filename = "_rest_clean_preprocessed_epochs.mat"; 
RestingClean.ID = "RestingClean"; 

Artefacts.Path = "../../01_Data/20180424EEGresting/Artifacts/"; 
Artefacts.Filename = "_rest_art_preprocessed_epochs.mat"; 
Artefacts.ID = "Artefacts"; 

Spikes.Path = "../../01_Data/20180424EEGspikes/"; 
Spikes.Filename = "pre-processed.mat"; 
Spikes.ID = "Spikes"; 

FileNameChan204 = '../../01_Data/EGI 204.Geneva Average 13.10-10.xyz';
chanlocs = readlocs(FileNameChan204);

subjects = ["03", "07", "11", "14", "16", "19"];

A = adj_generator([[chanlocs.X]', [chanlocs.Y]', [chanlocs.Z]'], R);

for iSub = 1:numel(subjects)
    
    data = {Spikes, RestingClean, Artefacts};
    
    for iD = 1:numel(data)
        % STEP 1: load the data
        % load spike data, clean resting data and artefact resting data
        load(strcat(data{iD}.Path, "sub-", subjects(iSub), data{iD}.Filename));
        if strcmp(data{iD}.ID, "Spikes")
            data{iD}.Signal = data_prep; 
            clear data_prep
        elseif strcmp(data{iD}.ID,"RestingClean")
            data{iD}.Signal = data_clean_epochs; 
            clear data_clean_epochs
        elseif strcmp(data{iD}.ID, "Artefacts")
            data{iD}.Signal = data_art_epochs; 
            clear data_art_epochs
        end
        % STEP 2: perform the GFT of the data
        for iEp = 1:size(data{iD}.Signal, 3)
            data{iD}.GFT(:,:,iEp) = gft(data{iD}.Signal(:, :, iEp), A); 
        end
        
        % STEP 3: perform a non-parametric t-test on the mean and the variance of
        % the data
        
        
        
        
        
        
        
    
    end
    
    figure, 
    plot(data{1}.Signal(:,:,5)'); 
    
    clear data



end




