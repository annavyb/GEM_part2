% 2017 08 10
% Processing the data that I got from HUG in Aug 2017 --> Namely Vallelian
% and Jungo

% I take the inspiration from the script mainSeizureDataProcessing for
% geenrating the matlab structure

% The principle goal of this code is to transform the .edf file into the
% .mat structure

clc
clear
close all

dbstop if error

set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaulttextinterpreter','none')

%% Adding paths and necessary parameters
cd('/Users/annavybornova/EPFL/Master_4/GEM/02MATLAB_Code');

addpath('/Users/annavybornova/EPFL/Master_4/eeglab14_0_0b')
addpath('/Users/annavybornova/EPFL/Master_4/GEM/02MATLAB_Code');

PathToData = '/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_Seizure_Data';
PathToGftCode = '/Users/annavybornova/EPFL/3eme_annee/Bachelor project/Matlab_code';

addpath(PathToData);
addpath(PathToGftCode);

OptionsFiltering.LowCf=1;
OptionsFiltering.HighCf=30;
OptionsFiltering.FiltOrder=8;

%%
% adding chanlocs
chanlocs = ...
    readlocs('/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_Seizure_Data/Electrode_loc/EGI204.GenevaAverage13.10-10.xyz');

%% Loading data
Subjects = {'Jungo'};

% for-loop to process each patient

for iSub = 1:length(Subjects)
    
    disp(['Now processing --> ', Subjects{iSub}]);
    
    % define the subject-specific file location
    
    if strcmp(Subjects{iSub}, 'Jesus')
        
        PathToFile = {'/Users/annavybornova/EPFL/Master_4/GEM/01Data/20170808DataHUGNewMarkers/Jesus'};
        FilesResting = {'JESUS_DA_CONCEI_0626_1056_resting_state.Export'};
        FilesSeizure = {'JESUS_DA_CONCEI_0626_1056_seizure_clip.Export'};
        
    elseif strcmp(Subjects{iSub}, 'Vallelian')
        PathToFile = {'/Users/annavybornova/EPFL/Master_4/GEM/01Data/20170808DataHUGNewMarkers/Vallelian'};
        FilesResting = {'VALLELIAN BENOI 0117 1540_resting state.Export'};
        FilesSeizure = {'VALLELIAN BENOI 0117 1540_seizure1.Export', ...
            'VALLELIAN BENOI 0117 1540_seizure2.Export', ...
            'VALLELIAN BENOI 0117 1540_seizure3.Export', ...
            'VALLELIAN BENOI 0117 1540_seizure4.Export'};
        
    elseif strcmp(Subjects{iSub}, 'Oguey')
        PathToFile = {'/Users/annavybornova/EPFL/Master_4/GEM/01Data/20170808DataHUGNewMarkers/Oguey'};
        FilesResting = {'OGUEY 0805 2113_resting state.Export'};
        FilesSeizure = {'OGUEY 0805 2113_seizure1.Export'};
        
    elseif strcmp(Subjects{iSub}, 'Jungo')
        PathToFile = {'/Users/annavybornova/EPFL/Master_4/GEM/01Data/20170808DataHUGNewMarkers/Jungo'};
        FilesResting = {'JUNGO NATHALIE 0617 2019_resting state_fil'};
        FilesSeizure = {'JUNGO NATHALIE new 0617 2019_seg_1.Export.i3DS2.To204.i3DS2.To204.Export'};
        
    end
    
    % opening the edf files in matlab
    
    ALLEEG = pop_biosig([PathToFile{1}, '/', FilesResting{1}, '.edf']);
    ALLEEG.chanlocs = chanlocs;
    
    save([PathToFile{1}, '/', FilesResting{1}, '.mat'], 'ALLEEG')
    clear ALLEEG
    
    for iSeiz = 1:length(FilesSeizure)
        ALLEEG = pop_biosig([PathToFile{1}, '/', FilesSeizure{iSeiz}, '.edf']);
        ALLEEG.chanlocs = chanlocs;
        
        save([PathToFile{1}, '/', FilesSeizure{iSeiz}, '.mat'], 'ALLEEG')
        clear ALLEEG
        
    end
    
    
end






