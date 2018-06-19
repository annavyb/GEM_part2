% 2017 04 04 
% In this script I do the pre-processing steps as in the BTOP paper and I
% do the eeplot of the data

% STEP1: Average Common Referencing
% STEP2: Filtering & plot the data scroll

clc
clear
close all

dbstop if error
%% Adding paths and necessary parameters
cd('/Users/annavybornova/EPFL/Master_4/GEM/02MATLAB_Code');

addpath('/Users/annavybornova/EPFL/Master_4/eeglab14_0_0b')
addpath('/Users/annavybornova/EPFL/Master_4/GEM/02MATLAB_Code');

PathToData = '/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_Seizure_Data';
PathToGftCode = '/Users/annavybornova/EPFL/3eme_annee/Bachelor project/Matlab_code';

addpath(PathToData);
addpath(PathToGftCode);

xlsRange = 1;
OptionsFiltering.LowCf=1;
OptionsFiltering.HighCf=30;
OptionsFiltering.FiltOrder=8;

ConsiderOnlyDifferentMarker = true;

%%
% adding chanlocs
chanlocs = ...
    readlocs('/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_Seizure_Data/Electrode_loc/EGI204.GenevaAverage13.10-10.xyz');

% Generate the adjacency matrix based on the chanlocs
R = 30;

disp('... Adjacency Matrix generation');
A = adj_generator([[chanlocs.X]', [chanlocs.Y]', [chanlocs.Z]'], R);
nNeighbors = sum(A);
avgNeighbors = mean(nNeighbors);

disp(['R = ', num2str(R)]);
disp(['avg number of neighbors = ', num2str(avgNeighbors)]);

%% Loading data

FileName = {'01EEGBarukcic','02EEGJesus','03EEGJungo', '04EEGOguey', '05EEGRizvi', ...
    '06EEGVallelian', '07EEGAlessi', '08EEGBersier', '09EEGVarela'};
cd('../01Data/HUG_Seizure_Data/');
SelectedFiles = [1:length(FileName)];
ind = 1;
ind1 = 1;

for iFile = SelectedFiles
    
    disp(['...Now processing file: ', FileName{iFile}]);
    
    dirName = [PathToData,'/', FileName{iFile}];
    files = dir( fullfile(dirName,'*.edf') );
    MarkerFile = dir( fullfile(dirName,'*.xlsx') );
    cd(FileName{iFile});
    
    ALLEEG = pop_biosig(files.name);
    
    
    
    ALLEEG.chanlocs = chanlocs; % Objectives --> some of the subjects
    % have less than 204 channels -->
    % identify them: Rizvi
    Markers = xlsread(MarkerFile.name);
    ind1 = ind1 + 1;
    
    if size(ALLEEG.data, 1) ~= 204
        warning('The number of channels is different from 204')
        DifferentNumberOfElectrodes(ind).Name = FileName{iFile};
        DifferentNumberOfElectrodes(ind).Number = size(ALLEEG.data, 1);
        ind = ind + 1;
    end
    
  
    
    
    
    %% STEP 1: AVERAGE COMMON REFERENCING
    
    ALLEEGreref = ALLEEG;
    ALLEEGreref.data = reref(ALLEEG.data);
    
    %% STEP 2: Filtering
    
    ALLEEGfilt = ALLEEGreref;
    
    OptionsFiltering.Srate = ALLEEGreref.srate;
    OptionsFiltering.PlotFreqz = false;
    OptionsFiltering.PlotSpectrum = false;
    OptionsFiltering.PlotSpectrumChannel = 4;
    
    DataF = filterEEG(ALLEEGreref.data, OptionsFiltering);
    ALLEEGfilt.data = DataF;
    
    % plot the filtered channels 
    eegplot(DataF, 'srate', ALLEEGfilt.srate, 'title', FileName{iFile},...
        'winlength', 10, 'dispchans', 102);
    cd ..
end 