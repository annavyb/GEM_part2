% 2017 08 14
% This script serves as a test for the function edfFilePreprocessingForSlepian.m

clc
clear all
close all

dbstop if error;

set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaulttextinterpreter','none')

% parameters setting
% Pipeline parameters initialization
Params = setParam(1); % see the function setParam for further information

SavePreprocessedDataFolder = '/Users/annavybornova/EPFL/Master_4/GEM/01Data/20170808DataHUGNewMarkers/PreprocessedData/';

%% Loading data
Subjects = {'Jesus','Vallelian', 'Oguey', 'Jungo'};

% load chanlocs file

FileNameChan257 = '/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_Seizure_Epochs_V04_05/EGI 257.Geneva Average 13.10-10.xyz';
FileNameChan204 = '/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_Seizure_Epochs_V04_05/EGI 204.Geneva Average 13.10-10.xyz';

chanlocs = readlocs(FileNameChan204);

% for-loop to process each patient

for iSub = 4 %1:length(Subjects)
    
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
        FilesSeizure = {'JungoSeizure1'};
        
    end
    
    %% Resting Data
    
    % load
    load([PathToFile{1}, '/', FilesResting{1}, '.mat']);
    
    ALLEEG = edfFilePreprocessingForSlepian(ALLEEG, chanlocs, Params);
    
    %     eegplot(ALLEEG.data, 'srate',ALLEEG.srate, ...
    %         'winlength', 60, 'dispchans', 50, 'eloc_file', chanlocs);
    %
    save([SavePreprocessedDataFolder, Subjects{iSub}, ...
        'RestingPreProcessedForSlepian.mat'], 'ALLEEG');
    %keyboard
    
    clear ALLEEG
    %% Seizure Data
    if ~strcmp(Subjects{iSub}, 'Jungo')
        for iSeiz = 1:length(FilesSeizure)
            load([PathToFile{1}, '/', FilesSeizure{iSeiz}, '.mat']);
            SeizureOnsetFileName = ['/Users/annavybornova/EPFL/Master_4/GEM/01Data/20170808DataHUGNewMarkers/', ...
                Subjects{iSub}, '/', Subjects{iSub},'Seizure', num2str(iSeiz), ...
                'Events.xlsx'];
            Params.SeizureOnsetFileName = SeizureOnsetFileName;
            ALLEEG = edfFilePreprocessingForSlepian(ALLEEG, chanlocs, Params);
            
            %         eegplot(ALLEEG.data, 'srate', ALLEEG.srate, ...
            %             'winlength', 60, 'dispchans', 50, 'eloc_file', chanlocs);
            
            save([SavePreprocessedDataFolder, Subjects{iSub}, ...
                'Seizure', num2str(iSeiz),'PreProcessedForSlepian.mat'], 'ALLEEG');
            %keyboard;
            
            clear ALLEEG
            Params = rmfield(Params, 'SeizureOnsetFileName');
        end
    else
        
        % STEP 2 : Filter the data between 1Hz and 30 Hz
        for iSeiz = 1:length(FilesSeizure)
            load([PathToFile{1}, '/', FilesSeizure{iSeiz}, '.mat']);
            SeizureOnsetFileName = ['/Users/annavybornova/EPFL/Master_4/GEM/01Data/20170808DataHUGNewMarkers/', ...
                Subjects{iSub}, '/', Subjects{iSub},'Seizure', num2str(iSeiz), ...
                'Events.xlsx'];
            Params.SeizureOnsetFileName = SeizureOnsetFileName;
            EEGSamplingRate = ALLEEG.srate;
            h1=fdesign.highpass('N,F3dB',Params.Filtering.FiltOrder,Params.Filtering.LowCf,...
                EEGSamplingRate);
            d1=design(h1,'Butter');
            h2=fdesign.lowpass('N,F3dB',Params.Filtering.FiltOrder,Params.Filtering.HighCf,...
                EEGSamplingRate);
            d2=design(h2,'Butter');
            
            DataEEGFilt = filterEEG(ALLEEG.data, d1, d2);
            
            ALLEEG.data = DataEEGFilt;
            
            ALLEEG.SeizureOnset = 600000000;
            ALLEEG.SeizureDuration = 1000;
            
            save([SavePreprocessedDataFolder, Subjects{iSub}, ...
                'Seizure', num2str(iSeiz),'PreProcessedForSlepian.mat'], 'ALLEEG');
            clear ALLEEG
            Params = rmfield(Params, 'SeizureOnsetFileName');
        end
        end
        
    end
