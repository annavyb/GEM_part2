% 2017 08 10
% Processing the data that I got from HUG in Aug 2017 --> Namely Vallelian
% and Jungo


% The principle goal of this code is to verify whether the events that we
% recover in the .mat structure from the edf file correspond to the actual
% seizure

% In order to achieve this goal I first filter the data and then I use
% eegplot function to visualize the signals

clc
clear
close all

dbstop if error

set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaulttextinterpreter','none')

% parameters setting
% Pipeline parameters initialization
Params = setParam(1); % see the function setParam for further information


%% Loading data
Subjects = {'Jesus','Vallelian'};

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
        
    end
    % load
    load([PathToFile{1}, '/', FilesResting{1}, '.mat']);
    
    % Go from 257 to 204 channel configuration
    DataAll = ALLEEG.data;
    
    FileNameChan257 = '/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_Seizure_Epochs_V04_05/EGI 257.Geneva Average 13.10-10.xyz';
    FileNameChan204 = '/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_Seizure_Epochs_V04_05/EGI 204.Geneva Average 13.10-10.xyz';
    
    [IndSameChannels, Labels] = findSameChannels(FileNameChan257, FileNameChan204, 257, 204);
    DataAll = DataAll(IndSameChannels, :);
    
    % readlocs function puts the channels in the different order; make sure
    % that the order in DataAll is exactly the same as in chanlocs
    
    chanlocs = readlocs(FileNameChan204);
    indSame = 1;
    for iChL = 1:length(chanlocs)
        for iChO = 1:length(Labels)
            if strcmp(chanlocs(iChL).labels, Labels(iChO))
                IndReorderToChanlocs(iChL) = iChO;
            end
        end
    end
    
    DataAll = DataAll(IndReorderToChanlocs, :);
    
    % filter
    EEGSamplingRate = ALLEEG.srate;
    h1=fdesign.highpass('N,F3dB',Params.Filtering.FiltOrder,Params.Filtering.LowCf,...
        EEGSamplingRate);
    d1=design(h1,'Butter');
    h2=fdesign.lowpass('N,F3dB',Params.Filtering.FiltOrder,Params.Filtering.HighCf,...
        EEGSamplingRate);
    d2=design(h2,'Butter');
    
    DataEEGFilt = filterEEG(DataAll, d1, d2);
    
    % plot
    eegplot(DataEEGFilt, 'srate', EEGSamplingRate, ...
        'winlength', 60, 'dispchans', 50, 'eloc_file', chanlocs);
    
    ALLEEG.data = DataEEGFilt;
    ALLEEG.actions = ['FILTERED_'];
    save(['/Users/annavybornova/EPFL/Master_4/GEM/01Data/20170808DataHUGNewMarkers/PreprocessedData/', ...
        Subjects{iSub}, 'Resting', ALLEEG.actions], 'ALLEEG');
    keyboard;
    
    clear ALLEEG DataAll DataEEGFilt
    
    
    for iSeiz = 1:length(FilesSeizure)
        % load
        load([PathToFile{1}, '/', FilesSeizure{iSeiz}, '.mat']);
        
        % Go from 257 to 204 channel configuration
        DataAll = ALLEEG.data;
        
        FileNameChan257 = '/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_Seizure_Epochs_V04_05/EGI 257.Geneva Average 13.10-10.xyz';
        FileNameChan204 = '/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_Seizure_Epochs_V04_05/EGI 204.Geneva Average 13.10-10.xyz';
        
        [IndSameChannels, Labels] = findSameChannels(FileNameChan257, FileNameChan204, 257, 204);
        DataAll = DataAll(IndSameChannels, :);
        
        % readlocs function puts the channels in the different order; make sure
        % that the order in DataAll is exactly the same as in chanlocs
        
        chanlocs = readlocs(FileNameChan204);
        indSame = 1;
        for iChL = 1:length(chanlocs)
            for iChO = 1:length(Labels)
                if strcmp(chanlocs(iChL).labels, Labels(iChO))
                    IndReorderToChanlocs(iChL) = iChO;
                end
            end
        end
        
        DataAll = DataAll(IndReorderToChanlocs, :);
        
        
        
        % filter
        
        DataEEGFilt = filterEEG(DataAll, d1, d2);
        
        
        % plot
        
        eegplot(DataEEGFilt, 'srate', EEGSamplingRate, ...
            'winlength', 60, 'dispchans', 50, 'eloc_file', chanlocs);
        
        ALLEEG.data = DataEEGFilt;
        ALLEEG.actions = ['FILTERED_'];
        save(['/Users/annavybornova/EPFL/Master_4/GEM/01Data/20170808DataHUGNewMarkers/PreprocessedData/', ...
            Subjects{iSub}, 'Seizure', num2str(iSeiz), ALLEEG.actions], 'ALLEEG');
        
        keyboard;
        
        clear ALLEEG DataAll DataEEGFilt
      
        
    end
    
    
end









