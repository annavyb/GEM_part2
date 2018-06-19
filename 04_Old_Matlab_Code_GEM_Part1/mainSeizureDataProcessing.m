% 2017 03 29

% General signal-processing pipeline for the Seizure Data (data from Pieter)
% Includes the following steps

% STEP1: Average Common Referencing
% STEP2: Filtering
% STEP3: Cut according to the markers (seizure vs non-seizure)
% STEP4: Cut into 1s epochs with 0.5s overlap
% STEP5: compute the GFT of the data
% STEP6: train a linear classifier to discriminate between
% seizure/non-seizure data
% STEP7: evaluate the performance in terms of
%       1) missclassification level (is GFT better than simple EEG?)
%       2) in terms of Betas (which channels are the most discriminative?)

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
SelectedFiles = [3, 6:length(FileName)];
ind = 1;
ind1 = 1;

for iFile = 9
    
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
    OptionsFiltering.PlotFreqz = true;
    OptionsFiltering.PlotSpectrum = true;
    OptionsFiltering.PlotSpectrumChannel = 4;
    
    DataF = filterEEG(ALLEEGreref.data, OptionsFiltering);
    ALLEEGfilt.data = DataF;
    
    % plot the filtered channels 
    eegplot(DataF, 'srate', ALLEEGfilt.srate, 'title', FileName{iFile},...
        'winlength', 10, 'dispchans', 102); 
    keyboard
    
    %% STEP 3: Cut according to the markers (seizure vs non-seizure)
    range = reshape(Markers', 1, 2*size(Markers, 1));
    range = [1, range, ALLEEG.pnts];
    count1 = 1;
    count2 = 1;
    for iRange = 1:(length(range)-1)
        
        s1 = range(iRange) ;
        s2 = range(iRange + 1);
        
        
        if mod(iRange,2) == 0 % the number is 2,4, etc
            
            Ictal(count1).ALLEEG = ALLEEGfilt;
            
            if s1~=s2
                Ictal(count1).ALLEEG.data = ALLEEGfilt.data(:, s1:s2);
                Ictal(count1).ALLEEG.times = ALLEEGfilt.times(s1:s2);
                Ictal(count1).ALLEEG.xmin = s1/ALLEEGfilt.srate;
                Ictal(count1).ALLEEG.xmax = s2/ALLEEGfilt.srate;
                Ictal(count1).ALLEEG.pnts = size(Ictal(count1).ALLEEG.data,2);
            else
%                 Ictal(count1).ALLEEG.data = ALLEEGfilt.data(:, s1:range(end));
%                 Ictal(count1).ALLEEG.times = ALLEEGfilt.times(s1:range(end));
%                 Ictal(count1).ALLEEG.xmin = s1/ALLEEGfilt.srate;
%                 Ictal(count1).ALLEEG.xmax = range(end)/ALLEEGfilt.srate;
%                 Ictal(count1).ALLEEG.pnts = size(Ictal(count1).ALLEEG.data,2);
            end
            
            count1 = count1+1;
            
        else
            InterIctal(count2).ALLEEG = ALLEEGfilt;
            
            InterIctal(count2).ALLEEG.data = ALLEEGfilt.data(:, s1:s2);
            InterIctal(count2).ALLEEG.times = ALLEEGfilt.times(s1:s2);
            InterIctal(count2).ALLEEG.xmin = s1/ALLEEGfilt.srate;
            InterIctal(count2).ALLEEG.xmax = s2/ALLEEGfilt.srate;
            InterIctal(count2).ALLEEG.pnts = size(InterIctal(count2).ALLEEG.data, 2);
            
            
            count2 = count2+1;
        end
        
        
    end
    for i = 1: length(Ictal)
        figure, 
        eegplot(Ictal(i).ALLEEG.data,'srate',Ictal(i).ALLEEG.srate ); 
    end
    
      for j = 1: length(InterIctal)
        figure, 
        eegplot(InterIctal(i).ALLEEG.data,'srate',InterIctal(i).ALLEEG.srate ); 
    end
    
    %% STEP5: Cut into 1s epochs with 0.5s overlap
    WINDOW_SIZE = ALLEEG.srate;
    WINDOW_OVERLAP = round(WINDOW_SIZE/2);
    
    OptionsWin.WinSize = WINDOW_SIZE;
    OptionsWin.WinOverlap = WINDOW_OVERLAP;
    
    for iIctal = 1:length(Ictal)
        cd Ictal
        WinIctal(iIctal).Windows = windowEEG(Ictal(iIctal).ALLEEG, OptionsWin);
        for iWindows = 1:length(WinIctal(iIctal).Windows)
            if ~isempty(WinIctal(iIctal).Windows(iWindows)) 
            k = WinIctal(iIctal).Windows(iWindows).ALLEEG; 
            save(['Ictal',...
                num2str(WinIctal(iIctal).Windows(iWindows).ALLEEG.times(1)), '_',...
                num2str(WinIctal(iIctal).Windows(iWindows).ALLEEG.times(end)), '.mat'],...
                'k');
            clear k
            end
        end
        cd ..
    end
    
    for iInterIctal = 1:length(InterIctal)
        WinInterIctal(iInterIctal).Windows = windowEEG(InterIctal(iInterIctal).ALLEEG, OptionsWin);
        cd InterIctal
        for iWindows = 1:length(WinInterIctal(iInterIctal).Windows)
           if ~isempty(WinInterIctal(iInterIctal).Windows(iWindows))
            k.ALLEEG = WinInterIctal(iInterIctal).Windows(iWindows).ALLEEG; 
            save(['InterIctal',...
                num2str(WinInterIctal(iInterIctal).Windows(iWindows).ALLEEG.times(1)), '_',...
                num2str(WinInterIctal(iInterIctal).Windows(iWindows).ALLEEG.times(end)), '.mat'],...
                'k');
            clear k
           end
          
        end
        cd ..
    end
    cd ..
    clear ALLEEG ALLEEGreref ALLEEGfilt Markers DataF Ictal InterIctal M range
end

% come back to the initial folder
cd('/Users/annavybornova/EPFL/Master_4/GEM/02MATLAB_Code');






