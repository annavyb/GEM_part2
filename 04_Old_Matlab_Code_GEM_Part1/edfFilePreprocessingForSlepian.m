function [ OutALLEEG ] = edfFilePreprocessingForSlepian(ALLEEG, chanlocs, Params)
%EDFFILEPREPROCESSINGFORSLEPIAN performs the preprocessing steps to
%generate a .mat file containing ALLEEG data structure that will be further
%used for the Slepian analysis.

% We proceed in the following steps:
% STEP 1 : reduce the number of electrodes to 204 instead of 257 electrodes
% and put them into the order that corresponds to the chanlocs file
% STEP 2 : Filter the data between 1Hz and 30 Hz
% STEP 3 : Include the seizure onset markers


%   Input

%   [ALLEEG] - [] - []
%   [chanlocs] - [] - []
%   [Params] - [] - []

%   Output
%   [OutALLEEG] - [] - []

% Author: Anna Vybornova
% 2017 08 14
% Copyright EPFL MIP::Lab 2017
% _______

% STEP 1 : reduce the number of electrodes to 204 instead of 257 electrodes
% and put them into the order that corresponds to the chanlocs file

DataAll = ALLEEG.data;

FileNameChan257 = '/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_Seizure_Epochs_V04_05/EGI 257.Geneva Average 13.10-10.xyz';
FileNameChan204 = '/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_Seizure_Epochs_V04_05/EGI 204.Geneva Average 13.10-10.xyz';

[IndSameChannels, Labels] = findSameChannels(FileNameChan257, FileNameChan204, 257, 204);
DataAll = DataAll(IndSameChannels, :);

% readlocs function puts the channels in the different order; make sure
% that the order in DataAll is exactly the same as in chanlocs

indSame = 1;
for iChL = 1:length(chanlocs)
    for iChO = 1:length(Labels)
        if strcmp(chanlocs(iChL).labels, Labels(iChO))
            IndReorderToChanlocs(iChL) = iChO;
        end
    end
end

DataAll = DataAll(IndReorderToChanlocs, :);

% STEP 2 : Filter the data between 1Hz and 30 Hz

EEGSamplingRate = ALLEEG.srate;
h1=fdesign.highpass('N,F3dB',Params.Filtering.FiltOrder,Params.Filtering.LowCf,...
    EEGSamplingRate);
d1=design(h1,'Butter');
h2=fdesign.lowpass('N,F3dB',Params.Filtering.FiltOrder,Params.Filtering.HighCf,...
    EEGSamplingRate);
d2=design(h2,'Butter');

DataEEGFilt = filterEEG(DataAll, d1, d2);

ALLEEG.data = DataEEGFilt;

% eegplot(DataEEGFilt, 'srate', EEGSamplingRate, ...
%             'winlength', 60, 'dispchans', 50, 'eloc_file', chanlocs);
%

% STEP 3 : Include the seizure onset markers

if isfield(Params, 'SeizureOnsetFileName')
    disp('Including seizure onset markers ...')
    [num] = xlsread(Params.SeizureOnsetFileName);
    ALLEEG.SeizureOnset = num(1);
    ALLEEG.SeizureDuration = num(2);
end

OutALLEEG = ALLEEG;

end

