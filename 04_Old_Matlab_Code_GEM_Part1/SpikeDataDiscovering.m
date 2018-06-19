% 2017 03 28
% Discovering spike data from Ana

clc;
clear all;
close all;

Patients = {'CB', 'CC', 'CJ', 'JG', 'JH', 'KR', 'MR', 'NA', 'NS', 'PJ', 'RS',...
    'TZ'};
EpSite = {'L','L', 'R','L','L','L', 'R', 'R', 'R', 'L','L','L'}; 
Folders = {'Avg_B', 'Avg_T'};

for iPatients = 1:length(Patients)
    PathToData = ['/Volumes/TOSHIBA EXT/02_EEG_EpilepsyData_Ana/SpikeData/', ...
        Patients{iPatients}, '/EEG/'];
    EpochInfo(iPatients).Name = Patients{iPatients}; 
    for iFolders = 1:length(Folders)
        switch iFolders 
            case 1
                
                PathAbove =  [PathToData, Folders{iFolders}];
                filesAbove = dir( fullfile(PathAbove,'*To*.ep') );
                for iFabove = 1: length(filesAbove)
                   Epoch = load([PathAbove, '/', filesAbove(iFabove).name]);
                   EpochInfo(iPatients).BA(iFabove).Length = size(Epoch,1); 
                   clear Epoch
                end
                clear filesAbove iFabove PathAbove
                
                Path = [PathToData, Folders{iFolders},'/Epochs']; 
                files = dir( fullfile(Path,'*To*.ep'));
                for iFiles = 1: length(files)
                    Epoch = load([Path, '/', files(iFiles).name]); 
                    EpochInfo(iPatients).B(iFiles).Length = size(Epoch,1);
                   
                    clear Epoch
                end
                clear files
                
            case 2
                
                PathAbove =  [PathToData, Folders{iFolders},EpSite{iPatients}];
                filesAbove = dir( fullfile(PathAbove,'*To*.ep') );
                for iFabove = 1: length(filesAbove)
                   Epoch = load([PathAbove, '/', filesAbove(iFabove).name]);
                   EpochInfo(iPatients).TA(iFabove).Length = size(Epoch,1); 
                   clear Epoch
                end
                clear filesAbove iFabove PathAbove
                
                Path = [PathToData, Folders{iFolders},EpSite{iPatients}, ...
                    '/Epochs'];
                files = dir( fullfile(Path,'*To*.ep') );
                for iFiles = 1: length(files)
                    Epoch = load([Path, '/', files(iFiles).name]);
                    EpochInfo(iPatients).T(iFiles).Length = size(Epoch,1);
                   
                    clear Epoch
                end
        end 
        
    end 
end