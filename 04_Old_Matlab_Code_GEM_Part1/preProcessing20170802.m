% 2017 08 02
% In this script I try to open the .mff files that I got on the 2017 08 26
% in HUG (from the NetStation)
addpath('/Users/annavybornova/EPFL/Master_4/GEM/02MATLAB_Code/mffimport');
Subjects = {'Jesus', 'Jungo', 'Oguey', 'Vallelian'};
for iFile = 1:length(Subjects)
    
    DataFolder = ['/Volumes/TOSHIBA EXT/Anna EEGs/', Subjects{iFile}, '/'];
    switch iFile
        case 1
            SeizureFileName = 'JESUS_DA_CONCEI_0626_1056_seizure clip.mff';
            RestingFileName = 'JESUS_DA_CONCEI_0626_1056_resting state.mff';
            
            PathSeizure = [DataFolder, SeizureFileName];
            PathResting = [DataFolder, RestingFileName];
            
        case 2
            %SeizureFileName =
            RestingFileName = 'JUNGO NATHALIE 0617 2019_resting state.mff';
            
            PathSeizure = [DataFolder, SeizureFileName];
            PathResting = [DataFolder, RestingFileName];
            
        case 3
            SeizureFileName = 'OGUEY 0805 2113_seizure1.mff';
            RestingFileName = 'OGUEY 0805 2113_resting state.mff';
            
            PathSeizure = [DataFolder, SeizureFileName];
            PathResting = [DataFolder, RestingFileName];
            
        case 4
            SeizureFileName{1} = 'VALLELIAN BENOI 0117 1540_seizure1.mff';
            SeizureFileName{2} = 'VALLELIAN BENOI 0117 1540_seizure2.mff';
            SeizureFileName{3} = 'VALLELIAN BENOI 0117 1540_seizure3.mff';
            SeizureFileNAme{4} = 'VALLELIAN BENOI 0117 1540_seizure4.mff';
            RestingFileName = 'VALLELIAN BENOI 0117 1540_resting state.mff';
            
            for iSeiz = 1: length(SeizureFileName)
                PathSeizure = [DataFolder, SeizureFileName{iSeiz}];
                PathResting = [DataFolder, RestingFileName];
            end
            
    end
    
%     %head = ft_read_data(PathSeizure, 'headerformat', 'egi_mff', ...
%     'fallback', [], ...
%     'checkmaxfilter', true, ...
%     'chanindx',[],...
%     'coordsys', 'head', ...
%     'coilaccuracy', []); 
    %function_header = read_mff_header(PathSeizure); 
    
   head  = read_mff_header(PathSeizure); 
    
end

