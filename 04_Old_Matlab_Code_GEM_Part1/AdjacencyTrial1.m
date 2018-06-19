% 2017 02 24
% In this script I explore the construction of the adjacency matrix for the
% 10-20 EEG cap with a bipolar configuration of the recording

% The thing is that I have the coordinates of the electrodes themselves
% However, in a bipolar montage we record the voltage difference between 2
% adjacent electrodes. Thus, to construct the graph, I choose that the
% nodes are located in the middle point between each pair of the electrode

clc
clear all
close all

dbstop if error
set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaulttextinterpreter','none')


% adding paths

PathToOldCode = ['/Users/annavybornova/EPFL/3eme_annee/Bachelor project/Matlab_code'];
addpath(PathToOldCode);


%% Step 1: Loading the EEG of interest

[Electrodes.Coord, Electrodes.Name, ~ ] = xlsread('/Users/annavybornova/EPFL/Master_4/SurrogateDataForEpilepsy/01Data/EegCap/Coord');
[File.TimeBeginEnd, File.Name, ~ ] = xlsread('/Users/annavybornova/EPFL/Master_4/SurrogateDataForEpilepsy/01Data/BostonDatabase/Summary.xlsx');

for iFile = 1:(length(File.Name)-1)
    
    FileName = File.Name{iFile};
    
    Onset = File.TimeBeginEnd(iFile, 1);
    EndSeizure = File.TimeBeginEnd(iFile,2)
    
    Path = ['/Users/annavybornova/EPFL/Master_4/SurrogateDataForEpilepsy/01Data/BostonDatabase/'];
    cd(Path);
    
    ALLEEG = pop_biosig(FileName);
    ChanName = {ALLEEG.chanlocs.labels};
    ALLEEG.chanlocs = [];
    
    %% STEP 2 : Recalculating channel locations and deleting some channels
    
    ind = 1;
    for iChan = 1:length(ChanName)
        
        name = ChanName{iChan};
        for iName = 1:length(name)
            if strcmp(name(iName), '-')
                SepInd = iName;
            end
        end
        Electrode1.Name = name(1:(SepInd-1));
        Electrode2.Name = name((SepInd+1):end);
        
        for iCoord = 1:(size(Electrodes.Coord,1))
            if strcmp(Electrode1.Name, Electrodes.Name(iCoord))
                Electrode1.Coord = Electrodes.Coord(iCoord, :);
            end
            if strcmp(Electrode2.Name, Electrodes.Name(iCoord))
                Electrode2.Coord = Electrodes.Coord(iCoord, :);
            end
        end
        if isfield(Electrode1, 'Coord') && isfield(Electrode2, 'Coord')
            ChanLoc(iChan, :) = (Electrode1.Coord + Electrode2.Coord)/2;
            ALLEEG.chanlocs(ind).labels = name;
            ALLEEG.chanlocs(ind).X = ChanLoc(iChan,1);
            ALLEEG.chanlocs(ind).Y = ChanLoc(iChan,2);
            ALLEEG.chanlocs(ind).Z = ChanLoc(iChan,3);
            SubstData(ind, :) = ALLEEG.data(iChan, :);
            ind = ind+1;
        else
            disp([name, ' is deleted ... '])
        end
        clear name SepInd Electrode1 Electrode2
        
    end
    ALLEEG.data = [];
    ALLEEG.data = SubstData;
    
    %% STEP 3: generation of the adjacency matrix with the new coordinates
    % Question what would be the optimal threshold distance?
    R = 80; % threshold distance for the definition of the adjacency matrix
    
    %figure,
    
    A = adj_generator([[ALLEEG.chanlocs.X]', [ALLEEG.chanlocs.Y]', [ALLEEG.chanlocs.Z]'], R);
    nNeighbors = sum(A);
    avgNeighbors = mean(nNeighbors);
%     imagesc(A);
%     title(['R = ', num2str(R), 'mm ', 'Avg Neighbors = ', num2str(avgNeighbors)],...
%         'FontSize', 16);
%     xticks(1:length(ChanName));
%     yticks(1:length(ChanName));
%     set(gca, 'XTickLabel', {ALLEEG.chanlocs.labels});
%     set(gca, 'YTickLabel', {ALLEEG.chanlocs.labels});
%     set(gca, 'XTickLabelRotation', 45);
%     set(gca, 'YTickLabelRotation', 45);
%     axis 'square'
    
    
    % I chose for now the adjacency matrix with an average of 4.3 neighbors
    
    %% STEP 4: Graph Fourier Transform
    
    L = Laplacian(A);
    %calculate eigen values of L
    L_eig = eig(L);
    L_eig_sorted = sort(L_eig);
    [V, D] = eig(L);
    F = V'*ALLEEG.data;
    
    %find the samples a little bit before the onset
    
    Indeces = find(((ALLEEG.times > (Onset*1000 - (5*6000))) + (ALLEEG.times < (EndSeizure*1000+(5*6000)))) == 2);
    
    figure,
    ax(1) = subplot(2,1,1)
    imagesc(F(:, Indeces))
    set(gca, 'XTick', 1:1000:length(Indeces));
    set(gca, 'XTickLabel', ALLEEG.times(Indeces(1:1000:length(Indeces))));
    set(gca, 'XTickLabelRotation', 45);
    title([FileName(1:8), ' Seizure: ', num2str(Onset), '-', num2str(EndSeizure)]);
    colorbar;
    ax(2) = subplot(2,1,2)
    plot(ALLEEG.times(Indeces), F(:,Indeces));
    xlabel('time [s]');
    
    savefig([FileName(1:8),'SpatialFrequencySeizure'])
    clear ALLEEG
end



