% 2017 03 20
% Here I would like to discover the HUG epilepsy database, check the
% markers etc

clc;
clear all;
close all;

dbstop if error

addpath('/Users/annavybornova/EPFL/Master_4/GEM/02MATLAB_Code')
PathToData = '/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_Seizure_Data';
PathToGftCode = ['/Users/annavybornova/EPFL/3eme_annee/Bachelor project/Matlab_code'];
addpath('PathToData');
addpath(PathToGftCode);
addpath('/Users/annavybornova/EPFL/Master_4/GEM/02MATLAB_Code'); 

FileName = {'02EEGJesus','03EEGJungo', ...
    '06EEGVallelian', '07EEGAlessi', '08EEGBersier'};
cd('../01Data/HUG_Seizure_Data/');

xlsRange = 1;
OptionsFiltering.LowCf=1;
OptionsFiltering.HighCf=30;
OptionsFiltering.FiltOrder=8;

%%
% adding chanlocs
chanlocs = ...
    readlocs('/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_Seizure_Data/Electrode_loc/EGI204.GenevaAverage13.10-10.xyz')

% visualize the electrode locations
% figure,
% topoplot(ones(1,204), chanlocs, 'electrodes', 'ptslabels');

% Generate the adjacency matrix based on the chanlocs
radia = [30];

for R = radia % threshold distance for the definition of the adjacency matrix
    
    A = adj_generator([[chanlocs.X]', [chanlocs.Y]', [chanlocs.Z]'], R);
    nNeighbors = sum(A);
    avgNeighbors = mean(nNeighbors);
    
    disp(['R = ', num2str(R)]);
    disp(['avg number of neighbors = ', num2str(avgNeighbors)]);
    
    % visualize the adjacency matrix
    figure,
    imagesc(A);
    set(gca, 'YTick', 1: 204);
    set(gca, 'XTick', 1:204);
    set(gca, 'YTickLabels', {chanlocs.labels}, 'FontSize', 6);
    set(gca, 'XTickLabels', {chanlocs.labels}, 'FontSize', 6);
    title(['Adjacency matrix. R = ', num2str(R),...
        'AvgNeigh = ', num2str(avgNeighbors)], 'FontSize', 12);
    set(gca, 'XTickLabelRotation', 50);
    set(gca, 'YTickLabelRotation', 50);
    
    
    %%
    for iFile = 3%:length(FileName)
        %% raw data
        dirName = [PathToData,'/', FileName{iFile}];
        files = dir( fullfile(dirName,'*.edf') );
        cd(FileName{iFile});
        ALLEEG = pop_biosig(files.name);
        
        %     figure,
        %     % I divide time by 1000 because it is in ms (at least I think so)
        %     plot(zscore(ALLEEG.data'),'DisplayName','ALLEEG.data');
        %     title(['ALLEEG Sub', FileName{iFile}, 'fs = ', num2str(ALLEEG.srate)]);
        %     xlabel('Samples');
        %     saveas(gcf,[ '/Users/annavybornova/EPFL/Master_4/GEM/05SavedData/20170321DiscoveringHugEpilepsyData/',...
        %         FileName{iFile},'zScoreChannels.png']);
        cd ..
        
       
        
        % Assess the spectrum of the channel
        Channel = 2;
        figure,
        dataFourier = abs(fft(ALLEEG.data'));
        f = ALLEEG.srate*(0:(size(dataFourier,1)/2))/(size(dataFourier,1));
        plot(f(2:end), dataFourier(2:(size(dataFourier,1))/2+1,Channel));
        title(['Fourier spectrum Sub ', FileName{iFile}, 'Ch', num2str(Channel)]);
        xlabel('frequency [Hz]');
        
        saveas(gcf,[ '/Users/annavybornova/EPFL/Master_4/GEM/05SavedData/20170321DiscoveringHugEpilepsyData/',...
            FileName{iFile},'FourierSpectrum.png']);
        %close all
        
        InfoToWrite(iFile).FileName = FileName{iFile};
        InfoToWrite(iFile).Fs = ALLEEG.srate;
        InfoToWrite(iFile).TotSamples =  length(ALLEEG.times);
        InfoToWrite(iFile).TotDuration = ALLEEG.times(end)/1000;
        
        
        %% pre-processing of the data
        % use the function from eeglab
        
        % average common referencing
        ALLEEGreref = ALLEEG;
        ALLEEGreref.data = reref(ALLEEG.data);
        
        %filtering
        ALLEEGfilt = ALLEEGreref;
        %     ALLEEGfilt.data = eegfilt(ALLEEGreref.data, ALLEEGreref.srate, 0, 1);
        %     ALLEEGfilt.data = eegfilt(ALLEEGfilt.data, ALLEEGfilt.srate,0,30);
        %     figure,
        %     dataFourierFilt = abs(fft(ALLEEGfilt.data'));
        %     f = ALLEEG.srate*(0:(size(dataFourier,1)/2))/(size(dataFourier,1));
        %     plot(f(2:end), dataFourierFilt(2:(size(dataFourier,1))/2+1,Channel));
        %     title(['Fourier spectrum Filt Sub ', FileName{iFile}, 'Ch', num2str(Channel)]);
        %     xlabel('frequency [Hz]');
        
        OptionsFiltering.Srate = ALLEEGreref.srate; 
        OptionsFiltering.PlotFreqz = true; 
        OptionsFiltering.PlotSpectrum = true; 
        OptionsFiltering.PlotSpectrumChannel = 2; 
        
        DataF = filterEEG(ALLEEGreref.data, OptionsFiltering); 
        ALLEEGfilt.data = DataF; 
        ALLEEGfiltz = ALLEEGfilt; 
        ALLEEGfiltz.data = (zscore(ALLEEGfilt.data'))'; 
        
        
         figure, 
        eegplot(ALLEEGfilt.data, 'srate', ALLEEG.srate); 
        %% Graph Fourier Transform
        L = Laplacian(A);
        %calculate eigen values of L
        L_eig = eig(L);
        L_eig_sorted = sort(L_eig);
        [V, D] = eig(L);
        
        F = V'*ALLEEGfilt.data;
        Fz = V'* ALLEEGfiltz.data; 
        
        figure,
        subplot(2,1,1)
        imagesc(F(:, 1:round(0.75*size(F,2))));
        title(['R = ', num2str(R), 'Sub ', FileName{iFile}, ' sr = ', num2str(ALLEEG.srate), 'GFT']);
        
        subplot(2,1,2)
        imagesc(ALLEEGfilt.data(:, 1:round(0.75*size(F,2))));
        
        figure, 
        subplot(2,1,1)
        imagesc(Fz(:, 1:round(0.75*size(Fz,2))));
        title(['R = ', num2str(R), 'Sub ', FileName{iFile}, ' sr = ', num2str(ALLEEG.srate), 'GFTz']);
        
        subplot(2,1,2)
        imagesc(ALLEEGfiltz.data(:, 1:round(0.75*size(F,2))));
        
        
        clear ALLEEG ALLEEGreref ALLEEGfilt D d1 d2 dataFourier dataFourierFilt ...
            f F h1 h2  L L_eig L_eig_sorted ...
            nn nNeighbors V
    end
    % dlmwrite('/Users/annavybornova/EPFL/Master_4/GEM/05SavedData/20170321DiscoveringHugEpilepsyData/HugSubjectInfo.txt', ...
    %     InfoToWrite);
end
