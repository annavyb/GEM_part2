% 2017 08 14
% In this script I perform the Slepian analysis of the new data from HUG.
% (For the moment I will start with 2 patients i.e. Jesus and Vallelian)
% These data were first exposed to the pre-processing steps with the script
% EdfFilePreprocessingForSlepianTest.m and the function
% edfFilePreprocessingForSlepian.m

% In this file I would like to proceed in the following way:

% STEP 0: Bad channel handling, i.e. identifying and removing the bad
% channels. Then further remplacement of the bad channels with the
% interpolation from the neighboring channels

% STEP 1:

clc;
clear all;
close all;

dbstop if error;
set(0,'defaultfigurecolor',[1 1 1]);
set(0,'defaulttextinterpreter','none');

Params = setParam(1);
Params

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add path to the GFT code
addpath('/Users/annavybornova/EPFL/3eme_annee/Bachelor project/Matlab_code');

% add path to the Slepian code
addpath('/Users/annavybornova/EPFL/Master_4/GEM/02MATLAB_Code/SPLdemo');

% add path to the ellipsoide code

addpath('/Users/annavybornova/EPFL/Master_4/GEM/02MATLAB_Code/fit_ellipse');


%% Settings

Files = {'Jesus', 'Vallelian', 'Oguey', 'Jungo'}; % do not change the order here, there is a switch coming later

iFile = 1; % select the file of interest

PlotEeg = 0; % do the eeg plot in order to visualize the bad channels

% define the File-specific parameters

switch iFile
    
    case 1% Case Jesus
        
        BadChannel.Resting = [164];
        BadChannel.Seizure = [];
        NumSeizureFiles = 1;
        CufOffSeizureTime(1) = 470; % seizure activity before the artifact (visual inspection)
        %T8 - 203, P8 - 177, Tp10 - 184
        CONST_SEIZURE_NODE = 184;
        
    case 2% Case Vallelian
        
        BadChannel.Resting = [16, 108, 183];
        BadChannel.Seizure = [108, 183];
        NumSeizureFiles = 4;
        CufOffSeizureTime(1) = 652; % only very small portion
        CufOffSeizureTime(2) = 565.5;
        CufOffSeizureTime(3) = 268;
        CufOffSeizureTime(4) = 24;
        
        %F3 - 55, P7 - 28, CP3 - 30, P5 - 32
        CONST_SEIZURE_NODE = 30;
        
    case 3 % Oguey
        
        BadChannel.Resting = [64, 104, 183];
        BadChannel.Seizure = [92, 101];
        NumSeizureFiles = 1;
        CufOffSeizureTime(1) = 624;
        %F4 - 150
        CONST_SEIZURE_NODE = 150;
    case 4 % Jungo
        BadChannel.Resting = [108, 182, 183];
        BadChannel.Seizure = [64, 104, 108, 182, 183];
        NumSeizureFiles = 1;
        CufOffSeizureTime(1) = 709;
        CONST_SEIZURE_NODE = 203 % ground truth seizure-generating node
        
        
        
end

FileNameResting = ['/Users/annavybornova/EPFL/Master_4/GEM/01Data/20170808DataHUGNewMarkers/PreprocessedDataForSlepian/',...
    Files{iFile}, 'RestingPreProcessedForSlepian.mat'];
DataResting = load(FileNameResting);

if PlotEeg
    eegplot(DataResting.ALLEEG.data, 'srate', DataResting.ALLEEG.srate, ...
        'winlength', 60, 'dispchans', 50); %, 'eloc_file', DataResting.ALLEEG.chanlocs);
end

% Bad channels handling
disp('-> Now Interpolating Bad Channels Resting')
DataResting.ALLEEG.nbchan = 204;
EEGOUT = pop_interp(DataResting.ALLEEG, BadChannel.Resting, 'spherical');
DataResting.ALLEEG.data = EEGOUT.data;


for iSeiz = 1:NumSeizureFiles
    FileNameSeizure = ['/Users/annavybornova/EPFL/Master_4/GEM/01Data/20170808DataHUGNewMarkers/PreprocessedDataForSlepian/',...
        Files{iFile},'Seizure', num2str(iSeiz),'PreProcessedForSlepian.mat'];
    DataSeizure{iSeiz} = load(FileNameSeizure);
    DataSeizure{iSeiz}.ALLEEG.nbchan = 204;
    
   
    
    % Bad channels handling
    disp('--> Now Interpolating Bad Channels Seizure')
    EEGOUT = pop_interp(DataSeizure{iSeiz}.ALLEEG, BadChannel.Resting, 'spherical');
    DataSeizure{iSeiz}.ALLEEG.data = EEGOUT.data;
    
     if PlotEeg
        eegplot(DataSeizure{iSeiz}.ALLEEG.data, 'srate',DataSeizure{iSeiz}.ALLEEG.srate, ...
            'winlength', 60, 'dispchans', 50, 'eloc_file', DataResting.ALLEEG.chanlocs);
    end
    
    SeizureOnset(iSeiz) = round(DataSeizure{iSeiz}.ALLEEG.srate...
        *DataSeizure{iSeiz}.ALLEEG.SeizureOnset/1000000); % seizure onset in samples
    CutOffSeizureSample(iSeiz) = CufOffSeizureTime(iSeiz)*DataSeizure{iSeiz}.ALLEEG.srate;
end

chanlocs = DataResting.ALLEEG.chanlocs;


%% The old code

%% Settings

TO_PLOT = 1; % plot figures or not
CONST_NORMALIZE= 1;     % normalization of the Adjacency matrix
CONST_OPERATOR=1;      % graph Laplacian
%CONST_W = 20;             % bandwidth
CONST_PRINT=0;
CONST_SORTED_GFT = 0; % select GFT vectors according to the most of the
% discriminative power (see SeizureDataDiscoveryV_04_29.m
% or the midterm presentation)
CONST_NORMALIZE_SLEPIAN_AMPLITUDE = 1;


%% Generate the adjacency matrix (something that does not need to be re-calculated at each iteration)

disp('... Adjacency Matrix generation');
Atemp = adj_generator([[chanlocs.X]', [chanlocs.Y]', [chanlocs.Z]'], Params.Gft.R);
nNeighbors = sum(Atemp);
avgNeighbors = mean(nNeighbors);

disp(['R = ', num2str(Params.Gft.R)]);
disp(['avg number of neighbors = ', num2str(avgNeighbors)]);

% converting A into a sparse matrix to make sure that it is conform with
% the SPL demo code

A = sparse(Atemp);
NODES_XY = [[chanlocs.X]' [chanlocs.Y]' [chanlocs.Z]'];

[A,D]=slepNormalize(A,CONST_NORMALIZE);

opts.issym=1;
opts.isreal=1;
opts.maxit=2500;
opts.disp=1;

[Utr,S1]=slepEigsLaplacian(A,D,204,opts); % Utr - eigen vectors, S1 - eigen values

%%
for CONST_W = 78 %20:80
    
    for CONST_SUBGRAPH_SIZE = 42 %20:80 % the size of the subgraph
        CONST_NUM_SLEP_TO_PLOT = 10;% number of Slepians to plot
        % what to output?
        CONST_GEN=1; % 0: global (Laplacian eigenvectors)
        % 1: Slepian with energy concentration
        % 2: Slepian with graph embedding distance
        
        % define the nodes that are includes into the selected region
        EuDist = sqrt((NODES_XY(:,1) - NODES_XY(CONST_SEIZURE_NODE, 1)).^2 ...
            + (NODES_XY(:,2) - NODES_XY(CONST_SEIZURE_NODE, 2)).^2 ...
            + (NODES_XY(:,3) - NODES_XY(CONST_SEIZURE_NODE, 3)).^2);
        
        [EuDistSorted, idx] = sort(EuDist);
        idx((CONST_SUBGRAPH_SIZE+1):end) = [];
        
        NODE_TYPE(idx)=1;
        CONST_NODES{1}=idx;
        
        % sanity check:
        if TO_PLOT
            figure,
            plot3(NODES_XY(:,1), NODES_XY(:,2), NODES_XY(:,3), 'O');
            hold on;
            plot3(NODES_XY(idx,1), NODES_XY(idx,2), NODES_XY(idx,3), 'ro');
            hold on;
            plot3(NODES_XY(CONST_SEIZURE_NODE,1), NODES_XY(CONST_SEIZURE_NODE,2),...
                NODES_XY(CONST_SEIZURE_NODE,3), 'k*');
            title (['Subgraph with N_s = ', num2str(CONST_SUBGRAPH_SIZE)] );
            legend('Non-included', 'Included in S', 'Seizure el');
        end
        
        CONST_U_W=CONST_W;
        
        %% Compute graph Slepians
        idx=CONST_NODES{1};
        if CONST_SORTED_GFT
            idxW = WSORTED(1:CONST_W);
        else
            idxW=1:CONST_W; % Band-limiting of the GFT spectrum
        end
        tmpS=(sum(diag(S1(idxW,idxW))));
        S1=S1/tmpS % some sort of normalization of the eigen values
        tmpUS=Utr(idx,idxW)*sqrt(S1(idxW,idxW)); % The product of the band.limited GFT and the another normalization of eigen values
        tmpU=Utr(idx,idxW);
        C=tmpU.'*tmpU; C=real((C+C.')/2);
        C2=tmpUS.'*tmpUS; C2=real((C2+C2.')/2);
        
        if CONST_GEN==0,
            SL0=Utr(:,idxW);
            eig0=diag(S1(idxW,idxW));
            conc=diag(C);
            cut=diag(C2);
        else
            lambda=CONST_GEN-1; % 0 favors Slepian, 1: favors graph embedding distance
            Cmix=(1-lambda)*C+lambda*C2;
            
            [myV,myD]=eig(Cmix);
            SL0=Utr(:,idxW)*myV;
            eig0=diag(SL0.'*(D-A)*SL0); eig0=eig0(1:CONST_W)/tmpS;
            conc=diag(myV.'*C*myV);
            cut=diag(myV.'*C2*myV);
        end;
        
        
        K = length(idxW)*CONST_SUBGRAPH_SIZE/size(NODES_XY,1);
        disp(['Shannon number = ', num2str(K)]);
        if TO_PLOT
            figure,
            clf;
            subplot(2,1,1)
            plotyy(1:CONST_W,conc,1:CONST_W,cut);
            title(['Shannon number = ', num2str(K)])
            legend({'concentration','graph distance'});
            
            for i = 1:CONST_NUM_SLEP_TO_PLOT
                
                subplot(2,CONST_NUM_SLEP_TO_PLOT,i+CONST_NUM_SLEP_TO_PLOT)
                if i == 1
                    topoplot(zeros(1, length(idx)), chanlocs(idx),  'electrodes', 'on', 'headrad', 0);
                    title(['ROI = ', num2str(CONST_SUBGRAPH_SIZE)]);
                else
                    topoplot(SL0(:,CONST_W-CONST_NUM_SLEP_TO_PLOT + i-1), chanlocs,'headrad', 0);
                    title(['SL =  ', num2str(CONST_W-CONST_NUM_SLEP_TO_PLOT + i-1), 'W = ', num2str(CONST_W)]);
                end
                
            end
            
            set(gcf, 'Position', [100, 100, 1049, 895]);
            %             saveas(gcf, [FOLDER_FIGURES, 'SL', 'Ns', num2str(CONST_SUBGRAPH_SIZE), 'W', ...
            %                 num2str(length(idxW)), num2str(CONST_SEIZURE_NODE), 'norm', ...
            %                 num2str(CONST_NORMALIZE_SLEPIAN_AMPLITUDE), '.png']);
        end
        
        %% Project the data into the space spanned by the first 2 Slepians
        
        SlResting = SL0'*DataResting.ALLEEG.data;
        %SlPreSeizure =
        for iSeiz = 1:length(DataSeizure)
            SlSeizure{iSeiz} = SL0'* DataSeizure{iSeiz}.ALLEEG.data(:, SeizureOnset(iSeiz):CutOffSeizureSample(iSeiz));
            SlPreSeizure{iSeiz} = SL0'* DataSeizure{iSeiz}.ALLEEG.data(:, 1:(SeizureOnset(iSeiz)-1));
        end
        
        
        normSl = @(sl) repmat(sqrt(sum(sl.^2,1)), size(sl,1), 1 );
        
        
        if CONST_NORMALIZE_SLEPIAN_AMPLITUDE
            
            disp('...Eliminating the effect of the global amplitude')
            SlResting = SlResting./normSl(SlResting);
            for iSeiz = 1:length(DataSeizure)
                SlSeizure{iSeiz} = SlSeizure{iSeiz}./normSl(SlSeizure{iSeiz});
                SlPreSeizure{iSeiz} = SlPreSeizure{iSeiz}./normSl(SlPreSeizure{iSeiz});
            end
            
        end
        
        
        
        
        %% Fit an ellipsoide to the before seizure data
        
        % draw the ellipse
        
        axis_handle = gca;
        ellipse_t = fit_ellipse( SlResting(1,:),SlResting(2,:),axis_handle )
        
        % What are the points outside the ellipse for the seizure and pre-seizure
        % data
        
        el = @(x,y,x0, y0, a, b, alpha) (((x-x0)*cos(-alpha) + (y-y0)*sin(-alpha)).^2/(a^2))...
            + (((x-x0)*sin(-alpha) - (y-y0)*cos(-alpha)).^2/(b^2));
        
        coeff = 1;
        step = 0.1;
        s = size(SlResting,2);
        
        while (s > 0.001*size(SlResting,2))
            try
                ind = (el(SlResting(1,:), SlResting(2,:), ellipse_t.X0, ...
                    ellipse_t.Y0, coeff*ellipse_t.a, coeff*ellipse_t.b, ellipse_t.phi )> 1);
                coeff = coeff + step
                s = sum(ind)
            catch
                ind = NaN;
                coeff = NaN;
                s = NaN;
            end
        end
        try
            ellipse_t.a = coeff*ellipse_t.a;
            ellipse_t.b= coeff*ellipse_t.b;
            
            % rotation matrix to rotate the axes with respect to an angle phi
            
            R = [ cos(ellipse_t.phi) sin(ellipse_t.phi); -sin(ellipse_t.phi) cos(ellipse_t.phi) ];
            
            
            % the axes
            ver_line        = [ [ellipse_t.X0 ellipse_t.X0]; ellipse_t.Y0+ellipse_t.b*[-1 1] ];
            horz_line       = [ ellipse_t.X0+ellipse_t.a*[-1 1]; [ellipse_t.Y0 ellipse_t.Y0] ];
            new_ver_line    = R*ver_line;
            new_horz_line   = R*horz_line;
            
            % the ellipse
            theta_r         = linspace(0,2*pi);
            ellipse_x_r     = ellipse_t.X0 + ellipse_t.a*cos( theta_r );
            ellipse_y_r     = ellipse_t.Y0 + ellipse_t.b*sin( theta_r );
            rotated_ellipse = R * [ellipse_x_r;ellipse_y_r];
            
            
            
            % What are the data points outside the resting state ellipse for the
            % pre-seizure and Seizure??
            
            
            for iSeiz = 1:length(DataSeizure)
                indSeiz{iSeiz} = (el(SlSeizure{iSeiz}(1,:), SlSeizure{iSeiz}(2,:), ellipse_t.X0, ...
                    ellipse_t.Y0, ellipse_t.a, ellipse_t.b, ellipse_t.phi )> 1);
                
                indPre{iSeiz} = (el(SlPreSeizure{iSeiz}(1,:), SlPreSeizure{iSeiz}(2,:), ellipse_t.X0, ...
                    ellipse_t.Y0, ellipse_t.a, ellipse_t.b, ellipse_t.phi )> 1);
            end
            %
            %             indAll = (el(SlAll(1,:), SlAll(2,:), ellipse_t.X0, ...
            %                 ellipse_t.Y0, ellipse_t.a, ellipse_t.b, ellipse_t.phi )> 1);
            %
            if TO_PLOT
                subplot(2,1,1)
                imagesc(indPre{1});
                ax = gca;
                set(gca, 'xticklabel', (ax.XTick)/DataSeizure{1}.ALLEEG.srate)
                title('Abnormal pattern Pre Seizure');
                xlabel('Time [s]');
                colormap('gray');
                %plot((1:length(indPre{1}))./DataSeizure{1}.ALLEEG.srate,indPre{1});
                subplot(2,1,2)
                %plot((1:length(indSeiz{1}))./DataSeizure{1}.ALLEEG.srate,indSeiz{1});
                imagesc(indSeiz{1});
                colormap('gray');
                ax = gca;
                set(gca, 'xticklabel', (SeizureOnset + ax.XTick)/DataSeizure{1}.ALLEEG.srate)
                xlabel('Time [s]');
                title('Abnormal pattern Seizure');
                
                 saveas(gcf, ['/Users/annavybornova/EPFL/Master_4/GEM/05SavedData/20170819NewDataHUGPerformance/',...
                    Files{iFile},num2str(iSeiz), 'Slepian1Slepian2.png']);
                
            end
            
            for iSeiz = 1:length(DataSeizure)
                perfSeiz{iSeiz}(CONST_W, CONST_SUBGRAPH_SIZE) = 100*sum(indSeiz{iSeiz})/length(indSeiz{iSeiz});
                perfPreSeiz{iSeiz}(CONST_W, CONST_SUBGRAPH_SIZE) = 100*sum(indPre{iSeiz})/length(indPre{iSeiz});
                abnormalPoints{iSeiz} = find(indPre{iSeiz} == 1)
                firstDetection{iSeiz}(CONST_W, CONST_SUBGRAPH_SIZE) = ...
                    abnormalPoints{iSeiz}(1)/DataSeizure{iSeiz}.ALLEEG.srate;
            end
            
            
            %%%%%%
            % PCA sanity check
            
            [coeff, ~] = pca(DataResting.ALLEEG.data', 'NumComponents',2);
            
            scoreResting = (coeff'* DataResting.ALLEEG.data)'; 
            scoreSeizure = (coeff'* DataSeizure{iSeiz}.ALLEEG.data)';
        
            
            PreSeizOutOfEllipsePCA = OutOfEllipse(scoreResting(:,1), scoreResting(:,2),...
                scoreSeizure(1:SeizureOnset, 1), scoreSeizure(1:SeizureOnset, 2));
            SeizureOutOfEllipsePCA = OutOfEllipse(scoreResting(:,1), scoreResting(:,2),...
                scoreSeizure(SeizureOnset:CutOffSeizureSample, 1), scoreSeizure(SeizureOnset:CutOffSeizureSample, 2));
          
            perfPreSeizPCA = 100*sum(PreSeizOutOfEllipsePCA)/length(PreSeizOutOfEllipsePCA); 
            perfSeizPCA = 100*sum(SeizureOutOfEllipsePCA)/length(SeizureOutOfEllipsePCA); 
            % visualize PC1 and PC2 
             
            if TO_PLOT
                
               figure
               subplot(1,2,1)
                topoplot(coeff(:,1), chanlocs,'headrad', 0);
                title('PC1'); 
               subplot(1,2,2)
               topoplot(coeff(:,2), chanlocs, 'headrad', 0); 
               title('PC2'); 
               
                saveas(gcf, ['/Users/annavybornova/EPFL/Master_4/GEM/05SavedData/20170819NewDataHUGPerformance/',...
             Files{iFile},num2str(iSeiz), 'PC1PC2.png']);
               
               
            % plot the time evolution of the abnormal pattern in PC1 and
            % PC2 space 
                figure, 
                subplot(2,1,1)
                imagesc(PreSeizOutOfEllipsePCA');
                ax = gca;
                set(gca, 'xticklabel', (ax.XTick)/DataSeizure{1}.ALLEEG.srate)
                title(['Abnormal pattern Pre Seizure', 'perf = ', num2str(perfPreSeizPCA)]);
                xlabel('Time [s]');
                colormap('gray');
                %plot((1:length(indPre{1}))./DataSeizure{1}.ALLEEG.srate,indPre{1});
                subplot(2,1,2)
                %plot((1:length(indSeiz{1}))./DataSeizure{1}.ALLEEG.srate,indSeiz{1});
                imagesc(SeizureOutOfEllipsePCA');
                colormap('gray');
                ax = gca;
                set(gca, 'xticklabel', (SeizureOnset + ax.XTick)/DataSeizure{1}.ALLEEG.srate)
                xlabel('Time [s]');
                title(['Abnormal pattern Seizure', 'perf = ', num2str(perfSeizPCA)]);
                 saveas(gcf, ['/Users/annavybornova/EPFL/Master_4/GEM/05SavedData/20170819NewDataHUGPerformance/',...
                    Files{iFile},num2str(iSeiz), 'PC1PC2timepoints.png']);
            end
            
            
            %%%%%%% GFT sanity check %%%%%%
           
            load('V.mat'); 
            gftCoeff = V(:,2:3); 
            
            gftScoreRest = (gftCoeff'* DataResting.ALLEEG.data)'; 
            gftScoreSeiz = (gftCoeff'* DataSeizure{iSeiz}.ALLEEG.data)';
             
             PreSeizOutOfEllipseGFT = OutOfEllipse(gftScoreRest(:,1), gftScoreRest(:,2),...
               gftScoreSeiz(1:SeizureOnset, 1), gftScoreSeiz(1:SeizureOnset, 2));
            SeizureOutOfEllipseGFT = OutOfEllipse(gftScoreRest(:,1), gftScoreRest(:,2),...
                gftScoreSeiz(SeizureOnset:CutOffSeizureSample, 1), gftScoreSeiz(SeizureOnset:CutOffSeizureSample, 2));
          
            perfPreSeizGFT = 100*sum(PreSeizOutOfEllipseGFT)/length(PreSeizOutOfEllipseGFT); 
            perfSeizGFT = 100*sum(SeizureOutOfEllipseGFT)/length(SeizureOutOfEllipseGFT); 
           
            % visualization 
               
            if TO_PLOT
                
               figure
               subplot(1,2,1)
                topoplot(gftCoeff(:,1), chanlocs,'headrad', 0);
                title('GFT1'); 
               subplot(1,2,2)
               topoplot(gftCoeff(:,2), chanlocs, 'headrad', 0); 
               title('GFT2'); 
                saveas(gcf, ['/Users/annavybornova/EPFL/Master_4/GEM/05SavedData/20170819NewDataHUGPerformance/',...
             Files{iFile},num2str(iSeiz), 'GFT1GFT2.png']);
               
            % plot the time evolution of the abnormal pattern in GFT1 and
            % GFT2 space 
                figure, 
                subplot(2,1,1)
                imagesc(PreSeizOutOfEllipseGFT');
                ax = gca;
                set(gca, 'xticklabel', (ax.XTick)/DataSeizure{1}.ALLEEG.srate)
                title(['Abnormal pattern Pre Seizure GFT', 'perf = ', num2str(perfPreSeizGFT)]);
                xlabel('Time [s]');
                colormap('gray');
                %plot((1:length(indPre{1}))./DataSeizure{1}.ALLEEG.srate,indPre{1});
                subplot(2,1,2)
                %plot((1:length(indSeiz{1}))./DataSeizure{1}.ALLEEG.srate,indSeiz{1});
                imagesc(SeizureOutOfEllipseGFT');
                colormap('gray');
                ax = gca;
                set(gca, 'xticklabel', (SeizureOnset + ax.XTick)/DataSeizure{1}.ALLEEG.srate)
                xlabel('Time [s]');
                title(['Abnormal pattern Seizure GFT', 'perf = ', num2str(perfSeizGFT)]);
                
                 saveas(gcf, ['/Users/annavybornova/EPFL/Master_4/GEM/05SavedData/20170819NewDataHUGPerformance/',...
                    Files{iFile},num2str(iSeiz), 'GFT1GFT2timepoints.png']);
            end
            
            clear indPre indSeiz indAll abnormalPoints
            
            
        catch
            
        end
        close all;
    end
    
end
%%
for iSeiz = 1:length(DataSeizure)
    figure,
    
    subplot(1,3,1)
    imagesc(perfSeiz{iSeiz}');
    title('Performance Seizure')
    xlabel('W');
    ylabel('Ns');
    axis square
    
    subplot(1,3,2)
    imagesc(perfPreSeiz{iSeiz}')
    title('Performance Pre Seizure')
    xlabel('W');
    ylabel('Ns');
    axis square
    
    subplot(1,3,3)
    imagesc(firstDetection{iSeiz})
    title('First detection')
    xlabel('W');
    ylabel('Ns');
    axis square
    %
    %     saveas(gcf, ['/Users/annavybornova/EPFL/Master_4/GEM/05SavedData/20170819NewDataHUGPerformance/',...
    %         Files{iFile},num2str(iSeiz), 'perf.fig']);
end

% save(['/Users/annavybornova/EPFL/Master_4/GEM/05SavedData/20170819NewDataHUGPerformance/',...
%         Files{iFile},num2str(iSeiz), 'workspace']);

