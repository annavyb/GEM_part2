% Plots for the optimal performances

clc
clear all
close all



TO_PLOT = 1;

CONST_NORM_AMPL_SLEP = 0; 

if CONST_NORM_AMPL_SLEP
   FOLDER_FIGURES = '/Users/annavybornova/EPFL/Master_4/GEM/05SavedData/20170717OptimalPerformanceJungo/NormAmpl'  
   load('perfPreSeizNorm.mat');
   load('perfSeizNorm.mat');
   CONST_PERF = 25;
else
    FOLDER_FIGURES = '/Users/annavybornova/EPFL/Master_4/GEM/05SavedData/20170717OptimalPerformanceJungo/';
    load('perfPreSeiz.mat');
    load('perfSeiz.mat');
    CONST_PERF = 40;
end
if TO_PLOT
    figure,
    
    subplot(1,2,1)
    imagesc(perfSeiz');
    title('Performance Seizure')
    xlabel('W');
    ylabel('Ns');
    colorbar;
    axis square;
    
    subplot(1,2,2)
    imagesc(perfPreSeiz')
    title('Performance Pre Seizure')
    xlabel('W');
    ylabel('Ns');
    colorbar;
    axis square
    
    set(gcf, 'Position', [100, 100, 1049, 895]);
    saveas(gcf, [FOLDER_FIGURES, 'AllPerformances.png']);
    
end

[W,Ns] = find(perfSeiz > CONST_PERF);

% Compute the Slepians corresponding to the performances > CONST_PERF

%% Settings
Files = {'JUNGO', 'OGUEY'}; % file number 1 - JUNGO, file 2 - OGUEY
iFile = 1;

% load the workspace for the pre-processed data from patient
load(['/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_',Files{iFile},...
    '_FOR_SLEPIAN/', Files{iFile},'preprocessed.mat'])


TO_PLOT = 1; % plot figures or not
CONST_NORMALIZE= 1;     % normalization of the Adjacency matrix
CONST_OPERATOR=1;      % graph Laplacian
CONST_PRINT=0;
CONST_SORTED_GFT = 0; % select GFT vectors according to the most of the
% discriminative power (see SeizureDataDiscoveryV_04_29.m
% or the midterm presentation)
CONST_NORMALIZE_SLEPIAN_AMPLITUDE = 0;
CONST_SEIZURE_NODE = 203 % ground truth seizure-generating node

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
for iP = 1:length(W)
    for CONST_W = W(iP)
        
        for CONST_SUBGRAPH_SIZE = Ns(iP) % the size of the subgraph
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
                saveas(gcf, [FOLDER_FIGURES, 'SL', 'Ns', num2str(CONST_SUBGRAPH_SIZE), 'W', ...
                    num2str(length(idxW)), num2str(CONST_SEIZURE_NODE), 'norm', ...
                    num2str(CONST_NORMALIZE_SLEPIAN_AMPLITUDE), '.png']);
            end
            
          % Store the Slepians of interest 
          
          Sl1(:, iP) = SL0(:,CONST_W); 
          SL2(:,iP) = SL0(:,CONST_W-1); 
            
        end
    end
end

%% Identify similar Channels

% Sl1 = zscore(Sl1); 
% Sl2 = zscore(SL2); 

% ImpChanSl1 = mean((Sl1)').^2./(std((Sl1)'.^2)); 
% ImpChanSl2 = mean((SL2)').^2./(std((SL2)'.^2)); 



% figure, 
% subplot(1,2,1); 
% topoplot(ImpChanSl1, chanlocs,'headrad', 0, 'electrodes', 'labels');
% 
% subplot(1,2,2); 
% topoplot(ImpChanSl2, chanlocs,'headrad', 0, 'electrodes', 'labels');






