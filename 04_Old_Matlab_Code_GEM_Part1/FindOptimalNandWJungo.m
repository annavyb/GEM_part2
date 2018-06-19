% 2017 07 14
% In this script I would like to investigate what would be the optimal
% combination of the number of nodes Ns and the frequency band W to obtain
% the best performances for the discrimination between the baseline,
% pre-seizure and seizure

% ( I take the inspiration from the script PhaseAnalysis0706 )


clc
clear
close all

dbstop if error;
set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaulttextinterpreter','none')

% add path to the GFT code
addpath('/Users/annavybornova/EPFL/3eme_annee/Bachelor project/Matlab_code');

% add path to the Slepian code
addpath('/Users/annavybornova/EPFL/Master_4/GEM/02MATLAB_Code/SPLdemo');

% add path to the ellipsoide code

addpath('/Users/annavybornova/EPFL/Master_4/GEM/02MATLAB_Code/fit_ellipse');



%% Settings
Files = {'JUNGO', 'OGUEY'}; % file number 1 - JUNGO, file 2 - OGUEY
iFile = 1;

% load the workspace for the pre-processed data from patient
load(['/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_',Files{iFile},...
    '_FOR_SLEPIAN/', Files{iFile},'preprocessed.mat'])


FOLDER_FIGURES = ['/Users/annavybornova/EPFL/Master_4/GEM/05SavedData/20170607PhaseAnalysis', Files{iFile},'/'];

TO_PLOT = 0; % plot figures or not
CONST_NORMALIZE= 1;     % normalization of the Adjacency matrix
CONST_OPERATOR=1;      % graph Laplacian
%CONST_W = 20;             % bandwidth
CONST_PRINT=0;
CONST_SORTED_GFT = 0; % select GFT vectors according to the most of the
% discriminative power (see SeizureDataDiscoveryV_04_29.m
% or the midterm presentation)
CONST_NORMALIZE_SLEPIAN_AMPLITUDE = 1;
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
for CONST_W = 10:100
    
    for CONST_SUBGRAPH_SIZE = 10:100 % the size of the subgraph
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
        
        %% Project the data into the space spanned by the first 2 Slepians
        
        SlBeforeSeizure = SL0'*DataBeforeSeizure;
        SlPreSeizure = SL0'*DataPreSeizure;
        SlSeizure = SL0'* DataSeizure;
        SlAll = SL0'*DataAll;
        
        normSl = @(sl) repmat(sqrt(sum(sl.^2,1)), size(sl,1), 1 );
        
        
            if CONST_NORMALIZE_SLEPIAN_AMPLITUDE
        
                disp('...Eliminating the effect of the global amplitude')
                SlBeforeSeizure = SlBeforeSeizure./normSl(SlBeforeSeizure);
                SlPreSeizure = SlPreSeizure./normSl(SlPreSeizure);
                SlSeizure = SlSeizure./normSl(SlSeizure);
                SlAll = SlAll./normSl(SlAll);
        
            end
        
        
        
        
        %% Fit an ellipsoide to the before seizure data
        
        % draw the ellipse
        
        axis_handle = gca;
        ellipse_t = fit_ellipse( SlBeforeSeizure(1,:),SlBeforeSeizure(2,:),axis_handle )
        
        % What are the points outside the ellipse for the seizure and pre-seizure
        % data
        
        el = @(x,y,x0, y0, a, b, alpha) (((x-x0)*cos(-alpha) + (y-y0)*sin(-alpha)).^2/(a^2))...
            + (((x-x0)*sin(-alpha) - (y-y0)*cos(-alpha)).^2/(b^2));
        
        coeff = 1;
        step = 0.1;
        s = size(SlBeforeSeizure,2);
        
        while (s > 0.001*size(SlBeforeSeizure,2))
            try
                ind = (el(SlBeforeSeizure(1,:), SlBeforeSeizure(2,:), ellipse_t.X0, ...
                    ellipse_t.Y0, coeff*ellipse_t.a, coeff*ellipse_t.b, ellipse_t.phi )> 1);
                coeff = coeff + step
                s = sum(ind)
            catch
                ind = NaN;
                coeff = NaN;
                s = NaN;
            end
        end
        
        ellipse_t.a = coeff*ellipse_t.a;
        ellipse_t.b= coeff*ellipse_t.b;
        
        % rotation matrix to rotate the axes with respect to an angle phi
        
        R = [ cos(ellipse_t.phi) sin(ellipse_t.phi); -sin(ellipse_t.phi) cos(ellipse_t.phi) ];
        
        try
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
            
            indPre = (el(SlPreSeizure(1,:), SlPreSeizure(2,:), ellipse_t.X0, ...
                ellipse_t.Y0, ellipse_t.a, ellipse_t.b, ellipse_t.phi )> 1);
            
            indSeiz = (el(SlSeizure(1,:), SlSeizure(2,:), ellipse_t.X0, ...
                ellipse_t.Y0, ellipse_t.a, ellipse_t.b, ellipse_t.phi )> 1);
            
            indAll = (el(SlAll(1,:), SlAll(2,:), ellipse_t.X0, ...
                ellipse_t.Y0, ellipse_t.a, ellipse_t.b, ellipse_t.phi )> 1);
        
             perfPreSeiz(CONST_W, CONST_SUBGRAPH_SIZE) = 100*sum(indPre)/length(indPre); 
             perfSeiz(CONST_W, CONST_SUBGRAPH_SIZE) = 100*sum(indSeiz)/length(indSeiz); 
             
             clear indPre indSeiz indAll
             
        catch
            
        end
        close all; 
    end
    
end

figure, 

subplot(1,2,1)
imagesc(perfSeiz'); 
title('Performance Seizure')
xlabel('W'); 
ylabel('Ns'); 

subplot(1,2,2)
imagesc(perfPreSeiz')
title('Performance Pre Seizure')
xlabel('W'); 
ylabel('Ns');
axis square