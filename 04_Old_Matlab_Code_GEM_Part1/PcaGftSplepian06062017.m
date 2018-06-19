% 2017 06 06 
% In this script I compare the projection into the Slepian space to the
% projection on the first 2nd and 3rd components of the Laplacian

% 2017 05 17

% GFT + SLEPIAN test for JUNGO patient (patient 4 from the paper of Willeke)

% /Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_JUNGO_FOR_SLEPIAN
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

FOLDER_FIGURES = ['/Users/annavybornova/EPFL/Master_4/GEM/05SavedData/20170606PcaGftSlepian', Files{iFile},'/'];

TO_PLOT = 1; % plot figures or not
CONST_NORMALIZE=1;     % normalization of the Adjacency matrix
CONST_OPERATOR=1;      % graph Laplacian
CONST_W = 82;             % bandwidth
CONST_PRINT=0;
CONST_SORTED_GFT = 0; % select GFT vectors according to the most of the
% discriminative power (see SeizureDataDiscoveryV_04_29.m
% or the midterm presentation)
CONST_NORMALIZE_SLEPIAN_AMPLITUDE = 0; 
ichan = 1;
for CONST_SEIZURE_NODE = 203 % ground truth seizure-generating node
    CONST_SUBGRAPH_SIZE = 20; % the size of the subgraph
    CONST_NUM_SLEP_TO_PLOT = 10; % number of lepians to plot
    % what to output?
    CONST_GEN=1; % 0: global (Laplacian eigenvectors)
    % 1: Slepian with energy concentration
    % 2: Slepian with graph embedding distance
    
    
    % Spatial frequencies sorted according to the maximum discriminative power
    %WSORTED = [68;39;2;50;41;64;81;83;5;93;47;12;53;42;34;69;43;78;7;62;80;56;57;77;72;54;136;24;38;48;79;108;9;3;6;36;70;60;37;51;74;63;20;87;107;109;35;71;141;16;163;28;85;45;171;30;52;75;103;97;58;159;31;23;100;67;176;92;26;137;44;61;91;140;196;59;19;123;4;101;76;166;84;49;73;40;88;65;86;46;193;144;172;158;90;175;15;148;55;11;167;189;143;113;174;162;199;8;145;161;105;147;152;98;96;33;99;17;131;138;184;139;185;82;18;13;177;142;89;160;121;128;153;112;104;125;66;150;14;132;102;198;119;116;111;180;194;120;146;95;124;22;94;183;117;29;27;154;10;149;106;129;156;21;157;135;195;126;134;200;115;32;169;114;170;168;186;133;178;179;192;164;181;130;201;203;173;204;187;155;165;25;118;197;191;127;151;110;202;122;182;188;190;1];
    %WSORTED = [39;2;41;5;47;12;42;34;43;24;38;48;9;3;6;36;70;60;37;51;74;63;20;87;107;109;35;71;141;16;163;28;85;45;171;30;52;75;103;97;58;159;31;23;100;67;176;92;26;137;44;61;91;140;196;59;19;123;4;101;76;166;84;49;73;40;88;65;86;46;193;144;172;158;90;175;15;148;55;11;167;189;143;113;174;162;199;8;145;161;105;147;152;98;96;33;99;17;131;138;184;139;185;82;18;13;177;142;89;160;121;128;153;112;104;125;66;150;14;132;102;198;119;116;111;180;194;120;146;95;124;22;94;183;117;29;27;154;10;149;106;129;156;21;157;135;195;126;134;200;115;32;169;114;170;168;186;133;178;179;192;164;181;130;201;203;173;204;187;155;165;25;118;197;191;127;151;110;202;122;182;188;190;1];
    
    %% Compute the adjacency matrix
    
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
    
    % define the nodes that are includes into the selected region
    
    EuDist = sqrt((NODES_XY(:,1) - NODES_XY(CONST_SEIZURE_NODE, 1)).^2 ...
        + (NODES_XY(:,2) - NODES_XY(CONST_SEIZURE_NODE, 2)).^2 ...
        + (NODES_XY(:,3) - NODES_XY(CONST_SEIZURE_NODE, 3)).^2);
    
    [EuDistSorted, idx] = sort(EuDist);
    idx((CONST_SUBGRAPH_SIZE+1):end) = [];
    
    NODE_TYPE(idx)=1;
    CONST_NODES{1}=idx;
    
    % sanity check:
    
    figure,
    plot3(NODES_XY(:,1), NODES_XY(:,2), NODES_XY(:,3), 'O');
    hold on;
    plot3(NODES_XY(idx,1), NODES_XY(idx,2), NODES_XY(idx,3), 'ro');
    hold on;
    plot3(NODES_XY(CONST_SEIZURE_NODE,1), NODES_XY(CONST_SEIZURE_NODE,2),...
        NODES_XY(CONST_SEIZURE_NODE,3), 'k*');
    title (['Subgraph with N_s = ', num2str(CONST_SUBGRAPH_SIZE)] );
    legend('Non-included', 'Included in S', 'Seizure el');
    
    %% Preprocess data, extract (truncated) graph spectrum
    
    [A,D]=slepNormalize(A,CONST_NORMALIZE);
    
    opts.issym=1;
    opts.isreal=1;
    opts.maxit=2500;
    opts.disp=1;
    
    [Utr,S1]=slepEigsLaplacian(A,D,204,opts); % Utr - eigen vectors, S1 - eigen values
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
    end
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
    if TO_PLOT
        set(gcf, 'Position', [100, 100, 1049, 895]);
        saveas(gcf, [FOLDER_FIGURES, 'SL', 'Ns', num2str(CONST_SUBGRAPH_SIZE), 'W', ...
            num2str(length(idxW)), num2str(CONST_SEIZURE_NODE), '.png']);
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
    
    
    %% Spatial PCA of the data (for the whole length of the recording)
    
    % STEP1: common average referencing of the data
    ComRef = @(x) x-repmat(mean(x,2), 1, size(x,2)); 
    DataAllCR = ComRef(DataAll); 
    DataBeforeSeizureCR = ComRef(DataBeforeSeizure); 
    DataPreSeizureCR = ComRef(DataPreSeizure); 
    DataSeizureCR = ComRef(DataSeizure); 
    
    [coeff, score] = pca(DataAllCR',  'NumComponents', 2); 
    
    % visualization of the first 2 components 
    
    figure, 
    subplot(1,2,1); 
    topoplot(coeff(:,1), chanlocs, 'electrodes', 'on', 'headrad', 0); 
    title('PC1')
    subplot(1,2,2); 
    topoplot(coeff(:,2), chanlocs,  'electrodes', 'on', 'headrad', 0); 
    title('PC2'); 
    
    set(gcf, 'Position', [100, 100, 1049, 895]);
    saveas(gcf, [FOLDER_FIGURES, 'S1S2', 'Ns', num2str(CONST_SUBGRAPH_SIZE), 'W', ...
        num2str(length(idxW)), 'PCcomponentsCR.png']);
    
    figure, 
    ax(1)=subplot(1,3,1)
    plot(score(1:300000,1), score(1:300000,2), 'ob'); 
    xlabel('PC1'); 
    ylabel('PC2'); 
    title('Before Seizure'); 
     axis square
    ax(2)=subplot(1,3,2)
    plot(score(480000:600000, 1), score(480000:600000, 2), 'og');
    xlabel('PC1'); 
    ylabel('PC2'); 
    title('Pre Seizure'); 
     axis square
    ax(3)=subplot(1,3,3)
    plot(score(600000:end, 1), score(600000:end, 2), 'or');
    xlabel('PC1'); 
    ylabel('PC2'); 
    title('Seizure');
     axis square
    linkaxes(ax, 'xy'); 
    
    set(gcf, 'Position', [100, 100, 1049, 895]);
    saveas(gcf, [FOLDER_FIGURES, 'S1S2', 'Ns', num2str(CONST_SUBGRAPH_SIZE), 'W', ...
        num2str(length(idxW)), 'DataInPCspaceCR.png']);
    
%     PreSeizOutOfEllipse = OutOfEllipse(score(1:300000,1), score(1:300000,2),...
%         score(480000:600000, 1), score(480000:600000, 2)); 
 %% Spatial PCA of the data (only for the Before Seizure data)  
 [coeff, ~] = pca(DataBeforeSeizureCR', 'NumComponents',2);
 score = (coeff'* DataAll)'; 
     % visualization of the first 2 components 
   PreSeizOutOfEllipse = OutOfEllipse(score(1:300000,1), score(1:300000,2),...
        score(480000:600000, 1), score(480000:600000, 2));  
   SeizureOutOfEllipse = OutOfEllipse(score(1:300000,1), score(1:300000,2),...
        score(600000:end, 1), score(600000:end, 2));
   AllOutOfEllipse = OutOfEllipse(score(1:300000,1), score(1:300000,2),...
        score(:, 1), score(:, 2));
    figure, 
    subplot(1,2,1); 
    topoplot(coeff(:,1), chanlocs, 'electrodes', 'on', 'headrad', 0); 
    title('PC1')
    subplot(1,2,2); 
    topoplot(coeff(:,2), chanlocs,  'electrodes', 'on', 'headrad', 0); 
    title('PC2'); 
    
    set(gcf, 'Position', [100, 100, 1049, 895]);
    saveas(gcf, [FOLDER_FIGURES, 'S1S2', 'Ns', num2str(CONST_SUBGRAPH_SIZE), 'W', ...
        num2str(length(idxW)), 'PCcomponentsCRBefore.png']);
    
    figure, 
    ax(1)=subplot(1,3,1)
    plot(score(1:300000,1), score(1:300000,2), 'ob'); 
    xlabel('PC1'); 
    ylabel('PC2'); 
    title('Before Seizure'); 
     axis square
    ax(2)=subplot(1,3,2)
    plot(score(480000:600000, 1), score(480000:600000, 2), 'og');
    xlabel('PC1'); 
    ylabel('PC2'); 
    title(['Pre Seizure. Score: ', num2str(100*sum(PreSeizOutOfEllipse)/...
        length(PreSeizOutOfEllipse)), ' %']); 
     axis square
    ax(3)=subplot(1,3,3)
    plot(score(600000:end, 1), score(600000:end, 2), 'or');
    xlabel('PC1'); 
    ylabel('PC2'); 
    title(['Seizure. Score: ', num2str(100*sum(SeizureOutOfEllipse)/...
        length(SeizureOutOfEllipse)), '%']);
     axis square
    linkaxes(ax, 'xy'); 
    
    set(gcf, 'Position', [100, 100, 1200, 400]);
    saveas(gcf, [FOLDER_FIGURES, 'S1S2', 'Ns', num2str(CONST_SUBGRAPH_SIZE), 'W', ...
        num2str(length(idxW)), 'DataInPCspaceCRBefore.png']);
    
    figure, 
    hax = axes; 
    plot((1:length(AllOutOfEllipse))/ALLEEG.srate, AllOutOfEllipse); 
    hold on; 
    line([600 600], get(hax,'YLim'), 'Color',[1 0 0], 'Linewidth', 1.5); 
    title('Out of Ellipse'); 
    xlabel('Time [s]'); 
    set(gcf, 'Position', [100, 100, 1200, 400]);
    saveas(gcf, [FOLDER_FIGURES, 'S1S2', 'Ns', num2str(CONST_SUBGRAPH_SIZE), 'W', ...
        num2str(length(idxW)), 'TimeCoursePCs.png']);
    
    clear score PreSeizOutOfEllipse SeizureOutOfEllipse AllOutOfEllipse
    %% Perform the decomposition in the space spanned by the Laplacian 2 
    % and Laplacian 3
    
    V = load('V.mat'); 
    V1 = V.V(:,2); 
    V2 = V.V(:,3); 
    
    score = ([V1 V2]' * DataAll)'; 
    
    % visualization of the first 2 components 
   PreSeizOutOfEllipse = OutOfEllipse(score(1:300000,1), score(1:300000,2),...
        score(480000:600000, 1), score(480000:600000, 2));  
   SeizureOutOfEllipse = OutOfEllipse(score(1:300000,1), score(1:300000,2),...
        score(600000:end, 1), score(600000:end, 2));
   AllOutOfEllipse = OutOfEllipse(score(1:300000,1), score(1:300000,2),...
        score(:, 1), score(:, 2));
    figure, 
    subplot(1,2,1); 
    topoplot(V1, chanlocs, 'electrodes', 'on', 'headrad', 0); 
    title('V2')
    subplot(1,2,2); 
    topoplot(V2, chanlocs,  'electrodes', 'on', 'headrad', 0); 
    title('V1'); 
    
    set(gcf, 'Position', [100, 100, 1049, 895]);
    saveas(gcf, [FOLDER_FIGURES, 'S1S2', 'Ns', num2str(CONST_SUBGRAPH_SIZE), 'W', ...
        num2str(length(idxW)), 'LaplacianEigenMaps.png']);
    
    figure, 
    ax(1) = subplot(1,3,1)
    plot(score(1:300000,1), score(1:300000,2), 'ob'); 
     xlabel('V2'); 
    ylabel('V3'); 
    title('Before Seizure'); 
     axis square
    ax(2)=subplot(1,3,2)
    plot(score(480000:600000, 1), score(480000:600000, 2), 'og');
    xlabel('V2'); 
    ylabel('V3'); 
    title(['Pre Seizure. Score: ', num2str(100*sum(PreSeizOutOfEllipse)/...
        length(PreSeizOutOfEllipse)), '%']); 
     axis square
    ax(3)=subplot(1,3,3)
    plot(score(600000:end, 1), score(600000:end, 2), 'or');
    xlabel('V2'); 
    ylabel('V3'); 
    title(['Seizure. Score: ', num2str(100*sum(SeizureOutOfEllipse)/...
        length(SeizureOutOfEllipse)), '%']);
     axis square
    linkaxes(ax, 'xy'); 
    
    set(gcf, 'Position', [100, 100, 1200, 400]);
    saveas(gcf, [FOLDER_FIGURES, 'S1S2', 'Ns', num2str(CONST_SUBGRAPH_SIZE), 'W', ...
        num2str(length(idxW)), 'DataInLaplacianSpace.png']);
    
    figure, 
    hax = axes; 
    plot((1:length(AllOutOfEllipse))/ALLEEG.srate, AllOutOfEllipse); 
    hold on; 
    line([600 600], get(hax,'YLim'), 'Color',[1 0 0], 'Linewidth', 1.5); 
    title('Out of Ellipse'); 
    xlabel('Time [s]'); 
    set(gcf, 'Position', [100, 100, 1200, 400]);
    saveas(gcf, [FOLDER_FIGURES, 'S1S2', 'Ns', num2str(CONST_SUBGRAPH_SIZE), 'W', ...
        num2str(length(idxW)), 'TimeCourseLaplacian.png']);
    
    
     clear score coeff PreSeizOutOfEllipse SeizureOutOfEllipse AllOutOfEllipse
 
  %% Spatial PCA of the data (non- CR only for the Before Seizure data)  
 [coeff, ~] = pca(DataBeforeSeizure', 'NumComponents',2);
 score = (coeff'* DataAll)'; 
     % visualization of the first 2 components 
   PreSeizOutOfEllipse = OutOfEllipse(score(1:300000,1), score(1:300000,2),...
        score(480000:600000, 1), score(480000:600000, 2));  
   SeizureOutOfEllipse = OutOfEllipse(score(1:300000,1), score(1:300000,2),...
        score(600000:end, 1), score(600000:end, 2));
   AllOutOfEllipse = OutOfEllipse(score(1:300000,1), score(1:300000,2),...
        score(:, 1), score(:, 2));
    figure, 
    subplot(1,2,1); 
    topoplot(coeff(:,1), chanlocs, 'electrodes', 'on', 'headrad', 0); 
    title('PC1')
    subplot(1,2,2); 
    topoplot(coeff(:,2), chanlocs,  'electrodes', 'on', 'headrad', 0); 
    title('PC2'); 
    
    set(gcf, 'Position', [100, 100, 1049, 895]);
    saveas(gcf, [FOLDER_FIGURES, 'S1S2', 'Ns', num2str(CONST_SUBGRAPH_SIZE), 'W', ...
        num2str(length(idxW)), 'PCcomponentsnonCRBefore.png']);
    
    figure, 
    ax(1)=subplot(1,3,1)
    plot(score(1:300000,1), score(1:300000,2), 'ob'); 
    xlabel('PC1'); 
    ylabel('PC2'); 
    title('Before Seizure'); 
     axis square
    ax(2)=subplot(1,3,2)
    plot(score(480000:600000, 1), score(480000:600000, 2), 'og');
    xlabel('PC1'); 
    ylabel('PC2'); 
    title(['Pre Seizure. Score: ', num2str(100*sum(PreSeizOutOfEllipse)/...
        length(PreSeizOutOfEllipse)), ' %']); 
     axis square
    ax(3)=subplot(1,3,3)
    plot(score(600000:end, 1), score(600000:end, 2), 'or');
    xlabel('PC1'); 
    ylabel('PC2'); 
    title(['Seizure. Score: ', num2str(100*sum(SeizureOutOfEllipse)/...
        length(SeizureOutOfEllipse)), '%']);
     axis square
    linkaxes(ax, 'xy'); 
    
    set(gcf, 'Position', [100, 100, 1200, 400]);
    saveas(gcf, [FOLDER_FIGURES, 'S1S2', 'Ns', num2str(CONST_SUBGRAPH_SIZE), 'W', ...
        num2str(length(idxW)), 'DataInPCspacenonCRBefore.png']);
    
    figure, 
    hax = axes; 
    plot((1:length(AllOutOfEllipse))/ALLEEG.srate, AllOutOfEllipse); 
    hold on; 
    line([600 600], get(hax,'YLim'), 'Color',[1 0 0], 'Linewidth', 1.5); 
    title('Out of Ellipse'); 
    xlabel('Time [s]'); 
    set(gcf, 'Position', [100, 100, 1200, 400]);
    saveas(gcf, [FOLDER_FIGURES, 'S1S2', 'Ns', num2str(CONST_SUBGRAPH_SIZE), 'W', ...
        num2str(length(idxW)), 'TimeCoursePCsnonCR.png']);
    
    clear score PreSeizOutOfEllipse SeizureOutOfEllipse AllOutOfEllipse
end