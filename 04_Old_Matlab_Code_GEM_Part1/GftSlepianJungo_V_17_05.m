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

FOLDER_FIGURES = ['/Users/annavybornova/EPFL/Master_4/GEM/05SavedData/20170517SlepianGft', Files{iFile},'/'];

TO_PLOT = 1; % plot figures or not
CONST_NORMALIZE=1;     % normalization of the Adjacency matrix
CONST_OPERATOR=1;      % graph Laplacian
CONST_W = 86;             % bandwidth
CONST_PRINT=0;
CONST_SORTED_GFT = 0; % select GFT vectors according to the most of the
% discriminative power (see SeizureDataDiscoveryV_04_29.m
% or the midterm presentation)
CONST_NORMALIZE_SLEPIAN_AMPLITUDE = 0; 
ichan = 1;
for CONST_SEIZURE_NODE = 203 % ground truth seizure-generating node
    CONST_SUBGRAPH_SIZE = 68; % the size of the subgraph
    CONST_NUM_SLEP_TO_PLOT = 10; % number of Slepians to plot
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
    plot3(NODES_XY(:,1), NODES_XY(:,2), NODES_XY(:,3), 'k');
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
    
    normSl = @(sl) repmat(sqrt(sum(sl.^2,1)), size(sl,1), 1 ); 
   
    
    SlBeforeSeizure = SL0'*DataBeforeSeizure;
    SlPreSeizure = SL0'*DataPreSeizure;
    SlSeizure = SL0'* DataSeizure;
    SlAll = SL0'*DataAll;
    
    
%     if CONST_NORMALIZE_SLEPIAN_AMPLITUDE
%         disp('...Eliminating the effect of the global amplitude')
%         SlBeforeSeizure = SlBeforeSeizure./normSl(SlBeforeSeizure);
%         SlPreSeizure = SlPreSeizure./normSl(SlPreSeizure);
%         SlSeizure = SlSeizure./normSl(SlSeizure);
%         SlAll = SlAll./normSl(SlAll);
%         
%     end
    
    % figure,
    % h(1) = subplot(2,2,1);
    % plot(SlBeforeSeizure(1,:), SlBeforeSeizure(5,:), 'bo');
    % xlabel('S1');
    % ylabel('S5');
    % h(2) = subplot(2,2,2);
    % plot(SlPreSeizure(1,:), SlPreSeizure(5,:), 'go');
    % xlabel('S1');
    % ylabel('S5');
    % h(3) = subplot(2,2,3);
    % plot(SlSeizure(1,:), SlSeizure(5,:), 'ro');
    % xlabel('S1');
    % ylabel('S5');
    % subplot(2,2,4)
    % hold on;
    % plot(SlPreSeizure(1,:), SlPreSeizure(5,:), 'go');
    % plot(SlSeizure(1,:), SlSeizure(5,:), 'ro');
    % plot(SlBeforeSeizure(1,:), SlBeforeSeizure(5,:), 'bo');
    %
    % xlabel('S1');
    % ylabel('S5');
    %
    % legend('Pre', 'Seizure', 'BeforeSeizure');
    % linkaxes(h, 'xy');
    %
    % figure,
    % h(1) = subplot(2,2,1);
    % plot(SlBeforeSeizure(1,:), SlBeforeSeizure(2,:), 'bo');
    % xlabel('S1');
    % ylabel('S2');
    %
    % axis square
    %
    % h(2) = subplot(2,2,2);
    % plot(SlPreSeizure(1,:), SlPreSeizure(2,:), 'go');
    % xlabel('S1');
    % ylabel('S2');
    %
    % axis square
    %
    % h(3) = subplot(2,2,3);
    % plot(SlSeizure(1,:), SlSeizure(2,:), 'ro');
    % xlabel('S1');
    % ylabel('S2');
    %
    % axis square
    %
    % subplot(2,2,4)
    % hold on;
    % plot(SlPreSeizure(1,:), SlPreSeizure(2,:), 'go');
    % plot(SlSeizure(1,:), SlSeizure(2,:), 'ro');
    % plot(SlBeforeSeizure(1,:), SlBeforeSeizure(2,:), 'bo');
    % xlabel('S1');
    % ylabel('S2');
    %
    % axis square
    %
    % legend('Pre', 'Seizure', 'BeforeSeizure');
    % linkaxes(h, 'xy');
    %
    % set(gcf, 'Position', [100, 100, 1049, 895]);
    % saveas(gcf, [FOLDER_FIGURES, 'S1S2', 'Ns', num2str(CONST_SUBGRAPH_SIZE), 'W', ...
    %     num2str(length(idxW)), '.png']);
    
    
    % figure,
    % hold on;
    %
    % SlopeS1S2Before = (tan(SlBeforeSeizure(1,:)./ SlBeforeSeizure(2,:)))';
    % SlopeS1S2Pre =    (tan(SlPreSeizure(1,:)./ SlPreSeizure(2,:)))';
    % SlopeS1S2Seizure =    (tan(SlSeizure(1,:)./ SlSeizure(2,:)))';
    %
    % figure,
    % bar([mean(SlopeS1S2Before) mean(SlopeS1S2Pre) mean(SlopeS1S2Seizure)]');
    %
    % figure,
    % errorbar([mean(SlopeS1S2Before) mean(SlopeS1S2Pre) mean(SlopeS1S2Seizure)], ...
    %     [std(SlopeS1S2Before) std(SlopeS1S2Pre) std(SlopeS1S2Seizure)]);
    %
    % xlim([-1 5])
    
    %[std(SlopeS1S2Before) std(SlopeS1S2Pre) std(SlopeS1S2Seizure)]');
    
    
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
        ind = (el(SlBeforeSeizure(1,:), SlBeforeSeizure(2,:), ellipse_t.X0, ...
            ellipse_t.Y0, coeff*ellipse_t.a, coeff*ellipse_t.b, ellipse_t.phi )> 1);
        coeff = coeff + step
        s = sum(ind)
    end
    
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
    
    indPre = (el(SlPreSeizure(1,:), SlPreSeizure(2,:), ellipse_t.X0, ...
        ellipse_t.Y0, ellipse_t.a, ellipse_t.b, ellipse_t.phi )> 1);
    
    indSeiz = (el(SlSeizure(1,:), SlSeizure(2,:), ellipse_t.X0, ...
        ellipse_t.Y0, ellipse_t.a, ellipse_t.b, ellipse_t.phi )> 1);
    
    indAll = (el(SlAll(1,:), SlAll(2,:), ellipse_t.X0, ...
        ellipse_t.Y0, ellipse_t.a, ellipse_t.b, ellipse_t.phi )> 1);
    
    
    % @ 2017 05 28
    % What are the channels outside the normal zone for a time window 30s with
    % the overlap of 15 s
    
    WIN = 30*ALLEEG.srate;
    OVERLAP = WIN/2;
    
    indPreWin = buffer(indPre, WIN, OVERLAP, 'nodelay');
    perfPreWin(ichan, :) = 100*sum(indPreWin)/WIN;
    
    indSeizWin = buffer(indSeiz, WIN, OVERLAP, 'nodelay');
    perfSeizWin(ichan, :) = 100*sum(indSeizWin)/WIN;
    
    indAllWin = buffer(indAll, WIN, OVERLAP, 'nodelay');
    idxAllWin = buffer(1:length(indAll), WIN, OVERLAP, 'nodelay');
    [~,WinSeizOnset] = find(idxAllWin == SeizureOnset);
    
    perfAllWin(ichan, :) = 100*sum(indAllWin)/WIN;
    
    if TO_PLOT
        figure,
        h(1) = subplot(1,3,1);
        hold on;
        plot(SlBeforeSeizure(1,:), SlBeforeSeizure(2,:), 'k.');
        
        hold_state = get( axis_handle,'NextPlot' );
        set( axis_handle,'NextPlot','add' );
        plot( new_ver_line(1,:),new_ver_line(2,:),'r' );
        plot( new_horz_line(1,:),new_horz_line(2,:),'r' );
        plot( rotated_ellipse(1,:),rotated_ellipse(2,:),'r' );
        set( axis_handle,'NextPlot',hold_state );
        
        plot(SlBeforeSeizure(1,ind), SlBeforeSeizure(2,ind), 'k.');
        title('Before Seizure');
        axis square
        xlabel('SL 1'); 
        ylabel('SL2'); 
        set(gca, 'FontSize', 14); 
        
        h(2) = subplot(1,3,2);
        
        hold on;
        plot(SlPreSeizure(1,:), SlPreSeizure(2,:), 'k.');
        
        hold_state = get( axis_handle,'NextPlot' );
        set( axis_handle,'NextPlot','add' );
        plot( new_ver_line(1,:),new_ver_line(2,:),'r' );
        plot( new_horz_line(1,:),new_horz_line(2,:),'r' );
        plot( rotated_ellipse(1,:),rotated_ellipse(2,:),'r' );
        set( axis_handle,'NextPlot',hold_state );
        perfPre(ichan) = 100*sum(indPre)/size(DataPreSeizure,2);
        plot(SlPreSeizure(1,indPre), SlPreSeizure(2,indPre), 'c.');
        title(['Pre Seizure, ', 'outEllipse/All = ', num2str(perfPre(ichan))]);
        axis square
         xlabel('SL 1'); 
        ylabel('SL2'); 
        set(gca, 'FontSize', 14); 
        
        h(3) = subplot(1,3,3);
        
        hold on;
        plot(SlSeizure(1,:), SlSeizure(2,:), 'k.');
        
        hold_state = get( axis_handle,'NextPlot' );
        set( axis_handle,'NextPlot','add' );
        plot( new_ver_line(1,:),new_ver_line(2,:),'r' );
        plot( new_horz_line(1,:),new_horz_line(2,:),'r' );
        plot( rotated_ellipse(1,:),rotated_ellipse(2,:),'r' );
        set( axis_handle,'NextPlot',hold_state );
        
        perfSeiz(ichan) = 100*sum(indSeiz)/size(DataSeizure,2);
        plot(SlSeizure(1,indSeiz), SlSeizure(2,indSeiz), 'c.');
        title(['Seizure', 'outEllipse/All = ', num2str(perfSeiz(ichan))]);
        
        axis square
        xlabel('SL 1'); 
        ylabel('SL2'); 
        set(gca, 'FontSize', 14); 
        
        linkaxes(h, 'xy');
        
        set(gcf, 'Position', [100, 100, 1200, 400]);
        saveas(gcf, [FOLDER_FIGURES, 'SL', 'Ns', num2str(CONST_SUBGRAPH_SIZE), 'W', ...
            num2str(length(idxW)),num2str(CONST_SEIZURE_NODE), 'Ellipse.png']);
        
        
        figure,
        hax = axes;
        plot((1:length(indAll))/1000, indAll);
        hold on;
        line([600 600], get(hax,'YLim'), 'Color',[1 0 0], 'Linewidth', 1.5);
        title(['Out-Of-Ellipse EEG', 'S1S2', 'Ns', num2str(CONST_SUBGRAPH_SIZE), 'W', ...
            num2str(length(idxW))]);
        xlabel('Time [s]');
        set(gcf, 'Position', [100, 100, 1200, 400]);
        saveas(gcf, [FOLDER_FIGURES, 'S1S2', 'Ns', num2str(CONST_SUBGRAPH_SIZE), 'W', ...
            num2str(length(idxW)),num2str(CONST_SEIZURE_NODE), 'TimeCourse.png']);
    end
    
    figure, 
    imagesc(indAll)
 
        xlabel('Time [s]');
        map = [0 0 0; 0 1 1]; 
        colormap(map); 
        set(gca, 'FontSize', 14);
        set(gca, 'XTickLabel', {100 200 300 400 500 600 700})
    % @2017 06 06
    % plot the data in the slepian space with a datapoints colormapped to time
    
    % for all the points in the slepian space
    Sl = {SlAll, SlBeforeSeizure, SlPreSeizure, SlSeizure};
    SlName = {'All', 'Before', 'Pre', 'Seizure'};
    
    figure,
    for iSub = 1:4
        ax(iSub) = subplot(2,2,iSub)
        scatter(Sl{iSub}(1,:), Sl{iSub}(2,:), [], (1:size(Sl{iSub},2))/size(Sl{iSub},2));
        xlabel('Sl1');
        ylabel('Sl2');
        title(SlName{iSub});
        hold on;
        
        hold_state = get( axis_handle,'NextPlot' );
        set( axis_handle,'NextPlot','add' );
        plot( new_ver_line(1,:),new_ver_line(2,:),'r' );
        plot( new_horz_line(1,:),new_horz_line(2,:),'r' );
        plot( rotated_ellipse(1,:),rotated_ellipse(2,:),'r' );
        set( axis_handle,'NextPlot',hold_state );
        colorbar,
    end
    
    linkaxes(ax, 'xy');
    saveas(gcf, [FOLDER_FIGURES, 'S1S2', 'Ns', num2str(CONST_SUBGRAPH_SIZE), 'W', ...
        num2str(length(idxW)),num2str(CONST_SEIZURE_NODE), 'MapTimeToSl.png']);
    
    saveas(gcf, [FOLDER_FIGURES, 'S1S2', 'Ns', num2str(CONST_SUBGRAPH_SIZE), 'W', ...
        num2str(length(idxW)),num2str(CONST_SEIZURE_NODE), 'MapTimeToSl.fig']);
    
    % only for the point that are outside of the ellipse
    id = {indAll, indPre, indSeiz};
    idname = {'All', 'Pre', 'Seizure'};
    
    clear Sl ax
    Sl = {SlAll, SlPreSeizure, SlSeizure};
    figure,
    for iSub = 1:length(id)
        ax(iSub) = subplot(2,2,iSub)
        scatter(Sl{iSub}(1,id{iSub}), Sl{iSub}(2,id{iSub}), [], (1:sum(id{iSub}))/sum(id{iSub}));
        xlabel('Sl1');
        ylabel('Sl2');
        title(idname{iSub});
        hold on;
        
        hold_state = get( axis_handle,'NextPlot' );
        set( axis_handle,'NextPlot','add' );
        plot( new_ver_line(1,:),new_ver_line(2,:),'r' );
        plot( new_horz_line(1,:),new_horz_line(2,:),'r' );
        plot( rotated_ellipse(1,:),rotated_ellipse(2,:),'r' );
        set( axis_handle,'NextPlot',hold_state );
        colorbar,
    end
    linkaxes(ax, 'xy');
    saveas(gcf, [FOLDER_FIGURES, 'S1S2', 'Ns', num2str(CONST_SUBGRAPH_SIZE), 'W', ...
        num2str(length(idxW)),num2str(CONST_SEIZURE_NODE), 'MapTimeToSlOutOfEllipse.png']);
    
    
    % map the phase of the points from the slepian space to the time
    
    phaseSl1 = atan(SlAll(1,:)./SlAll(2,:));
    phaseSl2 = atan2(SlAll(2,:), SlAll(1,:));
    
    R = sqrt(SlAll(2,:).^2 + SlAll(1,:).^2);
    
    % plot what happens in time
    if TO_PLOT
        figure,
        imagesc(phaseSl1);
        title('Phase in the Slepian space [-pi/2 pi/2] ');
        saveas(gcf, [FOLDER_FIGURES, 'S1S2', 'Ns', num2str(CONST_SUBGRAPH_SIZE), 'W', ...
            num2str(length(idxW)),num2str(CONST_SEIZURE_NODE), 'Phase1.png']);
    end
    
    if TO_PLOT
        figure,
        imagesc(phaseSl2);
        title('Phase in the Slepian Space [-pi pi]');
        saveas(gcf, [FOLDER_FIGURES, 'S1S2', 'Ns', num2str(CONST_SUBGRAPH_SIZE), 'W', ...
            num2str(length(idxW)),num2str(CONST_SEIZURE_NODE), 'Phase2.png']);
        
    end
    
    
    
    % plot what happens in terms of phase and amplitude + time
    if TO_PLOT
        figure, 
        scatter(phaseSl2, R, [], 1:length(phaseSl2)); 
        xlabel('Phase in Slepian space'); 
        ylabel('Amplitude in the Slepian Space'); 
        
    end
    
    % see what happens in a temporal window 
    W = 30*ALLEEG.srate; 
    O = W/2; 
    I = buffer(1:size(SlAll, 2), W, O, 'nodelay'); 
    
    if TO_PLOT
        for iSub = 1 : (size(I,2)-1)
            h(iSub) = subplot(6, 8, iSub); 
           scatter(phaseSl2(I(:,iSub)), R(I(:,iSub)), [], 1:size(I,1)); 
           title(['Win = ', num2str(iSub)]);  
        end
        linkaxes(h, 'xy'); 
        saveas(gcf, [FOLDER_FIGURES, 'S1S2', 'Ns', num2str(CONST_SUBGRAPH_SIZE), 'W', ...
            num2str(length(idxW)),num2str(CONST_SEIZURE_NODE), 'PhaseAmplTimeWin.png']);
    end
    
    
    % @2017 06 06
    % Spatial PCA 
    
    [coeff, score] = pca(DataAll'); 
    
   % plot the 2 principle maps 
   figure, 
   subplot(1,2,1); 
   topoplot(coeff(:,1), chanlocs, 'electrodes', 'on', 'headrad',0); 
   title('PC1'); 
   subplot(1,2,2); 
   topoplot(coeff(:,2), chanlocs, 'electrodes', 'on', 'headrad',0); 
   title('PC2');
   
   % plot the projection on the data into that together with the time 
   
   figure, 
   
   ax(1) = subplot(1,3,1); 
   plot(score(1:300000,1), score(1:300000,2), 'ob'); 
   xlabel('PC1'); 
   ylabel('PC2'); 
   title('Before Seizure: PC space'); 
   
   ax(2) = subplot(1,3,2); 
   plot(score(480000:600000,1), score(480000:600000,2), 'og'); 
   xlabel('PC1'); 
   ylabel('PC2'); 
   title('pre Seizure: PC space'); 
   
   ax(3) = subplot(1,3,3); 
   plot(score(600000:end,1), score(600000:end,2), 'or'); 
   xlabel('seizure: PC space '); 
   ylabel('PC2');
   
   linkaxes(ax, 'xy'); 
   set(gcf, 'Position', [100, 100, 1200, 400]);
   
    %eegplot(DataAll, 'srate', EEGSamplingRate,'winlength', 60, 'dispchans', 50);
    ichan = ichan+1;
    close all
end
%
% figure,
% subplot(1,2,1)
% topoplot(perfPre, chanlocs,  'electrodes', 'on', 'headrad', 0);
% title('3 min Preseizure: fraction of timepoints out of ellipse');
%
% subplot(1,2,2)
% topoplot(perfSeiz, chanlocs,  'electrodes', 'on', 'headrad', 0);
% title('2 min seizure: fraction of timepoints out of ellipse');

n = 1;
figure,
for iWin = 1 : size(perfPreWin,2)
    subplot(3, size(perfPreWin,2)/3, n);
    
    topoplot(perfPreWin(:,iWin), chanlocs,  'electrodes', 'on', 'headrad', 0);
    title(['Win = ', num2str(n)]);
    colorbar
    n = n +1;
end

n = 1;
figure,
for iWin = 1 : size(perfSeizWin,2)
    subplot(2, size(perfSeizWin,2)/2, n);
    
    topoplot(perfSeizWin(:,iWin), chanlocs,  'electrodes', 'on', 'headrad', 0);
    title(['Win = ', num2str(n)]);
    colorbar
    n = n +1;
end

n = 1;
figure,
for iWin = 1 : size(perfAllWin,2)
    subplot(6, round(size(perfAllWin,2)/6), n);
    
    topoplot(perfAllWin(:,iWin), chanlocs,  'electrodes', 'on', 'headrad', 0);
    title(['Win = ', num2str(n)]);
    colorbar
    n = n + 1;
end



