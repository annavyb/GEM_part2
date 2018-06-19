% 2017 05 28 
% The goal of this script is to see whether the 2 most dicriminant GFT
% allow to discriminate the seizure/preseizure/non-seizure temporal
% dynamics like the Slepian does. 

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

% load the workspace for the pre-processed data from JUNGO
%load('/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_JUNGO_FOR_SLEPIAN/JUNGOpreprocessed.mat')
load('/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_OGUEY_FOR_SLEPIAN/OGUEYpreprocessed.mat')
%% Compute the GFT 


disp('... Adjacency Matrix generation');
AL = adj_generator([[chanlocs.X]', [chanlocs.Y]', [chanlocs.Z]'], Params.Gft.R);
nNeighbors = sum(AL);
avgNeighbors = mean(nNeighbors);

disp(['R = ', num2str(Params.Gft.R)]);
disp(['avg number of neighbors = ', num2str(avgNeighbors)]);

GftBeforeSeizure = gft(DataBeforeSeizure, AL);
GftSeizure = gft(DataSeizure, AL);
GftPreSeizure = gft(DataPreSeizure, AL);
GftAll = gft(DataAll, AL);
% 
% % DATA EPOCHING AND LABELING
% 
% Params.Epochs.WinEpoch = 2*EEGSamplingRate;
% Params.Epochs.Overlap = Params.Epochs.WinEpoch/2;
% 
% Params.Epochs.Label.Seizure = 1;
% Params.Epochs.Label.BeforeSeizure = 0;
% Params.Epochs.Label.PreSeizure = 2;
% 
% disp(['Epoch length: ', num2str(Params.Epochs.WinEpoch/EEGSamplingRate), ' [s]' ])
% disp(['Epoch overlap: ', num2str(Params.Epochs.Overlap/EEGSamplingRate), ' [s]' ])
% % generating the data frame
% 
% % EEG
% 
% EpochSeizure.Eeg = epochIndecesGeneration(DataSeizure, Params.Epochs.WinEpoch,...
%     Params.Epochs.Overlap, Params.Epochs.Label.Seizure );
% 
% EpochBeforeSeizure.Eeg = epochIndecesGeneration(DataBeforeSeizure, ...
%     Params.Epochs.WinEpoch,...
%     Params.Epochs.Overlap, Params.Epochs.Label.BeforeSeizure );
% 
% EpochPreSeizure.Eeg = epochIndecesGeneration(DataPreSeizure, ...
%     Params.Epochs.WinEpoch,...
%     Params.Epochs.Overlap, Params.Epochs.Label.PreSeizure );
% 
% % GFT
% EpochSeizure.Gft = epochIndecesGeneration(GftSeizure, Params.Epochs.WinEpoch,...
%     Params.Epochs.Overlap, Params.Epochs.Label.Seizure );
% 
% EpochBeforeSeizure.Gft = epochIndecesGeneration(GftBeforeSeizure, ...
%     Params.Epochs.WinEpoch,...
%     Params.Epochs.Overlap, Params.Epochs.Label.BeforeSeizure );
% 
% EpochPreSeizure.Gft = epochIndecesGeneration(GftPreSeizure, ...
%     Params.Epochs.WinEpoch,...
%     Params.Epochs.Overlap, Params.Epochs.Label.PreSeizure );
% Epoch.Gft = [EpochSeizure.Gft EpochBeforeSeizure.Gft EpochPreSeizure.Gft];
% Epoch.Eeg = [EpochSeizure.Eeg EpochBeforeSeizure.Eeg EpochPreSeizure.Eeg];
% 
% 
%   
%     for iEp = 1:length(Epoch.Gft)
%         Epoch.Gft(iEp).Mean = mean(Epoch.Gft(iEp).Signal)';
%         Epoch.Eeg(iEp).Mean = mean(Epoch.Eeg(iEp).Signal)';
%         
%         Epoch.Gft(iEp).MeanAbs = mean(abs(Epoch.Gft(iEp).Signal))';
%         Epoch.Eeg(iEp).MeanAbs = mean(abs(Epoch.Eeg(iEp).Signal))';
%         
%         Epoch.Gft(iEp).Variance = std(Epoch.Gft(iEp).Signal)';
%         Epoch.Eeg(iEp).Variance = std(Epoch.Eeg(iEp).Signal)';
%     end
%     
%     LabelsTot = [Epoch.Eeg.Label];
%     
%     MeanGft = [Epoch.Gft(:).Mean]';
%     MeanEeg = [Epoch.Eeg(:).Mean]';
%     MeanAbsGft =[Epoch.Gft(:).MeanAbs]';
%     MeanAbsEeg = [Epoch.Eeg(:).MeanAbs]';
%     VarianceGft = [Epoch.Gft(:).Variance]';
%     VarianceEeg = [Epoch.Eeg(:).Variance]';
%     
% RANK FEATURES

% % According to the variance 
%     Ind1 = find(LabelsTot == 0);%non-seizure indeces
%     Ind2 = find(LabelsTot == 1);%seizure indeces
%     L = LabelsTot([Ind1 Ind2]);%corresponding labels
%     [IDX(:,1), Z(:,1)] = rankfeatures(  VarianceGft([Ind1, Ind2],:)',...
%                     L, 'Criterion', 'wilcoxon');
%     [IDX(:,2), Z(:,2)] = rankfeatures(  VarianceGft([Ind1, Ind2],:)',...
%                     L, 'Criterion', 'bhattacharyya');
% According to the coefficients 
%GFTcoeff = [GftBeforeSeizure GftSeizure]; 
GFTcoeff = [DataBeforeSeizure DataSeizure]; 
labels = [zeros(1,size(GftBeforeSeizure,2)), ones(1, size(GftSeizure,2))]; 
[IDXcoef(:,1), Zcoef(:,1)] = rankfeatures(  GFTcoeff,...
                    labels, 'Criterion', 'wilcoxon'); 
 
[IDXcoef(:,2), Zcoef(:,2)] = rankfeatures(  GFTcoeff,...
                    labels, 'Criterion', 'bhattacharyya');       

% According to the EEG coefficients 

                
 % Visualization 
%  figure, 
%  ax(1) = subplot(1,3,1); 
%  plot(GftBeforeSeizure(IDX(1,1), :), GftBeforeSeizure(IDX(2,1), :), 'bo'); 
%  title('wilcoxon before seizure'); 
%  xlabel(['GFT', num2str(IDX(1,1))]); 
%  ylabel(['GFT', num2str(IDX(2,1))]); 
%   axis square
%  ax(2) = subplot(1,3,2); 
%   plot(GftPreSeizure(IDX(1,1), :), GftPreSeizure(IDX(2,1), :), 'go');
%    title('wilcoxon pre seizure'); 
%  xlabel(['GFT', num2str(IDX(1,1))]); 
%  ylabel(['GFT', num2str(IDX(2,1))]); 
%   axis square
%  ax(3) = subplot(1,3,3); 
%   plot(GftSeizure(IDX(1,1), :), GftSeizure(IDX(2,1), :), 'ro');
%   title('wilcoxon seizure'); 
%  xlabel(['GFT', num2str(IDX(1,1))]); 
%  ylabel(['GFT', num2str(IDX(2,1))]); 
%   axis square
%  linkaxes(ax, 'xy')
 
  figure, 
 ax(1) = subplot(1,3,1); 
 plot(GftBeforeSeizure(IDXcoef(1,1), :), GftBeforeSeizure(IDXcoef(2,1), :), 'bo'); 
 title('wilcoxon before seizure'); 
 xlabel(['GFT', num2str(IDXcoef(1,1))]); 
 ylabel(['GFT', num2str(IDXcoef(2,1))]); 
  axis square
 ax(2) = subplot(1,3,2); 
  plot(GftPreSeizure(IDXcoef(1,1), :), GftPreSeizure(IDXcoef(2,1), :), 'go');
   title('wilcoxon pre seizure'); 
 xlabel(['GFT', num2str(IDXcoef(1,1))]); 
 ylabel(['GFT', num2str(IDXcoef(2,1))]); 
  axis square
 ax(3) = subplot(1,3,3); 
  plot(GftSeizure(IDXcoef(1,1), :), GftSeizure(IDXcoef(2,1), :), 'ro');
  title('wilcoxon seizure'); 
 xlabel(['GFT', num2str(IDXcoef(1,1))]); 
 ylabel(['GFT', num2str(IDXcoef(2,1))]); 
  axis square
 linkaxes(ax, 'xy')
 
 % Evaluate the time course 
 % draw the ellipse

axis_handle = gca;
ellipse_t = fit_ellipse( GftBeforeSeizure(IDXcoef(1,1),:),GftBeforeSeizure(IDXcoef(2,1),:),axis_handle )

 % What are the points outside the ellipse for the seizure and pre-seizure
% data 

el = @(x,y,x0, y0, a, b, alpha) (((x-x0)*cos(-alpha) + (y-y0)*sin(-alpha)).^2/(a^2))...
    + (((x-x0)*sin(-alpha) - (y-y0)*cos(-alpha)).^2/(b^2)); 

coeff = 1; 
step = 0.1;
s = size(GftBeforeSeizure,2); 

while (s > 0.001*size(GftBeforeSeizure,2))
ind = (el(GftBeforeSeizure(IDXcoef(1,1),:), GftBeforeSeizure(IDXcoef(2,1),:), ellipse_t.X0, ...
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

indPre = (el(GftPreSeizure(IDXcoef(1,1),:), GftPreSeizure(IDXcoef(2,1),:), ellipse_t.X0, ...
    ellipse_t.Y0, ellipse_t.a, ellipse_t.b, ellipse_t.phi )> 1); 

indSeiz = (el(GftSeizure(IDXcoef(1,1),:), GftSeizure(IDXcoef(2,1),:), ellipse_t.X0, ...
    ellipse_t.Y0, ellipse_t.a, ellipse_t.b, ellipse_t.phi )> 1); 

indAll = (el(GftAll(IDXcoef(1,1),:), GftAll(IDXcoef(2,1),:), ellipse_t.X0, ...
    ellipse_t.Y0, ellipse_t.a, ellipse_t.b, ellipse_t.phi )> 1); 

ichan = 1; 
figure, 
h(1) = subplot(1,3,1);
hold on; 
plot(GftBeforeSeizure(IDXcoef(1,1),:), GftBeforeSeizure(IDXcoef(2,1),:), 'bo'); 

hold_state = get( axis_handle,'NextPlot' );
set( axis_handle,'NextPlot','add' );
plot( new_ver_line(1,:),new_ver_line(2,:),'r' );
plot( new_horz_line(1,:),new_horz_line(2,:),'r' );
plot( rotated_ellipse(1,:),rotated_ellipse(2,:),'r' );
set( axis_handle,'NextPlot',hold_state );

plot(GftBeforeSeizure(IDXcoef(1,1),ind), GftBeforeSeizure(IDXcoef(2,1),ind), 'ro'); 
title('Before Seizure'); 
axis square 

h(2) = subplot(1,3,2); 

hold on; 
plot(GftPreSeizure(IDXcoef(1,1),:), GftPreSeizure(IDXcoef(2,1),:), 'bo'); 

hold_state = get( axis_handle,'NextPlot' );
set( axis_handle,'NextPlot','add' );
plot( new_ver_line(1,:),new_ver_line(2,:),'r' );
plot( new_horz_line(1,:),new_horz_line(2,:),'r' );
plot( rotated_ellipse(1,:),rotated_ellipse(2,:),'r' );
set( axis_handle,'NextPlot',hold_state );
perfPre(ichan) = 100*sum(indPre)/size(DataPreSeizure,2); 

plot(GftPreSeizure(IDXcoef(1,1),indPre), GftPreSeizure(IDXcoef(2,1),indPre), 'ro'); 
title(['Pre Seizure, ', 'outEllipse/All = ', num2str(perfPre(ichan))]); 
axis square 


h(3) = subplot(1,3,3); 

hold on; 
plot(GftSeizure(IDXcoef(1,1),:), GftSeizure(IDXcoef(2,1),:), 'bo'); 

hold_state = get( axis_handle,'NextPlot' );
set( axis_handle,'NextPlot','add' );
plot( new_ver_line(1,:),new_ver_line(2,:),'r' );
plot( new_horz_line(1,:),new_horz_line(2,:),'r' );
plot( rotated_ellipse(1,:),rotated_ellipse(2,:),'r' );
set( axis_handle,'NextPlot',hold_state );

perfSeiz(ichan) = 100*sum(indSeiz)/size(DataSeizure,2); 
plot(GftSeizure(IDXcoef(1,1),indSeiz), GftSeizure(IDXcoef(2,1),indSeiz), 'ro'); 
title(['Seizure', 'outEllipse/All = ', num2str(perfSeiz(ichan))]); 

axis square 

linkaxes(h, 'xy');


figure, 
hax = axes; 
plot((1:length(indAll))/1000, indAll); 
%imagesc(indAll); 
hold on; 
line([600 600], get(hax,'YLim'), 'Color',[1 0 0], 'Linewidth', 1.5); 
title(['Out-Of-Ellipse EEG']); 
xlabel('Time [s]'); 
set(gcf, 'Position', [100, 100, 1200, 400]);
                
figure, 
load('V.mat'); 

subplot(1,2,1); 
    topoplot(V(:,75), chanlocs, 'electrodes', 'on', 'headrad', 0); 
    title('V2')
    subplot(1,2,2); 
    topoplot(V(:,52), chanlocs,  'electrodes', 'on', 'headrad', 0); 
    title('V1'); 
    
    




