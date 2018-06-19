% in this script I am testing what are the values necessary for correct
% plotting of the topoplot
clc
clear
close all
%load('/Users/annavybornova/EPFL/3eme_annee/Bachelor project/Matlab_code/EEG_channels.mat')

ChanLocs = readlocs('/Users/annavybornova/EPFL/Master_4/GEM/01Data/EegCap/Standard-10-20-Cap19.ced');

Rmanual = [0.5111
0.5111
0.5111
0.5111
0.4
0.2556
0.2556
0.4
0.4
0.2556
0.2556
0.4
0.5111
0.5111
0.5111
0.5111
0.12
0.12
0.5111
0.5111]; 

Map = rand(1,64);
Map(1:30) = 1;
figure,
topoplot(Map, ChanLocs)

%% STEP1: understanding of which fields are important and which ones are not
% 1) the conclusion here is that the important fields are: 'theta', 'radius', 'X',
% 'Y', 'Z'
% 2) Understand how theta and radius are calculated
% 3) note that 'theta' is the angle in the X Y plane and not as it is
% defined in wikipedia
for i = 1: length(ChanLocs)
    r(i) = sqrt(ChanLocs(i).X.^2 + ChanLocs(i).Y.^2); % + ChanLocs(i).Z.^2);
    
    %r(i) = ChanLocs(i).X/sin(pi*ChanLocs(i).theta/180 - pi/2)
    
    X(i) = ChanLocs(i).radius*cos(ChanLocs(i).theta);
    Y(i) = ChanLocs(i).radius*sin(ChanLocs(i).theta);
    
    Theta(i) = -180*atan(ChanLocs(i).Y/ChanLocs(i).X)/pi;
   
    if ChanLocs(i).X < 0 && ChanLocs(i).Y > 0
        Theta(i) = Theta(i) - 180;
    elseif ChanLocs(i).X < 0 && ChanLocs(i).Y < 0
        Theta(i) = Theta(i) + 180;
    end
end


% fields = {'sph_theta', 'sph_phi', 'sph_radius'};
% for iF = 1:length(fields)
%     FieldToRemove = fields{iF};
%     TChanLocs = rmfield(ChanLocs,FieldToRemove);
%     figure,
%     topoplot(Map,TChanLocs);
% end

% see what happens to theta
figure,
plot(Theta);
hold on;
plot([ChanLocs(:).theta], 'r');

% see what happens to radius
figure,
%subplot(2,1,1)
plot(r, 'ro');
hold on;
%subplot(2,1,2)
plot([ChanLocs(:).radius], 'bo');
set(gca, 'XTickLabel', {ChanLocs.labels});
set(gca, 'XTickLabelRotation', 45);
legend('r', 'Correct r')

% see what happens to X and Y

figure,
subplot(2,1,1)
plot([ChanLocs(:).X]);
hold on;
plot(X, 'r');
legend('Correct X', 'Calculated X');
subplot(2,1,2)
plot([ChanLocs(:).Y]);
hold on;
plot(Y, 'r');
legend('Correct Y', 'Calculated Y');

% Conclusion here for the moment is that X is not rcos(theta) and Y is not
% rsin(theta)


%% STEP 2 investigation of finding the midpoint between the electrodes
BipolarElectrodes = {'Fp1-F7';'F7-T7';'T7-P7';'P7-O1';'Fp1-F3';...
    'F3-C3';'C3-P3';'P3-O1';'Fp2-F4';'F4-C4';'C4-P4';'P4-O2';...
    'Fp2-F8';'F8-T8';'T8-P8';'P8-O2';'Fz-Cz';'Cz-Pz';'P7-T7';'T8-P8'}

BipolarCoord = struct('X', [], 'Y', [], 'Z', []);
for i = 1: length(BipolarElectrodes)
    % find the individual coordinates of each electrode
    
    for ichar = 1:length(BipolarElectrodes{i})
        name = BipolarElectrodes{i};
        if strcmp(name(ichar), '-')
            SepInd = ichar;
        end
    end
    
    Electrode1.Name = name(1:(SepInd-1));
    Electrode2.Name = name((SepInd+1):end);
    
    for iChan = 1:length({ChanLocs(:).labels})
        if strcmp(Electrode1.Name, ChanLocs(iChan).labels)
            Electrode1.X = ChanLocs(iChan).X;
            Electrode1.Y = ChanLocs(iChan).Y;
            Electrode1.Z = ChanLocs(iChan).Z;
            Electrode1.radius = ChanLocs(iChan).radius;
            Electrode1.theta = ChanLocs(iChan).theta; 
        end
        
        if strcmp(Electrode2.Name, ChanLocs(iChan).labels)
            Electrode2.X = ChanLocs(iChan).X;
            Electrode2.Y = ChanLocs(iChan).Y;
            Electrode2.Z = ChanLocs(iChan).Z;
            Electrode2.radius = ChanLocs(iChan).radius;
            Electrode2.theta = ChanLocs(iChan).theta; 
        end
        
    end
    BipolarCoord(i).labels = name;
    BipolarCoord(i).X  =  (Electrode1.X + Electrode2.X)/2;
    BipolarCoord(i).Y  =  (Electrode1.Y + Electrode2.Y)/2;
    BipolarCoord(i).Z  =  (Electrode1.Z + Electrode2.Z)/2;
    BipolarCoord(i).radius = sqrt(Electrode1.X.^2 + Electrode2.Y.^2);
    BipolarCoord(i).theta = -180*atan(BipolarCoord(i).Y/BipolarCoord(i).X)/pi;
    
    if BipolarCoord(i).X < 0 &&  BipolarCoord(i).Y  > 0
        BipolarCoord(i).theta = BipolarCoord(i).theta - 180;
    elseif BipolarCoord(i).X < 0 && BipolarCoord(i).Y < 0
        BipolarCoord(i).theta = BipolarCoord(i).theta + 180;
    end
end

% plot the cap with the new and the old coordinates 
 figure, 
 
 for iChan = 1:length({ChanLocs(:).labels})
     plot3(ChanLocs(iChan).X,ChanLocs(iChan).Y, ChanLocs(iChan).Z, 'ok');
     hold on; 
 end 
 hold on; 
 for iBipolar = 1:length([BipolarCoord(:).X])
    plot3(BipolarCoord(iBipolar).X, BipolarCoord(iBipolar).Y,...
        BipolarCoord(iBipolar).Z, 'or'); 
    hold on; 
 end

 figure, 

 % plot the cap with the new and the old coordinates 
 figure, 
 
 for iChan = 1:length({ChanLocs(:).labels})
     plot(ChanLocs(iChan).X,ChanLocs(iChan).Y,  'ok');
     hold on; 
 end 
 hold on; 
 for iBipolar = 1:length([BipolarCoord(:).X])
    plot(BipolarCoord(iBipolar).X, BipolarCoord(iBipolar).Y, 'or'); 
    hold on; 
 end
 
figure,
 for iChan = 1:length({ChanLocs(:).labels})
     
     subplot(5,5,iChan)
     plot(ChanLocs(iChan).X,ChanLocs(iChan).Y, 'ok');
     hold on; 
     plot([0,ChanLocs(iChan).X],[0,ChanLocs(iChan).Y], '-k')
     hold on; 
     plot([0,sin(pi*ChanLocs(iChan).theta/180)], [0, cos(pi*ChanLocs(iChan).theta/180)], '-r');
     xlim([-1,1]); 
     ylim = ([-1,1]);
     axis([-1, 1, -1, 1])
     hold on; 
     axis square
     title(ChanLocs(iChan).labels)
 end
 
 figure, 
 for iBipolar = 1:length([BipolarCoord(:).X])
    subplot(5,5,iBipolar)
    plot(BipolarCoord(iBipolar).X, BipolarCoord(iBipolar).Y, 'ok'); 
    hold on; 
    plot([0, BipolarCoord(iBipolar).X], [0, BipolarCoord(iBipolar).Y],'-k' )
    hold on; plot([0,sin(pi*BipolarCoord(iBipolar).theta/180)],...
        [0, cos(pi*BipolarCoord(iBipolar).theta/180)], '-r');
    title(BipolarElectrodes{iBipolar});
    axis([-1, 1, -1, 1])
    axis square
    
 end
%% investigation of ho the topoplot is influenced with the correct radius
% and the recalculated radius
a = rand(1,19)
figure, 
subplot(2,1,1)
topoplot(a, ChanLocs);

for i = 1:length(ChanLocs)
    ChanLocs(i).radius = r(i)/2;
end
subplot(2,1,2)
topoplot(a, ChanLocs);

for i = 1:length(BipolarCoord)
    BipolarCoord(i).radius = Rmanual(i); 
end

figure, 
topoplot(rand(1,20), BipolarCoord, 'electrodes','labels'); 





 
 






