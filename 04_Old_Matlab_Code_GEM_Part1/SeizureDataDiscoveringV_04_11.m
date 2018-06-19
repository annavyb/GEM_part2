% 2017 04 11
% In this script I compare the performance of various classifiers (start
% with linear) on the pre-seizure and seizure epochs

% STEP 1: load data
%
% STEP 2: data pre-processing:
%   i) Filtering between 1hz and 30 hz
%   ii) Plot the data using eeglab eegplot
%   iii) Check whether the seiture marker seems to be correct
%   iv) Identify the bad channels
%   v) Interpolate the bad channels
%
% STEP 3: Data epoching & and labeling
%   i) Divide the 5 first minutes of the recording into the epochs of 5s,
%       overlap 2.5 s
%   ii) Divide the first 20s after the seizure onset marker into the epochs
%   of 5s with 2.5s overlap
%   iii) Do not forget to provide the label for the epoch, i.e. 0 - no
%   seizure, 1 - seizure
%
% STEP 4: Obtain the GFT of the data
%
% STEP 5: Feature computation
% compute the temporal mean of each EEG channel/ GFT
% compute the variance of each EEG chanel/GFT
%

clc;
clear all;
close all;

dbstop if error;
set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaulttextinterpreter','none')

addpath('/Users/annavybornova/EPFL/3eme_annee/Bachelor project/Matlab_code')

%% STEP 0: Initialize the parameters

Params.FeatureRanking.Method = 'RANK_SUM';% Can take values :
% 1) 'RANK_SUM';
% 2) 'FISHER';
% 3) 'TTEST'; 
Params.FeatureRanking.Values = 'MEAN'; % Can take values :
% 1) 'ALL_CHANNELS';
% 2) 'MEAN'
% 3) 'VARIANCE'
% 4) 'MEAN_AND_VARIANCE'
Params.BadChannels = 'INTERPOLATE'; % Can take the values :
% 1) 'REMOVE'
% 2) 'INTERPOLATE'

disp('PARAMETERS:');
disp(['T-test method:',Params.FeatureRanking.Method]);
disp(['T-test computed on: ', Params.FeatureRanking.Values]);
disp(['Bad Channels handling: ', Params.BadChannels ]);



%% STEP 1: load data

FileName = {'JUNGO NATHALIE 0617 2019_seg_3.mat', 'OGUEY 0805 2113_seg_4.mat'};

iFile = 1; % the file of interest
FigNameSave = ['/Users/annavybornova/EPFL/Master_4/GEM/05SavedData/20170426SeizureDataJUNGO_OGUEY/', ...
    FileName{iFile}(1:3),'_', Params.BadChannels(1:3),'_',Params.FeatureRanking.Method(1:4),...
    '_', Params.FeatureRanking.Values(1:3)];

disp(['... Now processing ', FileName{iFile}]);

PathToFile = ['/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_Seizure_Epochs_V04_05/',...
    FileName{iFile}];
load(PathToFile);

% Define the patient-specific information

switch iFile
    case 1
        SeizureOnset = evt_Marks{2,2};
        DataAll = JUNGONATHALIE06172019_seg_3mff;
        %BadChannels = [257, 240, 238, 201, 178, 119, 37];
    case 2
        DataAll = OGUEY08052113_seg_4mff;
        SeizureOnset = evt_Events{2,2};
end

% load the electrode location file for 257 electrodes
chanlocs = ...
    readlocs('/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_Seizure_Epochs_V04_05/EGI 257.Geneva Average 13.10-10.xyz');

% load the electrode location file for 204 electrodes
chanlocs204 = ...
    readlocs('/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_Seizure_Epochs_V04_05/EGI 204.Geneva Average 13.10-10.xyz')

% find out what channels were removed when we pass from 257 configuration
% to the 204 configuration

ii = 1;
for iCh257 = 1:length(chanlocs)
    for iCh204 = 1:length(chanlocs204)
        if strcmp(chanlocs(iCh257).labels, chanlocs204(iCh204).labels)
            TrueIndCh257Ch204(ii) = iCh257; % the corresponding indeces
            ii = ii+1;
        end
    end
end


fid = fopen('EGI 257.Geneva Average 13.10-10.xyz');
%header = fscanf(fid,'%f',[2,1]);

for len = 1:257
    position(len,:) = fscanf(fid,'%f',[3])';
    channame{len} = fscanf(fid,'%s',[1]);
end

fclose(fid);

channame{1} = '1';

% data should be in the same order that the chanlocs
ind = 1;
for iChan = 1:size(DataAll,1)
    for jChan = 1:size(DataAll,1)
        if strcmp(channame{jChan},chanlocs(iChan).labels)
            Ind(ind).channame = jChan;
            Ind(ind).chanlocs = iChan;
            ind = ind+1;
        end
    end
end

% put the data in the same order as are the chanlocs
DataAll = DataAll([Ind.channame], :);

% only keep the data corresponding to the 204 chan configuration



DataAll = DataAll(TrueIndCh257Ch204, :);
chanlocs = chanlocs(TrueIndCh257Ch204);

figure, 
subplot(1,2,1)
topoplot(ones(1,204), chanlocs, 'electrodes', 'on', ...
        'style', 'both', 'headrad', 0, 'emarker', {'.','k',[],3} );
chanlocs257 = readlocs('/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_Seizure_Epochs_V04_05/EGI 257.Geneva Average 13.10-10.xyz'); 
subplot(1,2,2)
 topoplot(ones(1,257), chanlocs257, 'electrodes', 'on', ...
        'style', 'both', 'headrad', 0), 'emarker', {'.','k',[],3};

% what are the bad channels that are identified after the visual inspection
% with eegplot (flat/noisy channels)


switch iFile
    case 1
        BadChannels = [64 104 108 182 183]; % eye channels "32" - 89;
        %                                               "AF8" 163
        % what do these channels corre
    case 2
        BadChannels = [64 86 92 101 108 183]; %
        % what do this indeces correspond to? --> 'FP1'    '45'    '9'    '81'     '144'    '178'
        
end

%   v) Bad channels handling
if strcmp(Params.BadChannels, 'REMOVE')
    DataAll(BadChannels, : ) = [];
    chanlocs(BadChannels) = [];
elseif  strcmp(Params.BadChannels, 'INTERPOLATE')
    load('/Users/annavybornova/EPFL/Master_4/GEM/01Data/DataStructure.mat');
    ALLEEG.data = DataAll;
    ALLEEG.chanlocs = chanlocs;
    EEGOUT = pop_interp(ALLEEG, BadChannels, 'spherical');
    DataAll = EEGOUT.data;
end

% --> sanity check is ok: works

%% STEP 2: data pre-processing:
%   i) Filtering between 1hz and 30 hz

% --> filter design

FiltOrder = 8;
LowCf = 1;
HighCf = 30;

h1=fdesign.highpass('N,F3dB',FiltOrder,LowCf,EEGSamplingRate);
d1=design(h1,'Butter');
h2=fdesign.lowpass('N,F3dB',FiltOrder,HighCf,EEGSamplingRate);
d2=design(h2,'Butter');

% Time limits for the epoching (needed to reduce the time of the filtering)
% 5min at the beginning of the recording
% 30 s after the seizure onset

DurationInterictal = 5*60; % in seconds
SamplesInteriactal = round(DurationInterictal*EEGSamplingRate);

DurationIctal = 30; % in seconds
SamplesIctal = DurationIctal*EEGSamplingRate;

DataBeforeSeizure = DataAll(:, 1:SamplesInteriactal);
DataSeizure = DataAll(:, SeizureOnset:(SeizureOnset+SamplesIctal));

%filter the data
for iChan=1:size(DataSeizure,1)
    DataBeforeSeizure(iChan, :) = filtfilt(d1.sosMatrix,d1.ScaleValues,...
        double(DataBeforeSeizure(iChan,:)));
    DataBeforeSeizure(iChan,:) = filtfilt(d2.sosMatrix,d2.ScaleValues,...
        double(DataBeforeSeizure(iChan,:)));
    
    DataSeizure(iChan, :) = filtfilt(d1.sosMatrix,d1.ScaleValues,...
        double(DataSeizure(iChan,:)));
    DataSeizure(iChan,:)=filtfilt(d2.sosMatrix,d2.ScaleValues,...
        double(DataSeizure(iChan,:)));
    
end

%   ii) Plot the data using eeglab eegplot
eegplot([DataBeforeSeizure, DataSeizure], 'srate', EEGSamplingRate, ...
    'winlength', 60, 'dispchans', 50);

% iii) Visualization of the seizure map
figure,
for iMap = 1:30
    subplot(5,6,iMap)
    topoplot(DataSeizure(:,500*iMap), chanlocs, 'electrodes', 'on', ...
        'style', 'both', 'headrad', 0);
end

set(gcf, 'Position', [100, 100, 1049, 895]);
saveas(gcf,[FigNameSave,'SeizureMapsEEG.png']);

%% STEP 3: Data epoching & and labeling
%   i) Divide the 5 first minutes of the recording into the epochs of 5s,
%       overlap 2.5 s

WinEpoch = 2*EEGSamplingRate;
Overlap = WinEpoch/2;

% generating the data frame

IndecesEpochBeforeSeizure = buffer([1:size(DataBeforeSeizure,2)], WinEpoch,...
    Overlap, 'nodelay');
IndecesEpochBeforeSeizure(2:(end-1),:) = []; % leave only the info about the
%                                        beginning & the end of the window

IndecesEpochSeizure =  buffer(1:size(DataSeizure,2), WinEpoch,...
    Overlap, 'nodelay');
IndecesEpochSeizure(:, end) = []; % because of the zeros
IndecesEpochSeizure(2:(end-1),:) = []; % leave only the info about the
%                                        beginning & the end of the window



%   ii) Divide the first 20s after the seizure onset marker into the epochs
%   of 5s with 2.5s overlap
%   iii) Do not forget to provide the label for the epoch, i.e. 0 - no
%   seizure, 1 - seizure

%% STEP 4: Obtain the GFT of the data

% (i) Generate the adjacency matrix based on the chanlocs
R = 80;

disp('... Adjacency Matrix generation');
A = adj_generator([[chanlocs.X]', [chanlocs.Y]', [chanlocs.Z]'], R);
nNeighbors = sum(A);
avgNeighbors = mean(nNeighbors);

disp(['R = ', num2str(R)]);
disp(['avg number of neighbors = ', num2str(avgNeighbors)]);

% (ii) Compute GFT

L = Laplacian(A);
%calculate eigen values of L
L_eig = eig(L);
L_eig_sorted = sort(L_eig);
[V, D] = eig(L);

GftDataBeforeSeizure = V'*DataBeforeSeizure;
GftDataSeizure = V'* DataSeizure;

if strcmp(Params.FeatureRanking.Values, 'ALL_CHANNELS')
    figure,
    subplot(2,1,1)
    imagesc((zscore(GftDataSeizure, 0, 2)))
    title('GFT Seizure');
    subplot(2,1,2)
    imagesc(zscore(DataSeizure, 0, 2))
    title('EEG Seizure');
    set(gcf, 'Position', [100, 100, 1049, 895]);
    saveas(gcf,[FigNameSave,'GFT_SEIZURE.png']);
    
    figure,
    subplot(2,1,1)
    imagesc((zscore(GftDataBeforeSeizure(:,30000:60000), 0, 2)));
    title('GFT before seizure ' );
    
    subplot(2,1,2)
    imagesc((zscore(DataBeforeSeizure(:,30000:60000), 0, 2)));
    title('EEG before seizure ' );
    set(gcf, 'Position', [100, 100, 1049, 895]);
    saveas(gcf,[FigNameSave,'GFT_BEFORE_SEIZURE.png']);   
end

% eegplot(GftDataBeforeSeizure(:,180000:210000),...
%     'srate', EEGSamplingRate,  'winlength', 60, 'dispchans', 50);
%
% eegplot(DataBeforeSeizure(:,180000:210000),...
%     'srate', EEGSamplingRate,  'winlength', 60, 'dispchans', 50);
%% STEP 5: Epochs generation for 2 classes & Feature computation

% (i) Generation for the seizures
for iE = 1:size(IndecesEpochSeizure, 2)
    
    EegSeizure(iE).Signal = DataSeizure(:, IndecesEpochSeizure(1,iE):IndecesEpochSeizure(2,iE));
    EegSeizure(iE).Mean = mean(EegSeizure(iE).Signal, 2);
    EegSeizure(iE).Variance = var(EegSeizure(iE).Signal,0,2);
    EegSeizure(iE).Label = 1;
    
    GftSeizure(iE).Signal = GftDataSeizure(:, IndecesEpochSeizure(1,iE):IndecesEpochSeizure(2,iE));
    GftSeizure(iE).Mean = mean(GftSeizure(iE).Signal, 2);
    GftSeizure(iE).Variance = var(GftSeizure(iE).Signal,0,2);
    GftSeizure(iE).Label = 1;
    
end

for iE = 1:size(IndecesEpochBeforeSeizure, 2)
    
    EegBeforeSeizure(iE).Signal = DataBeforeSeizure(:, ...
        IndecesEpochBeforeSeizure(1,iE):IndecesEpochBeforeSeizure(2,iE));
    
    EegBeforeSeizure(iE).Mean = mean(EegBeforeSeizure(iE).Signal, 2);
    EegBeforeSeizure(iE).Variance = var(EegBeforeSeizure(iE).Signal,0,2);
    EegBeforeSeizure(iE).Label = 0;
    
    GftBeforeSeizure(iE).Signal = GftDataBeforeSeizure(:, ...
        IndecesEpochBeforeSeizure(1,iE):IndecesEpochBeforeSeizure(2,iE));
    
    GftBeforeSeizure(iE).Mean = mean(GftBeforeSeizure(iE).Signal, 2);
    GftBeforeSeizure(iE).Variance = var(GftBeforeSeizure(iE).Signal,0,2);
    GftBeforeSeizure(iE).Label = 0;
    
end

% (ii) Visualization of the mean and variance for each epoch for
% seizure/non-seizure in both GFT and EEG data

% SEIZURE

if strcmp(Params.FeatureRanking.Values, 'MEAN')
    figure,
    subplot(2,1,1)
    imagesc([EegSeizure.Mean]);
    title('EEG seizure Mean', 'FontSize', 16);
    
    subplot(2,1,2)
    imagesc([GftSeizure.Mean]);
    title('GFT seizure Mean', 'FontSize', 16);
    
    set(gcf, 'Position', [100, 100, 1049, 895]);
    saveas(gcf,[FigNameSave,'SEIZURE_MEAN.png']);
    
    % NON-SEIZURE
    
    figure,
    subplot(2,1,1)
    imagesc([EegBeforeSeizure.Mean]);
    title('EEG before seizure Mean',  'FontSize', 16);
    
    subplot(2,1,2)
    imagesc([GftBeforeSeizure.Mean]);
    title('GFT before seizure Mean',  'FontSize', 16);
    
    set(gcf, 'Position', [100, 100, 1049, 895]);
    saveas(gcf,[FigNameSave,'BEFORE_SEIZURE_MEAN.png']);
    
    % ALL
    
    figure,
    subplot(2,1,1)
    imagesc([EegBeforeSeizure.Mean EegSeizure.Mean]);
    title('EEG Mean', 'FontSize', 16);
    
    subplot(2,1,2)
    imagesc([GftBeforeSeizure.Mean GftSeizure.Mean]);
    title('GFT Mean', 'FontSize', 16);
    
elseif strcmp(Params.FeatureRanking.Values, 'VARIANCE')
    
    % Seizure 
    
    figure,
    subplot(2,1,1)
    imagesc([EegSeizure.Variance]);
    title('EEG seizure Variance', 'FontSize', 16);
    
    subplot(2,1,2)
    imagesc([GftSeizure.Variance]);
    title('GFT seizure Variance',  'FontSize', 16);
    set(gcf, 'Position', [100, 100, 1049, 895]);
    saveas(gcf,[FigNameSave,'SEIZURE_VARIANCE.png']);
    
    
    % Non Seizure
    
    figure,
    subplot(2,1,1)
    imagesc([EegBeforeSeizure.Variance]);
    title('EEG before seizure variance',  'FontSize', 16);
    
    subplot(2,1,2)
    imagesc([GftBeforeSeizure.Variance]);
    title('GFT before seizure variance',  'FontSize', 16);
    
    set(gcf, 'Position', [100, 100, 1049, 895]);
    saveas(gcf,[FigNameSave,'BEFORE_SEIZURE_VARIANCE.png']);
    
end


%% STEP 5 bis A: Feature ranking according to t-test

% 1) Use All channels for ranking
if strcmp(Params.FeatureRanking.Values, 'ALL_CHANNELS')
    
   
        
        % for the GFT
        for iSeiz = 1:length(GftSeizure)
            for iNonSeiz = 1:length(GftBeforeSeizure)
                for iFreq = 1:size(GftSeizure(1).Signal, 1)
                     if strcmp(Params.FeatureRanking.Method, 'RANK_SUM')
                    [p{iFreq}(iSeiz, iNonSeiz), h{iFreq}(iSeiz, iNonSeiz)] = ...
                        ranksum(GftSeizure(iSeiz).Signal(iFreq, :),...
                        GftBeforeSeizure(iNonSeiz).Signal(iFreq, :));
                     elseif strcmp(Params.FeatureRanking.Method, 'FISHER')
                         
                         
                         
                     end
                end
            end
        end
        
        % the same for the EEG channels
        
        for iSeiz = 1:length(EegSeizure)
            for iNonSeiz = 1:length(EegBeforeSeizure)
                for iFreq = 1:size(EegSeizure(1).Signal, 1)
                    [pEeg{iFreq}(iSeiz, iNonSeiz), hEeg{iFreq}(iSeiz, iNonSeiz)] = ...
                        ranksum(EegSeizure(iSeiz).Signal(iFreq, :),...
                        EegBeforeSeizure(iNonSeiz).Signal(iFreq, :));
                end
            end
        end
    
    
    for iFreq = 1:length(p)
        MedianPvalue(iFreq) = median(median(p{iFreq}));
        MedianPvalueEeg(iFreq) = median(median(pEeg{iFreq}));
    end
    
    [SortedPvalue, IndPvalueSorted] = sort(MedianPvalue);
    [SortedPvalueEeg, IndPvalueSortedEeg] = sort(MedianPvalueEeg);
    
    % 2) Use Mean values for the ranking
    
    
elseif  strcmp(Params.FeatureRanking.Values, 'MEAN')
    
    
    Concat.Gft.Seizure = [GftSeizure(:).Mean]';
    Concat.Gft.BeforeSeizure = [GftBeforeSeizure(:).Mean]';
    Concat.Eeg.Seizure = [EegSeizure(:).Mean]';
    Concat.Eeg.BeforeSeizure = [EegBeforeSeizure(:).Mean]';
elseif strcmp(Params.FeatureRanking.Values, 'VARIANCE')
    Concat.Gft.Seizure = [GftSeizure(:).Variance]';
    Concat.Gft.BeforeSeizure = [GftBeforeSeizure(:).Variance]';
    Concat.Eeg.Seizure = [EegSeizure(:).Variance]';
    Concat.Eeg.BeforeSeizure = [EegBeforeSeizure(:).Variance]';
end

    % for the GFT
    for iFreq = 1:size(Concat.Gft.Seizure, 2)
        if strcmp(Params.FeatureRanking.Method, 'RANK_SUM')
        [p(iFreq), h(iFreq)] = ...
            ranksum(Concat.Gft.Seizure(:,iFreq),...
            Concat.Gft.BeforeSeizure(:,iFreq));
        
        [pEeg(iFreq), hEeg(iFreq)] = ...
            ranksum(Concat.Eeg.Seizure(:,iFreq),...
            Concat.Eeg.BeforeSeizure(:,iFreq));
        
        [IDX, Z] = rankfeatures([Concat.Gft.Seizure; Concat.Gft.BeforeSeizure]',...
            [GftSeizure(:).Label GftBeforeSeizure(:).Label], 'Criterion', 'wilcoxon'); 
        elseif strcmp(Params.FeatureRanking.Method, 'FISHER')
            [IDX, Z] = rankfeat([Concat.Gft.Seizure; Concat.Gft.BeforeSeizure],...
            [GftSeizure(:).Label GftBeforeSeizure(:).Label], 'fisher'); 
        end
        
    end


[SortedPvalue, IndPvalueSorted] = sort(p);
[SortedPvalueEeg, IndPvalueSortedEeg] = sort(pEeg);



%% STEP 5 bis B : feature ranking according to the p-value

figure,
subplot(2,1,1)
plot(SortedPvalue);
xlabel('sorted frequency');
ylabel('p-value');
title('P-value ranksum test for each spatial frequency')

subplot(2,1,2)
plot(IndPvalueSorted, 'o');
xlabel('sorted frequency');
ylabel('frequency Value');

saveas(gcf,[FigNameSave,'FEATURE_RANKING_GFT.png']);

figure,
subplot(2,1,1)
plot(SortedPvalueEeg);
xlabel('sorted channels');
ylabel('p-value');
title('Median p-value ranksum test for each EEG channel')

subplot(2,1,2)
plot(IndPvalueSortedEeg, 'o');
xlabel('sorted channels');
ylabel('Channel Index');

saveas(gcf,[FigNameSave,'FEATURE_RANKING_EEG.png']);

figure,
plot(SortedPvalue, 'b');
hold on;
plot(SortedPvalueEeg, 'r');
title('Sorted p-values');
ylabel('p-value')
legend('GFT', 'EEG')
saveas(gcf,[FigNameSave,'FEATURE_RANKING_GFT_EEG.png']);

% Conclusion from here : I could use in principle 20 first features for the
% GFT

% Plot first 14 spatial frequencies

figure,
for f = 1:20
    subplot(4,5,f)
    topoplot(V(:, IndPvalueSorted(f)), chanlocs, 'headrad', 0);
    title(['f = ', num2str(IndPvalueSorted(f))]);
end
set(gcf, 'Position', [100, 100, 1049, 895]);
saveas(gcf,[FigNameSave,'SIGNIFICANT_SP_FREQ.png']);
% Plot the significant channels

ch = [30, 60, 90, 120];
figure,
for i = ch
    subplot(2,2,i/30)
    topoplot(ones(1,length(chanlocs(IndPvalueSortedEeg(1:i)))) ,...
        chanlocs(IndPvalueSortedEeg(1:i)),'electrodes', ...
        'labels','headrad', 0);
    title(['# Channels = ', num2str(i)]);
end

set(gcf, 'Position', [100, 100, 1049, 895]);
saveas(gcf,[FigNameSave,'SIGNIFICANT_EEG_CHAN.png']);


%% STEP 6: train the classifier using only the significant features

% 1) Let's do it for a mean for the GFT
if strcmp(Params.FeatureRanking.Values, 'MEAN') || ...
        strcmp(Params.FeatureRanking.Values, 'ALL_CHANNELS');
    ConcatGft.Features = [[GftBeforeSeizure(:).Mean],[GftSeizure(:).Mean]]';
    ConcatGft.Features = ConcatGft.Features(:,IndPvalueSorted); % Sort according to the p-value
    ConcatGft.Labels = [[GftBeforeSeizure(:).Label],[GftSeizure(:).Label]];
    
    ConcatEeg.Features = [[EegBeforeSeizure(:).Mean],[EegSeizure(:).Mean]]';
    ConcatEeg.Features = ConcatEeg.Features(:,IndPvalueSortedEeg); % Sort according to the p-value
    ConcatEeg.Labels = [[EegBeforeSeizure(:).Label],[EegSeizure(:).Label]];
    
elseif strcmp(Params.FeatureRanking.Values, 'VARIANCE')
    
    ConcatGft.Features = [[GftBeforeSeizure(:).Variance],[GftSeizure(:).Variance]]';
    ConcatGft.Features = ConcatGft.Features(:,IndPvalueSorted); % Sort according to the p-value
    ConcatGft.Labels = [[GftBeforeSeizure(:).Label],[GftSeizure(:).Label]];
    
    ConcatEeg.Features = [[EegBeforeSeizure(:).Variance],[EegSeizure(:).Variance]]';
    ConcatEeg.Features = ConcatEeg.Features(:,IndPvalueSortedEeg); % Sort according to the p-value
    ConcatEeg.Labels = [[EegBeforeSeizure(:).Label],[EegSeizure(:).Label]];
    
end

classifier = {'diaglinear', 'diagquadratic'};
nFolds = 3;
partGft = cvpartition(ConcatGft.Labels, 'Kfold', nFolds);

for k = 1:2
    for j = 1:nFolds
        repartition(partGft)
        indTrain = training(partGft, j);
        indTest = test(partGft,j);
        
        Train.set = ConcatGft.Features(indTrain, :);
        Train.labels = ConcatGft.Labels(indTrain);
        Testing.set = ConcatGft.Features(indTest, :);
        Testing.labels = ConcatGft.Labels(indTest);
        
        
        TrainEeg.set = ConcatEeg.Features(indTrain, :);
        TrainEeg.labels = ConcatEeg.Labels(indTrain);
        TestingEeg.set = ConcatEeg.Features(indTest, :);
        TestingEeg.labels = ConcatEeg.Labels(indTest);
        
        for i = 1:length(chanlocs)
            [classGft, errTrainGft(j,i,k)] = classify(Testing.set(:,1:i), Train.set(:,1:i), Train.labels, classifier{k});
            [classEeg, errTrainEeg(j,i,k)] = classify(TestingEeg.set(:,1:i), TrainEeg.set(:,1:i), TrainEeg.labels, classifier{k});
            
            [classGftTrain, errTrainGft(j,i,k)] = classify(Train.set(:,1:i), Train.set(:,1:i), Train.labels, classifier{k});
            [classEegTrain, errTrainEeg(j,i,k)] = classify(TrainEeg.set(:,1:i), TrainEeg.set(:,1:i), TrainEeg.labels, classifier{k});
            
            errTestGft(j,i,k) = classerr(classGft, Testing.labels);
            errTestEeg(j,i,k) = classerr(classEeg, TestingEeg.labels);
            
            
            
            errClassTestGft{k}(j,i,:) = class_Spec_err(classGft, Testing.labels); 
            errClassTestEeg{k}(j,i,:) = class_Spec_err(classEeg, TestingEeg.labels);
            
            errClassTrainGft{k}(j,i,:) = class_Spec_err(classGftTrain, Train.labels); 
            errClassTrainEeg{k}(j,i,:) = class_Spec_err(classEegTrain, TrainEeg.labels);
        end
    end
end

for i = 1:2
    figure,
    subplot(2,1,1);
    errorbar(mean(100*errTrainGft(:,:,i)), std(100*errTrainGft(:,:,i)));
    title([classifier{i}, 'Gft Train Error']);
    xlabel('Sorted features');
    ylabel('Misclassification');
    
    subplot(2,1,2);
    errorbar(mean(errTestGft(:,:,i)), std(errTestGft(:,:,i)));
    title([classifier{i}, 'Gft Test Error']);
    xlabel('Sorted features');
    ylabel('Misclassification');
    
    set(gcf, 'Position', [100, 100, 1049, 895]);
    saveas(gcf,[FigNameSave,'GFT',classifier{i},'.png']);
    
    figure,
    subplot(2,1,1)
    errorbar(mean(100*errTrainEeg(:,:,i)), std(100*errTrainEeg(:,:,i)));
    title([classifier{i}, 'Eeg Train Error']);
    xlabel('Sorted features');
    ylabel('Misclassification');
    
    subplot(2,1,2);
    errorbar(mean(errTestEeg(:,:,i)), std(errTestEeg(:,:,i)));
    title([classifier{i}, 'Eeg Test Error']);
    xlabel('Sorted features');
    ylabel('Misclassification');
    
    set(gcf, 'Position', [100, 100, 1049, 895]);
    saveas(gcf,[FigNameSave,'EEG',classifier{i},'.png']);
end

for i = 1:2 
   
    
    figure, 
    
    subplot(2,2,1)
    errorbar(mean(errClassTrainGft{i}(:,:, 1)), std(errClassTrainGft{i}(:,:,1))); 
    title([classifier{i}, 'Class specific error GFT Train: Before Seizure']);
    xlabel('Sorted features');
    ylabel('Misclassification');
    
    subplot(2,2,2)
    
    errorbar(mean(errClassTrainGft{i}(:,:, 2)), std(errClassTrainGft{i}(:,:,2))); 
    title([classifier{i}, 'Class specific error GFT Train : Seizure']);
    xlabel('Sorted features');
    ylabel('Misclassification');
    
    
    subplot(2,2,3)
    errorbar(mean(errClassTestGft{i}(:,:, 1)), std(errClassTestGft{i}(:,:,1))); 
    title([classifier{i}, 'Class specific error GFT Test: Before Seizure']);
    xlabel('Sorted features');
    ylabel('Misclassification');
    
    subplot(2,2,4)
    
    errorbar(mean(errClassTestGft{i}(:, :, 1)), std(errClassTestGft{i}(:,:,2))); 
    title([classifier{i}, 'Class specific error GFT Test: Seizure']);
    xlabel('Sorted features');
    ylabel('Misclassification');
    
     
    set(gcf, 'Position', [100, 100, 1049, 895]);
    saveas(gcf,[FigNameSave,'GFTclassSpecific',classifier{i},'.png']);
    
end

for i = 1:2 
   
    
    figure, 
    
    subplot(2,2,1)
    errorbar(mean(errClassTrainEeg{i}(:,:, 1)), std(errClassTrainEeg{i}(:,:,1))); 
    title([classifier{i}, 'Class specific error EEG Train: Before Seizure']);
    xlabel('Sorted features');
    ylabel('Misclassification');
    
    subplot(2,2,2)
    
    errorbar(mean(errClassTrainEeg{i}(:,:, 2)), std(errClassTrainEeg{i}(:,:,2))); 
    title([classifier{i}, 'Class specific error EEG Train : Seizure']);
    xlabel('Sorted features');
    ylabel('Misclassification');
    
    
    subplot(2,2,3)
    errorbar(mean(errClassTestEeg{i}(:,:, 1)), std(errClassTestEeg{i}(:,:,1))); 
    title([classifier{i}, 'Class specific error EEG Test: Before Seizure']);
    xlabel('Sorted features');
    ylabel('Misclassification');
    
    subplot(2,2,4)
    
    errorbar(mean(errClassTestEeg{i}(:, :, 1)), std(errClassTestEeg{i}(:,:,2))); 
    title([classifier{i}, 'Class specific error EEG Test: Seizure']);
    xlabel('Sorted features');
    ylabel('Misclassification');
     set(gcf, 'Position', [100, 100, 1049, 895]);
    saveas(gcf,[FigNameSave,'EEGclassSpecific',classifier{i},'.png']);
    
end

clear Train Testing TestingEeg TrainEeg ConcatGft partGft errTestGft errTrainGft

% %% STEP 6 : classification
% for iF = 1:20
%
%     TrainingEeg(:,iF) = EegBeforeSeizure(iF).Mean;
%     LabelTraining(iF) = EegBeforeSeizure(iF).Label;
%
%     TrainingEeg(:,iF+20) = EegSeizure(iF).Mean;
%     LabelTraining(iF+20) = EegSeizure(iF).Label;
%
%     TrainingGft(:,iF) = GftBeforeSeizure(iF).Mean;
%     LabelGftTraining(iF) = GftBeforeSeizure(iF).Label;
%
%     TrainingGft(:,iF+20) = GftSeizure(iF).Mean;
%     LabelGftTraining(iF+20) = GftSeizure(iF).Label;
%
% end
%
% for iF = 1:9
%
%     TestEeg(:,iF) = EegBeforeSeizure(iF+40).Mean;
%     LabelTest(iF) = EegBeforeSeizure(iF+40).Label;
%
%     TestEeg(:,iF+9) = EegSeizure(iF+20).Mean;
%     LabelTest (iF+9) = EegSeizure(iF+20).Label;
%
%     TestGft(:,iF) = GftBeforeSeizure(iF+40).Mean;
%     LabelGftTest(iF) = GftBeforeSeizure(iF+40).Label;
%
%     TestGft(:,iF+9) = GftSeizure(iF+20).Mean;
%     LabelGftTest (iF+9) = GftSeizure(iF+20).Label;
%
% end
%
% %[classGft, errGft] = classify(TestGft', TrainingGft', LabelGftTraining', 'diaglinear');
% [classGft, errGft] = classify(TestGft', TrainingGft', LabelGftTraining', 'diagquadratic');
% [classEeg, errEeg] = classify(TestEeg', TrainingEeg', LabelTraining', 'diagquadratic');
%





















