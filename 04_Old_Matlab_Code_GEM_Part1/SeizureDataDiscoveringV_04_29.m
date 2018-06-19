
% 2017 04 29
% Here I refine the script SeizureDataDiscoveryV_04_11.m to make it more
% compact, more modulable and in this way reduce the probability of making
% errors

% Here will go the pipeline description

clc;
clear all;
close all;

dbstop if error;
set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaulttextinterpreter','none')

addpath('/Users/annavybornova/EPFL/3eme_annee/Bachelor project/Matlab_code');

% Pipeline parameters initialization
Params = setParam(1); % see the function setParam for further information

%% STEP 1: LOAD DATA AND PUT IT IN THE CORRECT ORDER

FileName = {'JUNGO NATHALIE 0617 2019_seg_3.mat', 'OGUEY 0805 2113_seg_4.mat'};

iFile = 2; % the file of interest
FigNameSave = ['/Users/annavybornova/EPFL/Master_4/GEM/05SavedData/20170502SeizureDataPresentation/', ...
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
        BadChannels = [64 104 108 182 183];
    case 2
        DataAll = OGUEY08052113_seg_4mff;
        SeizureOnset = evt_Events{2,2};
        BadChannels = [64 86 92 101 108 183];
        Params.Duration.Ictal = 3*60;
        disp(['Changing seizure duration to ', num2str(Params.Duration.Ictal),...
            ' [s]']);
        
end

% Go from 257 to 204 channel configuration

FileNameChan257 = '/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_Seizure_Epochs_V04_05/EGI 257.Geneva Average 13.10-10.xyz';
FileNameChan204 = '/Users/annavybornova/EPFL/Master_4/GEM/01Data/HUG_Seizure_Epochs_V04_05/EGI 204.Geneva Average 13.10-10.xyz';

[IndSameChannels, Labels] = findSameChannels(FileNameChan257, FileNameChan204, 257, 204);
DataAll = DataAll(IndSameChannels, :);

% readlocs function puts the channels in the different order; make sure
% that the order in DataAll is exactly the same as in chanlocs

chanlocs = readlocs(FileNameChan204);
indSame = 1;
for iChL = 1:length(chanlocs)
    for iChO = 1:length(Labels)
        if strcmp(chanlocs(iChL).labels, Labels(iChO))
            IndReorderToChanlocs(iChL) = iChO;
        end
    end
end

DataAll = DataAll(IndReorderToChanlocs, :);

%@sanitycheck done! order in DataAll corresponds to the order in chanlocs
%and channels in DataAll correspond to the 204 channel configuration

%% STEP 2: BAD CHANNELS HANDLING

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

%% STEP 3: DIVIDING DATA INTO THE SEGMENTS OF INTEREST

% Time limits for the epoching (needed to reduce the time of the filtering)
% 5min at the beginning of the recording
% 2min right before the seizure
% 30 s after the seizure onset

SamplesInteriactal = Params.Duration.Interictal*EEGSamplingRate;
SamplesIctal = Params.Duration.Ictal*EEGSamplingRate;
SamplesPreIctal = Params.Duration.PreIctal*EEGSamplingRate;

DataBeforeSeizure = DataAll(:, 1:SamplesInteriactal);
DataSeizure = DataAll(:, SeizureOnset:(SeizureOnset+SamplesIctal));
DataPreSeizure = DataAll(:, (SeizureOnset - SamplesPreIctal -1): ...
    SeizureOnset-1);
DataAll(:, (SeizureOnset+SamplesIctal):end) = [];

%% STEP 4: DATA FILTERING

h1=fdesign.highpass('N,F3dB',Params.Filtering.FiltOrder,Params.Filtering.LowCf,...
    EEGSamplingRate);
d1=design(h1,'Butter');
h2=fdesign.lowpass('N,F3dB',Params.Filtering.FiltOrder,Params.Filtering.HighCf,...
    EEGSamplingRate);
d2=design(h2,'Butter');

DataBeforeSeizure = filterEEG(DataBeforeSeizure, d1, d2);
DataPreSeizure = filterEEG(DataPreSeizure, d1, d2);
DataSeizure = filterEEG(DataSeizure, d1, d2);
DataAll = filterEEG(DataAll,d1,d2); 

if Params.PlotAll
    eegplot([DataBeforeSeizure, DataPreSeizure, DataSeizure], 'srate', EEGSamplingRate, ...
        'winlength', 60, 'dispchans', 50);
end
%% STEP 5: COMPUTE THE GFT OF THE DATA

disp('... Adjacency Matrix generation');
A = adj_generator([[chanlocs.X]', [chanlocs.Y]', [chanlocs.Z]'], Params.Gft.R);
nNeighbors = sum(A);
avgNeighbors = mean(nNeighbors);

disp(['R = ', num2str(Params.Gft.R)]);
disp(['avg number of neighbors = ', num2str(avgNeighbors)]);

GftBeforeSeizure = gft(DataBeforeSeizure, A);
GftSeizure = gft(DataSeizure, A);
GftPreSeizure = gft(DataPreSeizure, A);


% STEP 6: DATA EPOCHING AND LABELING

Params.Epochs.WinEpoch = 2*EEGSamplingRate;
Params.Epochs.Overlap = Params.Epochs.WinEpoch/2;

Params.Epochs.Label.Seizure = 1;
Params.Epochs.Label.BeforeSeizure = 0;
Params.Epochs.Label.PreSeizure = 2;

disp(['Epoch length: ', num2str(Params.Epochs.WinEpoch/EEGSamplingRate), ' [s]' ])
disp(['Epoch overlap: ', num2str(Params.Epochs.Overlap/EEGSamplingRate), ' [s]' ])
% generating the data frame

% EEG

EpochSeizure.Eeg = epochIndecesGeneration(DataSeizure, Params.Epochs.WinEpoch,...
    Params.Epochs.Overlap, Params.Epochs.Label.Seizure );

EpochBeforeSeizure.Eeg = epochIndecesGeneration(DataBeforeSeizure, ...
    Params.Epochs.WinEpoch,...
    Params.Epochs.Overlap, Params.Epochs.Label.BeforeSeizure );

EpochPreSeizure.Eeg = epochIndecesGeneration(DataPreSeizure, ...
    Params.Epochs.WinEpoch,...
    Params.Epochs.Overlap, Params.Epochs.Label.PreSeizure );

% GFT
EpochSeizure.Gft = epochIndecesGeneration(GftSeizure, Params.Epochs.WinEpoch,...
    Params.Epochs.Overlap, Params.Epochs.Label.Seizure );

EpochBeforeSeizure.Gft = epochIndecesGeneration(GftBeforeSeizure, ...
    Params.Epochs.WinEpoch,...
    Params.Epochs.Overlap, Params.Epochs.Label.BeforeSeizure );

EpochPreSeizure.Gft = epochIndecesGeneration(GftPreSeizure, ...
    Params.Epochs.WinEpoch,...
    Params.Epochs.Overlap, Params.Epochs.Label.PreSeizure );
Epoch.Gft = [EpochSeizure.Gft EpochBeforeSeizure.Gft EpochPreSeizure.Gft];
Epoch.Eeg = [EpochSeizure.Eeg EpochBeforeSeizure.Eeg EpochPreSeizure.Eeg];

%% STEP 7: PRELIMINARY STATISTICAL TESTS
if Params.PreliminaryStatistics
    alpha = 0.05;
    disp('... Computing preliminary statistics')
    disp(['Alpha level: ', num2str(alpha)]);
    
    
    for iEp = 1:length(Epoch.Gft)
        Epoch.Gft(iEp).Mean = mean(Epoch.Gft(iEp).Signal)';
        Epoch.Eeg(iEp).Mean = mean(Epoch.Eeg(iEp).Signal)';
        
        Epoch.Gft(iEp).MeanAbs = mean(abs(Epoch.Gft(iEp).Signal))';
        Epoch.Eeg(iEp).MeanAbs = mean(abs(Epoch.Eeg(iEp).Signal))';
        
        Epoch.Gft(iEp).Variance = std(Epoch.Gft(iEp).Signal)';
        Epoch.Eeg(iEp).Variance = std(Epoch.Eeg(iEp).Signal)';
    end
    
    LabelsTot = [Epoch.Eeg.Label];
    
    MeanGft = [Epoch.Gft(:).Mean]';
    MeanEeg = [Epoch.Eeg(:).Mean]';
    MeanAbsGft =[Epoch.Gft(:).MeanAbs]';
    MeanAbsEeg = [Epoch.Eeg(:).MeanAbs]';
    VarianceGft = [Epoch.Gft(:).Variance]';
    VarianceEeg = [Epoch.Eeg(:).Variance]';
    
%     Metrics = { MeanGft, MeanEeg, MeanAbsGft, MeanAbsEeg, VarianceGft,...
%         VarianceEeg };

tests = {'ranksum', 'ttest2', 'signrank', 'ttest'};

Metrics = { MeanGft, MeanEeg, VarianceGft, VarianceEeg };

metricName = {'mean', 'variance'};
    
    testDataType = [1,0; 1,2; 2,0];
    %testDataType = [1,0];
    %testDataTypeLabel = {'Seiz vs non-seiz', 'Seiz vs pre-seiz', 'pre-seiz vs non-seiz' };
    testDataTypeLabel = {'Seiz vs non-seiz'};
    for iT = 1%1:length(testDataType)
        for iM = 1:length(Metrics)
            [p{iT}{iM}, indSorted{iT}{iM}, indSign{iT}{iM}] = calculateStatistics(Metrics{iM},...
                LabelsTot, alpha,  testDataType(iT, 1), testDataType(iT, 2));
            
        end
        
        figure,
        ind = 1;
        for i = 1:length(Metrics)/2
            for j = 2%1: size(p{iT}{1},1)
                subplot(2,2, ind)
                ind = ind + 1;
                
                plot(p{iT}{2*i - 1}(j,:), '-k');
                hold on;
                plot(p{iT}{2*i}(j,:), '--k');
                legend('Gft', 'Eeg');
                xlabel('Sorted frequency/channel', 'FontSize', 14);
                ylabel('pValue','FontSize', 14);
                %title([testDataTypeLabel{iT}, metricName{i}, tests{j}]);
                title([metricName{i},'', tests{j}]);
                set(gca, 'FontSize', 14);
                %ylim([0,alpha]);
                xlim([0, 210]);
                axis square
            end
        end
        
        set(gcf, 'Position', [100, 100, 1049, 895]);
        saveas(gcf,[FigNameSave,testDataTypeLabel{iT},...
            'PvaluePreliminaryStatisctics.png']);
        
%        figure, 
%        topoplot(ones(1, length(indSign{1,1}{1,2}{1,1})), chanlocs(indSign{1,1}{1,1}{1,1}), ...
%             'electrodes','labels','headrad', 0); 
       
        
        
        
%         figure,
%         ind = 1;
%         for i = 1:length(Metrics)/2
%             for j = 1: size(p{iT}{1},1)
%                 subplot(3, 4 , ind)
%                 ind = ind + 1;
%                 
%                 plot(indSign{iT}{2*i - 1}{j}, 'ok');
%                 hold on;
%                 %                 plot(indSorted{2*i}(j,:), '*k');
%                 legend('Gft');
%                 xlabel('Sorted frequency');
%                 ylabel('frequency value');
%                 title([testDataTypeLabel{iT},metricName{i}, tests{j}]);
%                 ylim([1, 204])
%                 xlim([0, 204]);
%             end
%         end
        
% 
%         set(gcf, 'Position', [100, 100, 1049, 895]);
%         saveas(gcf,[FigNameSave,testDataTypeLabel{iT},...
%             'GftSortedPreliminaryStatisctics.png']);
    end
end


% USING THE RANKFEATURE FUNCTION IN MATLAB

%RankFeatName = {'ttest', 'entropy', 'bhattacharyya', 'roc', 'wilcoxon', 'fisher'};
RankFeatName = {'ttest', 'fisher'};


for iT = 1%:length(testDataType)
    figure,
    Ind1 = find(LabelsTot == testDataType(iT, 1));
    Ind2 = find(LabelsTot == testDataType(iT, 2));
    L = LabelsTot([Ind1 Ind2]);
    for iM = 1:length(Metrics)
        
        for iR = 1 %:length(RankFeatName)
            if iR ~= length(RankFeatName)
                [IDX{iT}{iM}(iR,:), Z{iT}{iM}(iR,:)] = rankfeatures( Metrics{iM}([Ind1, Ind2],:)',...
                    L, 'Criterion', RankFeatName{iR});
            elseif iR == length(RankFeatName)
                [IDX{iT}{iM}(iR,:), Z{iT}{iM}(iR,:)] = rankfeat( Metrics{iM}([Ind1, Ind2],:),...
                    L, RankFeatName{iR});
            end
        end
        
    end
    clear Ind1 Ind2 L
    
    %@ 2017 06 11
    
    figure
    f = ones(1,204); 
    f(IDX{1}{3}(1:10)) = 0; 
    subplot(2,1,1)
    
    imagesc(f)
   
    title('Discriminative frequencies')
    colormap gray
    
    V = load('V.mat'); 
for i = 1:10 
    
    subplot(2,10,10+i)
        topoplot(V.V(:,IDX{1}{3}(i)), chanlocs,  'electrodes', 'on', 'headrad', 0);
        title(['f = ', num2str(IDX{1}{3}(i))]);
   
end
%     ind = 1;
%     for i = 1:length(Metrics)/2
%         for j = 1:size(IDX{iT}{1},1)
%             subplot(3,3,ind)
%             plot(IDX{iT}{2*i-1}(j,1:20), '*k');
%             ind = ind+1;
%             title([testDataTypeLabel{iT}, ' ', metricName{i}, ' ',...
%                 RankFeatName{j} ]);
%             ylim([0 210]);
%             xlabel('First 20 ranked features');
%             ylabel('GFT frequency');
%             
%         end
%     end
%     
%     set(gcf, 'Position', [100, 100, 1049, 895]);
%     saveas(gcf,[FigNameSave,testDataTypeLabel{iT},...
%         'First20impFeatures.png']);
    
end

%@ 2017 05 22 

% IDX{1}{3} - corresponds to the variance of the GFT 
numGft = 20; 
for iM = 1:size(IDX{1}{3},1 )
    subplot(1, size(IDX{1}{3},1 ), iM )
    vec = zeros(1, length(chanlocs)); 
    vec(IDX{1}{3}(iM, 1:numGft)) = ones(1, numGft) ; 
    imagesc(vec'); 
    colormap gray 
    title(RankFeatName{iM}); 
end 


%% Assess what are the most important EEG channels
for iT = 1%:length(testDataType)
    figure,
    ind = 1;
    for i = 1:length(Metrics)/2
        for j = 1:size(IDX{iT}{1},1)
            
            Channels = IDX{iT}{2*i}(j,1:20);
            subplot(2,2,ind)
            topoplot(Z{iT}{2*i}(j,1:20), chanlocs(Channels), ...
                'electrodes','labels','headrad', 0);
            ind = ind + 1;
            
            title([testDataTypeLabel{iT}, ' ', metricName{i}, ' ',...
                RankFeatName{j} ]);
            
        end
    end
    set(gcf, 'Position', [100, 100, 1049, 895]);
    saveas(gcf,[FigNameSave,testDataTypeLabel{iT},...
        'EEGFirst20impFeatures.png']);
end
close all
%% Assess the sign of the GFT

% GftPhase = sign([GftBeforeSeizure, GftPreSeizure, GftSeizure]);
% figure,
% imagesc(GftPhase);
% colormap gray
% 
% h(1) = subplot(2,1,2)
% imagesc([DataBeforeSeizure, DataPreSeizure, DataSeizure]);
% linkaxes(h, 'x');

%
% SumGftPhase = sum(GftPhase);
% figure,
% plot(crosscorr(abs(SumGftPhase)));

%% STEP 8: classifier implementation

% Metrics = { MeanGft, MeanEeg, MeanAbsGft, MeanAbsEeg, VarianceGft,...
%     VarianceEeg };

% Hyperparameter --> feature ranking; number of input features

classifier = {'diaglinear', 'diagquadratic'};

% classification seizure/ non-seizure [1,0]

nFolds = 3;
Ind1 = find(LabelsTot == testDataType(1, 1));

Ind2 = find(LabelsTot == testDataType(1, 2));

indEq = randperm(length(Ind2), length(Ind1));
Ind2 = Ind2(indEq);

L = LabelsTot([Ind1 Ind2]);
partition = cvpartition(L, 'Kfold', nFolds);

for iR = 1%:3
    for iM = 1:length(Metrics)
        DataSet = Metrics{iM}([Ind1 Ind2], :);
        DataSet = zscore(DataSet);
        
        
        for k = 1:length(classifier)
            for j = 1:nFolds
                repartition(partition)
                indTrain = training(partition, j);
                indTest = test(partition,j);
                
                TrainSet = DataSet(indTrain,:);
                TrainLabels = L(indTrain)';
                TestSet = DataSet(indTest,:);
                TestLabels = L(indTest);
                
                % feature ranking
                if ~strcmp(RankFeatName{iR}, 'fisher')
                    IdxRanking = rankfeatures(TrainSet', TrainLabels, 'Criterion', RankFeatName{iR});
                elseif strcmp(RankFeatName{iR}, 'fisher')
                    IdxRanking = rankfeat(TrainSet, TrainLabels, RankFeatName{iR});
                end
                for iF = 1:length(chanlocs)
                    if ~strcmp(classifier{k}, 'SVM')
                        [class, errTrain] = classify(TestSet(:,IdxRanking(1:iF)), ...
                            TrainSet(:,IdxRanking(1:iF)), TrainLabels, classifier{k});
                        perf = classperf(TestLabels, class);
                        Acc{k}{iM}(j,iF) = perf.CorrectRate;
                        AccTrain{k}{iM}(j,iF) = 1-errTrain;
                        
                    elseif strcmp(classifier{k}, 'SVM')
                        SVMStruct = svmtrain(TrainSet(:,IdxRanking(1:iF)), TrainLabels);
                        Group = svmclassify(SVMStruct,TestSet(:,IdxRanking(1:iF)));
                        GroupTrain = svmclassify(SVMStruct, TrainSet(:,IdxRanking(1:iF)));
                        perf = classperf(TestLabels, Group);
                        Acc{k}{iM}(j, iF) = perf.CorrectRate;
                        perfTrain = classperf(TrainLabels, GroupTrain);
                        AccTrain{k}{iM}(j, iF) = perfTrain.CorrectRate;
                    end
                    
                    
                end
                
            end
            clear TrainSet TestSet indTrain indTest IdxRanking class TrainLabels ...
                TestLabels SVMStruct Group GroupTrain perf perfTrain
        end
        clear DataSet
    end
    
  
    
    MName = {'MeanGft', 'MeanEeg', 'StdGft', 'StdEeg'}
    
    for iCl = 1:2
        figure,
        ind = 1
        for iM = 1: length(Acc{iCl})
            subplot(2,2,ind)
            errorbar(mean(AccTrain{iCl}{iM}), std(AccTrain{iCl}{iM}));
            ind = ind +1;
            xlabel('Feature number');
            ylabel('Overall Accuracy');
            title(['Train', classifier{iCl}, MName{iM}, RankFeatName{iR}]);
            ylim([0,1]);
        end
        
        set(gcf, 'Position', [100, 100, 1049, 895]);
        saveas(gcf,[FigNameSave,['Train', classifier{iCl},...
            RankFeatName{iR}],'.png']);
        
%         figure,
%         ind = 1
%         for iM = 1: length(Acc{iCl})
%             subplot(3,2,ind)
%             errorbar(mean(Acc{iCl}{iM}(1:100)), std(Acc{iCl}{iM}));
%             ind = ind +1;
%             xlabel('Feature number');
%             ylabel('Overall Accuracy');
%             title(['Test', classifier{iCl}, MName{iM}, RankFeatName{iR}]);
%             ylim([0,1]);
%         end
%         set(gcf, 'Position', [100, 100, 1049, 895]);
%         saveas(gcf,[FigNameSave,['Test', classifier{iCl},...
%             RankFeatName{iR}],'.png']);
        
        figure,
        ind = 1
        for iM = 1: length(Acc{iCl})
            subplot(2,2,ind)
            plot(mean(AccTrain{iCl}{iM}));
            hold on;
            plot(mean(Acc{iCl}{iM}));
            legend('Train', 'Test');
            ind = ind +1;
            xlabel('Feature number', 'FontSize', 14);
            ylabel('Overall Accuracy','FontSize', 14);
            title([classifier{iCl}, MName{iM}, RankFeatName{iR}], 'FontSize', 14);
            ylim([0,1]);
            xlim([0,100]); 
            axis square 
        end
        set(gcf, 'Position', [100, 100, 1049, 895]);
        saveas(gcf,[FigNameSave,['TrainAndTest', classifier{iCl},...
            RankFeatName{iR}],'.png']);
    end
end


%% STEP 8: CLASSIFICATION WITH A SELECTED CLASSIFIER

% Variance gft variance eeg; fisher; 10 features

Ind1 = find(LabelsTot == testDataType(1, 1));
Ind2 = find(LabelsTot == testDataType(1, 2));

indEq = randperm(length(Ind2), length(Ind1));
Ind2 = Ind2(indEq);

L = LabelsTot([Ind1 Ind2]);
partition = cvpartition(L, 'Kfold', nFolds);
repartition(partition)
for iR = 1%:3;
iF = 10;
k = 2;
j = 1; 
for iM = 1:4
    DataSet = Metrics{iM}([Ind1 Ind2], :);
    DataSet = zscore(DataSet);
       
        indTrain = training(partition, j);
        indTest = test(partition,j);
        
        TrainSet = DataSet(indTrain,:);
        TrainLabels = L(indTrain)';
        TestSet = DataSet(indTest,:);
        TestLabels = L(indTest);
        
        % feature ranking
        if ~strcmp(RankFeatName{iR}, 'fisher')
            IdxRanking = rankfeatures(TrainSet', TrainLabels, 'Criterion', RankFeatName{iR});
        elseif strcmp(RankFeatName{iR}, 'fisher')
            IdxRanking = rankfeat(TrainSet, TrainLabels, RankFeatName{iR});
        end
        
        [class, errTrain] = classify(TestSet(:,IdxRanking(1:iF)), ...
            TrainSet(:,IdxRanking(1:iF)), TrainLabels, classifier{k});
        perf = classperf(TestLabels, class);
        ATest(iM) = perf.CorrectRate;
        ATrain(iM) = 1-errTrain;
        ConfMat{iM}{iR} = perf.DiagnosticTable; 
        ClassError{iM}{iR} =  perf.SampleDistributionByClass;
        Features{iM}{iR} = IdxRanking(1:iF);
       IdxRankingTot{iM}{iR} = IdxRanking; 
        
    clear TrainSet TestSet indTrain indTest IdxRanking class TrainLabels ...
        TestLabels SVMStruct Group GroupTrain 
    clear DataSet
end
end

for iM = 1:4
   for iR = 1%:3 
       figure, 
       if ~(mod(iM,2) == 0)
           f = zeros(204,1); 
           f(Features{iM}{iR}) = 1; 
           imagesc(f); 
            set(gcf, 'Position', [100, 100, 200, 895]);
            colormap gray
           xlabel('Feature Number'); 
           ylabel('GFT Frequency'); 
           title(['10Feat', classifier{k}, MName{iM}, RankFeatName{iR}])
            saveas(gcf,[FigNameSave,'10Feat', classifier{k}, MName{iM}, RankFeatName{iR},'.png']);
       else
           topoplot(ones(1,length(Features{iM}{iR})), chanlocs(Features{iM}{iR}), ...
                'electrodes','labels','headrad', 0);
             title(['10Feat', classifier{k}, MName{iM}, RankFeatName{iR}]); 
             saveas(gcf,[FigNameSave,'10Feat', classifier{k}, MName{iM}, RankFeatName{iR},'.png']);
       end
   end
end
clear L


%@ 2017 05 03 

FeaturesGftVariance  = Features{3}{1}; 

% calculating the eigen decomposition of A

L = Laplacian(A);
%calculate eigen values of L
L_eig = eig(L);
L_eig_sorted = sort(L_eig);
[V, D] = eig(L);

Vred = V(:,FeaturesGftVariance); % select the columns of V corresponding to
% to the selected spatial frequecies 

OriginalMap = [DataBeforeSeizure, DataPreSeizure, DataSeizure]; 

s = zeros(1,204); 
s(FeaturesGftVariance) = 1; 
GftMapRed = diag(s) * V' * OriginalMap;
MapRed = V*GftMapRed;

% compare the original map and the reduced map 

ind = 1; 
%coeff = [(1:2:30) * 1000, (30:2:48)* 1000, (48:2:60) * 1000 ]; 

coeff = [(30:48)* 1000, (48:60) * 1000 ]; 

figure,
ind = 1; 
for i = 1:length(coeff)
subplot(8,4 ,ind)
topoplot(MapRed(:,coeff(i)), chanlocs,'headrad', 0);
            ind = ind+1; 
end

figure, 
ind = 1; 
for i = 1:length(coeff)
subplot(8,4 ,ind)
topoplot(OriginalMap(:,coeff(i)), chanlocs,'headrad', 0);
            ind = ind+1; 
end

%ind = ind +1;
figure, 
   eegplot(MapRed, 'srate', EEGSamplingRate, ...
        'winlength', 60, 'dispchans', 50);
    
   eegplot(OriginalMap, 'srate', EEGSamplingRate, ...
        'winlength', 60, 'dispchans', 50);





figure,
ind = 1
for iV = 1:size(Vred,2)
    subplot(2,5,ind)
    topoplot(Vred(:,iV), chanlocs,'headrad', 0);
            
    title(['f = ', num2str(FeaturesGftVariance(iV))]); 
    ind = ind+1; 
end


figure, 
imagesc(sign(GftMapRed(sort(FeaturesGftVariance), : )));
colormap gray

figure, 

imagesc(sign(GftMapRed(sort(indSign), : )));
colormap gray


%% Look at the Fourier Phase 
DataForFft = [DataBeforeSeizure, DataPreSeizure, DataSeizure]; 
DataFft = fft(DataForFft,[],2);
for i = 1: 204 
    for j = 1:size(DataFft,2)
    PhaseFft(i,j) = (angle(DataFft(i, j))); 
    end
end

figure, 
imagesc(PhaseFft);
















