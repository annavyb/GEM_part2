%  Seizure detection approach using the projection in the space
%  spanned by K optimal Slepians, where K is the Shannon number 
%
%  The normal/abnormal points are identified with the treshold malahanobis
%  calculated in the Slepian space (axes of the covariance matrix of the 
%  hyperellipsoide of the data) 
% -------------------------------------------------------------------------
%
%  Version: V1 2018 06 16
%  Author: VYBORNOVA Anna 
% 
% ------------------------- BEGIN CODE ------------------------------------


clc 
clear all
close all

dbstop if error;
set(0,'defaultfigurecolor',[1 1 1]);
set(0,'defaulttextinterpreter','none');

% set the parameters
Params = setParams(); 

% adding the paths to the Slepian DEMO code and to the GFT package 
addpath("../SPLdemo");
addpath("../GFTpackage"); 

%% Load the example of the data
% load the chanlocs --> to execute this line eeglab need to be launched, as
% readlocs is the eeglabfunction

FileNameChan204 = '../../01_Data/EGI 204.Geneva Average 13.10-10.xyz'; 
chanlocs = readlocs(FileNameChan204);

% resting clean data + spike data (pre-processed)

PathResting = '../../01_Data/20180424EEGresting/Clean'; 
PathSpike = '../../01_Data/20180424EEGspikes'; 

subjects = ["03", "07", "11", "14", "16"]; 
subID = 1; 

% for i = 1:length(subjects)
%    load(strcat(PathResting,"/sub-", subjects(i), "_rest_clean_preprocessed.mat")); 
%    load(strcat(PathSpike,"/sub-", subjects(i), "pre-processed.mat")); 
%    
%    
% end

load(strcat(PathResting,'/sub-', subjects(subID), '_rest_clean_preprocessed.mat')); 
load(strcat(PathSpike,'/sub-', subjects(subID), 'pre-processed.mat')); 

% load the datastructure
load('../../01_Data/DataStructure.mat');

%% Compute the Slepian basis 

disp('... Adjacency Matrix generation');
Atemp = adj_generator([[chanlocs.X]', [chanlocs.Y]', [chanlocs.Z]'], Params.Gft.R);
nNeighbors = sum(Atemp);
avgNeighbors = mean(nNeighbors);
[SL0, Params.Slepian] = computeSlepians(Atemp, chanlocs, Params.Slepian); 


%%  


