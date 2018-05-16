% In this script I visualize the new data 

clc 
clear
close all

% %% spikes 
Path = "/Users/annavybornova/EPFL/Master_4/GEM/GEM_part2/01_Data/20180424EEGspikes/";
Filename = "sub-19.mat";

load(strcat(Path, Filename))

figure, 
plot(data(:,channels_to_keep, 1))

%% resting state 
clear Path Filename

Path = "/Users/annavybornova/EPFL/Master_4/GEM/GEM_part2/01_Data/20180424EEGresting/Resting/";
Filename = "sub-19_rest.mat"; 

load(strcat(Path, Filename))

figure, 
plot(EEG)


% 








