% Separating the "useful" resting state data from the
% artifacts 


clc 
clear all
close all 


% The general parameters 

Path = "/Users/annavybornova/EPFL/Master_4/GEM/GEM_part2/01_Data/20180424EEGresting/Resting/"; 
Path_save = "/Users/annavybornova/EPFL/Master_4/GEM/GEM_part2/01_Data/20180424EEGresting/";  
subjects = ["03", "07", "11", "13", "14", "16", "19"]; 

for i = 1:length(subjects)
    % loads the resting-state data and separates the artifacts from the
    % clean data. Saves the structure
    
    filename = strcat(Path, "sub-", subjects(i), "_rest.mat"); 
    load(filename); 
    
    data_art = []; 
    
    data_clean = EEG(:, 1 : mrk_art(1, 1)-1); 
    count = 0; 
    for j = 1 : length(mrk_art)

        data_art = [data_art , EEG(:, mrk_art(j,1): mrk_art(j,2))]; 
        
        if (j == length(mrk_art))
            data_clean = [data_clean, EEG(:, (mrk_art(j,2)+1): end)]; 
        else
            data_clean = [data_clean, EEG(:, (mrk_art(j,2)+1): (mrk_art(j+1, 1)-1))]; 
        end
        
    end
    save(strcat(Path_save, "Clean/", "sub-", subjects(i), "_rest_clean.mat"), "data_clean")
    save(strcat(Path_save, "Artifacts/", "sub-", subjects(i), "_rest_art.mat"), "data_art")
    
end