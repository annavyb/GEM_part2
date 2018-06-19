function [Epoch] = epochIndecesGeneration(Data, WinEpoch, Overlap, Label)
% This function the generation the epochs of the data 
%
%--------------------------------------------------------------------------
%
%  Input:
%       'Data' -[double] - [EEG data in the format rows: channels x columns:
%                           timepoints]
%       'WinEpoch' - [double 1x1] - [length of the epoching window in samples]
%       'Overlap'  - [double 1x1] - [overlap between the epoching windows]
%       'Label'    - [double/string/bool 1x1] - [the label of the epoch]
%
%  Output:
%       'Epoch' - [structure] - [structure containing the information for 
%                                epoch generation]
%               'Epoch.Begin' - [double] - [the index of the beginning of a 
%                                           given epoching window]    
%               'Epoch.End' - [double] - [the index of the end of a 
%                                           given epoching window]
%               'Epoch.Label' - [double/string/bool] - [the label of a 
%                                           given epoching window]
% -------------------------------------------------------------------------
%
%  Version: V1 2017 04 29 
%  Author: VYBORNOVA Anna 
% 
% ------------------------- BEGIN CODE ------------------------------------
IndecesEpoch =  buffer(1:size(Data,2), WinEpoch,...
    Overlap, 'nodelay');
IndecesEpoch(:, end) = []; % because of the zeros
IndecesEpoch(2:(end-1),:) = []; % leave only the info about the
%                                        beginning & the end of the window
for iEp = 1:size(IndecesEpoch,2)
    Epoch(iEp).Begin = IndecesEpoch(1,iEp);
    Epoch(iEp).End = IndecesEpoch(2,iEp);
    Epoch(iEp).Label = Label;
    Epoch(iEp).Signal = Data(:, IndecesEpoch(1,iEp):IndecesEpoch(2,iEp))';
end

end

