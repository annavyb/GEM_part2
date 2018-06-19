function [DataFiltered] = filterEEG(Data, d1, d2)
% This function performs a band-pass filtering of the EEG data
%
%--------------------------------------------------------------------------
%
%  Input:
%       'Data' -[double] - [EEG data in the format rows: channels x columns:
%                           timepoints]
%       'd1' - 
%       'd2' -
%
%  Output:
%       'DataFiltered' - [double] - [filtered data matrix rows:channel x 
%                                   columns:timepoints]
% -------------------------------------------------------------------------
%
%  Version: V1 2017 04 29 
%  Author: VYBORNOVA Anna 
% 
% ------------------------- BEGIN CODE ------------------------------------


if size(Data, 1) ~= 204
    warning('The number of electrodes is different from 204')
end

disp('Filtering data...');

%filter the data
for nn=1:size(Data, 1)
    DataFiltered1(nn,:)=filtfilt(d1.sosMatrix,d1.ScaleValues,...
        double(Data(nn,:)));
    DataFiltered(nn,:)=filtfilt(d2.sosMatrix,d2.ScaleValues,...
        double(DataFiltered1(nn,:)));
end

end

