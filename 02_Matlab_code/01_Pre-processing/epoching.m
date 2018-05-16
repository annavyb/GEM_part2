function [data_out] = epoching(data_in, epoch_size)
% Performs the epoching of the input data  
%--------------------------------------------------------------------------
%
%  Input:
%       'data_in' -[2D array (channels x timepoints)] - [the input data]
%       'epoch_size' - [integer] - [number of points in the epoch, i.e. 1s 
%       epoch with fs = 1000 is 1000] 
%  Output:
%        'data_out' - [3D array (channels x epoch_size x number of epochs] -
%                       - the result of the epoching
% -------------------------------------------------------------------------
%
%  Version: V1 2018 05 14  
%  Author: VYBORNOVA Anna 
% 
% ------------------------- BEGIN CODE ------------------------------------

ind_end = fix(size(data_in,2)/epoch_size)*epoch_size;
data_out = reshape(data_in(:, 1:ind_end), size(data_in,1),...
    epoch_size, ind_end/epoch_size); 
    
    
end

