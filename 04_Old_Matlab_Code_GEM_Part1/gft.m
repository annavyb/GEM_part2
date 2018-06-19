function [ GftData] = gft(Data, A)
% Graph Fourier Transform (GFT) of the data 
%
%--------------------------------------------------------------------------
%
%  Input:
%       'Data' -[double] - [EEG data in the format rows: channels x columns:
%                           timepoints]
%       'A' - [double nChan x nChan] - [Adjacency matrix]
%
%  Output:
%       'GftData' - [structure] - [ GFT of the Data in the format rows: channels 
%                           x columns: timepoints]
%
% -------------------------------------------------------------------------
%
%  Version: V1 2017 04 29 
%  Author: VYBORNOVA Anna 
% 
% ------------------------- BEGIN CODE ------------------------------------


L = Laplacian(A);
%calculate eigen values of L
L_eig = eig(L);
L_eig_sorted = sort(L_eig);
[V, D] = eig(L);

GftData = V'* Data;

end

