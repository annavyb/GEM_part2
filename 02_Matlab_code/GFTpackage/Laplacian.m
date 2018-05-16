function [L] = Laplacian(A)
%Laplacian is a function that computes the Laplacian matrix from adjacency
%matrix
%   Input: Adjacency matrix. Output: Laplacian matrix
n = size(A);%determining the size of Matrix A 
d = sum(A,2); % a column vector, where each element represents a sum of each row of A 
D = []; % creation of an empty matrix for the degree matrix

%% COMPUTING THE DEGREE MATRIX
for i = 1:n(1)
    D(i,i) = d(i);
    
end
%% COMPUTING THE ADJACENCY MATRIX 
L = D-A;


end

