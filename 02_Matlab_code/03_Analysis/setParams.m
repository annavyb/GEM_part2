function [P] = setParams()
% setting the parameters for the Slepian analysis 
%
%--------------------------------------------------------------------------
%
%  Input: 
%
%  Output:
%       'P' - [structure] - [ the structure of the parameters]
%
%       'Gft' - graph Fourier Transform related parameters 
%       'P.Gft.R' - the radius for the ajacency matrix construction 
%
%       'Slepian' - the Slepian-related parameters 
%       'P.Slepian.CONST_NORMALIZE' - see the SPLdemo package
%
%       'P.Slepian.opt' - options needed for the function 
%                         slepEigsLaplacian.m in the SPLdemo folder
%  The structure 'P' has the following fields 
% -------------------------------------------------------------------------
%
%  Version: V1 2018 06 18
%  Author: VYBORNOVA Anna 
% 
% ------------------------- BEGIN CODE ------------------------------------

P.Gft.R = 80; % distance in mm

P.Slepian.CONST_NORMALIZE = 1; 
P.Slepian.CONST_OPERATOR = 1;
P.Slepian.CONST_NORMALIZE_SLEPIAN_AMPLITUDE = 1;
P.Slepian.CONST_W = 78;
P.Slepian.CONST_SUBGRAPH_SIZE = 42;
P.Slepian.CONST_NUM_SLEP_TO_PLOT = 10;
P.Slepian.CONST_GEN=1;
P.Slepian.CONST_SEIZURE_NODE = 10; 

P.Slepian.opts.issym=1;
P.Slepian.opts.isreal=1;
P.Slepian.opts.maxit=2500;
P.Slepian.opts.disp=1;

end

