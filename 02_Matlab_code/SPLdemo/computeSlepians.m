function [SL0, Options] = computeSlepians(Atemp, chanlocs, Options)
% COMPUTESLEPIANS - computes the Slepian basis functions 
%
%--------------------------------------------------------------------------
%
%  Input:
%
%  Output:
% -------------------------------------------------------------------------
%
%  Version: V1 2018 06 18
%  Author: VYBORNOVA Anna 
% 
% ------------------------- BEGIN CODE ------------------------------------

A = sparse(Atemp);
NODES_XY = [[chanlocs.X]' [chanlocs.Y]' [chanlocs.Z]'];
[A,D]=slepNormalize(A,Options.CONST_NORMALIZE);
[Utr,S1]=slepEigsLaplacian(A,D,204,Options.opts); % Utr - eigen vectors, S1 - eigen values

% define the nodes that are includes into the selected region
EuDist = sqrt((NODES_XY(:,1) - NODES_XY(Options.CONST_SEIZURE_NODE, 1)).^2 ...
    + (NODES_XY(:,2) - NODES_XY(Options.CONST_SEIZURE_NODE, 2)).^2 ...
    + (NODES_XY(:,3) - NODES_XY(Options.CONST_SEIZURE_NODE, 3)).^2);

[EuDistSorted, idx] = sort(EuDist);
idx((Options.CONST_SUBGRAPH_SIZE+1):end) = [];

NODE_TYPE(idx)=1;
Options.CONST_NODES{1}=idx;

Options.CONST_U_W = Options.CONST_W;

idx=Options.CONST_NODES{1};


idxW=1:Options.CONST_W; % Band-limiting of the GFT spectrum
tmpS=(sum(diag(S1(idxW,idxW))));
S1 = S1/tmpS; % some sort of normalization of the eigen values

tmpUS=Utr(idx,idxW)*sqrt(S1(idxW,idxW)); % The product of the band.limited GFT and the another normalization of eigen values
tmpU=Utr(idx,idxW);

C=tmpU.'*tmpU; C=real((C+C.')/2);
C2=tmpUS.'*tmpUS; C2=real((C2+C2.')/2);

if Options.CONST_GEN==0,
    SL0=Utr(:,idxW);
    eig0=diag(S1(idxW,idxW));
    conc=diag(C);
    cut=diag(C2);
else
    lambda=Options.CONST_GEN-1; % 0 favors Slepian, 1: favors graph embedding distance
    Cmix=(1-lambda)*C+lambda*C2;

    [myV,myD]=eig(Cmix);
    SL0=Utr(:,idxW)*myV;
    eig0=diag(SL0.'*(D-A)*SL0); eig0=eig0(1:Options.CONST_W)/tmpS;
    conc=diag(myV.'*C*myV);
    cut=diag(myV.'*C2*myV);
end;

end

