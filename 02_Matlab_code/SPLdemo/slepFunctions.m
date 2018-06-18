function [SL,SL_conc,SL_nodes]=slepFunctions(Utr0,CONST_NODES,CONST_NUM_SLEP)

% CONST_NUM_SLEP: number of Slepian functions to compute

fprintf('Computing Slepian functions: ');

iter_S=1;

for iter_M=1:length(CONST_NODES),
    
    if mod(iter_M-1,10)==0,
    fprintf('%04d/%04d',iter_M,length(CONST_NODES));
    end;
    
    CONST_M=CONST_NODES(iter_M);
    
    %mask=zeros(1,msize);
    %mask(CONST_M)=1;
    
    %M=diag(mask);
    %D=Utr0'*M*Utr0;
    D=Utr0(CONST_M,:)'*Utr0(CONST_M,:);
    
    [V,S2]=eig(D);
    
    conc=real(diag(S2));
    [conc_sorted,conc_idx]=sort(conc,'descend');
    
    UtrV=real(Utr0*V);
    
    idx=find(conc_sorted>0);
    idx=idx(1:min(length(idx),CONST_NUM_SLEP));
    
    for iter_S2=1:length(idx),
        SL(:,iter_S)=UtrV(:,conc_idx(idx(iter_S2)));
        SL_conc(iter_S)=conc_sorted(idx(iter_S2));
        SL_nodes(iter_S)=iter_M;
        iter_S=iter_S+1;
    end;
        
    if mod(iter_M-1,10)==0,
    fprintf('\b\b\b\b\b\b\b\b\b');
    end;
end;

fprintf('\n');

