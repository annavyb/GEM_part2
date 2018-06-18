function [SL,SL_conc,SL_nodes,SL_shannon,SL_V]=...
    slepFunctionsGroup(Utr0,CONST_NODES,CONST_THR)

% CONST_THR

fprintf('Computing Slepian functions: ');

iter_S=1;

for iter_M=1:length(CONST_NODES),
    
    if mod(iter_M-1,10)==0,
    fprintf('%04d/%04d',iter_M,length(CONST_NODES));
    end;
    
    CONST_M=CONST_NODES{iter_M};
    
    %mask=zeros(1,msize);
    %mask(CONST_M)=1;
    
    %M=diag(mask);
    %D=Utr0'*M*Utr0;
    D=Utr0(CONST_M,:)'*Utr0(CONST_M,:);
    
    [V,S2]=eig(D);
    
    conc=real(diag(S2));
    [conc_sorted,conc_idx]=sort(conc,'descend');
    
    UtrV=real(Utr0*V);
    
    % compute Shannon number
    SL_shannon(iter_M)=(length(CONST_M)/size(Utr0,1)*size(Utr0,2));
    SL_V{iter_M}=V(:,conc_idx);

    if CONST_THR==1,
        idx=1;
    elseif CONST_THR>1,
        idx=1:CONST_THR;
    else
        idx=find(conc_sorted>CONST_THR);
        if length(idx)==0,
            idx=1;
        end;
    end;
    
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

