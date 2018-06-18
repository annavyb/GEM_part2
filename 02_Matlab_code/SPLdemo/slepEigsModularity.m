function [Utr,S1]=slepEigsModularity(A,D,CONST_W,opts)

deg=diag(D);

if size(A,1)<500,

    [Utr,S1]=eig(full(A-deg*deg'/sum(deg))); 

    [S1,idx]=sort(diag(S1),'descend');
    S1=diag(S1);
    Utr=Utr(:,idx);

    Utr=Utr(:,1:max(CONST_W));
    S1=S1(1:max(CONST_W),1:max(CONST_W));
        
else
    
    [Utr,S1]=eigs(@Binline,size(A,1),max(CONST_W),'la',opts);

    [S1,idx]=sort(diag(S1),'descend');
    S1=diag(S1);
    Utr=Utr(:,idx);

end;

    function y=Binline(x)
        y=A*x-sum(deg.*x)/sum(deg)*deg;
    end

end
