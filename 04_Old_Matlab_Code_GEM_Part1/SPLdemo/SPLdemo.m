
addpath('/Users/annavybornova/EPFL/Master_4/GEM/02MATLAB_Code/SPLdemo/stlTools')
%% Settings
CONST_NUM_NODES=4567;
CONST_NORMALIZE=1;
CONST_OPERATOR=1;      % graph Laplacian
CONST_W=9;             % bandwidth
CONST_PRINT=0;

% what to output?
CONST_GEN=1; % 0: global (Laplacian eigenvectors)
             % 1: Slepian with energy concentration
             % 2: Slepian with graph embedding distance

%% Load data
[vertices,faces,normals,name] = stlRead('STL-files/oncapintada.stl');
i=[ faces(:,1); faces(:,1); faces(:,2); faces(:,2); faces(:,3); faces(:,3)];
j=[ faces(:,2); faces(:,3); faces(:,1); faces(:,3); faces(:,1); faces(:,2)];
A=sparse(i,j,ones(size(i,1),1),CONST_NUM_NODES,CONST_NUM_NODES); % the adjacency matrix
NODES_XY=vertices; % the 3d coordinates of my EEG cap will go here 

NODE_TYPE=ones(CONST_NUM_NODES,1)*2;
idx=find(NODES_XY(:,2)<-.4); % criteria to select a subgraph --> just the head of the animal
NODE_TYPE(idx)=1;
CONST_NODES{1}=idx;

%% Preprocess data, extract (truncated) graph spectrum

[A,D]=slepNormalize(A,CONST_NORMALIZE);

opts.issym=1;
opts.isreal=1;
opts.maxit=2500;
opts.disp=1;

[Utr,S1]=slepEigsLaplacian(A,D,CONST_W,opts); % Utr - eigen vectors, S1 - eigen values  
CONST_U_W=CONST_W;

%% Compute graph Slepians
idx=CONST_NODES{1};
idxW=1:CONST_W; % Band-limiting of the GFT spectrum
tmpS=(sum(diag(S1(idxW,idxW))));
S1=S1/tmpS % some sort of normalization of the eigen values 
tmpUS=Utr(idx,idxW)*sqrt(S1(idxW,idxW)); % The product of the band.limited GFT and the another normalization of eigen values
tmpU=Utr(idx,idxW);
C=tmpU.'*tmpU; C=real((C+C.')/2);
C2=tmpUS.'*tmpUS; C2=real((C2+C2.')/2);

if CONST_GEN==0,
    SL0=Utr(:,idxW);
    eig0=diag(S1(idxW,idxW)); 
    conc=diag(C);
    cut=diag(C2);
else
    lambda=CONST_GEN-1; % 0 favors Slepian, 1: favors graph embedding distance
    Cmix=(1-lambda)*C+lambda*C2;
    
    [myV,myD]=eig(Cmix);
    SL0=Utr(:,idxW)*myV;
    eig0=diag(SL0.'*(D-A)*SL0); eig0=eig0(1:CONST_W)/tmpS;
    conc=diag(myV.'*C*myV);
    cut=diag(myV.'*C2*myV);
end;
figure(1);
clf;
plotyy(1:CONST_W,conc,1:CONST_W,cut);
legend({'concentration','graph distance'});

%%
figure(2);
clf;
set(gcf,'Position',[590 528 530 400]);
for iter=1:9, 
    %subplot(3,3,iter);
    order=iter;
    row=2-floor((order-1)/3);
    col=mod(order-1,3);
    axes('Position',[1/3*col+.03 1/3*row 1/3 1/3]);
    name='';
    tmp=SL0(:,iter);
    [~,idx]=max(abs(tmp));
    if sign(tmp(idx(1)))<0,
        tmp=-tmp;
    end;
    stlPlotColor(vertices,faces,tmp,name); 
    axis off;
    axes('Position',[1/3*col 1/3*row 1/3 1/3]);
    axis off;
    text(0.02,.9,sprintf('%d',iter),'FontSize',14,'Interpreter','tex');
    %if CONST_GEN>0,
    %    text(0.02,.8,sprintf('\\mu =%4.3f',conc(iter)),'FontSize',12,'Interpreter','tex');
    %    text(0.02,.7,sprintf('\\mu_e=%4.3f',cut(iter)),'FontSize',12,'Interpreter','tex');
    %else
    PREFIX={'','',''};
    PREFIX{CONST_GEN+1}='\bf';
    text(0.02,.8,sprintf('%s\\lambda=%4.3f',PREFIX{1},eig0(iter)),'FontSize',12,'Interpreter','tex');
    text(0.02,.7,sprintf('%s\\mu=%4.3f',PREFIX{2},conc(iter)),'FontSize',12,'Interpreter','tex');
    text(0.02,.6,sprintf('%s\\xi=%4.3f',PREFIX{3},cut(iter)),'FontSize',12,'Interpreter','tex');
    %end;
end;

if CONST_PRINT==1,
    drawnow;
    set(gcf, 'color', 'none');
    set(gca, 'color', 'none');
    cd export_fig
    switch CONST_GEN,
        case 0,
            fname='fig-mesh-lapl.png';
        case 1,
            fname='fig-mesh-slep-conc.png';
        case 2,
            fname='fig-mesh-slep-emb.png';
    end;
    export_fig(sprintf('../%s',fname),'-m2');
    cd ..
end;
