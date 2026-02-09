lambda=1;
mu=0.001;
K=100;
Outdir = '/pbmc';
Indir  = '/pbmc';
num_smooth=50; %optional
%barcode_filtered='barcodes.txt'; 
%optional input, if use, please remove the first '%'. if do not use, add % before the line.
%%%%%%%%%%%%%%%%%%do not edit below.

[X PeakO Symbol PeakName barcode Symbol_location Peak_location chr]=read_ATAC_GEX(Indir);
if exist('barcode_filtered','var')==1
barcode_use=importdata(barcode_filtered);
[~,idx]=ismember(barcode_use,barcode);
X=X(:,idx);
PeakO=PeakO(:,idx);
barcode=barcode_use;
end
%%%input used X, PeakO, Symbol, PeakName, Symbol_location, Peak_location. X, and PeakO are log2(1+X) transformed.
%%filter
if size(X,2)<30000
X=full(X);
PeakO=full(PeakO);
end

a2=sum(X>0,2);
top2=maxk(a2,5000);
cut2=top2(end);
X=X(a2>cut2,:);
Symbol_location=Symbol_location(a2>cut2,:);
Symbol=Symbol(a2>cut2,:);

%%filter top 5% coverage bins
a1=sum(PeakO>0,2);
top1=maxk(a1,floor(0.05*length(a1)));
cut1=top1(end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
keepIdx = find(a1 <= cut1); % Find valid row indices
if ~issparse(PeakO)
    PeakO = sparse(PeakO);
end

PeakO = PeakO(keepIdx, :);
Peak_location=Peak_location(keepIdx,:);
PeakName=PeakName(keepIdx,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PeakO=PeakO(a1<=cut1,:);
%%Peak_location=Peak_location(a1<=cut1,:);
%%PeakName=PeakName(a1<=cut1,:);
PeakO=full(PeakO);
%%filter based on fold cahnge compared to local neighborhood
O_sum=sum(PeakO,2);
O_exp=bin_background(O_sum,Peak_location,chr);
O_fold=O_sum./(O_exp+0.00001);
top1=maxk(O_fold,50000);
cut1=top1(end);
PeakO=PeakO(O_fold>cut1,:);
Peak_location=Peak_location(O_fold>cut1,:);
PeakName=PeakName(O_fold>cut1,:);
D=RE_TG_disBased(Symbol_location,Peak_location);
[PeakO,tf,idf] = tfidf_trans(PeakO);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('num_smooth','var')==0
num_smooth=floor(sqrt(size(X,2)+size(PeakO,2)));
num_smooth1=floor(sqrt(size(PeakO,2)));
num_smooth2=floor(sqrt(size(X,2)));
else
num_smooth1=num_smooth;
num_smooth2=num_smooth;
end

[W1,H1,W2,H2]=coupledNMF_con(PeakO,X,D,K,lambda,mu);
H1=H1./sum(H1,2)*size(H1,2);
H2=H2./sum(H2,2)*size(H2,2);
n1=size(PeakO,2);
[S,H10,H20]=Fold_sim_based_clustering3_1([H1,H2],num_smooth1,num_smooth2,n1);
S1=S(1:n1);
S2=S(1+n1:end);
save([Outdir,'/Paramaters_',num2str(lambda),'_',num2str(K),'.mat'],'W1','W2','H1','H2','H10','H20','lambda','K','S1','S2')
dlmwrite([Outdir,'/scATAC_label.txt'],S1,'\t')
dlmwrite([Outdir,'/scRNA_label.txt'],S2,'\t')
dlmwrite([Outdir,'/scATAC_embedding.txt'],H10,'\t')
dlmwrite([Outdir,'/scRNA_embedding.txt'],H20,'\t')