function [S,H1,H2]=Fold_sim_based_clustering3_1(O,k1,k2,m)
%m number of cell in ATAC
%%%construct similarity matrix
%sim=O'*O;
gamma = 1;
n=size(O,2);
sim=squareform(pdist(O','cosine'));
sim=max(max(sim))-sim;
%sim=sim.*(eye(size(sim,1))--0);
KK = full(sum(sim));
twom = sum(KK);
sim_norm=sim - KK'*KK/twom;
[d1 f1]=maxk(sim_norm(1:m,:),k1);
[d2 f2]=maxk(sim_norm(1+m:end,:),k2);
f2=f2+m;
KNN_P1=[repelem([1:n]',k1,1) reshape(f1,k1*n,1)];
KNN_P2=[repelem([1:n]',k2,1) reshape(f2,k2*n,1)];
KNN_P=[KNN_P1;KNN_P2];
KNN=sparse(KNN_P(:,2),KNN_P(:,1),ones(size(KNN_P,1),1),size(sim_norm,2),size(sim_norm,2));
%O_smoothed=O*KNN/(k1+k2);
%H1=O_smoothed(:,1:m);
%H2=O_smoothed(:,1+m:end);
%%module detection
%addpath D:\duren\Matlab\GenLouvain-master
%A=1*(KNN+KNN'>0);
dis=pdist(KNN','jaccard');
A=1-squareform(dis);
%A=sim;
KK = full(sum(A));
twom = sum(KK);
B = A - gamma*KK'*KK/twom;
%%impute H
[d1 f1]=maxk(B(1:m,:),k1);
[d2 f2]=maxk(B(1+m:end,:),k2);
f2=f2+m;
KNN_P1=[repelem([1:n]',k1,1) reshape(f1,k1*n,1)];
KNN_P2=[repelem([1:n]',k2,1) reshape(f2,k2*n,1)];
KNN_P=[KNN_P1;KNN_P2];
KNN=sparse(KNN_P(:,2),KNN_P(:,1),ones(size(KNN_P,1),1),size(B,2),size(B,2));
O_smoothed=O*KNN/(k1+k2);
H1=O_smoothed(:,1:m);
H2=O_smoothed(:,1+m:end);
clear sim sim_norm dis
[S,Q] = genlouvain(B);
[num2str(length(unique(S))),' clusters detected']