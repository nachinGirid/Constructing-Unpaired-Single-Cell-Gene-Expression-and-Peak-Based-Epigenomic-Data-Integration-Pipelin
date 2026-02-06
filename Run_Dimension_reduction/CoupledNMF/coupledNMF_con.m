function [W1,H1,W2,H2,A]=coupledNMF_con(PeakO,X,D,K,lambda,mu,seed,maxiter)
if nargin < 4
	K=15;
end
if nargin < 5
	lambda=1;
end
if nargin < 6
	mu=0.001;
end
if nargin < 7
seed = 1;
end
if nargin < 8
maxiter=50;
end
opt = statset('MaxIter',100,'Display','final');
rng(seed)
w0=rand(size(X,1),K);
rng(seed)
h0=rand(K,size(X,2));
[w,h] = nmf2(X,K,'Options',opt,'Algorithm','mult','w0',w0,'h0',h0);
rng(seed)
W00=rand(size(PeakO,1),K);
rng(seed)
H100=rand(K,size(PeakO,2));
[W0,H10,H20,A0,dnorm]=CoupledNMF_3Fac_con(PeakO,X,D,K,maxiter,lambda,W00,H100,h,D,mu);
%[W1,W2,H1,H2,dnorm]=CoupledNMF_3Fac_fineTune_con(PeakO,X,K,5*maxiter,W0,A*W0,H10,H20,mu);
cost_diff = 0;
attempt = 0;
corr_diff=1;
mu=mu*2;
while ((corr_diff>1e-4)&(attempt <=20))
    attempt = attempt + 1
    W = W0;
    H1 = H10;
    H2 = H20;
    A = A0;
    [W0,H10,H20,A0,dnorm]=CoupledNMF_3Fac_con(PeakO,X,D,K,10,lambda,W,H1,H2,A,mu);
    cost_diff = (dnorm(1,1)+dnorm(1,2))-(dnorm(end,1)+dnorm(end,2))
    corr_diff = dnorm(end,4) - dnorm(1,4)
    if cost_diff > 1e-4
        mu=mu*2;
    end
    if cost_diff < -2e-3
        break;
    end
end
W2=A*W;
W1=W;
end