function [O1,tf1,idf]=tfidf_trans(O)
%tf=O./(sum(O)+1);
%tf1=log1p(tf*100000);
O=1*(O>0);
tf1=O./log(1+sum(O));
idf=log(1+size(O,2)./(1+sum(O>0,2)));
O1=tf1.*idf;
O1(isnan(O1))=0;