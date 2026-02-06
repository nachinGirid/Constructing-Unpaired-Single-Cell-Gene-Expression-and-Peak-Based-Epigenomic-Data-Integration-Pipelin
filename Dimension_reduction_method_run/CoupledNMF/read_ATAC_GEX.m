function [E O Symbol PeakName barcode Symbol_location PeakName_location chr]=read_ATAC_GEX(foldername,cell_feature_order)
if nargin < 2
    cell_feature_order = 1;
end
a=dlmread([foldername,'/matrix.mtx'],' ',3,0);
if cell_feature_order ~=1
    a=a(:,[2,1,3]);
end
fileID = fopen([foldername,'/features.tsv']);
C = textscan(fileID,'%s %s %s %s %f %f','Delimiter','\t');
fclose(fileID);
chr=unique(C{1,4});
chr=chr(~cellfun(@isempty,regexp(chr,'^chr')));
isATAC=ismember(C{1,3},'Peaks');
PeakName=C{1,2}(isATAC);
Symbol=C{1,2}(~isATAC);
features=C{1,2};
[~,chr_idx]=ismember(C{1,4},chr);
Symbol_location=[chr_idx(~isATAC) C{1,5}(~isATAC)];
PeakName_location=[chr_idx(isATAC) C{1,5}(isATAC)];
barcode=importdata([foldername,'/barcodes.tsv']);
if isa(barcode, 'struct')
    barcode = barcode.textdata;
end
%%rna
[d f]=ismember(features,Symbol);
a_rna=a(ismember(a(:,1),find(d==1)),:);
a_rna(:,1)=f(a_rna(:,1));
E=sparse(a_rna(:,1),a_rna(:,2),log2(1+a_rna(:,3)),length(Symbol),length(barcode));
%%atac
[d f]=ismember(features,PeakName);
a_atac=a(ismember(a(:,1),find(d==1)),:);
a_atac(:,1)=f(a_atac(:,1));
O=sparse(a_atac(:,1),a_atac(:,2),log10(1+a_atac(:,3)),length(PeakName),length(barcode));
%%
idx=Symbol_location(:,1)>0;
Symbol=Symbol(idx);
Symbol_location=Symbol_location(idx,:);
E=E(idx,:);
idx=PeakName_location(:,1)>0;
PeakName=PeakName(idx);
PeakName_location=PeakName_location(idx,:);
O=O(idx,:);
%%non-zero 10 cell
idx=sum(E'>0)'>0;
Symbol=Symbol(idx);
Symbol_location=Symbol_location(idx,:);
E=E(idx,:);
idx=sum(O'>0)'>10;
PeakName=PeakName(idx);
PeakName_location=PeakName_location(idx,:);
O=O(idx,:);