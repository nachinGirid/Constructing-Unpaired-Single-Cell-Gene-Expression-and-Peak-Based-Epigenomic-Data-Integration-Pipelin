function D=RE_TG_disBased(Symbol_location,Peak_location)
d0=2*10^(5);
RE_TG=[];
for i=1:size(Symbol_location,1)
    id=find((Peak_location(:,1)==Symbol_location(i,1)).*(abs(Peak_location(:,2)-Symbol_location(i,2))<10^6)==1);
    dis=abs(Peak_location(id,2)-Symbol_location(i,2));
    RE_TG=[RE_TG;[i*ones(length(id),1) id exp(-dis/d0)]];
end
D=sparse(RE_TG(:,1),RE_TG(:,2),RE_TG(:,3),size(Symbol_location,1),size(Peak_location,1));