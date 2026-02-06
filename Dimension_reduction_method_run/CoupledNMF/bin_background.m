function O_exp=bin_background(O,Peak_location,chr)

O_exp=O;
for j=1:size(chr,1)
    idx=find(Peak_location(:,1)==j);
    if length(idx)>0
    length(idx)
    O1=O(idx,:);
    Peak_location1=Peak_location(idx,:);
    temp=[];
    for i=1:size(Peak_location1,1)
        id=find((abs(Peak_location1(:,2)-Peak_location1(i,2))<10^6)==1);
        temp=[temp;[i*ones(length(id),1),id,1/length(id)*ones(length(id),1)]];
        length(id)
    end
    G=sparse(temp(:,1),temp(:,2),temp(:,3),size(Peak_location1,1),size(Peak_location1,1));
    O_exp(idx,:)=G*O(idx,:);
    end
end