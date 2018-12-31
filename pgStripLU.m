function [pgon] = pgStripLU(idxs,T)
%��Strip�����к�idxs,��ȡLU����ж�Ӧ��pgon
idxlArray = find(T.LU_Strip(:,1)==idxs);
p=[];
for i=1:length(idxlArray)
    idxl=idxlArray(i);
    w=T.LWH(idxl,1) + T.margin(idxl,1 ) + T.margin(idxl,2 );
    h=T.LWH(idxl,2) +  T.margin(idxl,3 ) + T.margin(idxl,4 );
%         LU.LWH(1,:) =  LU.LWH(1,:) +  T.margin(1,: ) + T.margin(2,: ); %��ȣ����ң�
%     LU.LWH(2,:) =  LU.LWH(2,:) +  T.margin(3,: ) + T.margin(4,: ); %���ȣ����£�
    
    x=T.CoordLUBin(idxl,1)-T.margin(idxl,1 );
    y=T.CoordLUBin(idxl,2)-T.margin(idxl,3 );
    p=[p;pgRectangle(x,y,w,h);[NaN,NaN]];
end
pgon = polyshape(p);


