function [pgon] = pgStripLU(idxs,T)
%由Strip的序列号idxs,获取LU表格中对应的pgon
idxlArray = find(T.LU_Strip(:,1)==idxs);
p=[];
for i=1:length(idxlArray)
    idxl=idxlArray(i);
    w=T.LWH(idxl,1) + T.margin(idxl,1 ) + T.margin(idxl,2 );
    h=T.LWH(idxl,2) +  T.margin(idxl,3 ) + T.margin(idxl,4 );
%         LU.LWH(1,:) =  LU.LWH(1,:) +  T.margin(1,: ) + T.margin(2,: ); %宽度（左右）
%     LU.LWH(2,:) =  LU.LWH(2,:) +  T.margin(3,: ) + T.margin(4,: ); %长度（上下）
    
    x=T.CoordLUBin(idxl,1)-T.margin(idxl,1 );
    y=T.CoordLUBin(idxl,2)-T.margin(idxl,3 );
    p=[p;pgRectangle(x,y,w,h);[NaN,NaN]];
end
pgon = polyshape(p);


