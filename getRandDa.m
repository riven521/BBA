function [da] = getRandDa(n)
     rng(1);
    ID = randi([11,25],1,n); %ID¿‡–Õ ˝
    LWH = zeros(n,3);
    [uniID] = unique(ID);
    for i=1:length(uniID)
    idx = find(ID(:)==uniID(i));
    LWH(idx,1) = randi([5,15]);
    LWH(idx,2) = randi([5,15]);
        for j=1:length(idx)
            LWH(idx(j),3) = randi([1,5]);
        end            
    end
    da.LUArray.ID = ID;
    da.LUArray.LWHREAL = LWH';
    da.BinArray.LWHREAL = [53; 55; 33];
    da.BinArray.BUFF = [3;3;3];
    da.LUArray.BUFF = [2;2;0];
    da.BinArray.LWH = da.BinArray.LWHREAL - da.BinArray.BUFF;
    da.LUArray.LWH = da.LUArray.LWHREAL + da.LUArray.BUFF;
end