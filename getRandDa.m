function [da] = getRandDa(n)
% 普通函数:随机产生输入数据
     rng(1);
    ID = randi([11,15],1,n); %ID类型数
    LWH = zeros(n,3);
    [uniID] = unique(ID);
    for i=1:length(uniID)
    idx = find(ID(:)==uniID(i));
    LWH(idx,1) = randi([5,35]);
    LWH(idx,2) = randi([5,35]);
        for j=1:length(idx)
            LWH(idx(j),3) = randi([1,8]);
        end            
    end
    da.LUArray.ID = ID;
    da.LUArray.LWH = LWH';
    da.BinArray.LWH = [53; 55; 10];
    da.BinArray.BUFF = [0;0;0];
    da.LUArray.BUFF = [0;0];
%     da.BinArray.LWH = da.BinArray.LWHREAL - da.BinArray.BUFF;
%     da.LUArray.LWH = da.LUArray.LWHREAL + da.LUArray.BUFF;
    da.BinArray.Weight = 1000;
    da.LUArray.Weight = da.LUArray.LWH(1,:);
end