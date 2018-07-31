function [da] = getRandDa(n)
% 普通函数:随机产生输入数据
      rng(1);
%      rng('default');
                                                ID = randi([11,18],1,n); %ID类型数
    LWH = zeros(n,3);
    isRota = zeros(1,n);
    [uniID] = unique(ID);
    nUniId =  length(uniID);
    isRotaId = randi([0,1], 1, nUniId);
    for i=1: nUniId
    idx = find(ID(:)==uniID(i));
                                                LWH(idx,1) = randi([100,1000]);
                                                LWH(idx,2) = randi([100,1000]);
                                                isRota(1,idx) = isRotaId(i);
        for j=1:length(idx)
                                                LWH(idx(j),3) = randi([500,1000]);%randi([1000,1500]);
        end            
    end
    
    da.LUArray.ID = ID;
    da.LUArray.isRota = isRota; %isRota: 判断某类型LU是否可以旋转,默认均可旋转,0表示不可以旋转
    
    da.LUArray.LWH = LWH';
                                                da.BinArray.LWH = [2400   5000   2400];%[2400   12750   2400];
    da.BinArray.BUFF = [50;50;0];
    da.LUArray.BUFF = [100;100];%[2;2];
                    %     da.BinArray.LWH = da.BinArray.LWHREAL - da.BinArray.BUFF;
                    %     da.LUArray.LWH = da.LUArray.LWHREAL + da.LUArray.BUFF;
    da.BinArray.Weight = 1000;
    da.LUArray.Weight = da.LUArray.LWH(1,:);
end


% function [da] = getRandDa(n)
% % 普通函数:随机产生输入数据
%      rng(222121);
% %      rng('default');
%                                                 ID = randi([11,13],1,n); %ID类型数
%     LWH = zeros(n,3);
%     [uniID] = unique(ID);
%     for i=1:length(uniID)
%     idx = find(ID(:)==uniID(i));
%                                                 LWH(idx,1) = randi([5,15]);
%                                                 LWH(idx,2) = randi([5,15]);
%         for j=1:length(idx)
%                                                 LWH(idx(j),3) = randi([1,8]);
%         end            
%     end
%     da.LUArray.ID = ID;
%     da.LUArray.LWH = LWH';
%                                                 da.BinArray.LWH = [53; 55; 10];
%     da.BinArray.BUFF = [0;0;0];
%     da.LUArray.BUFF = [0.1;0.1];%[2;2];
% %     da.BinArray.LWH = da.BinArray.LWHREAL - da.BinArray.BUFF;
% %     da.LUArray.LWH = da.LUArray.LWHREAL + da.LUArray.BUFF;
%     da.BinArray.Weight = 1000;
%     da.LUArray.Weight = da.LUArray.LWH(1,:);
% end