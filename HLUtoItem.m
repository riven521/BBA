function [da] = HLUtoItem(da,ParaArray)
% 输入: LU: ID LWH
% 输出: 
% da.LUArray (1 LWH (已知)
% da.LUArray (2 LUBeItemArray (dim1:序号LU在某个item dim2:lu进入顺序(底层-高层)) 
% da.LUArray (3 LUorder 函数内LU排序顺序)
% da.ItemArray (1 LWH ) 

%% 初始化
% nDim LU维度 nLU LU数量 nItem Item数量 nLUid LU种类 tmpheightBin Bin最大高度
nDim = size(da.LUArray.LWH,1);
nLU = size(da.LUArray.LWH,2);
nItem = nLU;
nLUid = size(unique(da.LUArray.ID),2);

tmpUniqueBin = unique(da.BinArray.LWH(1:nDim,:)','rows')';
heightBin = tmpUniqueBin(3);
clear tmpUniqueBin;

%% LU排序并初始化
% 获取LUorder
% 获取LWHLuSort
% 获取IDLuSort
% sort 获取先id后高度的排序 %

tmpLUMatrix = [da.LUArray.ID; da.LUArray.LWH(1:nDim,:)]; %1:ID 4:Hight
[~,tepLUorder] = sortrows(tmpLUMatrix',[1 4],{'ascend','descend'});
% tepLUorder = 1:length(da.LUArray.ID)';
% tepLUorder = [2 3 4 1 5]';
LUorder = tepLUorder'; %获取 LU排序(先ID后高度)后的顺序
LWHLuSort = da.LUArray.LWH(:,LUorder);
IDLuSort = da.LUArray.ID(:,LUorder);
clear tmpLUMatrix tepLUorder;

%% 55 LU->Item转换 
% 获取LUBeItemArraySort : 每个排序后LU在哪个Item内  以及顺序
% 获取LWHItem:  新生成的Item的长宽高
LWHItem = zeros(nDim,nItem);   %Item的宽长高
LUBeItemArraySort = zeros(2,nLU); %dim1:属于第几个Item dim2:属于该Item第几个排放 555
itemBeLUArray = zeros(1,nItem);  % 每个Item内堆垛的LU数量 后期不用
iItem = 1;

for iLUid=1:nLUid
    heightLeft = heightBin;
    for iLU=1:nLU
        if IDLuSort(iLU) == iLUid %仅对当前LU对应Luid在该iLUid内的进行操作
            if heightLeft < LWHLuSort(nDim,iLU) %如当前LU高度满足(高度在第nDim行)
                iItem =  iItem + 1;
                heightLeft = heightBin;%  
            end              
            heightLeft = heightLeft - LWHLuSort(nDim,iLU); %更新剩余高度
            LWHItem(1:2,iItem) = LWHLuSort(1:2,iLU); %更新item长宽
            LWHItem(3,iItem) = LWHItem(3,iItem) + LWHLuSort(nDim,iLU); %更新item高度
            itemBeLUArray(iItem) = itemBeLUArray(iItem) + 1;
            LUBeItemArraySort(1,iLU) = iItem;
            LUBeItemArraySort(2,iLU) = itemBeLUArray(iItem);
        end
    end
    iItem =  iItem + 1;
end

%% 后处理 并赋值到da
% 获取LUBeItemArray : 每个排序后LU在哪个Item内  以及顺序
% 获取LWHItem:  新生成的Item的长宽高
% 获取LUorder: LU的排序
da.ItemArray.LWH = LWHItem(:,LWHItem(1,:)>0); % 去除未使用的
LUBeItemArray=LUBeItemArraySort;
LUBeItemArray(:,LUorder) = LUBeItemArraySort;
da.LUArray.LUorder = LUorder;
da.LUArray.LUBeItemArray = LUBeItemArray;

% 测试代码
% da.ItemArray.LWH 
% LUBeItemArray
% da.LUArray.LWH

% 输出主要结果:获得每个item包含的 原始 LU序号
% % for iItem = 1:max(da.LUArray.LUBeItemArray(1,:))
% %     [~,idx] = find(da.LUArray.LUBeItemArray(1,:)==iItem);
% %     fprintf('item %d 的长宽高为:  ',iItem);
% %     fprintf('( %d ) ',da.ItemArray.LWH(:,iItem));
% %     fprintf('\n');
% %     fprintf('item %d 包含 original LU 索引号(长宽高)为  \n  ',iItem);
% %     fprintf('%d ',idx);
% %     fprintf('( %d ) ', da.LUArray.LWH(:,idx));
% %     fprintf('\n');
% % end

end