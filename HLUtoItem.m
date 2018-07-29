function [da] = HLUtoItem(da,ParaArray)
% 重要函数:LU堆垛后形成Item %  行数:长宽高(row);  列数:托盘数量(coloum);
% Input ---  LUArray: ID LWH Weight
% Output --- LUArray: order LUBeItemArray
% Output --- ItemArray: ID LWH Weight
% da.LUArray.LUBeItemArray (2行*n列: 行1: LU在第几个item 行2:LU进入该item的顺序(底-高)) 
% da.LUArray.order (1行*n列: LU排序顺序)
% da.ItemArray.ID （1行*n列: ITEM的类型 等同内部LU类型)
% da.ItemArray.LWH (3行*n列: ITEM的长宽高-长宽与LU相同,高度为堆垛后高度)
% da.ItemArray.Weight (1行*n列: ITEM的重量)
% 嵌套函数(先排序 后转换)
% getLUorder
% getSortedLU
% getItem

    %% 初始化
% nDim LU维度 nLU LU数量 nItem Item数量 nLUid LU种类 
% heightBin Bin最大高度
nDim = size(da.LUArray.LWH,1);
nLU = size(da.LUArray.LWH,2);
nItem = nLU;
nLUid = size(unique(da.LUArray.ID),2);

tmpUniqueBin = unique(da.BinArray.LWH(1:nDim,:)','rows')';
heightBin = tmpUniqueBin(3);
clear tmpUniqueBin;

    %% LU排序
    % getLUorder - 获取LU的顺序(重点是高度递减排序)
    % getSortedLU - 获取按order排序后的LU:sortedLUArray
    da.LUArray.order = getLUorder(); %获取 LU排序(先ID递增,后高度递减)
    sortedLUArray = getSortedLU(da.LUArray.order);

    %% 55 LU->Item转换
    IDItem = zeros(1,nItem);   %Item的ID类型
    WeightItem = zeros(1,nItem);   %Item的重量
    LWHItem = zeros(nDim,nItem);   %Item的宽长高
% %     RotaFlagItem = ones(1,nItem)*5;   %Item的旋转flag（0和1,其它为未使用） TOBE
% DELE
    LUBeItemArraySort = zeros(2,nLU); %dim1:属于第几个Item dim2:属于该Item第几个排放 555
    
    getItem(); %555 转换
    
    da.ItemArray.ID = IDItem(:,IDItem(1,:)>0); % 去除未使用的     
    da.ItemArray.Weight = WeightItem(:,WeightItem(1,:)>0); % 去除未使用的 
% %     da.ItemArray.itemRotaFlag = RotaFlagItem(:,RotaFlagItem(1,:)<=1); %
% TO BE DELE

    da.ItemArray.LWH = LWHItem(:,LWHItem(1,:)>0); % 去除未使用的 
    da.LUArray.LUBeItemArray(:,da.LUArray.order) = LUBeItemArraySort; % da.LUArray.LUBeItemArray : 每个排序后LU在哪个Item内  以及顺序

    %% 额外变量+测试script
    % 额外变量
    getITEMIDArray();
    % 输出主要结果:获得每个item包含的 原始 LU序号
    printscript();
    
    %% 嵌套函数
    function order = getLUorder()
        tmpLUMatrix = [da.LUArray.ID; da.LUArray.LWH(1:nDim,:)];
        [~,tepLUorder] = sortrows(tmpLUMatrix',[1 4],{'ascend','descend'}); %1:ID 4:Hight
        % tepLUorder = 1:length(da.LUArray.ID)'; %直接赋值1:n
        % tepLUorder = [2 3 4 1 5]';
        order = tepLUorder';
    end

    function lu = getSortedLU(LUorder)
        lu = structfun(@(x) x(:,LUorder),da.LUArray,'UniformOutput',false);
    end

    function getItem()
        iItem = 1;
        itemBeLUArray = zeros(1,nItem);  % 每个Item内堆垛的LU数量 后期不用
        for iLUid=1:nLUid
            heightLeft = heightBin;
            for iLU=1:nLU
                if sortedLUArray.ID(iLU) == iLUid %仅对当前LU对应Luid在该iLUid内的进行操作
                    if heightLeft < sortedLUArray.LWH(nDim,iLU) %如当前LU高度满足(高度在第nDim行)
                        iItem =  iItem + 1;
                        heightLeft = heightBin;%
                    end
                    heightLeft = heightLeft - sortedLUArray.LWH(nDim,iLU); %更新剩余高度
                    LWHItem(1:2,iItem) = sortedLUArray.LWH(1:2,iLU); %更新item长宽
                    LWHItem(3,iItem) = LWHItem(3,iItem) + sortedLUArray.LWH(nDim,iLU); %更新item高度
                    WeightItem(1,iItem) = WeightItem(1,iItem) + sortedLUArray.Weight(1,iLU); %更新item重量
                    IDItem(1,iItem) = sortedLUArray.ID(1,iLU); %更新ID类型
%                     RotaFlagItem(1,iItem) = sortedLUArray.RotaFlag(1,iLU); %更新RotaFlag类型                    
                    itemBeLUArray(iItem) = itemBeLUArray(iItem) + 1;
                    LUBeItemArraySort(1,iLU) = iItem;
                    LUBeItemArraySort(2,iLU) = itemBeLUArray(iItem);
                end
            end
            iItem =  iItem + 1;
        end
    end

    %%  获取ITEMID类型相关数据(同类型ID的体积，面积，重量)
    function getITEMIDArray()
%          printstruct(da);
        da.ItemIDArray.ID = unique(da.ItemArray.ID);
        nItemID = numel(da.ItemIDArray.ID);        

        for iID = 1:nItemID
        da.ItemIDArray.Weight(iID) = sum(da.ItemArray.Weight .* (da.ItemArray.ID == da.ItemIDArray.ID(iID)) );
        da.ItemIDArray.Volume(iID) = sum(prod(da.ItemArray.LWH) .* (da.ItemArray.ID == da.ItemIDArray.ID(iID)) );
        da.ItemIDArray.Area(iID) = sum(prod(da.ItemArray.LWH(1:2,:)) .* (da.ItemArray.ID == da.ItemIDArray.ID(iID)) );
        end
        
    end

    function printscript()
        for iItem = 1:max(da.LUArray.LUBeItemArray(1,:))
            [~,idx] = find(da.LUArray.LUBeItemArray(1,:)==iItem);
            fprintf('item %d 的长宽高为:  ',iItem);
            fprintf('( %d ) ',da.ItemArray.LWH(:,iItem));
            fprintf('\n');
            fprintf('item %d 包含 original LU 索引号(长宽高)为  \n  ',iItem);
            fprintf('%d ',idx);
            fprintf('( %d ) ', da.LUArray.LWH(:,idx));
            fprintf('\n');
        end
    end
end




%% 后处理 并赋值到da
% 获取LUBeItemArray : 每个排序后LU在哪个Item内  以及顺序
% 获取LWHItem:  新生成的Item的长宽高
% 获取LUorder: LU的排序
% da.ItemArray.LWH = LWHItem(:,LWHItem(1,:)>0); % 去除未使用的
% LUBeItemArray=LUBeItemArraySort;
% LUBeItemArray(:,LUorder) = LUBeItemArraySort;
% da.LUArray.LUorder = LUorder;
% da.LUArray.LUBeItemArray = LUBeItemArray;

% % function outStructArray = SortArrayofStruct( structArray, fieldName )
% %     %UNTITLED2 Summary of this function goes here
% %     %   Detailed explanation goes here
% %     if ( ~isempty(structArray) &&  length (structArray) > 0)
% %       [~,I] = sort(arrayfun (@(x) x.(fieldName), structArray)) ;
% %       outStructArray = structArray(I) ;        
% %     else 
% %         disp ('Array of struct is empty');
% %     end      
% % end
