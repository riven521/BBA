function [d] = HLUtoItem(d,ParaArray)
% 重要函数:LU堆垛后形成Item %  行数:长宽高(row);  列数:托盘数量(coloum);
% Input ---  LUArray: ID LWH Weight （LU: 保持原有顺序）
% Output --- LUArray: order LUBeItemArray （LU: 保持原有顺序）(ORDER是进入Item算法的LU顺序)
% Output --- ItemArray: ID LWH Weight（ITEM:没有顺序，算法计算后的顺序）

% d.LUArray.LUBeItemArray (2行*n列: 行1: LU在第几个item 行2:LU进入该item的顺序(底-高)) 
% d.LUArray.order (1行*n列: LU排序顺序)
% d.ItemArray.ID （1行*n列: ITEM的类型 等同内部LU类型)
% d.ItemArray.LWH (3行*n列: ITEM的长宽高-长宽与LU相同,高度为堆垛后高度)
% d.ItemArray.Weight (1行*n列: ITEM的重量)
% 嵌套函数(先排序 后转换)
% getLUorder
% getSortedLU
% getItem
% 

    %% 初始化
% nDim LU维度 nLU LU数量 nItem Item数量 nLUid LU种类 
% heightBin Bin最大高度
nDim = size(d.LUArray.LWH,1);
nLU = size(d.LUArray.LWH,2);
nItem = nLU;
nLUid = size(unique(d.LUArray.ID),2);

tmpUniqueBin = unique(d.BinArray.LWH(1:nDim,:)','rows')';
heightBin = tmpUniqueBin(3);
clear tmpUniqueBin;

    %% LU排序
    % getLUorder - 获取LU的顺序(重点是高度递减排序)
    % getSortedLU - 获取按order排序后的LU:sortedLUArray
    d.LUArray.order = getLUorder(); %获取 LU排序(先ID递增,后高度递减)
    sortedLUArray = getSortedLU(d.LUArray.order);

    %% 55 LU->Item转换
    ItemID = zeros(1,nItem);   %Item的ID类型
    ItemisRota = ones(1,nItem)*2;   %Item的可旋转类型(初始为2)
    
    ItemWeight = zeros(1,nItem);   %Item的重量
    ItemLWH = zeros(nDim,nItem);   %Item的宽长高
    LUBeItemSort = zeros(2,nLU);     %dim1:属于第几个Item dim2:属于该Item第几个排放 555
                % RotaFlagItem = ones(1,nItem)*5;   %Item的旋转flag（0和1,其它为未使用） TOBE DELE    
    
     getItem();     %555 LU到Item的转换
    
    d.ItemArray.ID = ItemID(:,ItemID(1,:)>0);                                    % 去除未使用的  
    d.ItemArray.isRota = ItemisRota(:,ItemisRota(1,:)<=1);               % 去除未使用的     
    d.ItemArray.Weight = ItemWeight(:,ItemWeight(1,:)>0);            % 去除未使用的 
    d.ItemArray.LWH = ItemLWH(:,ItemLWH(1,:)>0);                        % 去除未使用的 
    d.LUArray.LUBeItemArray(:,d.LUArray.order) = LUBeItemSort; % d.LUArray.LUBeItemArray : 每个排序后LU在哪个Item内  以及顺序
                %  d.ItemArray.itemRotaFlag = RotaFlagItem(:,RotaFlagItem(1,:)<=1); %TOBE DELE


    %% 额外变量+测试script
    % 额外变量
    getITEMIDArray();
    % 输出主要结果:获得每个item包含的 原始 LU序号
     printscript();
    
    %% 嵌套函数
    function order = getLUorder()
        tmpLUMatrix = [d.LUArray.ID; d.LUArray.LWH(1:nDim,:)];
        [~,tepLUorder] = sortrows(tmpLUMatrix',[1 4],{'ascend','descend'}); %1:ID 4:Hight
%         tepLUorder = 1:length(d.LUArray.ID)'; %直接赋值1:n
        % tepLUorder = [2 3 4 1 5]';
        order = tepLUorder';
    end

    function lu = getSortedLU(LUorder)
        lu = structfun(@(x) x(:,LUorder),d.LUArray,'UniformOutput',false);
    end

    % 将LU转换为Item的重要函数
    function getItem()
        iItem = 1;
        tmpitemBeLU = zeros(1,nItem);  % 每个Item内堆垛的LU数量 后期不用
        for iLUid=1:nLUid
            heightLeft = heightBin;
            for iLU=1:nLU
                if sortedLUArray.ID(iLU) == iLUid %仅对当前LU对应Luid在该iLUid内的进行操作
                    if heightLeft < sortedLUArray.LWH(nDim,iLU) %如当前LU高度满足(高度在第nDim行)
                        iItem =  iItem + 1;
                        heightLeft = heightBin; %
                    end
                    
                    heightLeft = heightLeft - sortedLUArray.LWH(nDim,iLU); %更新剩余高度
                    ItemLWH(1:2,iItem) = sortedLUArray.LWH(1:2,iLU); %更新item长宽
                    ItemLWH(3,iItem) = ItemLWH(3,iItem) + sortedLUArray.LWH(nDim,iLU); %更新item高度
                    ItemWeight(1,iItem) = ItemWeight(1,iItem) + sortedLUArray.Weight(1,iLU); %更新item重量
                    
                    tmpitemBeLU(iItem) = tmpitemBeLU(iItem) + 1;
                    LUBeItemSort(1,iLU) = iItem;
                    LUBeItemSort(2,iLU) = tmpitemBeLU(iItem);
                    
                    ItemID(1,iItem) = sortedLUArray.ID(1,iLU); %更新ID类型
                    ItemisRota(1,iItem) = sortedLUArray.isRota(1,iLU);  %更新ID可旋转类型
                                        %  RotaFlagItem(1,iItem) = sortedLUArray.RotaFlag(1,iLU); %更新RotaFlag类型                    
                end
            end
            iItem =  iItem + 1;
        end
    end

    %%  获取ITEMID类型相关数据(同类型ID的体积，面积，重量，item是否可旋转)
    function getITEMIDArray()
                %  printstruct(d);  
        d.ItemIDArray.ID = unique(d.ItemArray.ID);
        nItemID = numel(d.ItemIDArray.ID);        

        for iID = 1:nItemID
        d.ItemIDArray.Weight(iID) = sum(d.ItemArray.Weight .* (d.ItemArray.ID == d.ItemIDArray.ID(iID)) );
        d.ItemIDArray.Volume(iID) = sum(prod(d.ItemArray.LWH) .* (d.ItemArray.ID == d.ItemIDArray.ID(iID)) );
        d.ItemIDArray.Area(iID) = sum(prod(d.ItemArray.LWH(1:2,:)) .* (d.ItemArray.ID == d.ItemIDArray.ID(iID)) );
        d.ItemIDArray.isRota(iID) =  unique(d.ItemArray.isRota(d.ItemArray.ID == d.ItemIDArray.ID(iID))); if ~isscalar(d.ItemIDArray.isRota(iID)), error('致命错误'); end
        end
        printstruct(d);
    end

    function printscript()
        for iItem = 1:max(d.LUArray.LUBeItemArray(1,:))
            [~,idx] = find(d.LUArray.LUBeItemArray(1,:)==iItem);
            fprintf('item %d 的长宽高为:  ',iItem);
            fprintf('( %d ) ',d.ItemArray.LWH(:,iItem));
            fprintf('\n');
            fprintf('item %d 包含 original LU 索引号(长宽高)为  \n  ',iItem);
            fprintf('%d ',idx);
            fprintf('( %d ) ', d.LUArray.LWH(:,idx));
            fprintf('\n');
        end
    end
end




%% 后处理 并赋值到da
% 获取LUBeItemArray : 每个排序后LU在哪个Item内  以及顺序
% 获取LWHItem:  新生成的Item的长宽高
% 获取LUorder: LU的排序
% d.ItemArray.LWH = LWHItem(:,LWHItem(1,:)>0); % 去除未使用的
% LUBeItemArray=LUBeItemArraySort;
% LUBeItemArray(:,LUorder) = LUBeItemArraySort;
% d.LUArray.LUorder = LUorder;
% d.LUArray.LUBeItemArray = LUBeItemArray;

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
