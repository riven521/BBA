%% GET ITEM 相关属性
% 1: Item.isHeightFull ( Item内最大LU高度 * 0.95 >= 顶层间隙, 视为满层)  computeisHeightFullItem
% 4: Item.isNonMixed  计算每个Item是否为不需要混拼的可能
% 5: Item.nbItem 计算每个Item包含同LU/ID的数量
%% 函数
function   [Item] = cpuItem(Item,LU,Veh)

    %% 初始化
    sz = size(Item.Weight);
    % hVeh  = Veh.LWH(3,1);  % 暂时为用

    Item.isHeightFull = zeros(sz);    %Item的是否高度满层(初始为0)  %  Item.isWeightFine = ones(sz)*-1;    %Item的是否上轻下重(初始为-1) 修改为直接修复每个Item
    Item.MixOrder = ones(sz)*-1;    %1:ITEM可以完美放到所有strip,无需和其它ITEM混合; 0: 有必要混合; 将偶数个的Item提前进行Strip生成
        %Item.isMixedTile = zeros(sz);          % 0:当isNonMixed=1,全部为0;当isNonMixed=0时,部分余下必须混的尾垛ITEM赋值值为1. 
    Item.nbItem = zeros(sz);  
    
    %% V2 SECTION 0 计算ITEM的PID,LID,SID,由TABLE计算,方便知道什么是什么,不用1,2,3数字替换
    % 由混合的LU.DOC计算ITEM内包含的PID,LID,SID等数据 1808新增 计算Item.PID,LID,SID等使用
    T = struct2table(structfun(@(x) x',LU,'UniformOutput',false));
    
    nItem = size(Item.LWH,2);
    for iLuID=1:nItem
        f = T.LU_Item(:,1) == iLuID;                          %  f = logical(ones(height(t),1))
        Item.LID(:,iLuID) = {unique(T.ID(f))};             % NOTE: ITEM里的LID是LU的ID
        %         Item.LID(:,iItem) = {unique(t.LID(f))};
        Item.SID(:,iLuID) = {unique(T.SID(f))};
        Item.EID(:,iLuID) = {unique(T.EID(f))};
        Item.PID(:,iLuID) = {unique(T.PID(f))};
    end    
    %  t2 = struct2table(structfun(@(x) x',Item,'UniformOutput',false));


    %% SECTION 1 计算ITEM的isHeightFull
    Item.isHeightFull = computeisHeightFullItem(Item,LU);
    
    %% 	SECTION 2 Item.isNonMixed: Item排序依据
    [Item.MixOrder] = computeMixandTile(Item,Veh);  
    
    %% SECTION 3 计算ITEM的nbItem
    tmpItemLID = cell2mat(Item.LID);
    % V1
    %     for i=1:length(Item.Weight)
    %         Item.nbItem(i) = sum(tmpItemLID == tmpItemLID(i));
    %     end
    % V2
    Item.nbItem = sum(tmpItemLID == tmpItemLID');

end

%% 局部函数 %%
function TF = computeisHeightFullItem(Item,LU)
    global ISdiagItem
    
    n = length(Item.Weight);
    TF = deal(ones(1,n)*-1);    
    
    % (对角线>=顶层间隙, 视为满层; Item内最大LU高度 >= 顶层间隙, 视为满层)
    for iItem=1:n
        % ItemDiag : ITEM的对角线长度 
        ItemDiag = sqrt(Item.LWH(1,iItem)^2 + Item.LWH(2,iItem)^2);
        
        % maxHeightinLUofThisItem: Item内Lu的最高的高度
        flagLU = LU.LU_Item(1,:) == iItem;
        maxHeightinLUofThisItem = max(LU.LWH(3,flagLU));
        
        % hMargin: ITEM距离车顶的间隙 ************* V1  此方法可能造成pingpuall时,产生不必要的甩尾
        %         hMargin = hVeh - Item.LWH(3,iItem);
        
        % hMargin: 该ITEM距离所有车的所有Item堆垛的最高值的间隙(差距) ************* V2
        hMargin = max(Item.LWH(3,:)) - Item.LWH(3,iItem);
        
        %                 if abs(maxHeightinLUofThisItem - hMargin ) <=60
        %                 hMargin
        %                 maxHeightinLUofThisItem
        %                 end
            % V1: 相互冲突
                %         if ISdiagItem==1 && diagItem >= hMargin,  
                %             Item.isHeightFull(iItem) = 1;  else Item.isHeightFull(iItem) = 0; end
                %         if maxHeightinLUofThisItem >= hMargin,   Item.isHeightFull(iItem) = 1;  else  Item.isHeightFull(iItem) = 0; end    
            
        % V2: 任一满足, 均为满层
        % 条件1: 本堆垛内最大LU高度 > 本堆垛与最高堆垛的高度差 -> 满层 || 堆垛对角线 >  本堆垛与最高堆垛的高度差
        % 条件2: 本堆垛内最大LU高度*0.95 > 本堆垛与最高堆垛的高度差 %参数为2
        if ISdiagItem==1
            if maxHeightinLUofThisItem >= hMargin || ItemDiag >= hMargin, TF(iItem) = 1;  else  TF(iItem) = 0; end
        else
            if 0.95*maxHeightinLUofThisItem >= hMargin,   TF(iItem) = 1;  else  TF(iItem) = 0; end
        end
        
    end
    
end   


% v3：Item排序指标改为1个, nm：重要指标。
function [nm] = computeMixandTile(Item,Veh) 

    LuIDs = cellfun(@(x) x(1), Item.LID);   % arrayAllLID: 所有ITEM对应的LID值 向量形式
    uniLuLID = unique(LuIDs);
    
    n = length(LuIDs);
    [nm] = deal(zeros(1,n));
    
    % 循环: 托盘的ID/类型个数
    for iLuID=1:length(uniLuLID)
       
        % 1 计算nbmod（是否需要混合的重要指标）
        % 类型1 对于的LU flag标记
        flagItem = LuIDs ==  uniLuLID(iLuID);  % flagItem = ItemLID(:) == uniItemLID(iItem);  idExchange(LU.ID)必须有,否则改为后者，但影响巨大 
        
        ItemWidth = unique(Item.LWH(1,flagItem));     % Item宽度
        VehWidth = Veh.LWH(1,1);  % 车辆宽度        
        nbmaxItem= floor(VehWidth/ItemWidth);       % Item可放最大宽度层数
        
        nbItem = sum(flagItem);                                     % Item 数量
        
        % nbmod: 自身宽度数量 对 最大宽度数 取余 （为0，则自身无与其他类型堆垛的可能性；1 则必然 会出现混合strip）
        nbmod = mod(nbItem,nbmaxItem);        if nbItem ==0 || nbmod>nbItem, error('cpuItem种计算isNonMixed错误'); end
        
        % GapWidth = VehWidth - nbmaxItem*ItemWidth; % GapWidth 车辆安排最多个数后，剩余Item后的Gap宽度 （暂时不用）
        
        % 2 若正好为0，即刚好整数个条带
        if nbmod==0
            % 若余数为0，且只有一层Strip，降低其优先级，因为可能其宽度不满层，允许与其他混合（也可能通过gap调整改变）
            if nbItem==nbmaxItem                                        % && GapWidth > ItemWidth*0.5
                nm(flagItem)  = 2;
            else                                                                         % 否则指定为不允许和其它Item混合
                nm(flagItem)  = 1;
            end            
        elseif nbmod~=0            
            nm(flagItem)  = 1;  %先默认为1，继而对单独的必须混合的标记为3，排序放最后；没有2
            
            tmpSort=[Item.isHeightFull;Item.HLayer;Item.LWH(3,:)];
            [~, order]=sortrows(tmpSort(:,flagItem)', [1,2,3],{'ascend','ascend','ascend'});
            flagItemIdx = find(flagItem);
            flagmodIdx = flagItemIdx(order(1:nbmod));
            nm(flagmodIdx)=3;  %NOTE: flagmodIdx 部分Item标记
        end
        
    end
end

%% v2： v3：改2个参数为1个
% % function [nm, mt] = computeMixandTile(Item,Veh) 
% % 
% %     
% %     LuIDs = cellfun(@(x) x(1), Item.LID);   % arrayAllLID: 所有ITEM对应的LID值 向量形式
% %     uniLuLID = unique(LuIDs);
% %     
% %     n = length(LuIDs);
% %     [nm, mt] = deal(zeros(1,n));
% %     
% %     % 循环: 托盘的ID/类型个数
% %     for iLuID=1:length(uniLuLID)
% %        
% %         % 类型1 对于的LU flag标记
% %         flagItem = LuIDs ==  uniLuLID(iLuID);  % flagItem = ItemLID(:) == uniItemLID(iItem);  idExchange(LU.ID)必须有,否则改为后者，但影响巨大 
% %         
% %         ItemWidth = unique(Item.LWH(1,flagItem));     % Item宽度
% %         VehWidth = Veh.LWH(1,1);  % 车辆宽度        
% %         nbmaxItem= floor(VehWidth/ItemWidth);       % Item可放最大宽度层数
% %         
% %         nbItem = sum(flagItem);                                     % Item 数量
% %         
% %         % nbmod: 自身宽度数量 对 最大宽度数 取余 （为0，则自身无与其他类型堆垛的可能性；1 则必然 会出现混合strip）
% %         nbmod = mod(nbItem,nbmaxItem);        if nbItem ==0 || nbmod>nbItem, error('cpuItem种计算isNonMixed错误'); end
% %         
% %         % GapWidth 车辆安排最多个数后，剩余Item后的Gap宽度 （暂时不用）
% %         GapWidth = VehWidth - nbmaxItem*ItemWidth;
% %         
% %         % 计算：flagNonMixedItem（1 ：该类型Item不允许混合的，优先排序；0：非优先排序，与其它混合堆垛可能性大）
% %         % 余数不为0，则一定表明可以有空，必须和其它Item混合；
% %         if nbmod~=0 
% %             flagNonMixedItem = 0;
% %         else
% %             % 若余数为0，且只有一层Strip，剩余Gap宽度又较大，则指定为允许和其它Item混合；
% %             if nbItem==nbmaxItem && GapWidth > ItemWidth*0.5
% %                 flagNonMixedItem = 0; % flagNonMixedItem = 0;
% %             else % 否则指定为不允许和其它Item混合
% %                 flagNonMixedItem = 1;
% %             end
% %         end
% %         
% %         % 赋值：依据flagNonMixedItem，赋值Item的isNonMixed/isMixedTile （Item排序的重要依据）
% %         if   flagNonMixedItem 
% %             nm(flagItem) = 1;  % 如果Item为无需混合的, 则混合甩尾的标记均为0
% %             mt(flagItem) = 0;
% %         else  
% %             nm(flagItem)= 0;
% %             
% %             % mt  计算Item的isMixedTile ，即找出不允许混合的Items，并予以标记
% %             if nbmod == 0 % 即 nbItem==nbmaxItem && GapWidth > ItemWidth*0.5 其实不允许mix，但Gap太大，且只有一层
% %                 mt(flagItem)=1;  %NOTE: flagItem所有的Item都标记
% %             else            
% %                 tmpSort=[Item.isHeightFull;Item.HLayer;Item.LWH(3,:)];
% %                 [~, order]=sortrows(tmpSort(:,flagItem)', [1,2,3],{'ascend','ascend','ascend'});
% %                 flagItemIdx = find(flagItem);
% %                 flagmodIdx = flagItemIdx(order(1:nbmod));
% %                 mt(flagmodIdx)=2;  %NOTE: flagmodIdx 部分Item标记
% %             end
% %         end
% %     end
% % end


    %% v1: computeMixandTile
% % % %     %% 	SECTION 2 Item.isNonMixed: and  Item.isMixedTile
% % % %     %     [Item.isNonMixed, Item.isMixedTile] = computeMixandTile(Item) todo 改到一个函数内
% % % %     % ****************** Iten内是否为不需要混拼计算 ************ 开放
% % % %     % GET Item.isNonMixed: 计算每个Item是否为不需要混拼的可能
% % % %     ItemLID = cellfun(@(x) x(1), Item.LID); % arrayAllLID: 所有ITEM对应的LID值 向量形式
% % % %     uniItemLID = unique(ItemLID);
% % % %     % 循环: LID个数
% % % %     for iItem=1:length(uniItemLID)
% % % %         % Item i 对于的LU flag标记
% % % %         flagItem = ItemLID(:) == iItem;  % flagItem = ItemLID(:) == uniItemLID(iItem);  idExchange(LU.ID)必须有,否则改为后者，但影响巨大 
% % % %         
% % % %         ItemWidth = unique(Item.LWH(1,flagItem));     % Item宽度
% % % %         VehWidth = Veh.LWH(1,1);  % 车辆宽度
% % % %         
% % % %         nbmaxItem= floor(VehWidth/ItemWidth);     % Item可放最大宽度层数
% % % %         nbItem = sum(flagItem);
% % % %         
% % % %         % GapWidth 车辆安排最多个数后，剩余Item后的Gap宽度
% % % %         GapWidth = VehWidth - nbmaxItem*ItemWidth;
% % % %         
% % % %         % nbmod: 自身宽度数量 对 最大宽度数 取余 
% % % %         nbmod = mod(nbItem,nbmaxItem);        if nbItem ==0 || nbmod>nbItem, error('cpuItem种计算isNonMixed错误'); end
% % % %         
% % % %         % 计算：flagNonMixedItem（1 ：该类型Item不允许混合的，优先排序；0：非优先排序，与其它混合堆垛可能性大
% % % %         % 余数不为0，则一定表明可以有空，必须和其它Item混合；
% % % %         if nbmod~=0 
% % % %             flagNonMixedItem = 0;
% % % %         else
% % % %             % 若余数为0，且只有一层Strip，剩余Gap宽度又较大，则指定为允许和其它Item混合；
% % % %             if nbItem==nbmaxItem && GapWidth > ItemWidth*0.5
% % % %                 flagNonMixedItem = 0;
% % % %             else % 否则指定为不允许和其它Item混合
% % % %                 flagNonMixedItem = 1;
% % % %             end
% % % %         end
% % % %         
% % % %         % 赋值：依据flagNonMixedItem，赋值Item的isNonMixed/isMixedTile （Item排序的重要依据）
% % % %         if   flagNonMixedItem 
% % % %             Item.isNonMixed(flagItem) = 1;  % 如果Item为无需混合的, 则混合甩尾的标记均为0
% % % %             Item.isMixedTile(flagItem) = 0;
% % % %         else  
% % % %             Item.isNonMixed(flagItem)= 0;
% % % %             
% % % %             % Item.isMixedTile  计算Item的isMixedTile ，即找出不允许混合的Items，并予以标记
% % % %             if nbmod == 0 % 即 nbItem==nbmaxItem && GapWidth > ItemWidth*0.5 其实不允许mix，但Gap太大，且只有一层
% % % %                 Item.isMixedTile(flagItem)=1;  %NOTE: flagItem所有的Item都标记
% % % %             else            
% % % %                 tmpSort=[Item.isHeightFull;Item.HLayer;Item.LWH(3,:)];
% % % %                 [~, order]=sortrows(tmpSort(:,flagItem)', [1,2,3],{'ascend','ascend','ascend'});
% % % %                 flagItemIdx = find(flagItem);
% % % %                 flagmodIdx = flagItemIdx(order(1:nbmod));
% % % %                 Item.isMixedTile(flagmodIdx)=1;  %NOTE: flagmodIdx 部分Item标记
% % % %             end
% % % %         end
% % % %     end


    %% V1 : 由LU.DOC计算
% %     %% SECTION 0 计算ITEM的PID,LID,SID
% %     % 由混合的LU.DOC计算ITEM内包含的PID,LID,SID等数据 1808新增 计算Item.PID,LID,SID等使用
% %     LU.DOC=[LU.PID;LU.ID;LU.SID;zeros(size(LU.ID));zeros(size(LU.ID));...
% %         LU.LU_Item;];
% %     nItem = size(Item.LWH,2);
% %     for iItem=1:nItem
% %         tmp = LU.DOC([1,2,3], LU.DOC(6,:) == iItem);
% %         Item.PID(:,iItem) = num2cell(unique(tmp(1,:))',1); %unique(tmp(1,:))';
% %         Item.LID(:,iItem) = num2cell(unique(tmp(2,:))',1);
% %         Item.SID(:,iItem) =num2cell(unique(tmp(3,:))',1);
% %     end
    