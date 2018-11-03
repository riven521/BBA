%% GET ITEM 相关属性
% 1: Item.isHeightFull (对角线>=顶层间隙, 视为满层; Item内最大LU高度 >= 顶层间隙, 视为满层)
% 2: Item.Layer 计算每个Item内堆垛的层数
% 3: Item.isWeightFine 每个Item内是否满足上轻下重
% 4: Item.isNonMixed  计算每个Item是否为不需要混拼的可能
% 5: Item.nbItem 计算每个Item包含同LU/ID的数量
%% 函数
function   [Item,LU] = cpuItem(Item,LU,Veh)
    global ISdiagItem
    %% 初始化
    sz = size(Item.isRota);
    hVeh  = Veh.LWH(3,1);  

    Item.isHeightFull = zeros(sz);    %Item的是否高度满层(初始为0)
%     Item.isWeightFine = ones(sz)*-1;    %Item的是否上轻下重(初始为-1) 修改为直接修复每个Item
    Item.isNonMixed = ones(sz)*-1;    %Item是否非需要混合判定,将偶数个的Item提前进行Strip生成
    Item.isMixedTile = zeros(sz);    %Item混合,找出奇数个混合Item的尾托赋值为1

    %% SECTION 0 计算ITEM的PID,LID,SID
    % 由混合的LU.DOC计算ITEM内包含的PID,LID,SID等数据 1808新增 计算Item.PID,LID,SID等使用
    LU.DOC=[LU.PID;LU.ID;LU.SID;zeros(size(LU.ID));zeros(size(LU.ID));...
        LU.LU_Item;];
    nItem = size(Item.LWH,2);
    for iItem=1:nItem
        tmp = LU.DOC([1,2,3], LU.DOC(6,:) == iItem);
        Item.PID(:,iItem) = num2cell(unique(tmp(1,:))',1); %unique(tmp(1,:))';
        Item.LID(:,iItem) = num2cell(unique(tmp(2,:))',1);
        Item.SID(:,iItem) =num2cell(unique(tmp(3,:))',1);
    end

    %% SECTION 1 计算ITEM的isHeightFull
    % (对角线>=顶层间隙, 视为满层; Item内最大LU高度 >= 顶层间隙, 视为满层)
    for iItem=1:length(Item.isHeightFull)
        % diagItem: ITEM的对角线长度 
        diagItem = sqrt(Item.LWH(1,iItem)^2 + Item.LWH(2,iItem)^2);
        % maxHeightinLUofThisItem: Item内Lu的最高的高度
        flagLU = LU.LU_Item(1,:) == iItem;
        maxHeightinLUofThisItem = max(LU.LWH(3,flagLU));
        % hMargin: ITEM距离车顶的间隙 ************* 1  此方法可能造成pingpuall时,产生不必要的甩尾
        %         hMargin = hVeh - Item.LWH(3,iItem);
        % hMargin: ITEM距离所有Item的最高值的间隙 ************* 2
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
        if ISdiagItem==1
            if maxHeightinLUofThisItem >= hMargin || diagItem >= hMargin, Item.isHeightFull(iItem) = 1;  else  Item.isHeightFull(iItem) = 0; end
        else
            if maxHeightinLUofThisItem >= hMargin,   Item.isHeightFull(iItem) = 1;  else  Item.isHeightFull(iItem) = 0; end
        end
    end

    %% SECTION 2 计算ITEM的isWeightFine并进行修复
    % ****************** 上轻下重的判断+修复 ************ 开放
    % isWeightUpDown: ITEM增加判断是否上轻下重的判断Item.isWeightFine

    % V2版 修复上轻下重——check过程就修复了
    LU  = repairItems(LU);

% V1版 修复上轻下重——
% %     Item = isWeightUpDown(Item,LU);
% %     Item.isWeightFine
% %     % repairItemWeight: 如果存在上轻下重的case, 进行修复
% %     if ~all(Item.isWeightFine)
% %         [~,order] = find(Item.isWeightFine == 0);
% %         for iItem=1:length(order)
% %             LU = repairItemWeight(LU,order(iItem));
% %         end
% %     end
% %     
% %     Item = isWeightUpDown(Item,LU);
% %     Item.isWeightFine
% %     checktLU(LU)
% %     if ~all(Item.isWeightFine),   error('仍有上轻下重casse, 错误'); end

    %% SECTION 3 计算ITEM的isNonMixed/isMixedTile是否为不需要混拼/混拼排序找甩尾计算
    % ****************** Iten内是否为不需要混拼计算 ************ 开放
    % GET Item.isNonMixed: 计算每个Item是否为不需要混拼的可能
    ItemLID = cellfun(@(x) x(1), Item.LID); % arrayAllLID: 所有ITEM对应的LID值 向量形式
    uniItemLID = unique(ItemLID);
    % 循环: LID个数
    for iItem=1:length(uniItemLID)
        % Item i 对于的LU flag标记
        flagItem = ItemLID(:) == iItem;  % flagItem = ItemLID(:) == uniItemLID(iItem);  idExchange(LU.ID)必须有,否则改为后者，但影响巨大 
        
        ItemWidth = unique(Item.LWH(1,flagItem));     % Item宽度    
        VehWidth = Veh.LWH(1,1);  % 车辆宽度    
        maxWidthLayer= floor(VehWidth/ItemWidth); %Item可放宽度层数
        nb = sum(flagItem);
        nbmod = mod(nb,maxWidthLayer);
        if nb ==0 || nbmod>nb, error('cpuItem种计算isNonMixed错误'); end
        if nbmod == 0 %mod为0表明 不需要混合 不混合的提前在order中提前
            Item.isNonMixed(flagItem) = 1;
        else
            Item.isNonMixed(flagItem)= 0;
            % 计算Item的isMixedTile
            tmpSort=[Item.isHeightFull;Item.HLayer;Item.LWH(3,:)];
            [~, order]=sortrows(tmpSort(:,flagItem)', [1,2,3],{'ascend','ascend','ascend'});
            flagItemIdx = find(flagItem);
            flagmodIdx = flagItemIdx(order(1:nbmod));
            Item.isMixedTile(flagmodIdx)=1;
        end
    end
    
    %% SECTION 4 计算ITEM的nbItem
    tmpItemLID = cell2mat(Item.LID);
    % V1
    %     for i=1:length(Item.Weight)
    %         Item.nbItem(i) = sum(tmpItemLID == tmpItemLID(i));
    %     end
    % V2
    Item.nbItem = sum(tmpItemLID == tmpItemLID');

end

%% 局部函数 %%

%% 函数1 : isWeightFine: 判断LU是否上轻下重构成
%% v2 repairItems
function LU = repairItems(t)  % t必定是托盘类型的table或struct
% 0 结构体转换为table 预处理
if isstruct(t)
    t = struct2table(structfun(@(x) x', t,'UniformOutput',false));  end

% 1.1 初步核查LU_Item 如不存在获取子ITEMID 
if any(strcmp('LU_Item', t.Properties.VariableNames))
        t.ITEMID = t.LU_Item(:,1);
        t.ITEMSEQ = t.LU_Item(:,2); end

if any(strcmp('CoordLUBin', t.Properties.VariableNames))
    t.X = t.CoordLUBin(:,1);    t.Y = t.CoordLUBin(:,2);   t.Z = t.CoordLUBin(:,3);   end
    
%% 2 相同ITEMID下的CHEK
uniItemID = unique(t.ITEMID(:));
% 2.1 如果坐标CoordLUBin存在,用坐标系判断与ITEMID的差异(1:XY值(坐标是否相同); 2:重量和Z值)
for iItem = 1:length(uniItemID) %对ITEM进行循环
    flagIdx = t.ITEMID==uniItemID(iItem); %对单一ItemId去逻辑值;
    if any(strcmp('Y', t.Properties.VariableNames))
        vX = t{flagIdx,'X'};
        vY = t{flagIdx,'Y'};
        if any(vX ~= vX(1)) || any(vY ~= vY(1))
            error('相同ITEM,但X或Y坐标错位'); end
        
        v = t(flagIdx,{'Z','Weight','ITEMSEQ'});
    else
        v = t(flagIdx,{'Weight','ITEMSEQ'});
    end
    
    v = sortrows(v,'ITEMSEQ');
    
    if any(strcmp('Y', t.Properties.VariableNames))
        % 如果重量不是递减或Z高度不是递增
        if ~issorted(v.Z,'ascend') || ~issorted(v.Weight,'descend')
            issorted(v.Z,'ascend');  issorted(v.Weight,'descend');
            error('相同ITEM,但重量不是递减或Z高度不是递增'); end
        else
        if ~issorted(v.Weight,'descend')
                                                            %             uniItemID(iItem)
                                                            %             v
            t = repairItemWeight(t,uniItemID(iItem)); % 调整t的LU_Item
                                                            %             v = t(flagIdx,{'Weight','LU_Item'});
                                                            %             v = sortrows(v,'LU_Item');
        end
    end
end
LU = t;
if istable(LU)
    LU = table2struct(LU,'ToScalar',true);
    LU = (structfun(@(x) x',LU,'UniformOutput',false));
end

end

%% V1版 isWeightUpDown 无法找全所有上轻下重的case
% % function Item = isWeightUpDown(Item,LU)
% % for iItem = 1:max(LU.LU_Item(1,:)) %对ITEM进行循环
% %     [~,idx] = find(LU.LU_Item(1,:)==iItem);
% %     if iItem==21
% %         1
% %     end
% %     if isempty(idx), 
% %         1
% %     end
% %     nbLUinItem = length(idx);    
% %     % 对ITME内含2个以上LU的进行判断
% %     if nbLUinItem > 1 %Item包含不只一个Item,需要判断是否有轻重的变化
% %         currLUWeight = zeros(1,nbLUinItem);
% %         for iIdx = 1:nbLUinItem
% %             currIdx = idx(LU.LU_Item(2,idx) == iIdx);
% %             currLUWeight(iIdx) = LU.Weight(:,currIdx);          %                 currLUHight(iIdx) = LU.LWH(3,currIdx);
% %         end
% %         diff(currLUWeight)
% %         all(diff(currLUWeight) > 0) 
% %         if all(diff(currLUWeight) > 0) % 代表下轻上重
% %             % 修改LU.LU_Item的值             1 5: idx  为 1 2: LU.LU_Item(2,idx) == iIdx改为 2 1
% %             Item.isWeightFine(1,iItem) = 0;
% %         else
% %             Item.isWeightFine(1,iItem) = 1;
% %         end
% %     else  %ITEM内只有1个LU, 必定满足条件
% %         Item.isWeightFine(1,iItem) = 1; 
% %     end
% % end
% % end

%% 函数2 : repairItemWeight: LU如不是上轻下重,进行修复
% V2: 仅仅修改LU.LU_Item的第二行的值 LU.LU_Item(2,flagLU)
function LU = repairItemWeight(oLU,itemIdx)
if istable(oLU)
    LU = table2struct(oLU,'ToScalar',true);
    LU = (structfun(@(x) x',LU,'UniformOutput',false));
else
    LU = oLU;
end

flagLU = (LU.LU_Item(1,:)==itemIdx);   %找出本item对应的lu的index
tmpWeight=LU.Weight(:,flagLU);
[~,b] = sort(tmpWeight,'descend');
[~,b] = sort(b);

LU.LU_Item(2,flagLU) = b;

if istable(oLU) %返回的也要是table
    LU = struct2table(structfun(@(x) x',LU,'UniformOutput',false));
end
end


% ****************** 对角线开关修复 ************ 关闭
% % % V1: 先在insert2Item中赋予初始值, 此处再进行修复
% % % repairItemFull: 如果存在非满层的case, 进行微调:为0的改为1,当剩余高度小于ITEM对角线的高度, 视为FULL 满层
% % if ~all(Item.isHeightFull)
% %     [~,b] = find(Item.isHeightFull == 0);
% %     for i=1:length(b)
% %            Item = repairItemFull(Item,hVeh,b(i)); %DONE 
% %     end
% % end