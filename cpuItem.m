%% GET ITEM 相关属性
% 1: Item.isHeightFull (对角线>=顶层间隙, 视为满层; Item内最大LU高度 >= 顶层间隙, 视为满层)
% 2: Item.Layer 计算每个Item内堆垛的层数
% 3: Item.isWeightFine 每个Item内是否满足上轻下重
% 4: Item.isNonMixed  计算每个Item是否为不需要混拼的可能

%% 函数
function   [Item,LU] = cpuItem(Item,LU,Veh)
    %% 初始化
    sz = size(Item.isRota);
    hVeh  = Veh.LWH(3,1);  

    Item.isHeightFull = ones(sz)*-1;    %Item的是否高度满层(初始为-1)
    Item.isWeightFine = ones(sz)*-1;    %Item的是否上轻下重(初始为-1)
    Item.isNonMixed = ones(sz)*-1;    %Item是否非需要混合判定,将偶数个的Item提前进行Strip生成

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
        % hMargin: ITEM距离车顶的间隙
        hMargin = hVeh - Item.LWH(3,iItem);
        
        if diagItem >= hMargin,  Item.isHeightFull(iItem) = 1;  else Item.isHeightFull(iItem) = 0; end
        if maxHeightinLUofThisItem >= hMargin,   Item.isHeightFull(iItem) = 1;  else  Item.isHeightFull(iItem) = 0; end    
    end
       
    %% SECTION 2 计算ITEM的isWeightFine并进行修复
    % ****************** 上轻下重的判断+修复 ************ 开放
    % isWeightUpDown: ITEM增加判断是否上轻下重的判断Item.isWeightFine
    Item = isWeightUpDown(Item,LU);
    % repairItemWeight: 如果存在上轻下重的case, 进行修复
    if ~all(Item.isWeightFine)
        [~,b] = find(Item.isWeightFine == 0);
        for iItem=1:length(b)
            LU = repairItemWeight(LU,b(iItem));
        end
    end
    Item = isWeightUpDown(Item,LU);
    if ~all(Item.isWeightFine),   error('仍有上轻下重casse, 错误'); end

    %% SECTION 3 计算ITEM的isNonMixed是否为不需要混拼计算
    % ****************** Iten内是否为不需要混拼计算 ************ 开放
    % GET Item.isNonMixed: 计算每个Item是否为不需要混拼的可能
    ItemLID = cellfun(@(x) x(1), Item.LID); % arrayAllLID: 所有ITEM对应的LID值 向量形式
    for iItem=1:length(unique(ItemLID))
        % Item i 对于的LU flag标记
        flagItem = ItemLID(:) == iItem;
        ItemWidth = unique(Item.LWH(1,flagItem));     % Item宽度    
        VehWidth = Veh.LWH(1,1);  % 车辆宽度    
        maxWidthLayer= floor(VehWidth/ItemWidth); %Item可放宽度层数
        nb = sum(flagItem);
        if mod(nb,maxWidthLayer) == 0 %mod为0表明 不需要混合 不混合的提前在order中提前
            Item.isNonMixed(flagItem) = 1;
        else
            Item.isNonMixed(flagItem)= 0;
        end
    end
end

%% 局部函数 %%

%% 函数1 : isWeightFine: 判断LU是否上轻下重构成
function Item = isWeightUpDown(Item,LU)
for iItem = 1:max(LU.LU_Item(1,:)) %对ITEM进行循环
    [~,idx] = find(LU.LU_Item(1,:)==iItem);
    nbLUinItem = length(idx);
    % 对ITME内含2个以上LU的进行判断
    if length(idx) > 1 %Item包含不只一个Item,需要判断是否有轻重的变化
        currLUWeight = zeros(1,nbLUinItem);
        for iIdx = 1:nbLUinItem
            currIdx = idx(LU.LU_Item(2,idx) == iIdx);
            currLUWeight(iIdx) = LU.Weight(:,currIdx);          %                 currLUHight(iIdx) = LU.LWH(3,currIdx);
        end
        if diff(currLUWeight) > 0 % 代表下轻上重
            % 修改LU.LU_Item的值             1 5: idx  为 1 2: LU.LU_Item(2,idx) == iIdx改为 2 1
            Item.isWeightFine(1,iItem) = 0;
        else
            Item.isWeightFine(1,iItem) = 1;
        end
    else  %ITEM内只有1个LU, 必定满足条件
        Item.isWeightFine(1,iItem) = 1; 
    end
end
end

%% 函数2 : repairItemWeight: LU如不是上轻下重,进行修复
function LU = repairItemWeight(LU,itemIdx)
    [~,LUidx] = find(LU.LU_Item(1,:)==itemIdx); %找出本item对应的lu的index
    nbLUinItem = length(LUidx);
    currLUWeight = zeros(1,nbLUinItem);
    for iIdx = 1:nbLUinItem
        currIdx = LUidx(LU.LU_Item(2,LUidx) == iIdx);
        currLUWeight(iIdx) = LU.Weight(:,currIdx);
    end
    % 对Item内的LU进行排序获得b; 将进入顺序LU_Item的顺序2进行修复
    [~,b] = sort(currLUWeight,'descend');
    tt = LU.LU_Item(2,LUidx);
    LU.LU_Item(2,LUidx) = tt(b);
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