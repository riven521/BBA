%% 函数1 : isWeightFine: 判断LU是否上轻下重构成
% v3 repairItems
function ITEMSEQ = repairItems(LU)  % t必定是托盘类型的table或struct
    
T = getTableLU(LU);

uniItemID = unique(T.ITEMID(:));
for iItem = 1:length(uniItemID)         %对ITEM进行循环
    
    flagLUIdx = T.ITEMID==uniItemID(iItem); %对单一ItemId去逻辑值;
    
    subT = T(flagLUIdx,{'Weight','ITEMSEQ'});
    
    if isUpDownWeight(subT) 
        T.ITEMSEQ(flagLUIdx) = repairItemWeight(T,flagLUIdx); % 调整t的LU_Item
    end
    
end

ITEMSEQ = T.ITEMSEQ';

% LU = getSturctT(T); %可用不返回结构体LU
end

%% 函数2 : repairItemWeight: LU如不是上轻下重,进行修复
% V3: 不用struct还是用table修改 仅仅修改LU.LU_Item的第二行的值 LU.LU_Item(2,flagLU)
function b = repairItemWeight(T,flagLUIdx)
    tmpWeight=T.Weight(flagLUIdx);
    [~,b] = sort(tmpWeight,'descend');
    [~,b] = sort(b);
    T.ITEMSEQ(flagLUIdx) = b;
end











%% v2 repairItems
% % function LU = repairItems(t)  % t必定是托盘类型的table或struct
% %     chktLU(t)
% % % % 0 结构体转换为table 预处理
% % T = getTableLU(t)
% % T.Properties.VariableNames
% %    t=T
% % 
% % %% 2 相同ITEMID下的CHEK
% % 
% % uniItemID = unique(t.ITEMID(:));
% % % 2.1 如果坐标CoordLUBin存在,用坐标系判断与ITEMID的差异(1:XY值(坐标是否相同); 2:重量和Z值)
% % for iItem = 1:length(uniItemID) %对ITEM进行循环
% %     flagIdx = t.ITEMID==uniItemID(iItem); %对单一ItemId去逻辑值;
% %     if any(strcmp('Y', t.Properties.VariableNames))
% %         vX = t{flagIdx,'X'};
% %         vY = t{flagIdx,'Y'};
% %         if any(vX ~= vX(1)) || any(vY ~= vY(1))
% %             error('相同ITEM,但X或Y坐标错位'); end
% %         
% %         v = t(flagIdx,{'Z','Weight','ITEMSEQ'});
% %     else
% %         v = t(flagIdx,{'Weight','ITEMSEQ'});
% %     end
% %     
% %     v = sortrows(v,'ITEMSEQ');
% %     
% %     if any(strcmp('Y', t.Properties.VariableNames))
% %         % 如果重量不是递减或Z高度不是递增
% %         if ~issorted(v.Z,'ascend') || ~issorted(v.Weight,'descend')
% %             issorted(v.Z,'ascend')  
% %             issorted(v.Weight,'descend')
% %             error('相同ITEM,但重量不是递减或Z高度不是递增'); end
% %         else
% %         if ~issorted(v.Weight,'descend')
% %                                                             %             uniItemID(iItem)
% %                                                             %             v
% %             t = repairItemWeight(t,uniItemID(iItem)); % 调整t的LU_Item
% %                                                             %             v = t(flagIdx,{'Weight','LU_Item'});
% %                                                             %             v = sortrows(v,'LU_Item');
% %         end
% %     end
% % end
% % LU = t;
% % 
% % % 删除初始化获取的列
% % if any(strcmp('ITEMID', LU.Properties.VariableNames))
% %         LU.ITEMID = []; end
% % if any(strcmp('ITEMSEQ', LU.Properties.VariableNames))
% %         LU.ITEMSEQ = []; end
% %     if any(strcmp('X', LU.Properties.VariableNames))
% %         LU.X = []; end
% %     if any(strcmp('Y', LU.Properties.VariableNames))
% %         LU.Y = []; end
% %     if any(strcmp('Z', LU.Properties.VariableNames))
% %         LU.Z = []; end
% %     
% % if istable(LU)
% %     LU = table2struct(LU,'ToScalar',true);
% %     LU = (structfun(@(x) x',LU,'UniformOutput',false));
% % end
% % 
% % end






%%
% % %% 函数2 : repairItemWeight: LU如不是上轻下重,进行修复
% % % V2: 仅仅修改LU.LU_Item的第二行的值 LU.LU_Item(2,flagLU)
% % function LU = repairItemWeight(oLU,itemIdx)
% %     if istable(oLU)
% %         LU = table2struct(oLU,'ToScalar',true);
% %         LU = (structfun(@(x) x',LU,'UniformOutput',false));
% %     else
% %         LU = oLU;
% %     end
% % 
% %     flagLU = (LU.LU_Item(1,:)==itemIdx);   %找出本item对应的lu的index
% %     tmpWeight=LU.Weight(:,flagLU);
% %     [~,b] = sort(tmpWeight,'descend');
% %     [~,b] = sort(b);
% % 
% %     LU.LU_Item(2,flagLU) = b;
% % 
% %     if istable(oLU) %返回的也要是table
% %         LU = struct2table(structfun(@(x) x',LU,'UniformOutput',false));
% %     end
% % end


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

%% V1版 repairItemWeight 先在insert2Item中赋予初始值, 此处再进行修复
% ****************** 对角线开关修复 ************ 关闭
% % % V1: 先在insert2Item中赋予初始值, 此处再进行修复
% % % repairItemFull: 如果存在非满层的case, 进行微调:为0的改为1,当剩余高度小于ITEM对角线的高度, 视为FULL 满层
% % if ~all(Item.isHeightFull)
% %     [~,b] = find(Item.isHeightFull == 0);
% %     for i=1:length(b)
% %            Item = repairItemFull(Item,hVeh,b(i)); %DONE 
% %     end
% % end