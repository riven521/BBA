%% GET STRIP 相关属性
% 5 Strip.nbItem % 整数：冗余值, 具体ITEM的堆垛个数 车头摆放依据 -1：混合strip

%% 函数
function   [Strip] = cpuStripnbItem(Strip,Item,Veh)
Strip.nbItem = ones(size(Strip.Weight))*-1;   %单STRIP内部ITEM类型个数, 混合型默认为-1

%% 5,6,7
%Strip.nbItem: 混合:-1; 单纯: 对应Strip内部该Item的nbLID类型个数,数值越大,即该LU类型越多
LIDinItemsArray = cellfun(@(x) x(1), Item.LID); % arrayAllLID: 所有ITEM对应的LID值 向量形式

uniItem = unique(Item.Item_Strip(1,:));
for i=1:length(Strip.nbItem)
    if ~Strip.isMixed(1,i) %如是单纯型        
        cellLID = Item.LID(Item.Item_Strip(1,:) == uniItem(i)); % cellLID: 本Strip内的ITEM对应的LID值
        LIDinThisItemArray = cellfun(@(x)x(1), cellLID);
        if isscalar(unique(LIDinThisItemArray))
            
            Strip.nbItem(1,i) = sum(LIDinItemsArray == unique(LIDinThisItemArray));
            
        else
             error('单纯型STRIP内的ITEM的类型不同'); %arrayLID            
        end
    end
end
end
