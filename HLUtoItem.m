function [LU,Item,ItemID] = HLUtoItem(LU,Veh)
% 重要函数:LU堆垛后形成Item %  行数:长宽高(row);  列数:托盘数量(coloum);
% Input ---  LU: ID LWH Weight （LU: 保持原有顺序）
% Output --- LU: order LU_Item （LU: 保持原有顺序）(ORDER是进入Item算法的LU顺序)
% Output --- Item: ID LWH Weight ...（ITEM:没有顺序，算法计算后的顺序）

% LU.LU_Item           (2,n) : 行1: LU在第几个item 行2:LU进入该item的顺序(底-高))
% LU.order                (1,n):  LU排序顺序)
% Item.LID                （1,n): ITEM的类型 等同内部LU类型)
% Item.LWH              (3,n): ITEM的长宽高-长宽与LU相同,高度为堆垛后高度)
% Item.Weight          (1,n): ITEM的重量
% tmpItem_LU         (1,n): 行1 ITEM的LU数量

%% LU排序
% 获取LU的顺序(重点是高度递减排序)
[LU.order]  = getLUorder(LU); %获取 LU排序(先ID递增,后高度递减)
% printstruct(LU)
% 获取按order排序后的LU:sLU
if isSameCol(LU)
    sLU = structfun(@(x) x(:,LU.order),LU,'UniformOutput',false);
else
    error('不能使用structfun');
end

% printstruct(sLU)

%% 55 LU->Item转换

% 排序后的sLU, 经过堆垛获取垛Item,以及sLU在Item内的顺序
sz = size(sLU.ID);
nLU = sz(2);
hVeh  = Veh.LWH(3,1);  % tmpUniqueBin = unique(Veh.LWH(1:3,:)','rows')'; % hVeh = tmpUniqueBin(3);

% 仅需初始化需要自增的fields
    Item.LID = zeros(sz);            %Item的ID类型
    Item.SID = zeros(sz);
    Item.UID = zeros(sz);
    Item.PID = zeros(numel(unique(LU.PID)),sz(2));
    
    Item.isRota = ones(sz)*-1;    %Item的可旋转类型(初始为2)
    Item.Rotaed = ones(sz)*-1;
Item.LWH = zeros(3,nLU); % Item.LWH(1,:) = wStrip;   %dim1-宽度剩余  % Item.LWH(3,:) = hVeh; % 
Item.Weight = zeros(1,nLU); %Item的重量
% 临时使用
tmpItem_LU = zeros(1,nLU);  % 行1：每个Strip内的Item数量 ； 行2：每个Strip内的不同LUID数量
% sLU新增
sLU.LU_Item = zeros(2,sz(2));     %dim1:属于第几个Item dim2:属于该Item第几个排放

iItem = 1; iLU = 1; %iStrip代表item实质
while 1
    if iLU > nLU, break; end
    [thisItem,iItem] = getThisItem(iItem);
    insertLUToItem(thisItem,iLU);       
    iLU = iLU + 1;
end

% Get ITEM 务必可以放
    function [thisItem,iItem] = getThisItem(iItem)
        % 同样SID/UID 同样LUID Item高度满足 未考虑Weight等
        isflagCurr =hVeh - Item.LWH(3,iItem) >= sLU.LWH(3,iLU); %判断是否current's item剩余宽度 >= 当前iLU高度
        isSameID = Item.LID(iItem) == sLU.ID(iLU); %判断Item内部ID是否=当前iLU的ID
        isNewItem = Item.LWH(3,iItem) == 0; % 判断是否 new Item 高度==0            
        
        if isNewItem
                thisItem = iItem;
        else
            if isflagCurr && isSameID %如果高度允许 且 ID相同
                 thisItem = iItem;
            else
                iItem = iItem + 1;
                [thisItem,iItem] = getThisItem(iItem);
            end
        end
    end

% Put LU into thisItem
    function insertLUToItem(thisItem,iLU)
        %更新Height
        Item.LWH(3,thisItem) = Item.LWH(3,thisItem)  + sLU.LWH(3,iLU);
        Item.LWH(1:2,thisItem) = sLU.LWH(1:2,iLU);  %更新item长宽
        Item.Weight(1,thisItem) = Item.Weight(1,thisItem) + sLU.Weight(1,iLU); %更新item重量
        
        tmpItem_LU(1,thisItem) = tmpItem_LU(1,thisItem) + 1;
        
        sLU.LU_Item(1,iLU) = thisItem;
        sLU.LU_Item(2,iLU) = tmpItem_LU(1,thisItem);
        
                    %         tmpLUThisItem = sLU.LU_Item(1,:) == thisItem;
                    %         tmpItem_LU(2,iItem) = numel(unique(sLU.PID(1,tmpLUThisItem)));
        
        Item.LID(1,thisItem) = sLU.ID(1,iLU);               %更新ID类型
        Item.SID(1,thisItem) = sLU.SID(1,iLU);
        Item.UID(1,thisItem) = sLU.UID(1,iLU);
        Item.isRota(1,thisItem) = sLU.isRota(1,iLU);  %更新ID可旋转类型
        Item.Rotaed(1,thisItem) = sLU.Rotaed(1,iLU);  %更新ID旋转标记        
        

                        %         Item.PID(sLU.PID(1,iLU),thisItem) = Item.PID(sLU.PID(1,iLU),thisItem) + 1;         % 555 更新多行PID - 数值为出现次数
         Item.PID(sLU.PID(1,iLU),thisItem) =     1;                             % 555 更新多行PID - 数值为出现与否
        
    end

% LU内部更新,sLU依据order变化回来
if isSameCol(sLU)
    LU = reorderStruct(LU.order, sLU);
else
    error('不能使用structfun');
end

% Item去除未使用 %     Item.Rotaed(:,Item.itemorder) = sLU.Rotaed;
% 如果ITEM的列数全部相同
if isSameCol(Item)
    Item = structfun(@(x) x( : , Item.LID(1,:)>0 ), Item, 'UniformOutput', false);
else
    error('不能使用structfun');
end

% 额外变量 ItemID
% ItemID = getITEMIDArray(Item);
 ItemID = [];
%% 测试script TO BE FIX
% 输出主要结果:获得每个item包含的 原始 LU序号
printscript(LU,Item);

end

function printscript(LU,Item)
    for iItem = 1:max(LU.LU_Item(1,:))
        [~,idx] = find(LU.LU_Item(1,:)==iItem);
        fprintf('item %d 的长宽高为:  ',iItem);
        fprintf('( %d ) ',Item.LWH(:,iItem));
        fprintf('\n');
        fprintf('item %d 包含 original LU 索引号(长宽高)为  \n  ',iItem);
        fprintf('%d ',idx);
        fprintf('( %d ) ', LU.LWH(:,idx));
        fprintf('\n');
    end
end

function [tepLUorder] = getLUorder(LU)
tmpLUMatrix = [LU.SID; LU.ID; LU.PID; LU.LWH];
% [~,tepLUorder] = sortrows(tmpLUMatrix',[1, 2, 3, 6],{'ascend','ascend','ascend','descend'}); 
[~,tepLUorder] = sortrows(tmpLUMatrix',[1, 5, 2, 3, 6],{'ascend','descend','ascend','ascend','descend'}); 

% tmpLUMatrix = [LU.ID; LU.LWH; LU.SID; LU.PID];
% [~,tepLUorder] = sortrows(tmpLUMatrix',[5, 1, 6, 4],{'ascend','ascend','ascend','descend'}); %5:SID; 1:ID 4:Hight
%         tepLUorder = 1:length(LU.ID)'; %直接赋值1:n % tepLUorder = [2 3 4 1 5]';
if ~isrow(tepLUorder)
    tepLUorder = tepLUorder';
end
end


% 将LU转换为Item的重要函数
function [Item,LU_Item] = getItem(sLU,Veh)
%% 初始化
% nDim LU维度 nLU LU数量 nItem Item数量 nLUid LU种类
% heightBin Bin最大高度

sz = size(sLU.ID);

nLU = sz(2);
nLUid = size(unique(sLU.ID),2);

hVeh  = Veh.LWH(3,1);  % tmpUniqueBin = unique(Veh.LWH(1:3,:)','rows')'; % hVeh = tmpUniqueBin(3);

% 仅需初始化需要自增的fields
    % Item.LID = zeros(sz);             %Item的ID类型
    % Item.SID = zeros(sz);           
    % Item.UID = zeros(sz);           
    % Item.isRota = ones(sz)*2;    %Item的可旋转类型(初始为2)
Item.Weight = zeros(sz);     %Item的重量
Item.LWH = zeros(3,sz(2));  %Item的宽长高

LU_Item = zeros(2,sz(2));     %dim1:属于第几个Item dim2:属于该Item第几个排放

iItem = 1;
Item_LU = zeros(sz);  % 每个Item内堆垛的LU数量 后期不用
for iLUid=1:nLUid
    hLeft = hVeh;
    for iLU=1:nLU
        if sLU.ID(iLU) == iLUid %仅对当前LU对应Luid在该iLUid内的进行操作
            if hLeft < sLU.LWH(3,iLU) %如当前LU高度不满足(高度在第nDim行)
                iItem =  iItem + 1;
                hLeft = hVeh; 
            end
            
            hLeft = hLeft - sLU.LWH(3,iLU);                 %更新剩余高度
            Item.LWH(1:2,iItem) = sLU.LWH(1:2,iLU);  %更新item长宽
            Item.LWH(3,iItem) = Item.LWH(3,iItem) + sLU.LWH(3,iLU); %更新item高度
            Item.Weight(1,iItem) = Item.Weight(1,iItem) + sLU.Weight(1,iLU); %更新item重量
            
            Item_LU(iItem) = Item_LU(iItem) + 1;
            LU_Item(1,iLU) = iItem;
            LU_Item(2,iLU) = Item_LU(iItem);
            
            Item.LID(1,iItem) = sLU.ID(1,iLU);               %更新ID类型
            Item.SID(1,iItem) = sLU.SID(1,iLU);               
            Item.UID(1,iItem) = sLU.UID(1,iLU);               
            Item.isRota(1,iItem) = sLU.isRota(1,iLU);  %更新ID可旋转类型
            Item.Roated(1,iItem) = sLU.Rotaed(1,iLU);  %更新ID旋转标记
        end
    end
    iItem =  iItem + 1;
end

end

%%  获取ITEMID类型相关数据(同类型ID的体积，面积，重量，item是否可旋转)
function ItemID = getITEMIDArray(Item)

ItemID.ID = unique(Item.LID);
nItemID = numel(ItemID.ID);

for iID = 1:nItemID
    ItemID.Weight(iID) = sum(Item.Weight .* (Item.LID == ItemID.ID(iID)) );
    ItemID.Volume(iID) = sum(prod(Item.LWH) .* (Item.LID == ItemID.ID(iID)) );
    ItemID.Area(iID) = sum(prod(Item.LWH(1:2,:)) .* (Item.LID == ItemID.ID(iID)) );
    ItemID.isRota(iID) =  unique(Item.isRota(Item.LID == ItemID.ID(iID))); if ~isscalar(ItemID.isRota(iID)), error('致命错误'); end
end

end

