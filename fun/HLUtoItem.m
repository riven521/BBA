function [LU,Item] = HLUtoItem(LU,Veh)
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
% global ISmaxLayer

%% LU排序
% 获取LU的顺序( 555 )
[LU.order]  = getLUorder(LU); %获取 LU排序(先ID递增,后高度递减)

% LU.order = 1:length(LU.order)

% 获取按order排序后的LU: sLU
if isSameCol(LU),     sLU = structfun(@(x) x(:,LU.order),LU,'UniformOutput',false);  end

% LU.order(:,LU.order)   % os = sLU.order

%% 55 LU->Item转换

% 排序后的sLU, 经过堆垛获取垛Item,以及sLU在Item内的顺序
sz = size(sLU.ID);
nLU = sz(2);
hVeh  = Veh.LWH(3,1);  % tmpUniqueBin = unique(Veh.LWH(1:3,:)','rows')'; % hVeh = tmpUniqueBin(3);

% 仅需初始化需要自增的fields
%     Item.LID = zeros(sz);             %Item的ID类型
%     Item.SID = zeros(sz);
%     Item.UID = zeros(sz);
%     Item.PID = zeros(numel(unique(LU.PID)),sz(2));
    
        Item.isRota = ones(sz)*-1;    %Item的可旋转类型(初始为-1)
        Item.Rotaed = ones(sz)*-1;
    
        Item.HLayer = zeros(sz);    %Item的是否高度满层(初始为-1)

Item.LWH = zeros(3,nLU);     % Item.LWH(1,:) = wStrip;   %dim1-宽度剩余  % Item.LWH(3,:) = hVeh; % 
Item.Weight = zeros(1,nLU); %Item的重量

% 临时使用
tmpItem_LU = zeros(1,nLU);  % 行1：每个Strip内的Item数量 ； 行2：每个Strip内的不同LUID数量
% sLU新增
sLU.LU_Item = zeros(2,sz(2));     %dim1:属于第几个Item dim2:属于该Item第几个排放
   
iItem = 1; 
iLU = 1;            %isrip代表item实质

% 固定LU(基于顺序), 选择ITEM; 
while 1
    
    if iLU > nLU, break; end
    
    % 获取当前LU所在Item位置
    [thisItem,iItem] = getThisItem(iItem);
    
    % 获取Item的长宽高,重量,内部托盘数量HLayer, LU_Item
    insertLUToItem(thisItem,iLU);
    
    iLU = iLU + 1;
    
end


%% Get ITEM 务必可以放 NEXT FIT 如果非空Item，满足高度/层数,放入; 否则,换新Item放入
% V2 getThisItem
    function [thisItem,iItem] = getThisItem(iItem)
    % isflagHeight :   是否ITEM高度满足
    % isNewItem2 ： 是否ITEM属于新
    % isSameID2 ：   是否ITEM内的ID相同
    % isflagLayer ：   是否ITEM高度层数满足
    
        % 同样SID/UID 同样LUID Item高度满足 未考虑Weight等
        isflagHeight =hVeh - Item.LWH(3,iItem) >= sLU.LWH(3,iLU);    %判断是否current's item剩余宽度 >= 当前iLU高度
        
        flagLUinItem = sLU.LU_Item(1,:) == iItem;
        
        if ~any(flagLUinItem) %如果本iItem内不存在任意LU,即空Item
            isNewItem2 = 1;
        else
            isNewItem2 = 0;
            % 1 计算isSameID2
            if ~isscalar(unique(sLU.ID(flagLUinItem))),    error('Item内LUID不同,超预期错误');     end            
            isSameID2 = unique(sLU.ID(flagLUinItem)) ==  sLU.ID(iLU);  %改用V2版本:判断iLU与Item内LU是否属于同一个ID
            % 2 计算isflagLayer
            isflagLayer =  Item.HLayer(iItem) <  sLU.maxHLayer(iLU);  % 非空Item内的高度Layer < 此LU规定的最高高度Layer
        end
        
            % 老版本V1
                %         isSameID = Item.LID(iItem) == sLU.ID(iLU); %判断Item内部ID是否=当前iLU的ID
                %         isNewItem = Item.LWH(3,iItem) == 0; % 判断是否 new Item 高度==0
        
       % 如果是新TIEM, 一定可放；否则：如果高度满足 且最高层数允许 且 与本ITEM内的ID相同，也可放;
        if isNewItem2
                thisItem = iItem;
        else
            
            if isSameID2 && (isflagHeight && isflagLayer) %如果高度允许 且最高层数允许 且LU ID相同 && isflagLayer
                thisItem = iItem;
            else
                iItem = iItem + 1;
                [thisItem,iItem] = getThisItem(iItem);
            end
            

% %             if ~ISmaxLayer
% %             %if (isflagHeight && isSameID2) || isflagLayer %如果高度允许 且最高层数允许 且LU ID相同 && isflagLayer
% %             if isSameID2 && (isflagHeight && isflagLayer) %如果高度允许 且最高层数允许 且LU ID相同 && isflagLayer
% %                  thisItem = iItem;
% %             else
% %                 iItem = iItem + 1;
% %                 [thisItem,iItem] = getThisItem(iItem);
% %             end
% %             
% %             else % 首次一定不考虑层数, 否则初次堆垛可能有问题
% %             if isSameID2 && isflagHeight % || isflagLayer) 
% %                  thisItem = iItem;
% %             else
% %                 iItem = iItem + 1;
% %                 [thisItem,iItem] = getThisItem(iItem);
% %             end
% %             end
% %             
            
        end
        
    end


%% v1 getThisItem v2: 改进删除注释 该Next Fit 为 First Fit 
% %     function [thisItem,iItem] = getThisItem(iItem)
% %     % isflagHeight :   是否ITEM高度满足
% %     % isNewItem2 ： 是否ITEM属于新
% %     % isSameID2 ：   是否ITEM内的ID相同
% %     % isflagLayer ：   是否ITEM高度层数满足
% %     
% %         % 同样SID/UID 同样LUID Item高度满足 未考虑Weight等
% %          isflagHeight =hVeh - Item.LWH(3,iItem) >= sLU.LWH(3,iLU); %判断是否current's item剩余宽度 >= 当前iLU高度
% %         
% %         flagLUinItem = sLU.LU_Item(1,:) == iItem;
% %         if ~any(flagLUinItem) %如果本iItem内不存在任意LU,即空Item
% %             isNewItem2 = 1;
% %         else
% %             isNewItem2 = 0;
% %             % 1 计算isSameID2
% %             if ~isscalar(unique(sLU.ID(flagLUinItem))),    error('Item内LUID不同,超预期错误');     end            
% %             isSameID2 = unique(sLU.ID(flagLUinItem)) ==  sLU.ID(iLU);  %改用V2版本:判断iLU与Item内LU是否属于同一个ID
% %             % 2 计算isflagLayer
% %             isflagLayer =  Item.HLayer(iItem) <  sLU.maxHLayer(iLU);  % 非空Item内的高度Layer < 此LU规定的最高高度Layer
% %         end
% %         
% %             % 老版本V1
% %                 %         isSameID = Item.LID(iItem) == sLU.ID(iLU); %判断Item内部ID是否=当前iLU的ID
% %                 %         isNewItem = Item.LWH(3,iItem) == 0; % 判断是否 new Item 高度==0
% %         
% %        % 如果是新TIEM, 一定可放；否则：如果高度满足 且最高层数允许 且 与本ITEM内的ID相同，也可放;
% %         if isNewItem2
% %                 thisItem = iItem;
% %         else
% %             if isflagHeight && isflagLayer && isSameID2 %如果高度允许 且最高层数允许 且LU ID相同
% %                  thisItem = iItem;
% %             else
% %                 iItem = iItem + 1;
% %                 [thisItem,iItem] = getThisItem(iItem);
% %             end
% % 
% %         end
% %         
% %     end

%% Put LU into thisItem
    function insertLUToItem(thisItem,iLU)
        %更新Item的长宽高,种类,Lu_Item
        Item.LWH(3,thisItem) = Item.LWH(3,thisItem)  + sLU.LWH(3,iLU);
        Item.LWH(1:2,thisItem) = sLU.LWH(1:2,iLU);  %更新item长宽
        Item.Weight(1,thisItem) = Item.Weight(1,thisItem) + sLU.Weight(1,iLU); %更新item重量
        
        tmpItem_LU(1,thisItem) = tmpItem_LU(1,thisItem) + 1;
        
        sLU.LU_Item(1,iLU) = thisItem;
        sLU.LU_Item(2,iLU) = tmpItem_LU(1,thisItem);
        
                    %         tmpLUThisItem = sLU.LU_Item(1,:) == thisItem;
                    %         tmpItem_LU(2,iItem) = numel(unique(sLU.PID(1,tmpLUThisItem)));

% 注释掉Item的是否可旋转或旋转状态, 想不到哪里用
        Item.isRota(1,thisItem) = sLU.isRota(1,iLU);      %更新ID可旋转类型
        Item.Rotaed(1,thisItem) = sLU.Rotaed(1,iLU);   %更新ID旋转标记
        
% 注释掉Item内LU个数(通过sLU.LU_Item计算最好)        
        flagLUinItem = sLU.LU_Item(1,:) == thisItem;
        Item.HLayer(thisItem) = sum(flagLUinItem);          %更新Item内已安置层数或内部托盘个数
        
%         Item.LID(1,thisItem) = sLU.ID(1,iLU); %更新ID类型        
%         Item.SID(1,thisItem) = sLU.SID(1,iLU);   % Item.UID(1,thisItem) = sLU.UID(1,iLU);
        
%         Item.PID(sLU.PID(1,iLU),thisItem) = Item.PID(sLU.PID(1,iLU),thisItem) + 1;    % 555 更新多行PID - 数值为出现次数
%         Item.PID(sLU.PID(1,iLU),thisItem) = 1;      % 555 更新多行PID - 数值为出现与否     

    end

% LU内部更新,sLU依据order变化回来
if isSameCol(sLU)
    LU = getReorderStruct(LU.order, sLU);
else
    error('不能使用structfun');
end

% Item去除未使用 %     Item.Rotaed(:,Item.itemorder) = sLU.Rotaed;
% 如果ITEM的列数全部相同
if isSameCol(Item)
    Item = structfun(@(x) x( : , Item.LWH(1,:)>0 ), Item, 'UniformOutput', false);
else
    error('不能使用structfun');
end

%% 测试script TO BE FIX
% 输出主要结果:获得每个item包含的 原始 LU序号
% printscript(LU,Item);
end

%% getLUorder V2
function [ord] = getLUorder(LU)

T = getTableLU(LU);

% TODO 顺序

% [~,ord] = sortrows(T,{'SID','EID',...   % 离散: 相同SID/EID 递增排序 (顺序给定)
%     'isNonMixed','isMixedTile',...      % 离散: 不会出现非满垛 会出现非满垛时,属于非满垛的托盘标记1
%     'W','ID','maxHLayer','LID','H','PID'},...
%     {'ascend','ascend','descend','ascend','descend','ascend','descend','ascend','descend','descend'});  %NEW EID NEW maxHLayer

[~,ord] = sortrows(T,{'SID','EID',...   % 离散: 相同SID/EID 递增排序 (顺序给定)
    'W','ID','maxHLayer','LID','H','PID'},...
    {'ascend','ascend','descend','ascend','ascend','ascend','descend','descend'});  %NEW EID NEW maxHLayer

% 简化版,效果同上
% [~,ord2] = sortrows(T,{'SID','EID',...   % 离散: 相同SID/EID 递增排序 (顺序给定)
%     'isNonMixed','isMixedTile',...      % 离散: 不会出现非满垛 会出现非满垛时,属于非满垛的托盘标记1
%     'W','ID','maxHLayer','LID','H','PID'},...
%     {'ascend','ascend','descend','ascend','descend','ascend','descend','ascend','descend','descend'});  %NEW EID NEW maxHLayer


if ~isrow(ord),    ord = ord';  end

end

%% getLUorder V1 V2:删除注释,改为表格格式排序
% % function [tepLUorder] = getLUorder(LU)
% % % V1: ********** 不考虑isNonMixed
% % % tmpLUMatrix = [LU.SID; LU.ID; LU.PID; LU.LWH; LU.Weight];
% % %tmpLUMatrix = [LU.SID; LU.ID; LU.PID; LU.LWH(1,:); LU.LWH(2,:); LU.LWH(3,:); LU.Weight];
% % %[~,tepLUorder] = sortrows(tmpLUMatrix',[1, 5, 2, 3, 6, 7],{'ascend','descend','ascend','ascend','descend','descend'}); 
% % 
% % % 1 SID; 2 长度；3 ID；4 PID；5 高度；6 重量；
% % %tmpLUMatrix = [LU.SID; LU.LWH(2,:); LU.ID; LU.LID; LU.PID; LU.LWH(3,:); LU.Weight];
% % %[~,tepLUorder] = sortrows(tmpLUMatrix',[1, 2, 4, 5, 6, 7 ],{'ascend','descend','ascend','ascend','descend','descend'}); 
% % 
% % % V2: ********** 考虑isNonMixed
% % global ISisNonMixedLU ISisMixTileLU % TODO 考虑不满托盘,同样LULID下
% % tmpLUMatrix = [LU.SID; LU.isNonMixed; LU.isMixedTile; ...
% %                            LU.LWH(2,:); LU.ID; LU.LID; ...
% %                            LU.PID; LU.LWH(3,:); LU.Weight; LU.EID; LU.maxHLayer ]; %NEW EID NEW maxHLayer
% %                          
% % if ISisNonMixedLU==1
% %     if ISisMixTileLU==1
% %         % V4: 修改增LU的最高堆垛层11 （对相同LU内部按照maxHLayer递减排序）
% %         [~,tepLUorder] = sortrows(tmpLUMatrix',[1, 10, 2, 3, ...
% %                                     4, 5, 11, 6, 8, 7 ],{'ascend','ascend','descend','ascend','descend','ascend','descend','ascend','descend','descend'}); 
% %         % V3: 修改增jLU的EID排序
% %         %[~,tepLUorder] = sortrows(tmpLUMatrix',[1, 10, 2, 3, ...
% %            %                         4, 5, 6, 8, 7 ],{'ascend','ascend','descend','ascend','descend','ascend','ascend','descend','descend'}); 
% %         % V2: 修改为PID优先放在相同LID高度优先之后
% %         % [~,tepLUorder] = sortrows(tmpLUMatrix',[1, 2, 3, 4, 5, 6, 8, 7 ],{'ascend','descend','ascend','descend','ascend','ascend','descend','descend'}); 
% %                 % V1 : 问题在于LU.PID 零件号排序意义不大,无论递增或递减
% %                 % [~,tepLUorder] = sortrows(tmpLUMatrix',[1, 2, 3, 4, 5, 6, 7 ],{'ascend','descend','ascend','descend','ascend','ascend','descend'}); 
% %     else
% %         [~,tepLUorder] = sortrows(tmpLUMatrix',[1, 2, 4, 5, 6, 7 ],{'ascend','descend','descend','ascend','ascend','descend'}); 
% %     end
% % else
% %         [~,tepLUorder] = sortrows(tmpLUMatrix',[1, 4, 5, 6, 7 ],{'ascend','ascend','ascend','descend','descend'}); 
% % end
% % 
% % 
% % % [~,tepLUorder] = sortrows(tmpLUMatrix',[1, 4, 5, 6, 7 ],{'ascend','ascend','ascend','descend','descend'}); 
% % % [~,tepLUorder] = sortrows(tmpLUMatrix',[1, 2, 3, 4, 5, 6, 7 ],{'ascend','descend','ascend','ascend','ascend','descend','descend'}); 
% % 
% % % tmpLUMatrix = [LU.LWH(2,:); LU.ID; LU.LID; LU.PID; LU.LWH(3,:); LU.Weight]
% % 
% % 
% % % [~,tepLUorder] = sortrows(tmpLUMatrix',[1, 2, 3, 6],{'ascend','ascend','ascend','descend'}); 
% % 
% % 
% % % [~,tepLUorder] = sortrows(tmpLUMatrix',[1, 5, 2, 3, 7, 6],{'ascend','descend','ascend','ascend','descend','descend'}); 
% % 
% % % tmpLUMatrix = [LU.ID; LU.LWH; LU.SID; LU.PID];
% % % [~,tepLUorder] = sortrows(tmpLUMatrix',[5, 1, 6, 4],{'ascend','ascend','ascend','descend'}); %5:SID; 1:ID 4:Hight
% % %         tepLUorder = 1:length(LU.ID)'; %直接赋值1:n % tepLUorder = [2 3 4 1 5]';
% % if ~isrow(tepLUorder),    tepLUorder = tepLUorder'; end
% % 
% % end

%% COMMENT
% function printscript(LU,Item)
%     for iItem = 1:max(LU.LU_Item(1,:))
% %         [~,idx] = find(LU.LU_Item(1,:)==iItem);
% %         fprintf('item %d 的长宽高为:  ',iItem);
% %         fprintf('( %d ) ',Item.LWH(:,iItem));
% %         fprintf('\n');
% %         fprintf('item %d 包含 original LU 索引号(长宽高)为  \n  ',iItem);
% %         fprintf('%d ',idx);
% %         fprintf('( %d ) ', LU.LWH(:,idx));
% %         fprintf('\n');
% %         fprintf('item %d 包含 original LU 索引号(高)为  \n  ',iItem);
% %         fprintf('%d ',idx);
% %         fprintf('( %d ) ', LU.LWH(3,idx)); 
% %         fprintf('\n');
% %         fprintf('item %d 包含 original LU 重量为  \n  ',iItem);
% %         fprintf('%d ',idx);
% %         fprintf('( %d ) ', LU.Weight(:,idx));
% %         fprintf('\n');
% %                fprintf('item %d 包含 original LU ***为  \n  ',iItem);
% %         fprintf('%d ',idx);
% %         fprintf('( %d ) ', LU.LU_Item(2,idx)); 
% %         fprintf('\n'); 
% %         isWeightUpDown
% %         if length(idx) > 1 %Item包含不只一个Item,需要判断是否有轻重的变化
% %             currLUWeight = zeros(1,length(idx));
% %             currLUHight = zeros(1,length(idx));
% %             for iIdx = 1:length(idx)
% %                currIdx = idx(LU.LU_Item(2,idx) == iIdx);
% %                currLUWeight(iIdx) = LU.Weight(:,currIdx);
% %                currLUHight(iIdx) = LU.LWH(3,currIdx);      
% %             end
% %             if diff(currLUWeight) > 0 % 代表下轻上重
% %                 currLUWeight
% %             end
% %             if diff(currLUHight) >0  
% %                 currLUHight
% %             end
% %         end
% 
%     end
% end

