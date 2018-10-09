%% GET STRIP 相关属性
% 5 Strip.nbItem % 整数：冗余值, 具体ITEM的堆垛个数 车头摆放依据 -1：混合strip
% TODO 增加Strip.nbLU

%% 函数
function   [Strip] = cpuStripnbItem(Strip,Item,LU)
Strip.nbItem = ones(size(Strip.Weight))*-1;    %单STRIP内部ITEM类型个数, 混合型默认为-1
Strip.nbLU = ones(size(Strip.Weight))*-1;       %单STRIP内部LU类型个数, 混合型默认为-1

%% 7 Strip.nbItem Strip.nbLU
% Strip.nbItem: 混合:-1; 单纯: 对应Strip内部该Item的nbLID类型个数,数值越大,即该LU类型越多
% tmpItemLID = cell2mat(Item.LID);
% uniStripIdx = unique(Item.Item_Strip(1,:));

% 此处Strip一共43个,少于Item_Strip中的Strip个数,因为此时Strip是部分某个bin内的Strip
for iStrip=1:length(Strip.Weight)
    
    if ~Strip.isMixed(1,iStrip)  %如是单纯型
        
        if isfield(Strip,'f') && Strip.f(iStrip)~=1 % 如果是re进来的,需要判断下
               continue;
        end
        
        tmpLID = cell2mat(Strip.LID(iStrip));  % tmpLID给的是LID号; 并非索引顺序; LID号只有1个,索引顺序可能多个
        if ~isscalar(tmpLID), error('单纯型STRIP内的LID的类型不同'); end
        
        % 1 GET Strip.nbLU
        flagLID = LU.ID == tmpLID;
        Strip.nbLU(1,iStrip) = unique(LU.nbLID(flagLID)); % 相同LID号 对应的 nbLID 一定相同
        
        % 2 GET Strip.nbItem
        flagItem = ismember(Item.Item_Strip(1,:),iStrip);
        if ~any(flagItem),
            error('ismember(Item.Item_Strip(1,:),iStrip)错误,Item_Strip中序号问题 '); 
        end
        Strip.nbItem(1,iStrip) = unique(Item.nbItem(flagItem));
         
        end
end

% LIDinItemsArray = cell2mat(Item.LID); %所有ITEM对应的LID值 向量形式
% 
% for iStrip=1:length(Strip.Weight)
%     if ~Strip.isMixed(1,iStrip) %如是单纯型        
%         LU.ID
%         cellLID = Item.LID(Item.Item_Strip(1,:) == uniStripIdx(iStrip)); % cellLID: 本Strip内的ITEM对应的LID值
%         LIDinThisItemArray = cellfun(@(x)x(1), cellLID);
%         if isscalar(unique(LIDinThisItemArray))
%             Strip.nbLU(1,iStrip) = LU.nbLID(LU.ID(unique(LIDinThisItemArray)))
%             Strip.nbItem(1,iStrip) = sum(LIDinItemsArray == unique(LIDinThisItemArray));         
%             if Strip.nbItem(1,iStrip) > Strip.nbLU(1,iStrip)
%                 error('strip.nbItem(1,iItem)< Strip.nbLU((1,iItem)'); 
%             end
%             
%         else
%              error('单纯型STRIP内的ITEM的类型不同'); %arrayLID            
%         end
%     end
% end
end
