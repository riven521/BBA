%% GET STRIP 相关属性
% 5 Strip.nbItem % 整数：冗余值, 具体ITEM的堆垛个数 车头摆放依据 -1：混合strip
% TODO 增加Strip.nbLU

%% 函数
% V2: 似乎是为了还让Strip保持原来的数量, 增加对剩余Strip的判断, 结合f值 
    % 此处Strip个数,少于Item_Strip中的Strip个数,因为此时Strip是部分某个bin内的Strip
    function   [nbItem, nbLU, nbLULID] = cpuStripnbItem(Strip,Item,LU)
    nbItem = ones(size(Strip.Weight))*-1;    %单STRIP内部ITEM类型个数, 混合型默认为-1
    nbLU = ones(size(Strip.Weight))*-1;       %单STRIP内部LU类型个数, 混合型默认为-1
    nbLULID = ones(size(Strip.Weight))*-1;       %单STRIP内部LULID类型个数, 混合型默认为-1   

    if isfield(Strip,'f')
        % UPDATE LU.nbID AND  Item.nbItem 特别重要,对于re调度的来说
        tmpItemLID = cell2mat(Item.LID);
            LU.nbLID = ones(size(LU.nbLID))*-1;   %如不更新,可能导致保留原数据
            LU.nbID = ones(size(LU.nbID))*-1;   %如不更新,可能导致保留原数据
            Item.nbItem= ones(size(Item.nbItem))*-1;  %如不更新,可能导致保留原数据
        Item.nbItem(Item.f) = sum(tmpItemLID(Item.f) == tmpItemLID(Item.f)');
        LU.nbID(LU.f) = sum(LU.ID(LU.f) == LU.ID(LU.f)');
        LU.nbLID(LU.f) = sum(LU.LID(LU.f) == LU.LID(LU.f)');
    end
    
    for iStrip=1:length(Strip.Weight)
        if ~Strip.isMixed(1,iStrip)  %如是单纯型
            if isfield(Strip,'f') && Strip.f(iStrip)~=1 % 如果是re进来的,即包含f字段, 需要判断该Strip是否未1
                continue; 
            end
            
            % 1 GET Strip.nbLU
            %  TODO : flagLID相关可以删除
                    iStripLID = cell2mat(Strip.LID(iStrip));  if ~isscalar(iStripLID), error('单纯型STRIP内的LID的类型不同'); end% tmpLID给的是LID号; 并非索引顺序; LID号只有1个,索引顺序可能多个
                    if isfield(Strip,'f'),        flagLID(LU.f) = LU.ID(LU.f) == iStripLID;
                    else  flagLID = LU.ID == iStripLID; end         
            flagLID2 = ismember(LU.LU_Strip(1,:),iStrip);
                    if ~any(flagLID2),    error('ismember(LU.LU_Strip(1,:),iStrip)错误,LU_Strip中序号问题 ');     end
                  if unique(LU.nbID(flagLID))~= unique(LU.nbID(flagLID2)),    error('错误,LU_Strip中序号问题 ');     end
            nbLU(1,iStrip) = unique(LU.nbID(flagLID2)); % 相同LID号 对应的 nbLID 一定相同 NOTE: 要更新LU.nbLID
            
            % 1.1 GET Strip.nbLLU
            flagLULID = ismember(LU.LU_Strip(1,:),iStrip);
                 if ~any(flagLULID),    error('ismember(LU.LU_Strip(1,:),iStrip)错误,LU_Strip中序号问题 ');     end
                 if length(unique(LU.nbLID(flagLULID)))>1
                     nbLULID(1,iStrip) = -1;  % 如多余1个,表明是LID混合但ID不混合的Strip,
                 else
                     nbLULID(1,iStrip) = unique(LU.nbLID(flagLULID)); % 相同LID号 对应的 nbLID 一定相同 NOTE: 要更新LU.nbLID
                 end
            
            
            % 2 GET Strip.nbItem
            flagItem = ismember(Item.Item_Strip(1,:),iStrip);
                  if ~any(flagItem),    error('ismember(Item.Item_Strip(1,:),iStrip)错误,Item_Strip中序号问题 ');     end
            nbItem(1,iStrip) = unique(Item.nbItem(flagItem)); % NOTE: 要更新Item.nbItem
        end
    end
            
    

% % for iStrip=1:length(Strip.Weight)    
% %     if ~Strip.isMixed(1,iStrip)  %如是单纯型
% %         
% %         if isfield(Strip,'f') && Strip.f(iStrip)~=1 % 如果是re进来的,即包含f字段, 需要判断该Strip是否未1
% %                continue;
% %         end
% %         
% %         iStripLID = cell2mat(Strip.LID(iStrip));  % tmpLID给的是LID号; 并非索引顺序; LID号只有1个,索引顺序可能多个
% %                 if ~isscalar(iStripLID), error('单纯型STRIP内的LID的类型不同'); end
% %         
% %         % 1 GET Strip.nbLU
% %         LU.nbLID = sum(LU.ID==LU.ID');
% %         flagLID = LU.ID == iStripLID;
% %         flagLID2 = ismember(LU.LU_Strip(1,:),iStrip);
% %             if ~any(flagLID2),    error('ismember(LU.LU_Strip(1,:),iStrip)错误,LU_Strip中序号问题 ');     end
% %         Strip.nbLU(1,iStrip) = unique(LU.nbLID(flagLID)); % 相同LID号 对应的 nbLID 一定相同 NOTE: 要更新LU.nbLID
% %         
% %         % 2 GET Strip.nbItem
% %         flagItem = ismember(Item.Item_Strip(1,:),iStrip);
% %             if ~any(flagItem),    error('ismember(Item.Item_Strip(1,:),iStrip)错误,Item_Strip中序号问题 ');     end
% %         Strip.nbItem(1,iStrip) = unique(Item.nbItem(flagItem)); % NOTE: 要更新Item.nbItem
% %         
% % 
% %          %Strip.nbLU(1,iStrip) = unique(LU.nbLID(  LU.ID(:) == LIDinThisItemArray))
% %          %Strip.nbItem(1,iStrip) = sum(LIDinItemsArray == LIDinThisItemArray);
% %         end
% % end

% 
% function   [nbItem, nbLU] = cpuStripnbItem(Strip,Item,LU)
% Strip.nbItem = ones(size(Strip.Weight))*-1;    %单STRIP内部ITEM类型个数, 混合型默认为-1
% Strip.nbLU = ones(size(Strip.Weight))*-1;       %单STRIP内部LU类型个数, 混合型默认为-1

% % %% 7 Strip.nbItem Strip.nbLU
% % % Strip.nbItem: 混合:-1; 单纯: 对应Strip内部该Item的nbLID类型个数,数值越大,即该LU类型越多
% % % tmpItemLID = cell2mat(Item.LID);
% % 
% % % V2: 似乎是为了还让Strip保持原来的数量, 增加对剩余Strip的判断, 结合f值 
% % % 此处Strip个数,少于Item_Strip中的Strip个数,因为此时Strip是部分某个bin内的Strip
% % for iStrip=1:length(Strip.Weight)
% %     
% %     if ~Strip.isMixed(1,iStrip)  %如是单纯型
% %         
% %         if isfield(Strip,'f') && Strip.f(iStrip)~=1 % 如果是re进来的,即包含f字段, 需要判断该Strip是否未1
% %                continue;
% %         end
% %         
% %         tmpLID = cell2mat(Strip.LID(iStrip));  % tmpLID给的是LID号; 并非索引顺序; LID号只有1个,索引顺序可能多个
% %                 if ~isscalar(tmpLID), error('单纯型STRIP内的LID的类型不同'); end
% %         
% %         % 1 GET Strip.nbLU
% %         flagLID = LU.ID == tmpLID;
% %         Strip.nbLU(1,iStrip) = unique(LU.nbLID(flagLID)); % 相同LID号 对应的 nbLID 一定相同
% %         
% %         % 2 GET Strip.nbItem
% %         flagItem = ismember(Item.Item_Strip(1,:),iStrip);
% %         if ~any(flagItem),
% %             error('ismember(Item.Item_Strip(1,:),iStrip)错误,Item_Strip中序号问题 '); 
% %         end
% %         Strip.nbItem(1,iStrip) = unique(Item.nbItem(flagItem));
% %          
% %         end
% % end


% V1: 原始结果, 可以是部分Strip, 结果似乎正确.
% % uniStripIdx = unique(Item.Item_Strip(1,:));
% % 
% % LIDinItemsArray = cell2mat(Item.LID); %所有ITEM对应的LID值 向量形式
% % 
% % unique(LIDinItemsArray)
% % for i=1:length(unique(LIDinItemsArray))
% %     fprintf('LID %1.0f 有 %1.0f 个 \n', i,sum(LIDinItemsArray(:) == i));    
% % end
% % 
% % for iStrip=1:length(Strip.Weight)
% %     if ~Strip.isMixed(1,iStrip) %如是单纯型        
% %         LU.ID;
% %         cellLID = Item.LID(Item.Item_Strip(1,:) == uniStripIdx(iStrip)); % cellLID: 本Strip内的ITEM对应的LID值
% %         LIDinThisItemArray = cellfun(@(x)x(1), cellLID);
% %         LIDinThisItemArray = unique(LIDinThisItemArray)
% %         if isscalar(LIDinThisItemArray)
% %            
% %             %             Strip.nbLU(1,iStrip) = LU.nbLID(LU.ID(LIDinThisItemArray))
% %             a = unique(LU.nbLID(  LU.ID(:) == LIDinThisItemArray))
% %             Strip.nbLU(1,iStrip) = a
% %             Strip.nbItem(1,iStrip) = sum(LIDinItemsArray == LIDinThisItemArray);
% %             if Strip.nbItem(1,iStrip) > Strip.nbLU(1,iStrip)
% %                 error('strip.nbItem(1,iItem)< Strip.nbLU((1,iItem)'); 
% %             end
% %             
% %         else
% %              error('单纯型STRIP内的ITEM的类型不同'); %arrayLID            
% %         end
% %     end
% % end

end
