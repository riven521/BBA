%% GET STRIP 相关属性
% LU.LU_Strip, LU.CoordLUStrip
% 1 Strip.isMixed % 1：混合层； 0：单纯层
% 2 Strip.isHeightFull % 1：全部都是满层； 0：包含非满层 依据Item.isHeightFull 判断 
% 3 Strip.isWidthFull % % 1：全部都是宽度满层； 0：宽度非满层 依据 Strip宽度剩余 是否 可容纳 最小的Item宽度判断 
% 4 Strip.maxHeight % Strip内所有Item的最大高度

% 5 Strip.nbItem % 整数：冗余值, 具体ITEM的堆垛个数 车头摆放依据 -1：混合strip
% 6 Strip.isAllPured % ：1：单纯且该ITEM没有混合型； 0：单纯但其它strip有混合的； -1：混合strip
% 7 Strip.isSingleItem %: 1:单纯且单个； 0：单纯且多个；-1：混合strip

% 8 Strip.loadingrate   % 每个strip的装载比率
% 9 Strip.loadingrateLimit   % 每个strip的有限装载比率

% Strip.seqSW
% Strip.LID

%% 函数
function   [Strip,LU] = cpuStrip(Strip,Item,LU,Veh)
%% 初始化
    Strip.isMixed = ones(size(Strip.Weight))*-1;   %是否为混合型,包含多个LID 
    Strip.isHeightFull = ones(size(Strip.Weight))*-1;   %是否包含非Full的Item
    Strip.isHeightBalance = ones(size(Strip.Weight))*-1;   %是否Item的高度差异不大
    Strip.nbItem = ones(size(Strip.Weight))*-1;   %单STRIP内部ITEM类型个数, 混合型默认为-1
    Strip.isAllPured = ones(size(Strip.Weight))*-1;   %单STRIP对应LID是否包含混合STRIP, 包含混合型默认为-1
    Strip.isSingleItem = ones(size(Strip.Weight))*-1;   %单Strip内对应只有1个ITEM
    Strip.isWidthFull = ones(size(Strip.Weight))*-1;     %是否为宽度非Full的Item
    Strip.maxHeight = ones(size(Strip.Weight))*-1;     %Strip的最高高度.
    Strip.lowestHeight = ones(size(Strip.Weight))*-1;     %Strip的最低高度.
    Strip.meanHeight = ones(size(Strip.Weight))*-1;     %Strip的最低高度.
%     Strip.seqSW = ones(size(Strip.Weight))*-1;     %Strip的??? 暂未用

    Strip.Stripvolume = ones(size(Strip.Weight))*-1;  %每个strip的可用体积 = 高度*宽度(车辆的宽度)
    Strip.StripvolumeLimit = ones(size(Strip.Weight))*-1; % %每个strip的有限可用体积 = 高度*宽度(strip使用宽度=车辆宽度-strip剩余宽度)
    Strip.Itemvolume = ones(size(Strip.Weight))*-1; %每个strip包含的Item装载体积
    Strip.loadingrate = ones(size(Strip.Weight))*-1; % 每个strip的装载比率
    Strip.loadingrateLimit = ones(size(Strip.Weight))*-1;     % 每个strip的有限装载比率

%% 0.0 计算LU.LU_Strip, LU.CoordLUStrip, Strip内的PID,LID,SID
% 由混合的LU.DOC新增LU_STRIP, 计算STRIP内包含的PID,LID,SID等数据 1808新增
nbLU = size(LU.LWH,2);
LU.LU_Strip = zeros(2,nbLU);
LU.CoordLUStrip = zeros(3,nbLU);

% 更新LU_Strip
for iLU=1:nbLU
    % 更新LU_Strip第一行
    iItem = LU.LU_Item(1,iLU);   %iLU属于第几个Item, Item属于第几个Strip,则Lu属于第几个Strip
    LU.LU_Strip(1,iLU)= Item.Item_Strip(1,iItem);
    % 更新LU_Strip第二行
    fiItem = find(Item.Item_Strip(1,:) == Item.Item_Strip(1,iItem) & Item.Item_Strip(2,:) < Item.Item_Strip(2,iItem));
    nbLUfiItem = sum(ismember(LU.LU_Item(1,:),fiItem));
    LU.LU_Strip(2,iLU) = nbLUfiItem+LU.LU_Item(2,iLU); % 进入Strip顺序: 同一Strip内先前进入个数nbLUfiItem + 本iLU在Item的顺序
    % 更新LU.CoordLUStrip
    LU.CoordLUStrip(1,iLU) = Item.CoordItemStrip(1,iItem);
    LU.CoordLUStrip(2,iLU) = Item.CoordItemStrip(2,iItem);
        % fLU: 与iLU同属iItem 且 顺序晚于本iLU; 可能为空, 不影响.
    fLU = LU.LU_Item(1,:) == iItem & LU.LU_Item(2,:) < LU.LU_Item(2,iLU);
    LU.CoordLUStrip(3,iLU) = sum(LU.LWH(3,fLU));
end


LU.DOC=[LU.DOC; LU.LU_Strip];
nStrip = size(Strip.LW,2);
for iStrip=1:nStrip
    tmp = LU.DOC([1,2,3], LU.DOC(8,:) == iStrip);
    Strip.PID(:,iStrip) = num2cell(unique(tmp(1,:))',1);
    Strip.LID(:,iStrip) = num2cell(unique(tmp(2,:))',1);
    Strip.SID(:,iStrip) = num2cell(unique(tmp(3,:))',1);
end
    
%% 0: 计算strip装载率
Strip = computeLoadingRateStrip(Strip,Item,Veh); 
    
%% 1: STRIP.isMixed: STRIP增加判断是否单纯型/混合型判断
Strip = isMixedStrip(Strip);

%% 2: STRIP.isHeightFull: STRIP增加判断是否包含非Full的Item. 
Strip = isHeightBalanceStrip(Strip,Item);
Strip = isFullStrip(Strip,Item);

%% 3: STRIP.isWidthFull : STRIP增加判断是否包含宽度width非Full的Item. 
Strip = isWidthFullStrip(Strip,Item);

%% 4: STRIP.maxHeight : 计算Strip的最大高度
for i=1:length(Strip.maxHeight)
    % 计算最大值
    % Item.Item_Strip(1,:) == i) : Strip i 内部的Item flag
    Strip.maxHeight(i) = max(Item.LWH(3, Item.Item_Strip(1,:) == i));
    Strip.lowestHeight(i) = min(Item.LWH(3, Item.Item_Strip(1,:) == i));
    Strip.meanHeight(i) = Strip.maxHeight(i) - Strip.lowestHeight(i); 
    if any(Strip.meanHeight(i)<0), error('负值不可能;'); end
end

%% 5,6
%Strip.isAllPured：混合:-1; 单纯: 1 (混合strip内没有改ID) ; 0 (混合strip内含有该单纯strip的ID)
%Strip.isSingleItem: 混合: -1; 单纯: Strip内仅有一个Item,必定是单纯的.
LIDinItemsArray = cellfun(@(x) x(1), Item.LID); % arrayAllLID: 所有ITEM对应的LID值 向量形式
mixedStrip = find(Strip.isMixed(1,:) == 1);
mixedLID = [];
for m=1:length(mixedStrip)
     f = Item.Item_Strip(1,:) == mixedStrip(m);
     mixedLID = [mixedLID, LIDinItemsArray(f)];
end
mixedLID = unique(mixedLID);

uniItem = unique(Item.Item_Strip(1,:));
for i=1:length(Strip.isAllPured)
    if ~Strip.isMixed(1,i) %如是单纯型        
        cellLID = Item.LID(Item.Item_Strip(1,:) == uniItem(i)); % cellLID: 本Strip内的ITEM对应的LID值
            %         cellLID = Item.LID(Item.Item_Strip(1,:) == i); % cellLID: 本Strip内的ITEM对应的LID值
        LIDinThisItemArray = cellfun(@(x)x(1), cellLID);
        if isscalar(unique(LIDinThisItemArray)) 
            
            if isscalar(LIDinThisItemArray)
                Strip.isSingleItem(1,i) = 1;
            else
                Strip.isSingleItem(1,i) = 0;
            end
            
            if ismember(unique(LIDinThisItemArray),mixedLID)
                Strip.isAllPured(1,i) = 0;
            else
                Strip.isAllPured(1,i) = 1;
            end
            
            % Strip.nbItem(1,i) = sum(LIDinItemsArray ==
            % unique(LIDinThisItemArray)); %单独放入函数计算
            
        else
             error('单纯型STRIP内的ITEM的类型不同'); %arrayLID            
        end
    else % strip如是混合型     
        Strip.isSingleItem(1,i) = 0;
        Strip.isAllPured(1,i) = 0;
    end
end

%% 7
%Strip.nbItem: 混合:-; 单纯: 对应Strip内部该Item的nbLID类型个数,数值越大,即该LU类型越多
[Strip.nbItem, Strip.nbLU] = cpuStripnbItem(Strip,Item,LU);

end

%% 局部函数 %%

%% 函数1: 判断STRIP是否包含非HeightFull的Item
function Strip = isHeightBalanceStrip(Strip,Item) 
    % 1 循环判断Strip是否包含Item为full的,如包含,则Strip为full
    uniStrip = unique(Item.Item_Strip(1,:));

    % 2 循环判断Strip内部Item之间的最大差值, 是否<= 最小的对角线或一共绝对值, 如是,均为Full; 
    for i=1:length(Strip.isHeightFull)
        fItem = Item.Item_Strip(1,:) == uniStrip(i);        
        % 高度间隙
        maxHeightDiff = max(abs(diff(Item.LWH(3,fItem))));
        % 对比2: 最高的Item的1/3
        oneThirdsHeightItem = max(Item.LWH(3,fItem))*1/3
        
        % 对比2: 最大差值, 是否<= 1/3最高Item (Item高度平行, 即使很低,  也认为是满层)
        if ~isempty(maxHeightDiff)
            if maxHeightDiff <= oneThirdsHeightItem                 
                     Strip.isHeightBalance(i) = 1;
            else
                     Strip.isHeightBalance(i) = 0;       
            end
        else %均为空值, 即只有Strip一个堆垛, 算是均衡的
            if sum(fItem)~=1, error('该strip包含不是一个堆垛!'); end
            Strip.isHeightBalance(i) = 1;
        end
    end
    
    % 防错语句
    if any(Strip.isHeightBalance==-1), error('存在Strip.isHeightBalance未分配!'); end
end

%% V1 不单独考虑isHeightBalance高度均衡, 完全
% function Strip = isFullStrip(Strip,Item)
%     % 1 循环判断Strip是否包含Item为full的,如包含,则Strip为full
%     uniStrip = unique(Item.Item_Strip(1,:));
% %     for i=1:length(Strip.isHeightFull)
% %         %          if all(Item.isHeightFull(Item.Item_Strip(1,:) == i)) %如果本STRIP对应ITEM的isFull均为1,则本STRIP也为full
% %          if all(Item.isHeightFull(Item.Item_Strip(1,:) == uniStrip(i))) %如果本STRIP对应ITEM的isFull均为1,则本STRIP也为full
% %              Strip.isHeightFull(i) = 1;
% %          else
% %              Strip.isHeightFull(i) = 0;
% %          end
% %     end
%     
% 
%     % 2 循环判断Strip内部Item之间的最大差值, 是否<= 最小的对角线或一共绝对值, 如是,均为Full; 
%     for i=1:length(Strip.isHeightFull)
%         fItem = Item.Item_Strip(1,:) == uniStrip(i);        
%         diagItem = sqrt(Item.LWH(1,fItem).^2 + Item.LWH(2,fItem).^2);
%         % 高度间隙
%         maxHeightDiff = max(abs(diff(Item.LWH(3,fItem))));
%         
%         % 对比1: 最小的Item对角线
%         minDiagItem = min(diagItem);
%         % 对比2: 最高的Item的1/3
%         oneThirdsHeightItem = max(Item.LWH(3,fItem))*1/3
%         % 对比3: 绝对值
%         absHeight = 300;
%         
%         % 对比2: 最大差值, 是否<= 1/3最高Item (Item高度平行, 即使很低,  也认为是满层)
%         if ~isempty(maxHeightDiff)
%             if maxHeightDiff <= oneThirdsHeightItem                 
%                      Strip.isHeightFull(i) = 1;
%                      % 3 增加即使maxHeightDiff很小, 若整体高度低, 也视为非HeightFull
%                      if ~all(Item.isHeightFull(Item.Item_Strip(1,:) == uniStrip(i))) %如果本STRIP对应ITEM的isFull不是均为1,则本STRIP非full
%                          Strip.isHeightFull(i) = 0;        
%                      end
%             else
%                      Strip.isHeightFull(i) = 0;       
%                      % 3 增加即使maxHeightDiff很大, 若整体高度高, 也视为HeightFull
%                      if all(Item.isHeightFull(Item.Item_Strip(1,:) == uniStrip(i))) %如果本STRIP对应ITEM的isFull均为1,则本STRIP也为full
%                          Strip.isHeightFull(i) = 1;
%                      end
%             end
%         else
%             if all(Item.isHeightFull(Item.Item_Strip(1,:) == uniStrip(i))) %如果本STRIP对应ITEM的isFull均为1,则本STRIP也为full
%                 Strip.isHeightFull(i) = 1;
%             else
%                 Strip.isHeightFull(i) = 0;
%             end
%         end
%     end
%     
%     % 防错语句
%     if any(Strip.isHeightFull==-1), error('存在Strip.isHeightFull未分配!'); end
%      
% end

%% V2 考虑isHeightBalance高度均衡
function Strip = isFullStrip(Strip,Item)
    % 1 循环判断Strip是否包含Item为full的,如包含,则Strip为full
    uniStrip = unique(Item.Item_Strip(1,:));

    % 2 循环判断Strip内部Item之间的最大差值, 是否<= 最小的对角线或一共绝对值, 如是,均为Full; 
    for i=1:length(Strip.isHeightFull)
        if Strip.isHeightBalance(i) == 0 %如果高度不均衡,一定是高度不满 (不考虑高度都很高,但不均衡情形)
            Strip.isHeightFull(i) = 0;
        else %如果高度均衡,看Item是否都是满层,如不是,Strip高度也不满
            Strip.isHeightFull(i) = 1; %Strip单一Item也是高度均衡, 如它的Item是不满, 该Strip也是不满
            if ~all(Item.isHeightFull(Item.Item_Strip(1,:) == uniStrip(i))) %如果本STRIP对应ITEM的isFull不是均为1,则本STRIP非full
                Strip.isHeightFull(i) = 0;
            end
        end
    end

    % 防错语句
    if any(Strip.isHeightFull==-1), error('存在Strip.isHeightFull未分配!'); end

end

%% 函数2: 判断STRIP是否包含非WidthFull的Item
% ****************** Strip内是否包含Width非Full的Item计算 ************ 开放
function Strip = isWidthFullStrip(Strip,Item) 
    % 循环判断Strip是否Widthfull, 宽度间隙 >= 本Strip包含的Item的宽度最小值
    uniItem = unique(Item.Item_Strip(1,:));
    for i=1:length(Strip.isWidthFull)
        flagItem = Item.Item_Strip(1,:) == uniItem(i); 
        max(Item.LWH(1, flagItem))
        if Strip.LW(1, i) >= min(Item.LWH(1, flagItem))
            Strip.isWidthFull(i) = 0;
        else
            Strip.isWidthFull(i) = 1;
        end
    end
end

%% 函数3: 判断STRIP是否混合型
function Strip = isMixedStrip(Strip)
    % 循环判断Strip是否为混合型
    for i=1:length(Strip.isMixed)
         if numel(Strip.LID{i}) > 1
             Strip.isMixed(i) = 1;
         else
             Strip.isMixed(i) = 0;
         end
    end
end

function Strip = computeLoadingRateStrip(Strip,Item,Veh)
    % 初始化
    nStrip = size(Strip.LW,2);

    % 计算每个strip的装载率
    %每个strip的可用体积 = 高度*宽度(车辆的宽度)
    Strip.Stripvolume = Strip.LW(2,:)*Veh.LWH(1,1);
    %每个strip的有限可用体积 = 高度*宽度(strip使用宽度=车辆宽度-strip剩余宽度)
    Strip.StripvolumeLimit = Strip.LW(2,:) .* (Veh.LWH(1,1) - Strip.LW(1,:));
    a = Item.LWH;
    b = Item.Item_Strip;
    uniItem = unique(Item.Item_Strip(1,:));
    for iStrip =1:nStrip
        %每个strip包含的Item装载体积
        Strip.Itemvolume(iStrip)= sum(a(1, (b(1,:)==uniItem(iStrip))) .* a(2, (b(1,:)==uniItem(iStrip))));
    end
    %每个strip的装载比率
    Strip.loadingrate =  Strip.Itemvolume ./ Strip.Stripvolume;
    %每个strip的有限装载比率
    Strip.loadingrateLimit =  Strip.Itemvolume ./ Strip.StripvolumeLimit;
end
