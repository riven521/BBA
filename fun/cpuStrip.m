%% GET STRIP 相关属性
% LU.LU_Strip, LU.CoordLUStrip
% 1 Strip.isMixed % 1：混合层； 0：单纯层
% 2 Strip.isHeightFull % 1：全部都是满层； 0：包含非满层 依据Item.isHeightFull 判断 
% 3 Strip.isWidthFull % % 1：全部都是宽度满层； 0：宽度非满层 依据 Strip宽度剩余 是否 可容纳 最小的Item宽度判断 
% 4 Strip.maxHeight % Strip内所有Item的最大高度

% 5 Strip.nbItem % 整数：冗余值, 具体ITEM的堆垛个数 车头摆放依据 -1：混合strip
% 6 Strip.isAllPured % ：1：单纯且该ITEM没有混合型； 0：单纯但其它strip有混合的； -1：混合strip
% 7 Strip.isSingleItem %: 1:单纯且单个； 0：单纯且多个；-1：混合strip

% computeLoadingRateStrip:
% 8 Strip.loadingrate           % 每个strip的装载比率
% 9 Strip.loadingrateLimit   % 每个strip的有限装载比率

function   [Strip] = cpuStrip(Strip,Item,LU,Veh)

%% 初始化
    
    Strip.isHeightFull = ones(size(Strip.Weight))*-1;   %是否包含非Full的Item
    Strip.isHeightBalance = ones(size(Strip.Weight))*-1;   %是否Item的高度差异不大
    Strip.isWidthFull = ones(size(Strip.Weight))*-1;     %是否为宽度非Full的Item  若 宽度间隙 >= 本Strip包含的Item的宽度最小值， 则宽度不满，取值0
        
    Strip.isMixed = ones(size(Strip.Weight))*-1;         %是否为混合型,包含多个LID 
    Strip.isMixedSID = ones(size(Strip.Weight))*-1;   %是否为混合型,包含多个SID   暂时用不到 注释
    Strip.isMixedEID = ones(size(Strip.Weight))*-1;   %是否为混合型,包含多个EID   暂时用不到 注释
    
%     Strip.nbLID = ones(size(Strip.Weight))*-1;   %STRIP内LID种类数量,非实际托盘个数
%     Strip.nbSID= ones(size(Strip.Weight))*-1;    %STRIP内SID种类数量
%     Strip.nbEID= ones(size(Strip.Weight))*-1;    %STRIP内EID种类数量
         
    Strip.nbItem = ones(size(Strip.Weight))*-1;            %单STRIP内部ITEM类型个数, 混合型默认为-1
%     Strip.isAllPured = ones(size(Strip.Weight))*-1;       %单STRIP对应LID是否包含混合STRIP, 包含混合型默认为-1
%     Strip.isSingleItem = ones(size(Strip.Weight))*-1;   %单Strip内对应只有1个ITEM
        
    Strip.maxHeight = ones(size(Strip.Weight))*-1;     %Strip的最高高度.
    Strip.lowestHeight = ones(size(Strip.Weight))*-1;     %Strip的最低高度.
    Strip.meanHeight = ones(size(Strip.Weight))*-1;     %Strip的最低高度.
 
%     Strip.Stripvolume = ones(size(Strip.Weight))*-1;  %每个strip的可用体积 = 高度*宽度(车辆的宽度)
%     Strip.StripvolumeLimit = ones(size(Strip.Weight))*-1; % %每个strip的有限可用体积 = 高度*宽度(strip使用宽度=车辆宽度-strip剩余宽度)
%     Strip.Itemvolume = ones(size(Strip.Weight))*-1; %每个strip包含的Item装载体积
    Strip.loadingrate = ones(size(Strip.Weight))*-1; % 每个strip的装载比率
    Strip.loadingrateLimit = ones(size(Strip.Weight))*-1;     % 每个strip的有限装载比率
    
   % isGapBalance: 条带的宽度W,沿车厢场边的差值 
%     Strip.GapValue = ones(size(Strip.Weight))*-1;         %混合STRIP（isMixed==1）及宽度不满STRIP（此时最小长度为0）内的GAP值
%     Strip.isGapBalance = ones(size(Strip.Weight))*-1;   %GAP值超过最长ITEM的1/3即为不balance的
      
%     Strip.maxLULength= ones(size(Strip.Weight))*-1;     %Strip的最高高度.
%     Strip.lowestLULength = ones(size(Strip.Weight))*-1;     %Strip的最低高度.
%     Strip.meanLUWidth = ones(size(Strip.Weight))*-1;     %Strip的最低高度.
%     Strip.seqSW = ones(size(Strip.Weight))*-1;     %Strip的??? 暂未用

%% 1 V2 计算STRIP的PID,LID,SID
    % 由混合的LU.DOC新增LU_STRIP, 计算STRIP内包含的PID,LID,SID等数据 1808新增
    % 由TABLE计算,方便知道什么是什么,不用1,2,3数字替换
    t = struct2table(structfun(@(x) x',LU,'UniformOutput',false));
    
    nStrip = size(Strip.LW,2);
    for iStrip=1:nStrip
        f = t.LU_Strip(:,1) == iStrip;   
        Strip.LID(:,iStrip) = {unique(t.ID(f))};           % NOTE: STRIP里的LID是LU的ID
        %         Item.LID(:,iItem) = {unique(t.LID(f))};
        Strip.SID(:,iStrip) = {unique(t.SID(f))};
        Strip.EID(:,iStrip) = {unique(t.EID(f))};
        Strip.PID(:,iStrip) = {unique(t.PID(f))};
    end    
    %  t2 = struct2table(structfun(@(x) x',Strip,'UniformOutput',false));

%% 2: 计算strip装载率Strip.loadingrate,Strip.loadingrateLimit]  cpuBin使用
[Strip.loadingrate,Strip.loadingrateLimit] = computeLoadingRateStrip(Strip,Item,Veh); 

%% 3: 计算STRIP.maxHeight等 : 计算Strip的高度的最大,最小,差值和 均值高度 HStripSW使用
[Strip.maxHeight,Strip.lowestHeight,Strip.meanHeight,Strip.diffHeight] = computeStripHeight(Strip,Item);

%% 4: STRIP.isMixed: STRIP增加判断是否单纯型/混合型判断 HStripSW使用
[Strip.isMixed,Strip.isMixedSID,Strip.isMixedEID] = isMixedStrip(Strip);
% [Strip.nbLID,Strip.nbSID,Strip.nbEID] = computeStripnbIDs(Strip);   % 暂时不用

%% 5: STRIP.isHeightBalance: STRIP增加判断是否isHeightBalance(高度均衡判定):  HStripBalnace使用 
Strip.isHeightBalance = isHeightBalanceStrip(Strip);

% 备用函数,前一个版本
tmpisHeightBalance = isHeightBalanceStrip1(Strip,Item);
if sum(tmpisHeightBalance ~= Strip.isHeightBalance) > 0 
    error('isHeightBalanceStrip计算可能错误'); end

%% 6: STRIP.isHeightFullStrip: STRIP增加判断是否Strip时高度满层/或不满(包含非Full的Item. )
Strip.isHeightFull = isHeightFullStrip(Strip,Item);

%  if sum(~Strip.isHeightBalance & Strip.isHeightFull)
%   plotSolutionT(LU,Veh,0,0,1,1,3,'量大车头后Bin'); % Bin排序后
%  end
 
%% 7: STRIP.isWidthFull : STRIP增加判断是否包含宽度width非Full的Item. 
Strip.isWidthFull = isWidthFullStrip(Strip,Item);

%% 8：cpuStripnbItem
%Strip.nbItem: 混合:-; 单纯: 对应Strip内部该Item的nbLID类型个数,数值越大,即该LU类型越多
[Strip.nbItem, Strip.nbLU, Strip.nbLULID] = cpuStripnbItem(Strip,Item,LU);

[Strip.nLUID, Strip.nLULID] = getStripnID(LU);%此处需要计算，在Strip2bin需要用 


%% 暂时无用 注释 4.1 STRIP.GapValue 和 isGAPBalance(条带的宽度W,沿车厢场边的差值)的计算

% fmix = Strip.isMixed==1;
% fwid = Strip.isWidthFull==0;

% StripCheck 可多次调用
% chkStrip(Strip);

% Calc GapValue
% Strip.GapValue(~fmix) =0;
% Strip.GapValue(fmix) = Strip.maxLULength(fmix) - Strip.lowestLULength(fmix);
% Strip.GapValue(fwid) = Strip.maxLULength(fwid) - 0;

% Calc isGapBalance
% Strip.isGapBalance = Strip.GapValue < 0.33*Strip.maxLULength;

% fgap=Strip.isGapBalance ==0;
% Strip.LID
% printstruct(Strip)

%% 5,6 暂时注释
%Strip.isAllPured：混合:-1; 单纯: 1 (混合strip内没有改ID) ; 0 (混合strip内含有该单纯strip的ID)
%Strip.isSingleItem: 混合: -1; 单纯: Strip内仅有一个Item,必定是单纯的.
% % % LIDinItemsArray = cellfun(@(x) x(1), Item.LID); % arrayAllLID: 所有ITEM对应的LID值 向量形式
% % % mixedStrip = find(Strip.isMixed(1,:) == 1);
% % % mixedLID = [];
% % % for m=1:length(mixedStrip)
% % %      f = Item.Item_Strip(1,:) == mixedStrip(m);
% % %      mixedLID = [mixedLID, LIDinItemsArray(f)];
% % % end
% % % mixedLID = unique(mixedLID);
% % % 
% % % uniItem = unique(Item.Item_Strip(1,:));
% % % for i=1:length(Strip.isAllPured)
% % %     if ~Strip.isMixed(1,i) %如是单纯型        
% % %         cellLID = Item.LID(Item.Item_Strip(1,:) == uniItem(i)); % cellLID: 本Strip内的ITEM对应的LID值
% % %             %         cellLID = Item.LID(Item.Item_Strip(1,:) == i); % cellLID: 本Strip内的ITEM对应的LID值
% % %         LIDinThisItemArray = cellfun(@(x)x(1), cellLID);
% % %         if isscalar(unique(LIDinThisItemArray)) 
% % %             
% % %             if isscalar(LIDinThisItemArray)
% % %                 Strip.isSingleItem(1,i) = 1;
% % %             else
% % %                 Strip.isSingleItem(1,i) = 0;
% % %             end
% % %             
% % %             if ismember(unique(LIDinThisItemArray),mixedLID)
% % %                 Strip.isAllPured(1,i) = 0;
% % %             else
% % %                 Strip.isAllPured(1,i) = 1;
% % %             end
% % %             
% % %             % Strip.nbItem(1,i) = sum(LIDinItemsArray ==
% % %             % unique(LIDinThisItemArray)); %单独放入函数计算
% % %             
% % %         else
% % %              error('单纯型STRIP内的ITEM的类型不同'); %arrayLID            
% % %         end
% % %     else % strip如是混合型     
% % %         Strip.isSingleItem(1,i) = 0;
% % %         Strip.isAllPured(1,i) = 0;
% % %     end
% % % end


end



%% 局部函数 %%

%% 函数1:  V4 isHeightBalanceStrip: 判断STRIP是否包含非HeightFull的Item
% V4: % 改为:Strip.diffHeight
function isHeightBalance = isHeightBalanceStrip(Strip) 
    global parBalance
    
    n = length(Strip.Weight);
    isHeightBalance = deal(ones(1,n)*-1);    
    
    % 2 循环判断Strip内部Item之间的最大差值, 是否<= 最小的对角线或一共绝对值, 如是,均为Full; 
    for i=1:n
  
        % 本strip内最高的Item的1/3
        oneThirdsHeightItem = Strip.maxHeight(i)*parBalance;
        
         % 本strip内的堆垛高度差: 高度间隙  < 本strip内最高的Item的1/3
        if Strip.diffHeight(i) <= oneThirdsHeightItem    
             isHeightBalance(i) = 1;
        else
             isHeightBalance(i) = 0;    
        end
        
    end
    
    % 防错语句
    if any(isHeightBalance==-1), error('存在Strip.isHeightBalance未分配!'); end
    
end

%% V3 isHeightFullStrip: 若高度均衡,通过内堆垛是否高度满层判定 Item.isHeightFull; 否则不满层
% getOrderofLID / cpuBIn / HStripSW 使用(用途广泛) 555
function isHeightFull = isHeightFullStrip(Strip,Item)
    
    n = length(Strip.Weight);
    [isHeightFull] = deal(ones(1,n)*-1);
    
    uniStrip = unique(Item.Item_Strip(1,:));

    % 2 循环判断Strip内部Item之间的最大差值, 是否<= 最小的对角线或一共绝对值, 如是,均为Full; 
    for i=1:n

             %v2: 独立出高度均衡: 如果本STRIP对应ITEM的isFull不是均为1,则本STRIP非full
             % i.e. 如有不满的堆垛属于本strip,认为是高度不满层,可能均衡也可能不均衡,均衡的应该极少
            if ~all(Item.isHeightFull(Item.Item_Strip(1,:) == uniStrip(i)))                 
                isHeightFull(i) = 0;
            else
                isHeightFull(i) = 1;
                if Strip.isHeightBalance(i) == 0
%                     LU.LU_Strip(1,:) == i ??
                    warning('高度满层的条带,居然也会高度不均衡'); % todo fixme 
                end
            end
        
        % V1: 若高度不均衡 , 认为是高度不满层, 独立该对顶
% %         if Strip.isHeightBalance(i) == 0 %如果高度不均衡,一定是高度不满 (不考虑高度都很高,但不均衡情形)
% %             isHeightFull(i) = 0;
% %             
% %         else %如果高度均衡,看Item是否都是满层,如不是,Strip高度也不满
% %             isHeightFull(i) = 1; %Strip单一Item也是高度均衡, 如它的Item是不满, 该Strip也是不满
% %             
% %             if ~all(Item.isHeightFull(Item.Item_Strip(1,:) == uniStrip(i))) %如果本STRIP对应ITEM的isFull不是均为1,则本STRIP非full
% %                 isHeightFull(i) = 0;
% %             end
% %         end
        
    end

    % 防错语句
    if any(isHeightFull==-1), error('存在Strip.isHeightFull未分配!'); end
end

%% 函数2: 判断STRIP是否包含非WidthFull的Item isWidthFullStrip
% ****************** Strip内是否包含Width非Full的Item计算 宽度间隙 >= 本Strip包含的Item的宽度最小值 ************ 开放
function TF = isWidthFullStrip(Strip,Item) 

    n = length(Strip.Weight);
    [TF] = deal(ones(1,n)*-1);   
        
    % 循环判断Strip是否Widthfull, 宽度间隙 >= 本Strip包含的Item的宽度最小值
    uniItem = unique(Item.Item_Strip(1,:));
    
    for i=1:n
        flagItem = Item.Item_Strip(1,:) == uniItem(i);  % max(Item.LWH(1, flagItem));
        
        if Strip.LW(1, i) >= min(Item.LWH(1, flagItem))
            TF(i) = 0;
        else
            TF(i) = 1;
        end
        
    end
    
end


% computeStripnbIDs
% function [nbLID,nbSID,nbEID] = computeStripnbIDs(Strip)
% 
%     n = length(Strip.Weight);
%     [nbLID,nbSID,nbEID] = deal(zeros(1,n));
%         
%     % 循环判断Strip是否为不同LU.ID的混合型
%     for i=1:n        
% 
%          nbLID(i) = numel(Strip.LID{i});
% 
%          nbSID(i) = numel(Strip.SID{i});
%          
%          nbEID(i) = numel(Strip.EID{i});
%     end
%     
% end

%% 函数4: 计算computeLoadingRateStrip
function [loadingrate,loadingrateLimit] = computeLoadingRateStrip(Strip,Item,Veh)
    % 初始化
    nStrip = length(Strip.Weight);

    % 计算每个strip的装载率:
    %每个strip的可用体积 = 高度*宽度(车辆的宽度)
    Stripvolume = Strip.LW(2,:)*Veh.LWH(1,1);
    %每个strip的有限可用体积 = 高度*宽度(strip使用宽度=车辆宽度-strip剩余宽度)
    StripvolumeLimit = Strip.LW(2,:) .* (Veh.LWH(1,1) - Strip.LW(1,:));
    
    a = Item.LWH;
    b = Item.Item_Strip;
    
    uniItem = unique(Item.Item_Strip(1,:));    
    for iStrip =1:nStrip
        %每个strip包含的Item装载体积
        Itemvolume(iStrip)= sum(a(1, (b(1,:)==uniItem(iStrip))) .* a(2, (b(1,:)==uniItem(iStrip))));
    end
    
    %每个strip的装载比率
    loadingrate =  Itemvolume ./ Stripvolume;
    %每个strip的有限装载比率
    loadingrateLimit =  Itemvolume ./ StripvolumeLimit;
end

%% 函数5:   computeStripHeight 计算高度
function [maxHeight,lowestHeight,meanHeight,diffHeight] = computeStripHeight(Strip,Item)

    n = length(Strip.Weight);
    [maxHeight,lowestHeight,meanHeight,diffHeight] = deal(zeros(1,n));
    
    for i=1:n

        % 计算Strip内部高度
        % Item.Item_Strip(1,:) == i) : Strip i 内部的Item flag
        maxHeight(i) = max(Item.LWH(3, Item.Item_Strip(1,:) == i));
        lowestHeight(i) = min(Item.LWH(3, Item.Item_Strip(1,:) == i));
        meanHeight(i) = mean(Item.LWH(3, Item.Item_Strip(1,:) == i)); %         meanHeight(i) = maxHeight(i) - lowestHeight(i);         meanHeight(i) = (maxHeight(i) + lowestHeight(i))/2;
        diffHeight(i) = maxHeight(i) - lowestHeight(i);
        
        if any(meanHeight(i)<0), error('负值不可能;'); end

        % 计算Strip内部 LU 宽度 Width -> 目的 计算Strip的宽度差值,目前无用
        %Strip.maxLULength(i) = max(Item.LWH(2, Item.Item_Strip(1,:) == i));     %似乎为Item长度,目前无用,因为有了混合gap调整
        %Strip.lowestLULength(i) = min(Item.LWH(2, Item.Item_Strip(1,:) == i));
        
    end
end




%% 函数1.2: V3: %  isHeightBalanceStrip1     Strip.diffHeight 备用
function isHeightBalance = isHeightBalanceStrip1(Strip,Item) 
    global parBalance
    
    n = length(Strip.Weight);
    isHeightBalance = deal(zeros(1,n));    
    
    % 1 获取所有的Strip向量
    StripArray = unique(Item.Item_Strip(1,:));

    % 2 循环判断Strip内部Item之间的最大差值, 是否<= 最小的对角线或一共绝对值, 如是,均为Full; 
    for i=1:n
        fItem = Item.Item_Strip(1,:) == StripArray(i);        
        % 本strip内的堆垛高度差: 高度间隙
        maxHeightDiff = max(abs(diff(Item.LWH(3,fItem))));

        % 对比2: 本strip内的最高的Item的1/3
        oneThirdsHeightItem = max(Item.LWH(3,fItem))*parBalance;
        
        % 对比2: 最大差值, 是否<= 1/3最高Item (Item高度平行, 即使很低,  也认为是满层)
        if ~isempty(maxHeightDiff)
            if maxHeightDiff <= oneThirdsHeightItem                 
                     isHeightBalance(i) = 1;
            else
                     isHeightBalance(i) = 0;       
            end
        else %均为空值, 即只有Strip一个堆垛, 算是均衡的
            if sum(fItem)~=1, error('该strip包含不是一个堆垛!'); end
            isHeightBalance(i) = 1;
        end
    end
    
    % 防错语句
    if any(isHeightBalance==-1), error('存在Strip.isHeightBalance未分配!'); end
end


%% 如下全部为注释

%% %% 1 计算LU.LU_Strip, LU.CoordLUStrip  期待拿到外面（或许已经拿到外面）
% % nbLU = size(LU.LWH,2);
% % LU.LU_Strip = zeros(2,nbLU);
% % LU.CoordLUStrip = zeros(3,nbLU);
% % 
% % % 更新LU_Strip
% % for iLU=1:nbLU
% %     % 更新LU_Strip第一行
% %     iItem = LU.LU_Item(1,iLU);   %iLU属于第几个Item, Item属于第几个Strip,则Lu属于第几个Strip
% %     LU.LU_Strip(1,iLU)= Item.Item_Strip(1,iItem);
% %     % 更新LU_Strip第二行
% %     fiItem = find(Item.Item_Strip(1,:) == Item.Item_Strip(1,iItem) & Item.Item_Strip(2,:) < Item.Item_Strip(2,iItem));
% %     nbLUfiItem = sum(ismember(LU.LU_Item(1,:),fiItem));
% %     LU.LU_Strip(2,iLU) = nbLUfiItem+LU.LU_Item(2,iLU); % 进入Strip顺序: 同一Strip内先前进入个数nbLUfiItem + 本iLU在Item的顺序
% %     % 更新LU.CoordLUStrip
% %     LU.CoordLUStrip(1,iLU) = Item.CoordItemStrip(1,iItem);
% %     LU.CoordLUStrip(2,iLU) = Item.CoordItemStrip(2,iItem);
% %         % fLU: 与iLU同属iItem 且 顺序晚于本iLU; 可能为空, 不影响.
% %     fLU = LU.LU_Item(1,:) == iItem & LU.LU_Item(2,:) < LU.LU_Item(2,iLU);
% %     LU.CoordLUStrip(3,iLU) = sum(LU.LWH(3,fLU));
% % end


%% % V2: 修改输出 isHeightBalanceStrip
% % function Strip = isHeightBalanceStrip(Strip,Item) 
% %     global parBalance
% %     
% %     % 1 循环判断Strip是否包含Item为full的,如包含,则Strip为full
% %     StripArray = unique(Item.Item_Strip(1,:));
% % 
% %     % 2 循环判断Strip内部Item之间的最大差值, 是否<= 最小的对角线或一共绝对值, 如是,均为Full; 
% %     for i=1:length(StripArray)
% %         fItem = Item.Item_Strip(1,:) == StripArray(i);        
% %         % 本strip内的堆垛高度差: 高度间隙
% %         maxHeightDiff = max(abs(diff(Item.LWH(3,fItem))));
% %         % 对比2: 本strip内的最高的Item的1/3
% %         oneThirdsHeightItem = max(Item.LWH(3,fItem))*parBalance;
% %         
% %         % 对比2: 最大差值, 是否<= 1/3最高Item (Item高度平行, 即使很低,  也认为是满层)
% %         if ~isempty(maxHeightDiff)
% %             if maxHeightDiff <= oneThirdsHeightItem                 
% %                      Strip.isHeightBalance(i) = 1;
% %             else
% %                      Strip.isHeightBalance(i) = 0;       
% %             end
% %         else %均为空值, 即只有Strip一个堆垛, 算是均衡的
% %             if sum(fItem)~=1, error('该strip包含不是一个堆垛!'); end
% %             Strip.isHeightBalance(i) = 1;
% %         end
% %     end
% %     
% %     % 防错语句
% %     if any(Strip.isHeightBalance==-1), error('存在Strip.isHeightBalance未分配!'); end
% % end

%% V1: Strip是否高度均衡 isHeightBalanceStrip
% % function Strip = isHeightBalanceStrip(Strip,Item) 
% % global parBalance
% %     % 1 循环判断Strip是否包含Item为full的,如包含,则Strip为full
% %     uniStrip = unique(Item.Item_Strip(1,:));
% % 
% %     % 2 循环判断Strip内部Item之间的最大差值, 是否<= 最小的对角线或一共绝对值, 如是,均为Full; 
% %     for i=1:length(Strip.isHeightFull)
% %         fItem = Item.Item_Strip(1,:) == uniStrip(i);        
% %         % 高度间隙
% %         maxHeightDiff = max(abs(diff(Item.LWH(3,fItem))));
% %         % 对比2: 最高的Item的1/3
% %         oneThirdsHeightItem = max(Item.LWH(3,fItem))*parBalance;
% %         
% %         % 对比2: 最大差值, 是否<= 1/3最高Item (Item高度平行, 即使很低,  也认为是满层)
% %         if ~isempty(maxHeightDiff)
% %             if maxHeightDiff <= oneThirdsHeightItem                 
% %                      Strip.isHeightBalance(i) = 1;
% %             else
% %                      Strip.isHeightBalance(i) = 0;       
% %             end
% %         else %均为空值, 即只有Strip一个堆垛, 算是均衡的
% %             if sum(fItem)~=1, error('该strip包含不是一个堆垛!'); end
% %             Strip.isHeightBalance(i) = 1;
% %         end
% %     end
% %     
% %     % 防错语句
% %     if any(Strip.isHeightBalance==-1), error('存在Strip.isHeightBalance未分配!'); end
% % end

%% V1 isHeightFullStrip:不单独考虑isHeightBalance高度均衡, 完全
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

%% V2 isHeightFullStrip:考虑isHeightBalance高度均衡的isHeightFullStrip
% % function Strip = isHeightFullStrip(Strip,Item)
% %     
% %     n = length(Strip.Weight);
% % %     [isMixed,isMixedSID,isMixedEID,nbLID,nbSID,nbEID] = deal(zeros(1,n));
% %     
% %     % 1 循环判断Strip是否包含Item为full的,如包含,则Strip为full
% %     
% %     uniStrip = unique(Item.Item_Strip(1,:));
% % 
% %     % 2 循环判断Strip内部Item之间的最大差值, 是否<= 最小的对角线或一共绝对值, 如是,均为Full; 
% %     for i=1:n
% %         if Strip.isHeightBalance(i) == 0 %如果高度不均衡,一定是高度不满 (不考虑高度都很高,但不均衡情形)
% %             Strip.isHeightFull(i) = 0;
% %         else %如果高度均衡,看Item是否都是满层,如不是,Strip高度也不满
% %             Strip.isHeightFull(i) = 1; %Strip单一Item也是高度均衡, 如它的Item是不满, 该Strip也是不满
% %             if ~all(Item.isHeightFull(Item.Item_Strip(1,:) == uniStrip(i))) %如果本STRIP对应ITEM的isFull不是均为1,则本STRIP非full
% %                 Strip.isHeightFull(i) = 0;
% %             end
% %         end
% %     end
% % 
% %     % 防错语句
% %     if any(Strip.isHeightFull==-1), error('存在Strip.isHeightFull未分配!'); end
% % 
% % end

%% V1 计算STRIP的PID,LID,SID
% % LU.DOC=[LU.DOC; LU.LU_Strip];
% % nStrip = size(Strip.LW,2);
% % for iStrip=1:nStrip
% %     tmp = LU.DOC([1,2,3], LU.DOC(8,:) == iStrip);
% %     Strip.PID(:,iStrip) = num2cell(unique(tmp(1,:))',1);
% %     Strip.LID(:,iStrip) = num2cell(unique(tmp(2,:))',1);
% %     Strip.SID(:,iStrip) = num2cell(unique(tmp(3,:))',1);
% % end
    