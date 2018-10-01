%% GET STRIP 相关属性
% 1 Strip.isMixed % 1：混合层； 0：单纯层
% 2 Strip.isHeightFull % 1：全部都是满层； 0：包含非满层 依据Item.isHeightFull 判断 
% 3 Strip.isWidthFull % % 1：全部都是宽度满层； 0：宽度非满层 依据 Strip宽度剩余 是否 可容纳 最小的Item宽度判断 
% 4 Strip.maxHeight % Strip内所有Item的最大高度

% 5 Strip.nbItem % 整数：冗余值, 具体ITEM的堆垛个数 车头摆放依据 -1：混合strip
% 6 Strip.isAllPured % ：1：单纯且该ITEM没有混合型； 0：单纯但也由混合的； -1：混合strip
% 7 Strip.isSingleItem %: 1:单纯且单个； 0：单纯且多个；-1：混合strip

% 8 Strip.loadingrate   % 每个strip的装载比率
% 9 Strip.loadingrateLimit   % 每个strip的有限装载比率

% Strip.seqSW
% Strip.LID

%% 函数
function   [Strip] = cpuStrip(Strip,Item,Veh)
%% 初始化
    Strip.isMixed = ones(size(Strip.Weight))*-1;   %是否为混合型,包含多个LID 
    Strip.isHeightFull = ones(size(Strip.Weight))*-1;   %是否包含非Full的Item
    Strip.nbItem = ones(size(Strip.Weight))*-1;   %单STRIP内部ITEM类型个数, 混合型默认为-1
    Strip.isAllPured = ones(size(Strip.Weight))*-1;   %单STRIP对应LID是否包含混合STRIP, 包含混合型默认为-1
    Strip.isSingleItem = ones(size(Strip.Weight))*-1;   %单Strip内对应只有1个ITEM
    Strip.isWidthFull = ones(size(Strip.Weight))*-1;     %是否为宽度非Full的Item
    Strip.maxHeight = ones(size(Strip.Weight))*-1;     %Strip的最高高度.
%     Strip.seqSW = ones(size(Strip.Weight))*-1;     %Strip的??? 暂未用

    Strip.Stripvolume = ones(size(Strip.Weight))*-1;  %每个strip的可用体积 = 高度*宽度(车辆的宽度)
    Strip.StripvolumeLimit = ones(size(Strip.Weight))*-1; % %每个strip的有限可用体积 = 高度*宽度(strip使用宽度=车辆宽度-strip剩余宽度)
    Strip.Itemvolume = ones(size(Strip.Weight))*-1; %每个strip包含的Item装载体积
    Strip.loadingrate = ones(size(Strip.Weight))*-1; % 每个strip的装载比率
    Strip.loadingrateLimit = ones(size(Strip.Weight))*-1;     % 每个strip的有限装载比率
    
%% 0: 计算strip装载率
Strip = computeLoadingRateStrip(Strip,Item,Veh); 
    
%% 1: STRIP.isMixed: STRIP增加判断是否单纯型/混合型判断
Strip = isMixedStrip(Strip);

%% 2: STRIP.isHeightFull: STRIP增加判断是否包含非Full的Item. 
Strip = isFullStrip(Strip,Item);

%% 3: STRIP.isWidthFull : STRIP增加判断是否包含宽度width非Full的Item. 
Strip = isWidthFullStrip(Strip,Item);

%% 4: STRIP.maxHeight : 计算Strip的最大高度
for i=1:length(Strip.maxHeight)
    % 计算最大值
    % Item.Item_Strip(1,:) == i) : Strip i 内部的Item flag
    Strip.maxHeight(i) = max(Item.LWH(3, Item.Item_Strip(1,:) == i));
end

%% 5,6,7
%Strip.isAllPured：单STRIP对应LID是否包含混合STRIP, 包含混合型默认为-1
%Strip.nbItem: STRIP增加内部该ITEM的nbLID类型个数,数值越大,即该LU类型越多
%Strip.isSingleItem: Strip
arrayAllLID = cellfun(@(x) x(1), Item.LID); % arrayAllLID: 所有ITEM对应的LID值 向量形式
mixedStrip = find(Strip.isMixed(1,:) == 1);
mixedLID = [];
for m=1:length(mixedStrip)
     flagitemIdx = Item.Item_Strip(1,:) == mixedStrip(m);
     mixedLID = [mixedLID, arrayAllLID(flagitemIdx)];
end
mixedLID = unique(mixedLID);
for i=1:length(Strip.isAllPured)
    if ~Strip.isMixed(1,i) %如是单纯型        
        cellLID = Item.LID(Item.Item_Strip(1,:) == i); % cellLID: 本Strip内的ITEM对应的LID值
        arrayLID = cellfun(@(x)x(1), cellLID);
        if isscalar(unique(arrayLID)) 
            if isscalar(arrayLID)
                Strip.isSingleItem(1,i) = 1;
            else
                Strip.isSingleItem(1,i) = 0;
            end
            if ismember(unique(arrayLID),mixedLID)
                Strip.isAllPured(1,i) = 0;
            else
                Strip.isAllPured(1,i) = 1;
            end
            Strip.nbItem(1,i) = sum(arrayAllLID == unique(arrayLID));
        else
             error('单纯型STRIP内的ITEM的类型不同'); %arrayLID            
        end
    end
end


end

%% 局部函数 %%

%% 函数1: 判断STRIP是否包含非HeightFull的Item
function Strip = isFullStrip(Strip,Item) 
    % 循环判断Strip是否full
    for i=1:length(Strip.isHeightFull)
         if all(Item.isHeightFull(Item.Item_Strip(1,:) == i)) %如果本STRIP对应ITEM的isFull均为1,则本STRIP也为full
             Strip.isHeightFull(i) = 1;
         else
             Strip.isHeightFull(i) = 0;
         end
    end    
end

%% 函数2: 判断STRIP是否包含非WightFull的Item
% ****************** Strip内是否包含Width非Full的Item计算 ************ 开放
function Strip = isWidthFullStrip(Strip,Item) 
    % 循环判断Strip是否Widthfull, 宽度间隙 >= 本Strip包含的Item的宽度最小值
    for i=1:length(Strip.isWidthFull)
        flagItem = Item.Item_Strip(1,:) == i;        
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
    for iStrip =1:nStrip
        %每个strip包含的Item装载体积
        Strip.Itemvolume(iStrip)= sum(a(1, (b(1,:)==iStrip)) .* a(2, (b(1,:)==iStrip)));
    end
    %每个strip的装载比率
    Strip.loadingrate =  Strip.Itemvolume ./ Strip.Stripvolume;
    %每个strip的有限装载比率
    Strip.loadingrateLimit =  Strip.Itemvolume ./ Strip.StripvolumeLimit;
end
