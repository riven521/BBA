function [da] = HItemToStrip(da,ParaArray)
% 重要函数:Item放入Strip中 %  行数:长宽高(row);  列数:托盘数量(coloum);
% Input ---  ITEM:  ID LWH Weight
% Output --- ITEM: itemorder itemBeStripMatrix itemRotaFlag CoordItemStrip
% Output --- StripArray: LW Weight
% da.ItemArray (1 LWH (已知)
% da.ItemArray (2 itemBeStripMatrix  (dim1:序号item在某个strip dim2:item进入顺序(左->右) 
% da.ItemArray (3 CoordItemStrip Item在strip的坐标) 
% da.ItemArray (4 itemorder 函数内Item排序顺序)
% da.StripArray (1 LW )

%% 初始化
% nDim Item维度(2) nItem Item数量 nStrip Strip数量 
% widthStrip Strip最大宽度
nDim = size(da.ItemArray.LWH,1);  if nDim ==3, nDim = nDim-1;end
nItem = size(da.ItemArray.LWH,2);
nStrip = nItem;

tmpUniqueBin = unique(da.BinArray.LWH(1:nDim,:)','rows')';
widthStrip = tmpUniqueBin(1);
clear tmpUniqueBin;


    %% ITEM排序 555 
    % getITEMorder - 获取ITEM的顺序(重点是高度递减排序) % ITEM两种排序方式 高度/最短边
    da.ItemArray.itemorder = getITEMorder(); 
    % getSortedITEM - 获取按order排序后的ITEM:sortedItemArray
    sortedItemArray = getSortedITEM(da.ItemArray.itemorder);
%     printstruct(da) ;printstruct(sortedItemArray)    
% da.ItemArray
%% 增对Rotation增加变量 
itemorder = da.ItemArray.itemorder;
% 获取1 LWHItemSort 排序后的ITEM的长宽高 2 itemRotaSort 排序后的FLAG 全部为0 不重要
if ParaArray.whichRotation == 0, LWHItemSort = sortedItemArray.LWH(1:nDim,:); itemRotaSort = zeros(1,size(da.ItemArray.LWH,2)); end %itemRotaSort = sortedItemArray.itemRotaFlag; 
% 获取1 LWHItemSortHori Hori orienting后LWH 2 itemRotaSortHori 原始item是否rotation  0 代表horizontal orientation 1 代表 vertical orientation
if ParaArray.whichRotation == 1,     [LWHItemSortHori,itemRotaSortHori] = horiOrient(sortedItemArray.LWH(1:nDim,:)); end


%% 55 LU->Item->Strip转换 
% 获取itemBeStripMatrixSort: 每个排序后Item在哪个Strip内  以及顺序
% 获取LWStrip:  新生成的Strip的长宽
% 获取CoordItemStripSort  Item在strip的坐标值
LWStrip = zeros(nDim,nItem);   %strip长宽 dim2-长度(以最高的计算) (高度仅做参考,current 高度)
LWStrip(1,:) = widthStrip;   %dim1-宽度剩余 
stripBeItemArray = zeros(1,nStrip);  % 每个Strip内的Item数量 后期不用
itemBeStripMatrixSort = zeros(2,nItem); %dim1:属于第几个level dim2:属于该level第几个排放 555
CoordItemStripSort = zeros(2,nItem); %Item在strip的坐标值

% 55 获取thisLevel - 当前item要放入的level序号
% 循环往strip中安置item,即固定item,变化选择不同level(thisLevel)
% 注释：获取 FLAG        可放下当前item(iStrip)的至少一个level的集合） 
% 注释：获取 thisLevel   从FLAG中找到按规则的那个thisLevel, 并执行 insert函数

iLevel = 1; iStrip = 1; %iStrip代表item实质
while 1
    if iStrip > nItem, break; end
    
    % 依据不同规则找到能放入当前item的strips/levels中的一个
    thisLevel = getThisLevel();     %iLevel会在次函数内不断递增，永远指示当前最新的level
    insertItemToStrip(thisLevel);
    
%     plot2DStrip(); %迭代画图
    
    iStrip = iStrip + 1;
end

% plot2DStrip(); %一次性画图

% 后处理 并赋值到da
        %Matalb code gerator use:
        itemBeStripMatrix=itemBeStripMatrixSort;CoordItemStrip=CoordItemStripSort;

% 获取itemBeStripMatrix : 每个Item在哪个Strip内  以及顺序
% 获取CoordItemStrip : 每个Item在Strip的坐标
% 获取itemRotaFlag : 每个item是否Rotation的标志
% 获取LWStrip:  新生成的strip的长宽
                
    itemBeStripMatrix(:,itemorder) = itemBeStripMatrixSort;
    CoordItemStrip(:,itemorder) = CoordItemStripSort;
    da.ItemArray.itemBeStripMatrix = itemBeStripMatrix;
    da.ItemArray.CoordItemStrip = CoordItemStrip;

    % % % 即使没有RotationFlage 也有该数组 判断是否Rotation
    if ParaArray.whichRotation == 0,   itemRotaFlag(:,itemorder) = itemRotaSort;  end
    if ParaArray.whichRotation == 1,   itemRotaFlag(:,itemorder) = itemRotaSortHori;   end
    da.ItemArray.itemRotaFlag = itemRotaFlag;

    LWStrip = LWStrip(:,LWStrip(2,:)>0);
    da.StripArray.LW = LWStrip; % 去除未使用的Strip
    
    %% 测试script
    % 输出主要结果:获得每个level包含的 
    printscript();
%     dap = rmfield(da, {'BinArray', 'LUArray'});  printstruct(dap);


    %% 嵌套函数        
    function order = getITEMorder()
        tmpLWHItem = da.ItemArray.LWH(1:nDim,:);
        tmpIDItem =     da.ItemArray.ID(1,:);
        if ParaArray.whichSortItemOrder == 1 %Descend of 长(高)
            tmpLWH = [tmpIDItem; tmpLWHItem]; %额外增加ITEM的ID到第一行形成临时变量
            [~,order] = sortrows(tmpLWH',[3 1 ],{'descend','descend'}); %按高度,ID(相同高度时)递减排序            
        end        
        if ParaArray.whichSortItemOrder == 2  %Descend of shortest最短边  -> 增对Rotation增加变量
            tmpLWH = [tmpLWHItem; min(tmpLWHItem(1:nDim,:))]; %额外增加最短边到第三行形成临时变量tmpLWH
            [~,order] = sortrows(tmpLWH',[3 2 1],{'descend','descend','descend'}); %way1 按最短边,高度,宽度递减排序
             %BACKUP  [~,itemorder] = sort(tmpLWH(nDim+1,:),'descend'); %获取临时变量排序后的顺序 way2
        end
        if ~isrow(order), order=order'; end
    end

    function item = getSortedITEM(order)
        item = structfun(@(x) x(:,order),da.ItemArray,'UniformOutput',false);
    end

    function thisLevel = getThisLevel()
        % 不同whichStripH下,获得共同的thisLevel
        if ParaArray.whichStripH == 1 % 1 bestfit 2 firstfit 3 nextfit
            % 增对Rotation增加变量
            if ParaArray.whichRotation == 1
                % 找到可以rotation下的level:任一摆放方向可放入该iItem的level
                flag = find(LWStrip(1,1:iLevel) >= LWHItemSortHori(1,iStrip) |  ...
                    LWStrip(1,1:iLevel) >= LWHItemSortHori(2,iStrip));
            else
                % 常规条件下的选择：find宽度足够的多个level,并安置在最小剩余水平宽度的
                flag = find(LWStrip(1,1:iLevel) >= LWHItemSort(1,iStrip));
            end
            if isempty(flag)
                iLevel = iLevel + 1;% 如果宽度不满足，则level升级
                thisLevel = getThisLevel();
            else
                if length(flag)>1
                    1
                end
                % 获取thisLevel: 唯一与FF区别从这到thisLevel的计算（选中满足条件且最小的）
                tepAvailableLevelArray = LWStrip(1,1:iLevel);   %获取所有已安排或新安排的level的剩余水平平宽度向量tepAvailableLevelArray
                tepMinLeftWdith = min(tepAvailableLevelArray(flag));                      %找出tepAvailableLevelArray中可容纳本iITem的剩余水平宽度的最小值，更新为tepMinLeftWdith
                thisLevel = find(tepAvailableLevelArray==tepMinLeftWdith);            %找出与最小值对应的那个/些level
                if ~all(ismember(thisLevel,flag)),           error('Not all thisLevel belongs to flag ');          end
                if length(thisLevel)>1
                    thisLevel = thisLevel(1);
                end
            end
        elseif ParaArray.whichStripH == 2 % firstfit  % firstfit下能不能直接套用bestfit的代码?
            % 增对Rotation增加变量
            if ParaArray.whichRotation == 1
                % 找到可以rotation下的level:任一摆放方向可放入该iItem的level
                flag = find(LWStrip(1,1:iLevel) >= LWHItemSortHori(1,iStrip) |  ...
                    LWStrip(1,1:iLevel) >= LWHItemSortHori(2,iStrip));
            else
                % 常规条件下的选择：find宽度足够的多个level,并安置在第一个遇到的 唯一区别是thisLevel的获取
                flag = find(LWStrip(1,1:iLevel) >= LWHItemSort(1,iStrip));
            end
            if isempty(flag)
                iLevel = iLevel + 1;% 如果宽度不满足，则level升级
                 thisLevel = getThisLevel();
            else
              if length(flag)>1
                    1
                end
                thisLevel = flag(1);
            end
        elseif ParaArray.whichStripH == 3 % nextfit
            % 增对Rotation增加变量
            if ParaArray.whichRotation == 1 % nextfit下不能直接套用bestfit的代码 FIXME
                % 找到可以rotation下的level:任一摆放方向可放入该iItem的level
                flag = find(LWStrip(1,iLevel) >= LWHItemSortHori(1,iStrip) |  ...
                    LWStrip(1,iLevel) >= LWHItemSortHori(2,iStrip));
            else
                % 不同条件下的选择：如果当前item的宽<=当前strip的当前level的宽
                flag = find(LWStrip(1,iLevel) >= LWHItemSort(1,iStrip) );
            end
            if  isempty(flag) %注意与之前~flag的区别 FIXME ~flag
                iLevel = iLevel + 1;% 如果宽度不满足，则level升级
                thisLevel = getThisLevel();
            else
                if length(flag)>1
                    1
                end
                flag
                iLevel
                if ~all(ismember(iLevel,flag)),        
                    warning('Not all thisLevel belongs to flag ');     
                end
                thisLevel = iLevel; % 即为
            end
        end
    end

    function insertItemToStrip(thisLevel)
        CoordItemStripSort=CoordItemStripSort;LWStrip=LWStrip;
        itemRotaSortHori=itemRotaSortHori;LWHItemSortHori=LWHItemSortHori;
        stripBeItemArray=stripBeItemArray;itemBeStripMatrixSort=itemBeStripMatrixSort;
        
        % 1 更新CoordItemStripSort
        CoordItemStripSort(1,iStrip) = widthStrip - LWStrip(1,thisLevel);  %更新x坐标
        CoordItemStripSort(2,iStrip) = sum(LWStrip(2,1:thisLevel-1));      %更新y坐标 %如果iLevel=1,长（高）坐标为0；否则为求和
        
        % 2 更新stripBeItemArray
        % 2 更新LWStrip
        if ParaArray.whichRotation == 1 && ParaArray.whichRotationHori == 0
            % 判断语句
            flagHori = LWStrip(1,thisLevel) >=  LWHItemSortHori(1,iStrip); %判断是否 wleft >= longest
            flagVert = LWStrip(1,thisLevel) >= LWHItemSortHori(2,iStrip);  %判断是否 wleft >= shortest
            if ~flagVert & ~flagVert
                LWStrip(1,thisLevel)
                LWHItemSortHori(1,iStrip)
                LWHItemSortHori(2,iStrip)
            end
            isNewLevel = LWStrip(1,thisLevel) == widthStrip; % 判断是否 new Level
            flagLargerLengthHori = LWStrip(2,thisLevel) <  LWHItemSortHori(2,iStrip); %判断是否 要更新level长/高度
            flagLargerLengthVert = LWStrip(2,thisLevel) <  LWHItemSortHori(1,iStrip); %判断是否 要更新level长/高度
            
            % 更新strip信息
            % 如果新level,优先按照horizontally方式(即默认方式安置)
            if isNewLevel
                if flagHori
                    LWStrip(1,thisLevel) = LWStrip(1,thisLevel) - LWHItemSortHori(1,iStrip); %更新wleft by Horizontally
                    if flagLargerLengthHori
                        LWStrip(2,thisLevel) = LWHItemSortHori(2,iStrip); %更新strip高度lleft
                    end
                elseif flagVert
                    LWStrip(1,thisLevel) = LWStrip(1,thisLevel) ...
                        - LWHItemSortHori(2,iStrip); %更新wleft by Vertically
                    if flagLargerLengthVert
                        LWStrip(2,thisLevel) = LWHItemSortHori(1,iStrip); %更新strip高度lleft
                    end
                    % 当flagVert==1 不仅标记 还要把物品真正的rotation(反)过去
                    itemRotaSortHori(iStrip) = ~itemRotaSortHori(iStrip);
                    tep = LWHItemSortHori(1,iStrip);
                    LWHItemSortHori(1,iStrip) = LWHItemSortHori(2,iStrip);
                    LWHItemSortHori(2,iStrip) = tep;
                else      %如果wleft <shorter edge
                    error('数据错误,无论横竖都放不下');
                end
            else  % 如果非level,优先按照vertically方式放
                if flagVert
                    LWStrip(1,thisLevel) = LWStrip(1,thisLevel) ...
                        - LWHItemSortHori(2,iStrip); %更新wleft by Vertically
                    if flagLargerLengthVert
                        LWStrip(2,thisLevel) = LWHItemSortHori(1,iStrip); %更新strip高度lleft
                    end
                    % 当flagVert==1 不仅标记 还要把物品真正的rotation(反)过去
                    itemRotaSortHori(iStrip) = ~itemRotaSortHori(iStrip);
                    tep = LWHItemSortHori(1,iStrip);
                    LWHItemSortHori(1,iStrip) = LWHItemSortHori(2,iStrip);
                    LWHItemSortHori(2,iStrip) = tep;
                elseif flagHori
                    LWStrip(1,thisLevel) = LWStrip(1,thisLevel) ...
                        - LWHItemSortHori(1,iStrip); %更新wleft by Horizontally
                    if flagLargerLengthHori
                        LWStrip(2,thisLevel) = LWHItemSortHori(2,iStrip); %更新strip高度lleft
                    end
                else
                    error('level选择错误,无论横竖都放不下');
                end               
            end
            stripBeItemArray(thisLevel) = stripBeItemArray(thisLevel) + 1; %只要该level安放一个item,数量就增加1
        end
        
        if ParaArray.whichRotation == 1 && ParaArray.whichRotationHori == 1
           % 判断语句
            flagHori = LWStrip(1,thisLevel) >=  LWHItemSortHori(1,iStrip); %判断是否 wleft >= longest
            flagVert = LWStrip(1,thisLevel) >= LWHItemSortHori(2,iStrip);  %判断是否 wleft >= shortest
            isNewLevel = LWStrip(1,thisLevel) == widthStrip; % 判断是否 new Level
            flagLargerLengthHori = LWStrip(2,thisLevel) <  LWHItemSortHori(2,iStrip); %判断是否 要更新level长/高度
            flagLargerLengthVert = LWStrip(2,thisLevel) <  LWHItemSortHori(1,iStrip); %判断是否 要更新level长/高度
            
            % 更新strip信息
            % 如果新level,优先按照horizontally方式(即默认方式安置)
            if isNewLevel
                if flagHori
                    LWStrip(1,thisLevel) = LWStrip(1,thisLevel) - LWHItemSortHori(1,iStrip); %更新wleft by Horizontally
                    if flagLargerLengthHori
                        LWStrip(2,thisLevel) = LWHItemSortHori(2,iStrip); %更新strip高度lleft
                    end
                elseif flagVert
                    LWStrip(1,thisLevel) = LWStrip(1,thisLevel) ...
                        - LWHItemSortHori(2,iStrip); %更新wleft by Vertically
                    if flagLargerLengthVert
                        LWStrip(2,thisLevel) = LWHItemSortHori(1,iStrip); %更新strip高度lleft
                    end
                    % 当flagVert==1 不仅标记 还要把物品真正的rotation(反)过去
                    itemRotaSortHori(iStrip) = ~itemRotaSortHori(iStrip);
                    tep = LWHItemSortHori(1,iStrip);
                    LWHItemSortHori(1,iStrip) = LWHItemSortHori(2,iStrip);
                    LWHItemSortHori(2,iStrip) = tep;
                else      %如果wleft <shorter edge
                    error('数据错误,无论横竖都放不下');
                end
            else  % 如果非level,优先按照vertically方式放
                if flagHori
                    LWStrip(1,thisLevel) = LWStrip(1,thisLevel) ...
                        - LWHItemSortHori(1,iStrip); %更新wleft by Horizontally
                    if flagLargerLengthHori
                        LWStrip(2,thisLevel) = LWHItemSortHori(2,iStrip); %更新strip高度lleft
                    end                    
                elseif flagVert
                    LWStrip(1,thisLevel) = LWStrip(1,thisLevel) ...
                        - LWHItemSortHori(2,iStrip); %更新wleft by Vertically
                    if flagLargerLengthVert
                        LWStrip(2,thisLevel) = LWHItemSortHori(1,iStrip); %更新strip高度lleft
                    end
                    % 当flagVert==1 不仅标记 还要把物品真正的rotation(反)过去
                    itemRotaSortHori(iStrip) = ~itemRotaSortHori(iStrip);
                    tep = LWHItemSortHori(1,iStrip);
                    LWHItemSortHori(1,iStrip) = LWHItemSortHori(2,iStrip);
                    LWHItemSortHori(2,iStrip) = tep;
                else
                    error('level选择错误,无论横竖都放不下');
                end
            end
            stripBeItemArray(thisLevel) = stripBeItemArray(thisLevel) + 1; %只要该level安放一个item,数量就增加1            
        end
        
        if ParaArray.whichRotation == 1 && ParaArray.whichRotationHori == 2
           % 判断语句
            flagHori = LWStrip(1,thisLevel) >=  LWHItemSortHori(1,iStrip); %判断是否 wleft >= longest
            flagVert = LWStrip(1,thisLevel) >= LWHItemSortHori(2,iStrip);  %判断是否 wleft >= shortest
            isNewLevel = LWStrip(1,thisLevel) == widthStrip; % 判断是否 new Level
            flagLargerLengthHori = LWStrip(2,thisLevel) <  LWHItemSortHori(2,iStrip); %判断是否 要更新level长/高度
            flagLargerLengthVert = LWStrip(2,thisLevel) <  LWHItemSortHori(1,iStrip); %判断是否 要更新level长/高度
            
            % 更新strip信息
            % 如果新level,优先按照horizontally方式(即默认方式安置)
            if isNewLevel
                if flagVert
                    LWStrip(1,thisLevel) = LWStrip(1,thisLevel) ...
                        - LWHItemSortHori(2,iStrip); %更新wleft by Vertically
                    if flagLargerLengthVert
                        LWStrip(2,thisLevel) = LWHItemSortHori(1,iStrip); %更新strip高度lleft
                    end
                    % 当flagVert==1 不仅标记 还要把物品真正的rotation(反)过去
                    itemRotaSortHori(iStrip) = ~itemRotaSortHori(iStrip);
                    tep = LWHItemSortHori(1,iStrip);
                    LWHItemSortHori(1,iStrip) = LWHItemSortHori(2,iStrip);
                    LWHItemSortHori(2,iStrip) = tep;
                elseif flagHori
                    LWStrip(1,thisLevel) = LWStrip(1,thisLevel) - LWHItemSortHori(1,iStrip); %更新wleft by Horizontally
                    if flagLargerLengthHori
                        LWStrip(2,thisLevel) = LWHItemSortHori(2,iStrip); %更新strip高度lleft
                    end
                else      %如果wleft <shorter edge
                    error('数据错误,无论横竖都放不下');
                end
            else  % 如果非level,优先按照vertically方式放
                if flagVert
                    LWStrip(1,thisLevel) = LWStrip(1,thisLevel) ...
                        - LWHItemSortHori(2,iStrip); %更新wleft by Vertically
                    if flagLargerLengthVert
                        LWStrip(2,thisLevel) = LWHItemSortHori(1,iStrip); %更新strip高度lleft
                    end
                    % 当flagVert==1 不仅标记 还要把物品真正的rotation(反)过去
                    itemRotaSortHori(iStrip) = ~itemRotaSortHori(iStrip);
                    tep = LWHItemSortHori(1,iStrip);
                    LWHItemSortHori(1,iStrip) = LWHItemSortHori(2,iStrip);
                    LWHItemSortHori(2,iStrip) = tep;
                elseif flagHori
                    LWStrip(1,thisLevel) = LWStrip(1,thisLevel) ...
                        - LWHItemSortHori(1,iStrip); %更新wleft by Horizontally
                    if flagLargerLengthHori
                        LWStrip(2,thisLevel) = LWHItemSortHori(2,iStrip); %更新strip高度lleft
                    end                    
                else
                    error('level选择错误,无论横竖都放不下');
                end
            end
            stripBeItemArray(thisLevel) = stripBeItemArray(thisLevel) + 1; %只要该level安放一个item,数量就增加1            
        end
        
        if ParaArray.whichRotation == 0
            % 555更新strip信息
            LWStrip(1,thisLevel) = LWStrip(1,thisLevel) - LWHItemSort(1,iStrip); %更新wleft
            if LWHItemSort(2,iStrip) > LWStrip(2,thisLevel)         % 如果物品高度>本strip高度
                LWStrip(2,thisLevel) = LWHItemSort(2,iStrip); %更新strip高度lleft
            end
            stripBeItemArray(thisLevel) = stripBeItemArray(thisLevel) + 1; %只要该level安放一个item,数量就增加1
        end

        % 3 更新item归属strip信息
        itemBeStripMatrixSort(1,iStrip) = thisLevel;    %第几个level
        itemBeStripMatrixSort(2,iStrip) = stripBeItemArray(thisLevel); %本level下第几次安置
    end
    
    function printscript()
        % 测试代码
        % % LWHStrip
        % % itemBeStripMatrixSort
        % % LWHItemSort
        % % itemCoordMatrixSort
        % % itemBeStripMatrix
        % % LWHItem
        % % itemCoordMatrix
        %  printstruct(da);
        
        % 输出主要结果:获得从1开始每个strip包含的数据
        for iStrip = 1:max(da.ItemArray.itemBeStripMatrix(1,:))
            [~,idx] = find(da.ItemArray.itemBeStripMatrix(1,:)==iStrip);
            fprintf('strip %d 的剩余宽+最大长为:  ',iStrip);
            fprintf('( %d ) ',da.StripArray.LW(:,iStrip));
            fprintf('\n');
            fprintf('strip %d 包含 original Item 索引号(长宽)[旋转标志]{坐标}为  \n  ',iStrip);
            fprintf('%d ',idx);
            fprintf('( %d ) ', da.ItemArray.LWH(1:nDim,idx));fprintf('\n');
            fprintf('[ %d ] ', da.ItemArray.itemRotaFlag(:,idx));fprintf('\n');
            fprintf('{ %d } ', da.ItemArray.CoordItemStrip(:,idx));fprintf('\n');
            fprintf('\n');
        end
    end

    function plot2DStrip()
        %% 初始化
        wStrip = widthStrip;        
        hStrip = sum(LWStrip(2,itemBeStripMatrixSort(2,:)>0));        
        nstrip = sum(itemBeStripMatrixSort(2,:)>0);

        nIDType = unique(sortedItemArray.ID);
        nColors = hsv(length(nIDType)); %不同类型LU赋予不同颜色
        
        %% 画图
        % 1 画图：画本次Strip
        DrawRectangle([wStrip/2 hStrip/2 wStrip hStrip 0],'--');
        hold on;
        % 2 画图：逐个strip/item 画图
        for istrip = 1:nstrip
            % 找出当前istrip的物品索引
            idxDrawItem = find(itemBeStripMatrixSort(1,:)==istrip);
            % 获取该索引下的变量
            drawItemCoordMatrix = CoordItemStripSort(:,idxDrawItem);
            drawItemLWH = LWHItemSortHori(:,idxDrawItem);
            drawItemId = sortedItemArray.ID(:,idxDrawItem);

            % 画图：逐个item
            nThisItem = size(drawItemLWH,2);
            for iItem = 1:nThisItem
                % 画图：画本次iItem
                itemWidth = drawItemLWH(1,iItem);
                itemLength = drawItemLWH(2,iItem);
                itemCenter = [drawItemCoordMatrix(1,iItem)+itemWidth/2 ...
                    drawItemCoordMatrix(2,iItem)+itemLength/2 ];
                
                % 增加对rotation的判断 内部嵌套函数不需要
                %             if ParaArray.whichRotation == 1
                %                 drawItemRotaMatrix = itemRotaSortHori(:,idxDrawItem); %增加rotation后增
                %             end
                
                % %         if ParaArray.whichRotation == 1 && drawItemRotaMatrix(iItem)
                % %             itemWidth = drawItemLWH(2,iItem);
                % %             itemLength = drawItemLWH(1,iItem);
                % %             itemCenter = [drawItemCoordMatrix(1,iItem)+itemWidth/2 ...
                % %                         drawItemCoordMatrix(2,iItem)+itemLength/2 ];
                % %         end
                
                % 增加对本次iItem的类型（颜色）判断
                itemID = drawItemId(iItem);
                itemColor = 0.8*nColors(nIDType==itemID, : );
                
                DrawRectangle([itemCenter itemWidth itemLength 0],  '-',itemColor);
                hold on;
            end
        end
        % hold off;
    end
end

