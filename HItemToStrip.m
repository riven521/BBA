% function [d] = HItemToStrip(d,p)
function [LU,Item,Strip]= HItemToStrip(LU,Item,Veh,p)
% 重要函数:Item放入Strip中 %  行数:长宽高(row);  列数:托盘数量(coloum);
% Input ---  ITEM:  ID LWH Weight
% Output --- ITEM: itemorder Item_Strip itemRotaFlag CoordItemStrip
% Output --- Strip: LW Weight
% Item (1 LWH (已知)
% Item (2 Item_Strip  (dim1:序号item在某个strip dim2:item进入顺序(左->右) 
% Item (3 CoordItemStrip Item在strip的坐标) 
% Item (4 itemo rder 函数内Item排序顺序)
% Strip (1 LW )
% tmpStrip_Item         (2,n): % 行1：每个Strip内的Item数量 ； 行2：每个Strip内的不同LUID数量
%% 初始化
% nDim Item维度(2) nItem Item数量 nStrip Strip数量 
% widthStrip Strip最大宽度
nDim = size(Item.LWH,1); if nDim ==3, nDim = nDim-1;end
sz = size(Item.LWH);
nItem = sz(2);
wStrip = Veh.LWH(1,1);

%% 先判断ITEM是以Horizontal/Vertical 方式摆放（连带是否旋转）；再判断进入算法的顺序
% 无论是否允许旋转, 只看是否需要以Horizontal/Vertical方式摆放

%  [ItemLWRota, ItemRotaed] = placeItemHori(Item,1);  %第二个参数：1: Hori; 0: Vert；其它: 原封不动
% 获得原封不动的返回值:赋初始值
% x = Item.Rotaed
% [Item.Rotaed] = placeItemHori(Item.LWH,Item.isRota,2);  %第二个参数：1: Hori; 0: Vert；其它: 原封不动
%         if any(x~=Item.Rotaed),             error('111111111111');         end
% Item.LWH = getRotaedLWH(Item.LWH, Item.Rotaed, LU.buff); 

    %% ITEM排序 555
    % 获取Item的顺序 % ITEM两种排序方式 高度/最短边
%     printstruct(Item)    
    [Item.itemorder] = getITEMorder(Item,p.whichSortItemOrder );
    
    % 获取按order排序后的ITEM: sItem
    if isSameCol(Item)
        sItem = structfun(@(x) x(:,Item.itemorder),Item,'UniformOutput',false);
    else
        error('不能使用structfun');
    end
%     printstruct(sItem)

    % 下面用途不大，主要原因在于在Gpreproc中以及做了H/V放置处理了
% %     % 1和2以内的, sortedItemArray对应的LWHRota和Rotaed更新了->需求返回到原矩阵ItemArry中
% %     if p.whichRotationHori == 1 % 无论哪个level,都按照horizontally方式摆放
% %          x = sItem.Rotaed;
% %          sItem.LWH
% %         [ sItem.Rotaed] = placeItemHori(sItem.LWH,sItem.isRota,1);  %第二个参数：1: Hori; 0: Vert；其它: 原封不动       
% % %          if any(x~=sItem.Rotaed),                 error('111111111111');         end    
% %     end
% %     if p.whichRotationHori == 2 % 无论哪个level,都按照vertical方式摆放
% %                     x = sItem.Rotaed
% %         [ sItem.Rotaed] = placeItemHori(sItem.LWH,sItem.isRota,0);  %第二个参数：1: Hori; 0: Vert；其它: 原封不动
% % %                     if any(x~=sItem.Rotaed),                                  error('111111111111');         end
% %     end
% %     sItem.LWH = getRotaedLWH(sItem.LWH, sItem.Rotaed, LU.buff); 
% %     sItem.LWH
%% 55 LU->Item->Strip转换 
% 提取变量 此处只使用LWHRota和Rotaed; 不使用LWH
% ItemLWRotaSort = sItem.LWH(1:2,:); %ItemLWSortHori
% ItemRotaedSort = sItem.Rotaed; % ItemRotaSortHori % itemRotaSort = zeros(1,size(Item.LWH,2));
% ItemisRotaSort = sItem.isRota;
% ItemWeightSort = sItem.Weight;

% Itemorder = Item.itemorder;

%%
% 1 Strip初始化
% 获取Strip.LW:  新生成的Strip的长宽等
Strip.LW = zeros(2,nItem);   %strip长宽 dim2-长度(以最高的计算) (高度仅做参考,current 高度)
Strip.LW(1,:) = wStrip;   %dim1-宽度剩余 
Strip.Weight = zeros(1,nItem); % 初始赋值
    
% 2  临时
tmpStrip_Item = zeros(2,nItem);  % 行1：每个Strip内的Item数量 ； 行2：每个Strip内的不同LUID数量
% 3 sItem新增
% 获取CoordItemStripSort  Item在strip的坐标值
% 获取Item_Strip: 每个排序后Item在哪个Strip内  以及顺序
sItem.Item_Strip = zeros(2,nItem); %dim1:属于第几个level dim2:属于该level第几个排放 555
sItem.CoordItemStrip = zeros(2,nItem); %Item在strip的坐标值

% 55 获取thisLevel - 当前item要放入的level序号
% 循环往strip中安置item,即固定item,变化选择不同level(thisLevel)
% 注释：获取 FLAG        可放下当前item(iItem)的至少一个level的集合） 
% 注释：获取 thisLevel   从FLAG中找到按规则的那个thisLevel, 并执行 insert函数

iLevel = 1; iItem = 1; %iStrip代表item实质

while 1
    if iItem > nItem, break; end

    % 依据不同规则找到能放入当前item的strips/levels中的一个    
    [thisLevel,iLevel] = getThisLevel(iItem,iLevel,sItem, Strip, p);     %iLevel会在次函数内不断递增，永远指示当前最新的level
    
    insertItemToStrip(thisLevel,iItem);
    %     plot2DStrip(); %迭代画图    
    iItem = iItem + 1;
end
    
% figure(randi(1000));
% plot2DStrip();  %可能有问题: 一次性画图
figure(randi(1000));
plot3DStrip();

% Item相关：更新的按顺序返回（无更新的不需返回）
% LU内部更新,sLU依据order变化回来(主要为了sLU中新增的几个变量,要按顺序转回来)
% 获取Rotaed : 每个item是否Rotation的标志
% 获取LWHRota：每个item结合Rotaed标志后获得的LWH
if isSameCol(sItem)
    Item = reorderStruct(Item.itemorder, sItem);
else
    error('不能使用structfun');
end

                    % 获取Item_Strip : 每个Item在哪个Strip内  以及顺序
                    % 获取CoordItemStrip : 每个Item在Strip的坐标
                    %     Item.Item_Strip(:,Item.itemorder) = sItem.Item_Strip;
                    %     Item.CoordItemStrip(:,Item.itemorder) = sItem.CoordItemStrip;    
                    % ItemArray旋转相关
                    %         Item.Rotaed(:,Item.itemorder) = sItem.Rotaed;
                    %     ItemLWRota(:,Item.itemorder) = sItem.LWH(1:2,:);
                    %     Item.LWH = [ItemLWRota; Item.LWH(3,:)];             % 返回原始顺序的旋转后的ItemArray


                    % LUArray旋转相关,及时更新    已嵌入到rotateItem
                %     nbItem=length(Item.Rotaed);
                %     % 循环每个item
                %     for idxItem=1:nbItem
                %         tmpflagThisItem = (LU.LU_Item(1,:)==idxItem );
                %         % 对应位置LU.Rotaed更新
                %         if Item.Rotaed(idxItem)
                %             LU.Rotaed(tmpflagThisItem) = ~LU.Rotaed(tmpflagThisItem);
                %             % 对应位置LU.LWH更新
                %             LU.LWH(1, tmpflagThisItem) = Item.LWH(1, idxItem);
                %             LU.LWH(2, tmpflagThisItem) = Item.LWH(2, idxItem);
                %         end
                %     end
    
% Strip相关: 无顺序概念(去除初始化的多余列)
% Strip.LW
% 获取LWStrip:  新生成的strip的长宽
% 获取StripWeight:  新生成的strip的重量
% 如果Strip的列数全部相同
if isSameCol(Strip)
    Strip = structfun(@(x) x( : , Strip.Weight(1,:)>0 ), Strip, 'UniformOutput', false);
else
    error('不能使用structfun');
end
                 

%% 测试script
    % 输出主要结果:获得每个level包含的 
%     printscript();
%     printstruct(d);
    
    %% 嵌套函数        
    function insertItemToStrip(thisLevel,iItem)
%         %为了matlab的coder 无用
%         sItem.CoordItemStrip=sItem.CoordItemStrip;Strip.LW=Strip.LW;
%         ItemRotaSort=ItemRotaSort;ItemLWSort=ItemLWSort;
%         Strip_Item=Strip_Item;sItem.Item_Strip=sItem.Item_Strip;
        
        % 1 更新Item相关Sort数据
        %  1.1 更新CoordItemStripSort
        sItem.CoordItemStrip(1,iItem) = wStrip - Strip.LW(1,thisLevel);        %更新x坐标
        sItem.CoordItemStrip(2,iItem) = sum(Strip.LW(2,1:thisLevel-1));      %更新y坐标 %如果iLevel=1,长（高）坐标为0；否则为求和
        
        % 2 更新Strip相关数据（非排序）
        %  2.1 更新LWStrip
        updateLWStrip(); %不在判断是否允许旋转；不再判断是否属于新Level；不再判断是否当前可放入
                % % %         if sItem.isRota(iItem) == 1 %此Item可以旋转
                % % %             % 判断语句: 2个
                % % %             isflagCurr = Strip.LW(1,thisLevel) >=  sItem.LWH(1,iItem); %判断是否current's strip剩余宽度 >= 当前高度（非旋转）
                % % %             isNewLevel = Strip.LW(1,thisLevel) == wStrip; % 判断是否 new Level            
                % % %             % 更新strip信息
                % % %             if isNewLevel %无论如何,均可以放入,无论何种摆放,因此:直接更新（Item已按Hori/Vert摆放过）
                % % %                     updateLWStrip();
                % % %             else % 如果非新level 如可以摆放,放入; 否则,调换长宽（旋转）后放入
                % % %                     if isflagCurr
                % % %                         updateLWStrip();
                % % %                     else
                % % %                         error('11111111');
                % % % %                         rotateItem(); %不能直接注释
                % % % %                         updateLWStrip();
                % % %                     end
                % % %              end            
                % % %         elseif sItem.isRota(iItem) == 0 %此Item不可以旋转
                % % %             % 更新strip信息
                % % %             updateLWStrip();
                % % %         end
        
        %  2.2 更新Strip.Strip_Item 行1 本strip内包含几个Item
        tmpStrip_Item(1,thisLevel) = tmpStrip_Item(1,thisLevel) + 1; %只要该level安放一个item,数量就增加1
        
        %  2.3 更新本level对应的StripWeight: 
        Strip.Weight(thisLevel) =  Strip.Weight(thisLevel) + sItem.Weight(iItem);
        
        %  1.3 更新item归属strip信息itemBeStripMatrixSort
        sItem.Item_Strip(1,iItem) = thisLevel;    %第几个level
        sItem.Item_Strip(2,iItem) = tmpStrip_Item(1,thisLevel); %本level下第几次安置
        
                        %  2.2 更新Strip.Strip_Item 行2 本strip内包含几种Item
                %         itemThisLevel = sItem.Item_Strip(1,:) == thisLevel;
                %         Strip.Strip_Item(2,thisLevel) = numel(unique(sItem.LID(1,itemThisLevel)));

        % 4 二级嵌套函数- 基本不再需要
% %         function rotateItem()
% %             %  1 不仅标记Rotaed变化 还要把ITEM真正的rotate(反)过去
% %             sItem.Rotaed(iItem) = ~sItem.Rotaed(iItem);
% %             tep = sItem.LWH(1,iItem);
% %             sItem.LWH(1,iItem) = sItem.LWH(2,iItem);
% %             sItem.LWH(2,iItem) = tep;
% %             sItem.LWH(3,iItem) = sItem.LWH(3,iItem);
% %             
% %             %  2 Item内部的Lu也要随之更新
% %             tmpflagThisItem = (LU.LU_Item(1,:)==iItem);
% %             % 对应位置LU.Rotaed更新
% %             if Item.Rotaed(iItem)
% %                 LU.Rotaed(tmpflagThisItem) = ~LU.Rotaed(tmpflagThisItem);
% %                 % 对应位置LU.LWH的长宽更新，高度不更新
% %                 LU.LWH(1, tmpflagThisItem) = Item.LWH(1, iItem);
% %                 LU.LWH(2, tmpflagThisItem) = Item.LWH(2, iItem);
% %             end
% %         end
        
        function updateLWStrip()
            Strip.LW(1,thisLevel) = Strip.LW(1,thisLevel) - sItem.LWH(1,iItem); %更新wleft (摆放方向前面一定)
            Strip.LW(2,thisLevel) = max(Strip.LW(2,thisLevel), sItem.LWH(2,iItem)); %更新strip高度lleft(取最大值)
            
            % 更新Strip中包含ID类与否
%             Strip.LID(sItem.LID(1,iItem),thisLevel) =  1;         
%             Strip.SID(sItem.SID(1,iItem),thisLevel) =  1;        
%             Strip.UID(sItem.UID(1,iItem),thisLevel) = 1;         % 数值为出现与否

%              Strip.PID(:,thisLevel) = Strip.PID(:,thisLevel) + sItem.PID(:,iItem); % 数值为出现次数
%              Strip.PID(Strip.PID>0) = 1; % 数值改为出现与否
            
%             Strip.SID(sItem.SID(1,iItem),thisLevel) = 1;         % 555 更新多行PID
%             Strip.UID(sItem.UID(1,iItem),thisLevel) = 1;         % 555 更新多行PID
%             Item.PID(sLU.PID(1,iLU),thisItem) = 1;         % 555 更新多行PID
            
        end
        
    end
    
    function printscript()
        % 测试代码
        % % LWHStrip
        % % sItem.Item_Strip
        % % ItemLWSort
        % % itemCoordMatrixSort
        % % Item_Strip
        % % LWHItem
        % % itemCoordMatrix
        %  printstruct(d);
        
        % 输出主要结果:获得从1开始每个strip包含的数据
        for iStrip = 1:max(Item.Item_Strip(1,:))
            [~,idx] = find(Item.Item_Strip(1,:)==iStrip);
            fprintf('strip %d 的剩余宽+最大长为:  ',iStrip);
            fprintf('( %d ) ',Strip.LW(:,iStrip));
            fprintf('\n');
            fprintf('strip %d 包含 original Item 索引号(长宽)[旋转标志]{坐标}为  \n  ',iStrip);
            fprintf('%d ',idx);
            fprintf('( %d ) ', Item.LWH(1:nDim,idx));fprintf('\n');
            fprintf('[ %d ] ', Item.Rotaed(:,idx));fprintf('\n');  %sItem.Rotaed
            fprintf('{ %d } ', Item.CoordItemStrip(:,idx));fprintf('\n');
            fprintf('\n');
        end
    end
    
    % sItem进行三维作图
    function plot3DStrip()
        nIDType = unique(LU.ID);
        nColors = hsv(length(nIDType)); %不同类型Item的LU赋予不同颜色
        
        yxz = sItem.LWH;
        coord = sItem.CoordItemStrip;   coord(3,:) = 0;
        
        for i=1:numel(sItem.Weight)
            j=sItem.LID(i);
            LUColor = 0.8*nColors(nIDType==j{1}, : );
            plotcube(yxz(:,i)',coord(:,i)',0.7,LUColor);
            
            % Set the lable and the font size
            axis equal;         grid on;        view(60,40);
            xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
            xlim([0 Veh.LWH(1,1)]);         zlim([0 Veh.LWH(3,1)]);
        end
    end

    function plot2DStrip()
        %% 初始化
        wStrip = wStrip;        
        hStrip = sum(Strip.LW(2,sItem.Item_Strip(2,:)>0));        
        nstrip = sum(sItem.Item_Strip(2,:)>0);

        tmpLID = cellfun(@(x)x(1), sItem.LID);        
        nIDType = unique(tmpLID);
        nColors = hsv(length(nIDType)); %不同类型LU赋予不同颜色
        
        %% 画图
        % 1 画图：画本次Strip
        DrawRectangle([wStrip/2 hStrip/2 wStrip hStrip 0],'--', [0.5 0.5 0.5]);
        hold on;
        % 2 画图：逐个strip/item 画图
        for istrip = 1:nstrip
            % 找出当前istrip的物品索引
            idxDrawItem = find(sItem.Item_Strip(1,:)==istrip);
            % 获取该索引下的变量
            drawItemCoordMatrix = sItem.CoordItemStrip(:,idxDrawItem);
            drawItemLWH = sItem.LWH(:,idxDrawItem);
            drawItemId = tmpLID(:,idxDrawItem);

            % 画图：逐个item
            nThisItem = size(drawItemLWH,2);
            for iplotItem = 1:nThisItem
                % 画图：画本次iItem
                itemWidth = drawItemLWH(1,iplotItem);
                itemLength = drawItemLWH(2,iplotItem);
                itemCenter = [drawItemCoordMatrix(1,iplotItem)+itemWidth/2 ...
                    drawItemCoordMatrix(2,iplotItem)+itemLength/2 ];
                
                % 增加对本次iItem的类型（颜色）判断
                itemID = drawItemId(iplotItem);
                itemColor = 0.8*nColors(nIDType==itemID, : );
                
                DrawRectangle([itemCenter itemWidth itemLength 0],  '-',itemColor);
                hold on;
            end
        end
        % hold off;
    end
end

%% 局部函数 %%

%% 函数1: getITEMorder
% 给定ITEM的顺序,按NEXT FIT的方式插入STRIP（先插入SID小的; 后续高度/宽度： 后续LID）
function order = getITEMorder(Item,whichSortItemOrder)

%对SID排序: SID按给定顺序排序,序号小的在前面
szRow = cellfun(@(x)size(x,1), Item.SID);
if (max(szRow)~=min(szRow)),  error('同一ITEM不应该有多个SID');  end %同一Item应该只有一个SID
SIDorder = cell2mat(Item.SID);   %直接cell2mat转换; %ITEM按SID 1-n的顺序返回 

%对LID排序: LID无指定顺序, 仅在SID长宽全部一致,再按LID由小到达排序,其实没有意义(无SID/LID属于同一ITEM),最后看高度 
szRow = cellfun(@(x)size(x,1), Item.LID);
if (max(szRow)~=min(szRow)),  error('同一ITEM不应该有多个SID');  end %同一Item应该只有一个LID
LIDorder = cell2mat(Item.LID);   %直接cell2mat转换; %ITEM按SID 1-n的顺序返回 

% LIDorder = ones(1,length(SIDorder));
% *********** 考虑isNonMixed
tmpItem = [SIDorder; LIDorder; Item.LWH; Item.isNonMixed];  % tmpItem = [Item.SID; Item.LID; Item.LWH];  % tmpItem = [ Item.LWH];
        % [~,order] = sortrows(tmpItem',[1, 4, 2, 5],{'ascend','descend','descend','descend'}); 

% 增加isNonMixed
    % [~,order] = sortrows(tmpItem',[1,6,  4, 3, 2, 5],{'ascend','descend','descend','descend','descend','descend'});  
    % 增加判定Item先长度,后高度排序. % 会出现nbcol2错误:
    % 1: SID ; 2: isNonMixed; 3: Longth/Height; 4: Height 5:Width; 6: LID; 
    % [~,order] = sortrows(tmpItem',[1,6, 4, 5, 3, 2],{'ascend','descend','descend','descend','descend','descend'});  
% 目前顺序 : 1: SID ; 2: isNonMixed; 3: Longth/Height; 4:Width; 5: LID; 6: Height
[~,order] = sortrows(tmpItem',[1,6, 4, 3, 2, 5 ],{'ascend','descend','descend','descend','descend','descend'});  
% [~,order] = sortrows(tmpItem',[1,6, 2, 4, 3, 5 ],{'ascend','descend','descend','descend','descend','descend'});  

% *********** 不考虑isNonMixed
% % tmpItem = [SIDorder; LIDorder; Item.LWH; ];  % tmpItem = [Item.SID; Item.LID; Item.LWH];  % tmpItem = [ Item.LWH];
% % % [~,order] = sortrows(tmpItem',[1, 4, 2, 5],{'ascend','descend','descend','descend'}); 
% % [~,order] = sortrows(tmpItem',[1, 4, 3, 2, 5],{'ascend','descend','descend','descend','descend'});  

%按ITEM长度/随后宽度/随后高度 排序有问题 可能相同IDLU被分开
% 增加LUID 2: 确保即使长宽完全相同 但LUID相同的 也必须放一起
if ~isrow(order), order=order'; end
end

%% 函数2: getThisLevel
function [thisLevel,iLevel] = getThisLevel( iItem, iLevel, sItem,  Strip, p)
% 不同whichStripH下,获得共同的thisLevel
if p.whichStripH == 1 % 1 bestfit 2 firstfit 3 nextfit
    % 增对Rotation增加变量
    if sItem.isRota(iItem) == 1 %此Item可以旋转
        %             if ParaArray.whichRotation == 1
        % 找到可以rotation下的level:任一摆放方向可放入该iItem的level
        flag = find(Strip.LW(1, 1 : iLevel) >= sItem.LWH(1,iItem) |  ...
            Strip.LW(1, 1 : iLevel) >= sItem.LWH(2,iItem));
    else %此Item不可以旋转
        % 常规条件下的选择：find宽度足够的多个level,并安置在最小剩余水平宽度的
        flag = find(Strip.LW(1, 1 : iLevel) >= sItem.LWH(1,iItem));
    end
    if isempty(flag)
        iLevel = iLevel + 1;% 如果宽度不满足，则level升级
        [thisLevel,iLevel] = getThisLevel( iItem, iLevel, sItem,  Strip, p);
    else
        % 获取thisLevel: 唯一与FF区别从这到thisLevel的计算（选中满足条件且最小的）
        tmpLevels = Strip.LW(1,1:iLevel);   %获取所有已安排或新安排的level的剩余水平平宽度向量tmpLevels
        tepMinLeftWdith = min(tmpLevels(flag));                      %找出tepAvailableLevelArray中可容纳本iITem的剩余水平宽度的最小值，更新为tepMinLeftWdith
        thisLevel = find(tmpLevels==tepMinLeftWdith);            %找出与最小值对应的那个/些level
        if ~all(ismember(thisLevel,flag)),     error('Not all thisLevel belongs to flag ');          end
        if length(thisLevel)>1
            thisLevel = thisLevel(1);
        end
    end
elseif p.whichStripH == 2 % firstfit  % firstfit下能不能直接套用bestfit的代码?
    % 增对Rotation增加变量
    if sItem.isRota(iItem) == 1 %此Item可以旋转
        % 找到可以rotation下的level:任一摆放方向可放入该iItem的level
        flag = find(Strip.LW(1, 1 : iLevel) >= sItem.LWH(1,iItem) |  ...
            Strip.LW(1, 1 : iLevel) >= sItem.LWH(2,iItem));
    else
        % 常规条件下的选择：find宽度足够的多个level,并安置在第一个遇到的 唯一区别是thisLevel的获取
        flag = find(Strip.LW(1, 1 : iLevel) >= sItem.LWH(1,iItem));
    end
    if isempty(flag)
        iLevel = iLevel + 1;% 如果宽度不满足，则level升级
        [thisLevel,iLevel] = getThisLevel( iItem, iLevel, sItem, Strip, p);
    else
        thisLevel = flag(1);
        if ~all(ismember(thisLevel,flag)),     error('Not all thisLevel belongs to flag ');       end
    end
elseif p.whichStripH == 3 % nextfit
    % 增对Rotation增加变量
    if sItem.isRota(iItem) == 1 %此Item可以旋转 % nextfit下不能直接套用bestfit的代码
        % 判定当前level是否可以在任一摆放方向可放入该iItem flaged: 有内容表示可以，否则不可以
        flaged = find(Strip.LW(1,iLevel) >= sItem.LWH(1,iItem) );
        %                 flaged = find(Strip.LW(1,iLevel) >= sItem.LWH(1,iItem) |  ...
        %                                       Strip.LW(1,iLevel) >= sItem.LWH(2,iItem));
    else
        % 不同条件下的选择：如果当前item的宽<=当前strip的当前level的宽
        flaged = find(Strip.LW(1,iLevel) >= sItem.LWH(1,iItem) );
    end
    if  isempty(flaged)  %注意与之前~flag的区别
        iLevel = iLevel + 1;% 如果宽度不满足，则level升级
        [thisLevel,iLevel] = getThisLevel( iItem, iLevel, sItem, Strip, p) ;
    else
        if  isempty(flaged) ,   error(' 不可能的错误 ');      end
        thisLevel = iLevel; % 当前level一定放的下
    end
    %             end
    
end
end
