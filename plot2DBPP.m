function [] = plot2DBPP(d,ParaArray)

% 作图函数:二维BPP
% 初始化
nDim = size(d.Item.LWH,1);  if nDim ==3, nDim = nDim-1;end
nThisItem = size(d.Item.LWH,2);
nIDType = unique(d.Item.ID);
nColors = hsv(length(nIDType)); %不同类型LU赋予不同颜色
% tmpUniqueBin = unique(d.Veh.LWH(1:nDim,:)','rows')';
widthBin = d.Veh.LWH(1,1);
lengthBin = d.Veh.LWH(2,1);
clear tmpUniqueBin;

%% 结构体提取
Item_Bin = d.Item.Item_Bin;
CoordItemBin = d.Item.CoordItemBin;
itemLWH = d.Item.LWH;
itemID = d.Item.ID;
itemSID = d.Item.SID;
% if ParaArray.whichRotation == 1    
    ItemRotaed = d.Item.Rotaed; %增加rotation后增 itemRotaFlag
% end

% sort 排序按bin的顺序 逐个bin画图
[~,binorder] = sort(Item_Bin(1,:),'ascend');
% 获取排序后的各个变量
itemBeBinMatrixSort = Item_Bin(:,binorder);
CoordItemBinSort = CoordItemBin(:,binorder);
itemLWHSort = itemLWH(:,binorder);
itemIDSort = itemID(:,binorder);
itemSIDSort = itemSID(:,binorder);
% if ParaArray.whichRotation == 1    
    ItemRotaedSort = ItemRotaed(:,binorder); %增加rotation后增
% end

% 1 画个画布 宽度为nBin+1个bin宽 长（高）度为bin高
nBin = max(Item_Bin(1,:));
DrawRectangle([widthBin*(nBin+1)/2 lengthBin/2 widthBin*(nBin+1) lengthBin 0],'--');
hold on;
% 2 逐个bin 画图
iterWidth=0;    %每个bin在前1个bin的右侧 此为增加变量
for iBin = 1:nBin
    % 找出当前iBin的物品索引
    idxDrawItem = find(itemBeBinMatrixSort(1,:)==iBin);
               %     drawBinMatrix = ppbelongItemBinMatrix(:,idxDrawItem);
    % 获取该索引下的变量
    drawItemCoordMatrix = CoordItemBinSort(:,idxDrawItem);
    drawItemLWH = itemLWHSort(:,idxDrawItem);
    drawItemId = itemIDSort(:,idxDrawItem);
%     if ParaArray.whichRotation == 1    
% %         drawItemRotaMatrix = ItemRotaedSort(:,idxDrawItem); %增加rotation后增
%     end
    % 画图：画本次iBin
    binCenter = [iterWidth+widthBin/2 lengthBin/2];
    DrawRectangle([binCenter widthBin lengthBin 0],'--')    ;
    hold on;
    % 画图：逐个item
    nThisItem = size(drawItemLWH,2);
    for iItem = 1:nThisItem
        % 画图：画本次iItem
        itemWidth = drawItemLWH(1,iItem);
        itemLength = drawItemLWH(2,iItem);
        itemCenter = [iterWidth+drawItemCoordMatrix(1,iItem)+itemWidth/2 ...
            drawItemCoordMatrix(2,iItem)+itemLength/2 ];

        % 增加对rotation的判断
%         if ParaArray.whichRotation == 1 && drawItemRotaMatrix(iItem)
%             itemWidth = drawItemLWH(2,iItem);
%             itemLength = drawItemLWH(1,iItem);            
%             itemCenter = [iterWidth+drawItemCoordMatrix(1,iItem)+itemWidth/2 ...
%                         drawItemCoordMatrix(2,iItem)+itemLength/2 ];
%         end
        
        % 增加对本次iItem的类型（颜色）判断
        itemID = drawItemId(iItem);
        itemColor = 0.8*nColors(nIDType==itemID, : );        
        
%         for iSid = 1:length(unique(itemSIDSort))
%             itemSIDSort(iSid)
        DrawRectangle([itemCenter itemWidth itemLength 0],  '-',  itemColor); 
%         itemSIDSort
        hold on;
%         end
    end
    % 递增本次bin的宽度
    iterWidth = iterWidth + widthBin;
    hold on;
end
hold off;
end

%%
% % function [] = plot2DBPP(d,ParaArray)
% % Veh=d.Bin;
% % Item=d.Item;
% % 
% % nDim = size(Item.LWH,1);  if nDim ==3, nDim = nDim-1;end
% % nItem = size(Item.LWH,2);
% % nStrip = nItem;
% % 
% % % itemDataMatrix = d.Item.LWH(1:nDim,:);
% % tmpbinDataMatrix = d.Veh.LWH(1:nDim,:); tmpUniqueBin = unique(tmpbinDataMatrix','rows')';
% % widthBin = tmpUniqueBin(1);
% % lengthBin = tmpUniqueBin(2);
% % 
% % %%
% % % printstruct(d);
% % pbelongItemBinMatrix = Item.itemBeBinMatrixSort;
% % pcoordItemBinMatrix = Item.itemCoordMatrixSort;
% % pitemMatrix = Item.LWHSort;
% % if ParaArray.whichRotation == 1    
% %     pitemRotaMatrix = Item.itemRotaSortHori; %增加rotation后增
% % end
% % % sort 排序按bin的顺序 逐个bin画图
% % [~,binord] = sort(pbelongItemBinMatrix(1,:),'ascend');
% % % 获取排序后的各个变量
% % ppbelongItemBinMatrix = pbelongItemBinMatrix(:,binord);
% % ppcoordItemBinMatrix = pcoordItemBinMatrix(:,binord);
% % ppitemMatrix = pitemMatrix(:,binord);
% % if ParaArray.whichRotation == 1    
% %     ppitemRotaMatrix = pitemRotaMatrix(:,binord); %增加rotation后增
% % end
% % % 画个画布 宽度为nBin+1个bin宽 长（高）度为bin高
% % nBin = max(pbelongItemBinMatrix(1,:));
% % DrawRectangle([widthBin*(nBin+1)/2 lengthBin/2 widthBin*(nBin+1) lengthBin 0],'--');
% % hold on;
% % % 画图：逐个bin
% % iterWidth=0;    %每个bin在前1个bin的右侧 此为增加变量
% % for iBin = 1:nBin
% %     % 找出当前iBin的物品索引
% %     idxDrawItem = find(ppbelongItemBinMatrix(1,:)==iBin);
% % %     drawBinMatrix = ppbelongItemBinMatrix(:,idxDrawItem);
% %     % 获取该索引下的变量
% %     drawItemCoordMatrix = ppcoordItemBinMatrix(:,idxDrawItem);
% %     drawItemMatrix = ppitemMatrix(:,idxDrawItem);
% %     if ParaArray.whichRotation == 1    
% %         drawItemRotaMatrix = ppitemRotaMatrix(:,idxDrawItem); %增加rotation后增
% %     end
% %     % 画图：画本次iBin
% % %     binWidth = Veh.LWH(1,1);
% % %     binLength = Veh.LWH(2,1);
% %     binCenter = [iterWidth+widthBin/2 lengthBin/2];
% %     DrawRectangle([binCenter widthBin lengthBin 0],'--')    ;
% %     hold on;
% %     % 画图：逐个item
% %     nItem = size(drawItemMatrix,2);
% %     for iItem = 1:nItem
% %         % 画图：画本次iItem
% %         itemWidth = drawItemMatrix(1,iItem);
% %         itemLength = drawItemMatrix(2,iItem);
% %         itemCenter = [iterWidth+drawItemCoordMatrix(1,iItem)+itemWidth/2 ...
% %             drawItemCoordMatrix(2,iItem)+itemLength/2 ];
% %         % 增加对rotation的判断
% %         if ParaArray.whichRotation == 1 && drawItemRotaMatrix(iItem)
% %             itemWidth = drawItemMatrix(2,iItem);
% %             itemLength = drawItemMatrix(1,iItem);            
% %             itemCenter = [iterWidth+drawItemCoordMatrix(1,iItem)+itemWidth/2 ...
% %                         drawItemCoordMatrix(2,iItem)+itemLength/2 ];
% %         end
% %         DrawRectangle([itemCenter itemWidth itemLength 0],'r-');
% %         hold on;
% %     end
% %     % 递增本次bin的宽度
% %     iterWidth = iterWidth + widthBin;
% %     hold on;
% % end
% % hold off;
% % end


%% 非结构体画图
% % function [] = plot2DBPP(pbelongItemBinMatrix,pcoordItemBinMatrix,pitemMatrix,Veh)
% % % sort 排序按bin的顺序 逐个bin画图
% % [~,binord] = sort(pbelongItemBinMatrix(1,:),'ascend');
% % % 获取排序后的各个变量
% % ppbelongItemBinMatrix = pbelongItemBinMatrix(:,binord);
% % ppcoordItemBinMatrix = pcoordItemBinMatrix(:,binord);
% % ppitemMatrix = pitemMatrix(:,binord);
% % % 画个画布 宽度为nBin+1个bin宽 长（高）度为bin高
% % nBin = max(pbelongItemBinMatrix(1,:));
% % DrawRectangle([Veh.LWH(1,1)*(nBin+1)/2 Veh.LWH(2,1)/2 Veh.LWH(1,1)*(nBin+1) Veh.LWH(2,1) 0],'--');
% % hold on;
% % % 画图：逐个bin
% % iterWidth=0;    %每个bin在前1个bin的右侧 此为增加变量
% % for iBin = 1:nBin
% %     % 找出当前iBin的物品索引
% %     idxDrawItem = find(ppbelongItemBinMatrix(1,:)==iBin);
% % %     drawBinMatrix = ppbelongItemBinMatrix(:,idxDrawItem);
% %     % 获取该索引下的变量
% %     drawItemCoordMatrix = ppcoordItemBinMatrix(:,idxDrawItem);
% %     drawItemMatrix = ppitemMatrix(:,idxDrawItem);
% %     % 画图：画本次iBin
% %     binWidth = Veh.LWH(1,1);
% %     binLength = Veh.LWH(2,1);
% %     binCenter = [iterWidth+binWidth/2 binLength/2];
% %     DrawRectangle([binCenter binWidth binLength 0],'--')    ;
% %     hold on;
% %     % 画图：逐个item
% %     nItem = size(drawItemMatrix,2);
% %     for iItem = 1:nItem
% %         % 画图：画本次iItem
% %         itemWidth = drawItemMatrix(1,iItem);
% %         itemLength = drawItemMatrix(2,iItem);
% %         itemCenter = [iterWidth+drawItemCoordMatrix(1,iItem)+itemWidth/2 drawItemCoordMatrix(2,iItem)+itemLength/2 ];
% %         DrawRectangle([itemCenter itemWidth itemLength 0],'r-');
% %         hold on;
% %     end
% %     % 递增本次bin的宽度
% %     iterWidth = iterWidth + binWidth;
% %     hold on;
% % end
% % hold off;
% % end