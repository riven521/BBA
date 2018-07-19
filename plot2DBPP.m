function [] = plot2DBPP(da,ParaArray)
% 初始化
nDim = size(da.ItemArray.LWH,1);  if nDim ==3, nDim = nDim-1;end
nThisItem = size(da.ItemArray.LWH,2);
tmpUniqueBin = unique(da.BinArray.LWH(1:nDim,:)','rows')';
widthBin = tmpUniqueBin(1);
lengthBin = tmpUniqueBin(2);
clear tmpUniqueBin;

%% 结构体提取
itemBeBinMatrix = da.ItemArray.itemBeBinMatrix;
CoordItemBin = da.ItemArray.CoordItemBin;
itemLWH = da.ItemArray.LWH;
if ParaArray.whichRotation == 1    
    itemRotaFlag = da.ItemArray.itemRotaFlag; %增加rotation后增
end

% sort 排序按bin的顺序 逐个bin画图
[~,binorder] = sort(itemBeBinMatrix(1,:),'ascend');
% 获取排序后的各个变量
itemBeBinMatrixSort = itemBeBinMatrix(:,binorder);
CoordItemBinSort = CoordItemBin(:,binorder);
itemLWHSort = itemLWH(:,binorder);
if ParaArray.whichRotation == 1    
    itemRotaFlagSort = itemRotaFlag(:,binorder); %增加rotation后增
end

%% 1 画个画布 宽度为nBin+1个bin宽 长（高）度为bin高
nBin = max(itemBeBinMatrix(1,:));
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
    drawItemMatrix = itemLWHSort(:,idxDrawItem);
    if ParaArray.whichRotation == 1    
        drawItemRotaMatrix = itemRotaFlagSort(:,idxDrawItem); %增加rotation后增
    end
    % 画图：画本次iBin
    binCenter = [iterWidth+widthBin/2 lengthBin/2];
    DrawRectangle([binCenter widthBin lengthBin 0],'--')    ;
    hold on;
    % 画图：逐个item
    nThisItem = size(drawItemMatrix,2);
    for iItem = 1:nThisItem
        % 画图：画本次iItem
        itemWidth = drawItemMatrix(1,iItem);
        itemLength = drawItemMatrix(2,iItem);
        itemCenter = [iterWidth+drawItemCoordMatrix(1,iItem)+itemWidth/2 ...
            drawItemCoordMatrix(2,iItem)+itemLength/2 ];
        % 增加对rotation的判断
        if ParaArray.whichRotation == 1 && drawItemRotaMatrix(iItem)
            itemWidth = drawItemMatrix(2,iItem);
            itemLength = drawItemMatrix(1,iItem);            
            itemCenter = [iterWidth+drawItemCoordMatrix(1,iItem)+itemWidth/2 ...
                        drawItemCoordMatrix(2,iItem)+itemLength/2 ];
        end
        DrawRectangle([itemCenter itemWidth itemLength 0],'r-');
        hold on;
    end
    % 递增本次bin的宽度
    iterWidth = iterWidth + widthBin;
    hold on;
end
hold off;
end

%%
% % function [] = plot2DBPP(da,ParaArray)
% % BinArray=da.BinSArray;
% % ItemArray=da.ItemArray;
% % 
% % nDim = size(ItemArray.LWH,1);  if nDim ==3, nDim = nDim-1;end
% % nItem = size(ItemArray.LWH,2);
% % nStrip = nItem;
% % 
% % % itemDataMatrix = da.ItemArray.LWH(1:nDim,:);
% % tmpbinDataMatrix = da.BinArray.LWH(1:nDim,:); tmpUniqueBin = unique(tmpbinDataMatrix','rows')';
% % widthBin = tmpUniqueBin(1);
% % lengthBin = tmpUniqueBin(2);
% % 
% % %%
% % % printstruct(da);
% % pbelongItemBinMatrix = ItemArray.itemBeBinMatrixSort;
% % pcoordItemBinMatrix = ItemArray.itemCoordMatrixSort;
% % pitemMatrix = ItemArray.LWHSort;
% % if ParaArray.whichRotation == 1    
% %     pitemRotaMatrix = ItemArray.itemRotaSortHori; %增加rotation后增
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
% % %     binWidth = BinArray.LWH(1,1);
% % %     binLength = BinArray.LWH(2,1);
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
% % function [] = plot2DBPP(pbelongItemBinMatrix,pcoordItemBinMatrix,pitemMatrix,BinArray)
% % % sort 排序按bin的顺序 逐个bin画图
% % [~,binord] = sort(pbelongItemBinMatrix(1,:),'ascend');
% % % 获取排序后的各个变量
% % ppbelongItemBinMatrix = pbelongItemBinMatrix(:,binord);
% % ppcoordItemBinMatrix = pcoordItemBinMatrix(:,binord);
% % ppitemMatrix = pitemMatrix(:,binord);
% % % 画个画布 宽度为nBin+1个bin宽 长（高）度为bin高
% % nBin = max(pbelongItemBinMatrix(1,:));
% % DrawRectangle([BinArray.LWH(1,1)*(nBin+1)/2 BinArray.LWH(2,1)/2 BinArray.LWH(1,1)*(nBin+1) BinArray.LWH(2,1) 0],'--');
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
% %     binWidth = BinArray.LWH(1,1);
% %     binLength = BinArray.LWH(2,1);
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