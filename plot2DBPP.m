function [] = plot2DBPP(da,ParaArray)
% ��ʼ��
nDim = size(da.ItemArray.LWH,1);  if nDim ==3, nDim = nDim-1;end
nThisItem = size(da.ItemArray.LWH,2);
tmpUniqueBin = unique(da.BinArray.LWH(1:nDim,:)','rows')';
widthBin = tmpUniqueBin(1);
lengthBin = tmpUniqueBin(2);
clear tmpUniqueBin;

%% �ṹ����ȡ
itemBeBinMatrix = da.ItemArray.itemBeBinMatrix;
CoordItemBin = da.ItemArray.CoordItemBin;
itemLWH = da.ItemArray.LWH;
if ParaArray.whichRotation == 1    
    itemRotaFlag = da.ItemArray.itemRotaFlag; %����rotation����
end

% sort ����bin��˳�� ���bin��ͼ
[~,binorder] = sort(itemBeBinMatrix(1,:),'ascend');
% ��ȡ�����ĸ�������
itemBeBinMatrixSort = itemBeBinMatrix(:,binorder);
CoordItemBinSort = CoordItemBin(:,binorder);
itemLWHSort = itemLWH(:,binorder);
if ParaArray.whichRotation == 1    
    itemRotaFlagSort = itemRotaFlag(:,binorder); %����rotation����
end

%% 1 �������� ���ΪnBin+1��bin�� �����ߣ���Ϊbin��
nBin = max(itemBeBinMatrix(1,:));
DrawRectangle([widthBin*(nBin+1)/2 lengthBin/2 widthBin*(nBin+1) lengthBin 0],'--');
hold on;
% 2 ���bin ��ͼ
iterWidth=0;    %ÿ��bin��ǰ1��bin���Ҳ� ��Ϊ���ӱ���
for iBin = 1:nBin
    % �ҳ���ǰiBin����Ʒ����
    idxDrawItem = find(itemBeBinMatrixSort(1,:)==iBin);
               %     drawBinMatrix = ppbelongItemBinMatrix(:,idxDrawItem);
    % ��ȡ�������µı���
    drawItemCoordMatrix = CoordItemBinSort(:,idxDrawItem);
    drawItemMatrix = itemLWHSort(:,idxDrawItem);
    if ParaArray.whichRotation == 1    
        drawItemRotaMatrix = itemRotaFlagSort(:,idxDrawItem); %����rotation����
    end
    % ��ͼ��������iBin
    binCenter = [iterWidth+widthBin/2 lengthBin/2];
    DrawRectangle([binCenter widthBin lengthBin 0],'--')    ;
    hold on;
    % ��ͼ�����item
    nThisItem = size(drawItemMatrix,2);
    for iItem = 1:nThisItem
        % ��ͼ��������iItem
        itemWidth = drawItemMatrix(1,iItem);
        itemLength = drawItemMatrix(2,iItem);
        itemCenter = [iterWidth+drawItemCoordMatrix(1,iItem)+itemWidth/2 ...
            drawItemCoordMatrix(2,iItem)+itemLength/2 ];
        % ���Ӷ�rotation���ж�
        if ParaArray.whichRotation == 1 && drawItemRotaMatrix(iItem)
            itemWidth = drawItemMatrix(2,iItem);
            itemLength = drawItemMatrix(1,iItem);            
            itemCenter = [iterWidth+drawItemCoordMatrix(1,iItem)+itemWidth/2 ...
                        drawItemCoordMatrix(2,iItem)+itemLength/2 ];
        end
        DrawRectangle([itemCenter itemWidth itemLength 0],'r-');
        hold on;
    end
    % ��������bin�Ŀ��
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
% %     pitemRotaMatrix = ItemArray.itemRotaSortHori; %����rotation����
% % end
% % % sort ����bin��˳�� ���bin��ͼ
% % [~,binord] = sort(pbelongItemBinMatrix(1,:),'ascend');
% % % ��ȡ�����ĸ�������
% % ppbelongItemBinMatrix = pbelongItemBinMatrix(:,binord);
% % ppcoordItemBinMatrix = pcoordItemBinMatrix(:,binord);
% % ppitemMatrix = pitemMatrix(:,binord);
% % if ParaArray.whichRotation == 1    
% %     ppitemRotaMatrix = pitemRotaMatrix(:,binord); %����rotation����
% % end
% % % �������� ���ΪnBin+1��bin�� �����ߣ���Ϊbin��
% % nBin = max(pbelongItemBinMatrix(1,:));
% % DrawRectangle([widthBin*(nBin+1)/2 lengthBin/2 widthBin*(nBin+1) lengthBin 0],'--');
% % hold on;
% % % ��ͼ�����bin
% % iterWidth=0;    %ÿ��bin��ǰ1��bin���Ҳ� ��Ϊ���ӱ���
% % for iBin = 1:nBin
% %     % �ҳ���ǰiBin����Ʒ����
% %     idxDrawItem = find(ppbelongItemBinMatrix(1,:)==iBin);
% % %     drawBinMatrix = ppbelongItemBinMatrix(:,idxDrawItem);
% %     % ��ȡ�������µı���
% %     drawItemCoordMatrix = ppcoordItemBinMatrix(:,idxDrawItem);
% %     drawItemMatrix = ppitemMatrix(:,idxDrawItem);
% %     if ParaArray.whichRotation == 1    
% %         drawItemRotaMatrix = ppitemRotaMatrix(:,idxDrawItem); %����rotation����
% %     end
% %     % ��ͼ��������iBin
% % %     binWidth = BinArray.LWH(1,1);
% % %     binLength = BinArray.LWH(2,1);
% %     binCenter = [iterWidth+widthBin/2 lengthBin/2];
% %     DrawRectangle([binCenter widthBin lengthBin 0],'--')    ;
% %     hold on;
% %     % ��ͼ�����item
% %     nItem = size(drawItemMatrix,2);
% %     for iItem = 1:nItem
% %         % ��ͼ��������iItem
% %         itemWidth = drawItemMatrix(1,iItem);
% %         itemLength = drawItemMatrix(2,iItem);
% %         itemCenter = [iterWidth+drawItemCoordMatrix(1,iItem)+itemWidth/2 ...
% %             drawItemCoordMatrix(2,iItem)+itemLength/2 ];
% %         % ���Ӷ�rotation���ж�
% %         if ParaArray.whichRotation == 1 && drawItemRotaMatrix(iItem)
% %             itemWidth = drawItemMatrix(2,iItem);
% %             itemLength = drawItemMatrix(1,iItem);            
% %             itemCenter = [iterWidth+drawItemCoordMatrix(1,iItem)+itemWidth/2 ...
% %                         drawItemCoordMatrix(2,iItem)+itemLength/2 ];
% %         end
% %         DrawRectangle([itemCenter itemWidth itemLength 0],'r-');
% %         hold on;
% %     end
% %     % ��������bin�Ŀ��
% %     iterWidth = iterWidth + widthBin;
% %     hold on;
% % end
% % hold off;
% % end


%% �ǽṹ�廭ͼ
% % function [] = plot2DBPP(pbelongItemBinMatrix,pcoordItemBinMatrix,pitemMatrix,BinArray)
% % % sort ����bin��˳�� ���bin��ͼ
% % [~,binord] = sort(pbelongItemBinMatrix(1,:),'ascend');
% % % ��ȡ�����ĸ�������
% % ppbelongItemBinMatrix = pbelongItemBinMatrix(:,binord);
% % ppcoordItemBinMatrix = pcoordItemBinMatrix(:,binord);
% % ppitemMatrix = pitemMatrix(:,binord);
% % % �������� ���ΪnBin+1��bin�� �����ߣ���Ϊbin��
% % nBin = max(pbelongItemBinMatrix(1,:));
% % DrawRectangle([BinArray.LWH(1,1)*(nBin+1)/2 BinArray.LWH(2,1)/2 BinArray.LWH(1,1)*(nBin+1) BinArray.LWH(2,1) 0],'--');
% % hold on;
% % % ��ͼ�����bin
% % iterWidth=0;    %ÿ��bin��ǰ1��bin���Ҳ� ��Ϊ���ӱ���
% % for iBin = 1:nBin
% %     % �ҳ���ǰiBin����Ʒ����
% %     idxDrawItem = find(ppbelongItemBinMatrix(1,:)==iBin);
% % %     drawBinMatrix = ppbelongItemBinMatrix(:,idxDrawItem);
% %     % ��ȡ�������µı���
% %     drawItemCoordMatrix = ppcoordItemBinMatrix(:,idxDrawItem);
% %     drawItemMatrix = ppitemMatrix(:,idxDrawItem);
% %     % ��ͼ��������iBin
% %     binWidth = BinArray.LWH(1,1);
% %     binLength = BinArray.LWH(2,1);
% %     binCenter = [iterWidth+binWidth/2 binLength/2];
% %     DrawRectangle([binCenter binWidth binLength 0],'--')    ;
% %     hold on;
% %     % ��ͼ�����item
% %     nItem = size(drawItemMatrix,2);
% %     for iItem = 1:nItem
% %         % ��ͼ��������iItem
% %         itemWidth = drawItemMatrix(1,iItem);
% %         itemLength = drawItemMatrix(2,iItem);
% %         itemCenter = [iterWidth+drawItemCoordMatrix(1,iItem)+itemWidth/2 drawItemCoordMatrix(2,iItem)+itemLength/2 ];
% %         DrawRectangle([itemCenter itemWidth itemLength 0],'r-');
% %         hold on;
% %     end
% %     % ��������bin�Ŀ��
% %     iterWidth = iterWidth + binWidth;
% %     hold on;
% % end
% % hold off;
% % end