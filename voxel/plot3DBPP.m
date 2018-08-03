function [] = plot3DBPP(d,ParaArray)
% ��ͼ����:��άBPP
% ��ʼ��
nDim = size(d.ItemArray.LWH,1);
nThisItem = size(d.ItemArray.LWH,2);
tmpUniqueBin = unique(d.BinArray.LWH(1:nDim,:)','rows')';
widthBin = tmpUniqueBin(1);
lengthBin = tmpUniqueBin(2);
heightBin = tmpUniqueBin(3);
clear tmpUniqueBin;

%% Item�ṹ����ȡ
Item_Bin = d.ItemArray.Item_Bin;
CoordItemBin = d.ItemArray.CoordItemBin;
itemLWH = d.ItemArray.LWH;
if ParaArray.whichRotation == 1    
    itemRotaFlag = d.ItemArray.itemRotaFlag; %����rotation����
end

% sort ����bin��˳�� ���bin��ͼ
[~,binorder] = sort(Item_Bin(1,:),'ascend');
% ��ȡ�����ĸ�������
itemBeBinMatrixSort = Item_Bin(:,binorder);
CoordItemBinSort = CoordItemBin(:,binorder);
itemLWHSort = itemLWH(:,binorder);
if ParaArray.whichRotation == 1    
    itemRotaFlagSort = itemRotaFlag(:,binorder); %����rotation����
end

%% LU�ṹ����ȡ
LU_Item = d.LUArray.LU_Item;
CoordLUBin = d.LUArray.CoordLUBin;
LULWH = d.LUArray.LWH;
if ParaArray.whichRotation == 1    
    LURotaFlag = d.LUArray.LURotaFlag; %����rotation����
end

% % sort ����bin��˳�� ���bin��ͼ
% [~,itemorder] = sort(LU_Item(1,:),'ascend');
% % ��ȡ�����ĸ�������
% LUBeItemArraySort = LU_Item(:,itemorder);
% CoordLUBinSort = CoordLUBin(:,itemorder);
% LULWHSort = LULWH(:,itemorder);
% if ParaArray.whichRotation == 1    
%     LURotaFlagSort = LURotaFlag(:,itemorder); %����rotation����
% end


%% 1 �������� ���ΪnBin+1��bin�� �����ߣ���Ϊbin��
nThisBin = max(Item_Bin(1,:));
axis([0 lengthBin 0 widthBin 0 heightBin]);
hold on;
% 2 ���bin ��ͼ
% iterWidth=0;    %ÿ��bin��ǰ1��bin���Ҳ� ��Ϊ���ӱ���
for iBin = 1:nThisBin
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
%     binCenter = [iterWidth+widthBin/2 lengthBin/2];
%     DrawRectangle([binCenter widthBin lengthBin 0],'--');
%     hold on;
    % ��ͼ�����item
    nThisItem = size(drawItemMatrix,2);
    for iItem = 1:nThisItem
        % �ҳ���ǰiItme����������
        idxDrawLU = find(LU_Item(1,:)==iItem);
         % ��ȡ�������µı���
        drawLUCoordMatrix = CoordLUBin(:,idxDrawLU); %CoordLUBinSort no use
        drawLUMatrix = LULWH(:,idxDrawLU); %LULWHSort
        if ParaArray.whichRotation == 1
            drawLURotaMatrix = LURotaFlag(:,idxDrawLU); %����rotation���� LURotaFlagSort
        end    
        % ���ݣ���ȡ����item�ĳ����
        itemWidth = drawItemMatrix(1,iItem);
        itemLength = drawItemMatrix(2,iItem);
        itemHeight = 300;
%         itemCenter = [iterWidth+drawItemCoordMatrix(1,iItem)+itemWidth/2 ...
%             drawItemCoordMatrix(2,iItem)+itemLength/2 ];
        % ���Ӷ�rotation���ж�
        if ParaArray.whichRotation == 1 && drawItemRotaMatrix(iItem)
            itemWidth = drawItemMatrix(2,iItem);
            itemLength = drawItemMatrix(1,iItem);      
            itemHeight = 300;
        end        
        
        nThisLU = size(drawLUMatrix,2);
        for iLU = 1:nThisLU
            % ��ͼ��������LU
            LUWidth = itemWidth;
            LULength = itemLength;
            LUHeight = drawLUMatrix(3,iLU);    %to modify
            voxel([drawLUCoordMatrix(2,iLU) drawLUCoordMatrix(1,iLU) drawLUCoordMatrix(3,iLU)],...
                [LUWidth LULength LUHeight],'b',0.95);
            hold on;
        end
    end
    % ��������bin�Ŀ��
%     iterWidth = iterWidth + widthBin;
%     hold on;
end
hold off;
end

%%
% % function [] = plot2DBPP(d,ParaArray)
% % BinArray=d.BinSArray;
% % ItemArray=d.ItemArray;
% % 
% % nDim = size(ItemArray.LWH,1);  if nDim ==3, nDim = nDim-1;end
% % nItem = size(ItemArray.LWH,2);
% % nStrip = nItem;
% % 
% % % itemDataMatrix = d.ItemArray.LWH(1:nDim,:);
% % tmpbinDataMatrix = d.BinArray.LWH(1:nDim,:); tmpUniqueBin = unique(tmpbinDataMatrix','rows')';
% % widthBin = tmpUniqueBin(1);
% % lengthBin = tmpUniqueBin(2);
% % 
% % %%
% % % printstruct(d);
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