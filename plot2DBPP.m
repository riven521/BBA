function [] = plot2DBPP(d,ParaArray)

% ��ͼ����:��άBPP
% ��ʼ��
nDim = size(d.Item.LWH,1);  if nDim ==3, nDim = nDim-1;end
nThisItem = size(d.Item.LWH,2);
nIDType = unique(d.Item.ID);
nColors = hsv(length(nIDType)); %��ͬ����LU���費ͬ��ɫ
% tmpUniqueBin = unique(d.Veh.LWH(1:nDim,:)','rows')';
widthBin = d.Veh.LWH(1,1);
lengthBin = d.Veh.LWH(2,1);
clear tmpUniqueBin;

%% �ṹ����ȡ
Item_Bin = d.Item.Item_Bin;
CoordItemBin = d.Item.CoordItemBin;
itemLWH = d.Item.LWH;
itemID = d.Item.ID;
itemSID = d.Item.SID;
% if ParaArray.whichRotation == 1    
    ItemRotaed = d.Item.Rotaed; %����rotation���� itemRotaFlag
% end

% sort ����bin��˳�� ���bin��ͼ
[~,binorder] = sort(Item_Bin(1,:),'ascend');
% ��ȡ�����ĸ�������
itemBeBinMatrixSort = Item_Bin(:,binorder);
CoordItemBinSort = CoordItemBin(:,binorder);
itemLWHSort = itemLWH(:,binorder);
itemIDSort = itemID(:,binorder);
itemSIDSort = itemSID(:,binorder);
% if ParaArray.whichRotation == 1    
    ItemRotaedSort = ItemRotaed(:,binorder); %����rotation����
% end

% 1 �������� ���ΪnBin+1��bin�� �����ߣ���Ϊbin��
nBin = max(Item_Bin(1,:));
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
    drawItemLWH = itemLWHSort(:,idxDrawItem);
    drawItemId = itemIDSort(:,idxDrawItem);
%     if ParaArray.whichRotation == 1    
% %         drawItemRotaMatrix = ItemRotaedSort(:,idxDrawItem); %����rotation����
%     end
    % ��ͼ��������iBin
    binCenter = [iterWidth+widthBin/2 lengthBin/2];
    DrawRectangle([binCenter widthBin lengthBin 0],'--')    ;
    hold on;
    % ��ͼ�����item
    nThisItem = size(drawItemLWH,2);
    for iItem = 1:nThisItem
        % ��ͼ��������iItem
        itemWidth = drawItemLWH(1,iItem);
        itemLength = drawItemLWH(2,iItem);
        itemCenter = [iterWidth+drawItemCoordMatrix(1,iItem)+itemWidth/2 ...
            drawItemCoordMatrix(2,iItem)+itemLength/2 ];

        % ���Ӷ�rotation���ж�
%         if ParaArray.whichRotation == 1 && drawItemRotaMatrix(iItem)
%             itemWidth = drawItemLWH(2,iItem);
%             itemLength = drawItemLWH(1,iItem);            
%             itemCenter = [iterWidth+drawItemCoordMatrix(1,iItem)+itemWidth/2 ...
%                         drawItemCoordMatrix(2,iItem)+itemLength/2 ];
%         end
        
        % ���ӶԱ���iItem�����ͣ���ɫ���ж�
        itemID = drawItemId(iItem);
        itemColor = 0.8*nColors(nIDType==itemID, : );        
        
%         for iSid = 1:length(unique(itemSIDSort))
%             itemSIDSort(iSid)
        DrawRectangle([itemCenter itemWidth itemLength 0],  '-',  itemColor); 
%         itemSIDSort
        hold on;
%         end
    end
    % ��������bin�Ŀ��
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
% %     pitemRotaMatrix = Item.itemRotaSortHori; %����rotation����
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
% % %     binWidth = Veh.LWH(1,1);
% % %     binLength = Veh.LWH(2,1);
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
% % function [] = plot2DBPP(pbelongItemBinMatrix,pcoordItemBinMatrix,pitemMatrix,Veh)
% % % sort ����bin��˳�� ���bin��ͼ
% % [~,binord] = sort(pbelongItemBinMatrix(1,:),'ascend');
% % % ��ȡ�����ĸ�������
% % ppbelongItemBinMatrix = pbelongItemBinMatrix(:,binord);
% % ppcoordItemBinMatrix = pcoordItemBinMatrix(:,binord);
% % ppitemMatrix = pitemMatrix(:,binord);
% % % �������� ���ΪnBin+1��bin�� �����ߣ���Ϊbin��
% % nBin = max(pbelongItemBinMatrix(1,:));
% % DrawRectangle([Veh.LWH(1,1)*(nBin+1)/2 Veh.LWH(2,1)/2 Veh.LWH(1,1)*(nBin+1) Veh.LWH(2,1) 0],'--');
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
% %     binWidth = Veh.LWH(1,1);
% %     binLength = Veh.LWH(2,1);
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