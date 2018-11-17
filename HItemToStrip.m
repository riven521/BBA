% function [d] = HItemToStrip(d,p)
function [LU,Item,Strip]= HItemToStrip(LU,Item,Veh,p)
% ��Ҫ����:Item����Strip�� %  ����:�����(row);  ����:��������(coloum);
% Input ---  ITEM:  ID LWH Weight
% Output --- ITEM: itemorder Item_Strip itemRotaFlag CoordItemStrip
% Output --- Strip: LW Weight
% Item (1 LWH (��֪)
% Item (2 Item_Strip  (dim1:���item��ĳ��strip dim2:item����˳��(��->��) 
% Item (3 CoordItemStrip Item��strip������) 
% Item (4 itemo rder ������Item����˳��)
% Strip (1 LW )
% tmpStrip_Item         (2,n): % ��1��ÿ��Strip�ڵ�Item���� �� ��2��ÿ��Strip�ڵĲ�ͬLUID����
%% ��ʼ��
% nDim Itemά��(2) nItem Item���� nStrip Strip���� 
% widthStrip Strip�����
nDim = size(Item.LWH,1); if nDim ==3, nDim = nDim-1;end
sz = size(Item.LWH);
nItem = sz(2);
wStrip = Veh.LWH(1,1);

%% ���ж�ITEM����Horizontal/Vertical ��ʽ�ڷţ������Ƿ���ת�������жϽ����㷨��˳��
% �����Ƿ�������ת, ֻ���Ƿ���Ҫ��Horizontal/Vertical��ʽ�ڷ�

% ���ԭ�ⲻ���ķ���ֵ:����ʼֵ
%  [ItemLWRota, ItemRotaed] = placeItemHori(Item,1);  %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��
% [Item.Rotaed] = placeItemHori(Item.LWH,Item.isRota,2);  %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��
% Item.LWH = getRotaedLWH(Item.LWH, Item.Rotaed, LU.buff); 

    %% ITEM���� 555
    % ��ȡItem��˳�� % ITEM��������ʽ �߶�/��̱�
%     printstruct(Item)    
    [Item.itemorder] = getITEMorder(Item,p.whichSortItemOrder );

    % ��ȡ��order������ITEM: sItem
    if isSameCol(Item)
        sItem = structfun(@(x) x(:,Item.itemorder),Item,'UniformOutput',false);
    else
        error('����ʹ��structfun');
    end
%     printstruct(sItem)

    % ������;������Ҫԭ��������Gpreproc���Լ�����H/V���ô�����
% %     % 1��2���ڵ�, sortedItemArray��Ӧ��LWHRota��Rotaed������->���󷵻ص�ԭ����ItemArry��
% %     if p.whichRotationHori == 1 % �����ĸ�level,������horizontally��ʽ�ڷ�
% %          x = sItem.Rotaed;
% %          sItem.LWH
% %         [ sItem.Rotaed] = placeItemHori(sItem.LWH,sItem.isRota,1);  %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��       
% % %          if any(x~=sItem.Rotaed),                 error('111111111111');         end    
% %     end
% %     if p.whichRotationHori == 2 % �����ĸ�level,������vertical��ʽ�ڷ�
% %                     x = sItem.Rotaed
% %         [ sItem.Rotaed] = placeItemHori(sItem.LWH,sItem.isRota,0);  %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��
% % %                     if any(x~=sItem.Rotaed),                                  error('111111111111');         end
% %     end
% %     sItem.LWH = getRotaedLWH(sItem.LWH, sItem.Rotaed, LU.buff); 
% %     sItem.LWH
%% 55 LU->Item->Stripת�� 
% ��ȡ���� �˴�ֻʹ��LWHRota��Rotaed; ��ʹ��LWH
% ItemLWRotaSort = sItem.LWH(1:2,:); %ItemLWSortHori
% ItemRotaedSort = sItem.Rotaed; % ItemRotaSortHori % itemRotaSort = zeros(1,size(Item.LWH,2));
% ItemisRotaSort = sItem.isRota;
% ItemWeightSort = sItem.Weight;

% Itemorder = Item.itemorder;

%%
% 1 Strip��ʼ��
% ��ȡStrip.LW:  �����ɵ�Strip�ĳ����
Strip.LW = zeros(2,nItem);   %strip���� dim2-����(����ߵļ���) (�߶Ƚ����ο�,current �߶�)
Strip.LW(1,:) = wStrip;   %dim1-���ʣ�� 
Strip.Weight = zeros(1,nItem); % ��ʼ��ֵ
    
% 2  ��ʱ
tmpStrip_Item = zeros(2,nItem);  % ��1��ÿ��Strip�ڵ�Item���� �� ��2��ÿ��Strip�ڵĲ�ͬLUID����

% 3 sItem����
% ��ȡCoordItemStripSort  Item��strip������ֵ
% ��ȡItem_Strip: ÿ�������Item���ĸ�Strip��  �Լ�˳��
sItem.Item_Strip = zeros(2,nItem); %dim1:���ڵڼ���level dim2:���ڸ�level�ڼ����ŷ� 555
sItem.CoordItemStrip = zeros(2,nItem); %Item��strip������ֵ


iLevel = 1; iItem = 1; %iStrip����itemʵ��
while 1
    if iItem > nItem, break; end

    % ���ݲ�ͬ�����ҵ��ܷ��뵱ǰitem��strips/levels�е�һ��    
    [thisLevel,iLevel,sItem,LU] = getThisLevel(iItem,iLevel,sItem, Strip, LU, p);     %iLevel���ڴκ����ڲ��ϵ�������Զָʾ��ǰ���µ�level

    insertItemToStrip(thisLevel,iItem);

%          plot2DStrip(); %������ͼ    
    iItem = iItem + 1;
end


%%  LU.CoordLUStri LU.LU_Strip �ļ��� Ӧ��cpustrip��
% % nbLU = size(LU.LWH,2);
% % LU.LU_Strip = zeros(2,nbLU);
% % LU.CoordLUStrip = zeros(3,nbLU);
% % 
% % % ����LU_Strip
% % for iLU=1:nbLU
% %     % ����LU_Strip��һ��
% %     iItem = LU.LU_Item(1,iLU);   %iLU���ڵڼ���Item, Item���ڵڼ���Strip,��Lu���ڵڼ���Strip
% %     LU.LU_Strip(1,iLU)= sItem.Item_Strip(1,iItem);
% %     % ����LU_Strip�ڶ���
% %     fiItem = find(sItem.Item_Strip(1,:) == sItem.Item_Strip(1,iItem) & sItem.Item_Strip(2,:) < sItem.Item_Strip(2,iItem));
% %     nbLUfiItem = sum(ismember(LU.LU_Item(1,:),fiItem));
% %     LU.LU_Strip(2,iLU) = nbLUfiItem+LU.LU_Item(2,iLU); % ����Strip˳��: ͬһStrip����ǰ�������nbLUfiItem + ��iLU��Item��˳��
% %     % ����LU.CoordLUStrip
% %     LU.CoordLUStrip(1,iLU) = sItem.CoordItemStrip(1,iItem);
% %     LU.CoordLUStrip(2,iLU) = sItem.CoordItemStrip(2,iItem);
% %         % fLU: ��iLUͬ��iItem �� ˳�����ڱ�iLU; ����Ϊ��, ��Ӱ��.
% %     fLU = LU.LU_Item(1,:) == iItem & LU.LU_Item(2,:) < LU.LU_Item(2,iLU);
% %     LU.CoordLUStrip(3,iLU) = sum(LU.LWH(3,fLU));
% % end

%%
% Item��أ����µİ�˳�򷵻أ��޸��µĲ��践�أ�ʵ�����ջ�ȡ��Strip˳�򲢲���Ҫ
% LU�ڲ�����,sLU����order�仯����(��ҪΪ��sLU�������ļ�������,Ҫ��˳��ת����)
if isSameCol(sItem)
    Item = reorderStruct(Item.itemorder, sItem);
else
    error('����ʹ��structfun');
end

                    % ��ȡItem_Strip : ÿ��Item���ĸ�Strip��  �Լ�˳��
                    % ��ȡCoordItemStrip : ÿ��Item��Strip������
                    %     Item.Item_Strip(:,Item.itemorder) = sItem.Item_Strip;
                    %     Item.CoordItemStrip(:,Item.itemorder) = sItem.CoordItemStrip;    
                    % ItemArray��ת���
                    %         Item.Rotaed(:,Item.itemorder) = sItem.Rotaed;
                    %     ItemLWRota(:,Item.itemorder) = sItem.LWH(1:2,:);
                    %     Item.LWH = [ItemLWRota; Item.LWH(3,:)];             % ����ԭʼ˳�����ת���ItemArray


                    % LUArray��ת���,��ʱ����    ��Ƕ�뵽rotateItem
                %     nbItem=length(Item.Rotaed);
                %     % ѭ��ÿ��item
                %     for idxItem=1:nbItem
                %         tmpflagThisItem = (LU.LU_Item(1,:)==idxItem );
                %         % ��Ӧλ��LU.Rotaed����
                %         if Item.Rotaed(idxItem)
                %             LU.Rotaed(tmpflagThisItem) = ~LU.Rotaed(tmpflagThisItem);
                %             % ��Ӧλ��LU.LWH����
                %             LU.LWH(1, tmpflagThisItem) = Item.LWH(1, idxItem);
                %             LU.LWH(2, tmpflagThisItem) = Item.LWH(2, idxItem);
                %         end
                %     end
    
% Strip���: ��˳�����(ȥ����ʼ���Ķ�����)
% Strip.LW
% ��ȡLWStrip:  �����ɵ�strip�ĳ���
% ��ȡStripWeight:  �����ɵ�strip������
% ���Strip������ȫ����ͬ
if isSameCol(Strip)
    Strip = structfun(@(x) x( : , Strip.Weight(1,:)>0 ), Strip, 'UniformOutput', false);
else
    error('����ʹ��structfun');
end
                 

%% ����script
%     printscript();

    
    %% Ƕ�׺���        
    function insertItemToStrip(thisLevel,iItem)       
        % 1 ����Item���Sort����
        %  1.1 ����CoordItemStripSort
        sItem.CoordItemStrip(1,iItem) = wStrip - Strip.LW(1,thisLevel);        %����x����
        sItem.CoordItemStrip(2,iItem) = sum(Strip.LW(2,1:thisLevel-1));      %����y���� %���iLevel=1,�����ߣ�����Ϊ0������Ϊ���
        
        % 2 ����Strip������ݣ�������
        %  2.1 ����LWStrip
        Strip.LW(1,thisLevel) = Strip.LW(1,thisLevel) - sItem.LWH(1,iItem); %����wleft (�ڷŷ���ǰ��һ��)
        Strip.LW(2,thisLevel) = max(Strip.LW(2,thisLevel), sItem.LWH(2,iItem)); %����strip�߶�lleft(ȡ���ֵ)
        
        %  2.2 ����Strip.Strip_Item ��1 ��strip�ڰ�������Item
        tmpStrip_Item(1,thisLevel) = tmpStrip_Item(1,thisLevel) + 1; %ֻҪ��level����һ��item,����������1
        
        %  2.3 ���±�level��Ӧ��StripWeight: 
        Strip.Weight(thisLevel) =  Strip.Weight(thisLevel) + sItem.Weight(iItem);
        
        %  1.3 ����item����strip��ϢitemBeStripMatrixSort
        sItem.Item_Strip(1,iItem) = thisLevel;    %�ڼ���level
        sItem.Item_Strip(2,iItem) = tmpStrip_Item(1,thisLevel); %��level�µڼ��ΰ���
        
        
        
        %         updateLWStrip(); %�����ж��Ƿ�������ת�������ж��Ƿ�������Level�������ж��Ƿ�ǰ�ɷ���
        
                % % %         if sItem.isRota(iItem) == 1 %��Item������ת
                % % %             % �ж����: 2��
                % % %             isflagCurr = Strip.LW(1,thisLevel) >=  sItem.LWH(1,iItem); %�ж��Ƿ�current's stripʣ���� >= ��ǰ�߶ȣ�����ת��
                % % %             isNewLevel = Strip.LW(1,thisLevel) == wStrip; % �ж��Ƿ� new Level            
                % % %             % ����strip��Ϣ
                % % %             if isNewLevel %�������,�����Է���,���ۺ��ְڷ�,���:ֱ�Ӹ��£�Item�Ѱ�Hori/Vert�ڷŹ���
                % % %                     updateLWStrip();
                % % %             else % �������level ����԰ڷ�,����; ����,����������ת�������
                % % %                     if isflagCurr
                % % %                         updateLWStrip();
                % % %                     else
                % % %                         error('11111111');
                % % % %                         rotateItem(); %����ֱ��ע��
                % % % %                         updateLWStrip();
                % % %                     end
                % % %              end            
                % % %         elseif sItem.isRota(iItem) == 0 %��Item��������ת
                % % %             % ����strip��Ϣ
                % % %             updateLWStrip();
                % % %         end
        

        
                        %  2.2 ����Strip.Strip_Item ��2 ��strip�ڰ�������Item
                %         itemThisLevel = sItem.Item_Strip(1,:) == thisLevel;
                %         Strip.Strip_Item(2,thisLevel) = numel(unique(sItem.LID(1,itemThisLevel)));

        % 4 ����Ƕ�׺���- ����������Ҫ
% %         function rotateItem()
% %             %  1 �������Rotaed�仯 ��Ҫ��ITEM������rotate(��)��ȥ
% %             sItem.Rotaed(iItem) = ~sItem.Rotaed(iItem);
% %             tep = sItem.LWH(1,iItem);
% %             sItem.LWH(1,iItem) = sItem.LWH(2,iItem);
% %             sItem.LWH(2,iItem) = tep;
% %             sItem.LWH(3,iItem) = sItem.LWH(3,iItem);
% %             
% %             %  2 Item�ڲ���LuҲҪ��֮����
% %             tmpflagThisItem = (LU.LU_Item(1,:)==iItem);
% %             % ��Ӧλ��LU.Rotaed����
% %             if Item.Rotaed(iItem)
% %                 LU.Rotaed(tmpflagThisItem) = ~LU.Rotaed(tmpflagThisItem);
% %                 % ��Ӧλ��LU.LWH�ĳ�����£��߶Ȳ�����
% %                 LU.LWH(1, tmpflagThisItem) = Item.LWH(1, iItem);
% %                 LU.LWH(2, tmpflagThisItem) = Item.LWH(2, iItem);
% %             end
% %         end
        
%         function updateLWStrip()
% %             Strip.LW(1,thisLevel) = Strip.LW(1,thisLevel) - sItem.LWH(1,iItem); %����wleft (�ڷŷ���ǰ��һ��)
% %             Strip.LW(2,thisLevel) = max(Strip.LW(2,thisLevel), sItem.LWH(2,iItem)); %����strip�߶�lleft(ȡ���ֵ)
            
            % ����Strip�а���ID�����
%             Strip.LID(sItem.LID(1,iItem),thisLevel) =  1;         
%             Strip.SID(sItem.SID(1,iItem),thisLevel) =  1;        
%             Strip.UID(sItem.UID(1,iItem),thisLevel) = 1;         % ��ֵΪ�������

%              Strip.PID(:,thisLevel) = Strip.PID(:,thisLevel) + sItem.PID(:,iItem); % ��ֵΪ���ִ���
%              Strip.PID(Strip.PID>0) = 1; % ��ֵ��Ϊ�������
            
%             Strip.SID(sItem.SID(1,iItem),thisLevel) = 1;         % 555 ���¶���PID
%             Strip.UID(sItem.UID(1,iItem),thisLevel) = 1;         % 555 ���¶���PID
%             Item.PID(sLU.PID(1,iLU),thisItem) = 1;         % 555 ���¶���PID
            
%         end
        
    end
    
    function printscript()
        % ���Դ���
        % % LWHStrip
        % % sItem.Item_Strip
        % % ItemLWSort
        % % itemCoordMatrixSort
        % % Item_Strip
        % % LWHItem
        % % itemCoordMatrix
        %  printstruct(d);
        
        % �����Ҫ���:��ô�1��ʼÿ��strip����������
        for iStrip = 1:max(Item.Item_Strip(1,:))
            [~,idx] = find(Item.Item_Strip(1,:)==iStrip);
            fprintf('strip %d ��ʣ���+���Ϊ:  ',iStrip);
            fprintf('( %d ) ',Strip.LW(:,iStrip));
            fprintf('\n');
            fprintf('strip %d ���� original Item ������(����)[��ת��־]{����}Ϊ  \n  ',iStrip);
            fprintf('%d ',idx);
            fprintf('( %d ) ', Item.LWH(1:nDim,idx));fprintf('\n');
            fprintf('[ %d ] ', Item.Rotaed(:,idx));fprintf('\n');  %sItem.Rotaed
            fprintf('{ %d } ', Item.CoordItemStrip(:,idx));fprintf('\n');
            fprintf('\n');
        end
    end

    function plot2DStrip()
        %% ��ʼ��
        wStrip = wStrip;        
        hStrip = sum(Strip.LW(2,sItem.Item_Strip(2,:)>0));        
        nstrip = sum(sItem.Item_Strip(2,:)>0);

        tmpLID = cellfun(@(x)x(1), sItem.LID);        
        nIDType = unique(tmpLID);
        nColors = hsv(length(nIDType)); %��ͬ����LU���費ͬ��ɫ
        
        %% ��ͼ
        % 1 ��ͼ��������Strip
        DrawRectangle([wStrip/2 hStrip/2 wStrip hStrip 0],'--', [0.5 0.5 0.5]);
        hold on;
        % 2 ��ͼ�����strip/item ��ͼ
        for istrip = 1:nstrip
            % �ҳ���ǰistrip����Ʒ����
            idxDrawItem = find(sItem.Item_Strip(1,:)==istrip);
            % ��ȡ�������µı���
            drawItemCoordMatrix = sItem.CoordItemStrip(:,idxDrawItem);
            drawItemLWH = sItem.LWH(:,idxDrawItem);
            drawItemId = tmpLID(:,idxDrawItem);

            % ��ͼ�����item
            nThisItem = size(drawItemLWH,2);
            for iplotItem = 1:nThisItem
                % ��ͼ��������iItem
                itemWidth = drawItemLWH(1,iplotItem);
                itemLength = drawItemLWH(2,iplotItem);
                itemCenter = [drawItemCoordMatrix(1,iplotItem)+itemWidth/2 ...
                    drawItemCoordMatrix(2,iplotItem)+itemLength/2 ];
                
                % ���ӶԱ���iItem�����ͣ���ɫ���ж�
                itemID = drawItemId(iplotItem);
                itemColor = 0.8*nColors(nIDType==itemID, : );
                
                DrawRectangle([itemCenter itemWidth itemLength 0],  '-',itemColor);
                hold on;
            end
        end
        % hold off;
    end
end

%% �ֲ����� %%

%% ����1: getITEMorder
% ����ITEM��˳��,��NEXT FIT�ķ�ʽ����STRIP���Ȳ���SIDС��; �����߶�/��ȣ� ����LID��
function order = getITEMorder(Item,whichSortItemOrder)

%��SID����: SID������˳������,���С����ǰ��
szRow = cellfun(@(x)size(x,1), Item.SID);
if (max(szRow)~=min(szRow)),  error('ͬһITEM��Ӧ���ж��SID');  end %ͬһItemӦ��ֻ��һ��SID
SIDorder = cell2mat(Item.SID);   %ֱ��cell2matת��; %ITEM��SID 1-n��˳�򷵻� 

%��LID����: LID��ָ��˳��, ����SID����ȫ��һ��,�ٰ�LID��С��������,��ʵû������(��SID/LID����ͬһITEM),��󿴸߶� 
szRow = cellfun(@(x)size(x,1), Item.LID);
if (max(szRow)~=min(szRow)),  error('ͬһITEM��Ӧ���ж��SID');  end %ͬһItemӦ��ֻ��һ��LID
LIDorder = cell2mat(Item.LID);   %ֱ��cell2matת��; %ITEM��SID 1-n��˳�򷵻� 

% V2: ********** ����isNonMixed
global ISisNonMixed ISisMixTile
% Ŀǰ˳�� : 1: SID ; 2: isNonMixed;(��ͬLID��) һ�����濪ʼ: 3: Longth/Height; 4:Width; 5: LID; (3,4,5,����һ��) 6: Height
tmpItem = [SIDorder; Item.isNonMixed; Item.isMixedTile; Item.LWH(2,:); Item.LWH(1,:); LIDorder; Item.LWH(3,:); ];
if ISisNonMixed==1    
    if ISisMixTile==1
        [~,order] = sortrows(tmpItem',[1, 2, 3, 4, 5, 6, 7 ],{'ascend','descend','ascend','descend','descend','descend','descend'});
    else
        [~,order] = sortrows(tmpItem',[1, 2, 4, 5, 6, 7 ],{'ascend','descend','descend','descend','descend','descend'});
    end
else
    [~,order] = sortrows(tmpItem',[1, 4, 5, 6, 7 ],{'ascend','descend','descend','descend','descend'});
end

%%
% % % clc
% % % s = [1 1 1 2 2;
% % %        10 10 10 12 12;
% % %        8 8  8 7 7;
% % %        10 9 6 6 8]
% % %    d=(diff(s(4,:)))
% % %     s=[s;[100,d]]
% % % %     [x,order] = sortrows(s',[1, 2, 3, 4,5 ],{'ascend','ascend','descend','descend','ascend'});
% % %     [x,order] = sortrows(s',[1, 2, 3, 5 ],{'ascend','ascend','descend','ascend'});
% % %     x'
% % %     order'
% % %    1
%%
% itemPriority = getPriorityofItem(SIDorder,Item.isNonMixed, Item.isHeightFull, Item.LWH)
% tmpSort = [SIDorder; itemPriority];
% [~,order] = sortrows(tmpSort',[1,2],{'ascend','ascend'});

if ~isrow(order), order=order'; end
end


% % function itemPriority = getPriorityofItem(SIDorder, isNonMixed, isHeightFull, LWH)
% % 
% % itemPriority = zeros(1,size(isNonMixed,2));
% % uniOrd = unique(SIDorder);
% % 
% % % ��ͬSID��ȡ��ͬ��˳��, ��SID�����1��ʼ
% % for i=1:length(uniOrd)    
% %     idxSID = SIDord==uniOrd(i);
% %     
% %     % UNUSED ��ȡ��SID�ڵ�STRIP��Ӧ��ֵtLW;tLL;tL;tID;celltID.
% %     tLW = SLW(:,idxSID);   tLL = SLoadingRateLimit(:,idxSID);
% %     tID = SLID(:,idxSID);    celltID = tID;
% %     tPured = Spured(:,idxSID);
% %     tSingle = Ssingle(:,idxSID);   %����δ֪
% %     
% %     % FINALLY USED ���ȳ�ͷ�ڷŵ�˳��
% %     tnbItem = SnbItem(:,idxSID);
% %     tMixed = SisMixed(:,idxSID);
% %     tFull = SHeightfullfull(:,idxSID);
% %     tLR = SLoadingrate(:,idxSID);        % max(tL)  %  min(tL)  % mean(tL)
% %     
% %     priority = 1;
% %     SIDorder = zeros(1,size(tID,2));
% %     szRow = cellfun(@(x)size(x,1), tID);
% %     if isscalar(szRow) %�����STRIP��ֻ��1STRIP����Ϊ������STRIP����ֵΪ1
% %         SIDorder = 1;
% %     else
% %         [tID, ~] = padcat(tID{:});   if iscolumn(tID), tID = tID'; end
% %         
% % 
% %         tmpM = [tnbItem; tMixed; tFull; tLL; tLR; tLW ];  
% %         [~,order] = sortrows(tmpM',[1,2,3,5],{'descend','ascend','descend','descend'}); 
% % 
% %         if ~isrow(order), order=order'; end
% %         
% %         % 2 ����1�İڷ�˳��, ��������STRIP˳��torder
% %         while any(SIDorder==0)
% %             [~,o]=find(SIDorder(order) == 0,1,'first');
% %             SIDorder(order(o)) = priority;  
% %             priority=priority+1;
% %             
% %             % 2.2 �ҳ���order(o)λ��tLID��Ӧ������Strip.
% %             tnbItem = celltID{:,order(o)};
% %             [SIDorder,priority] = getAdjPriority(priority,order,SIDorder,tID,tnbItem);
% %         end
% %     end
% %     itemPriority(idxSID) = SIDorder;
% % end
% % end


% % % ����ITEM��˳��,��NEXT FIT�ķ�ʽ����STRIP���Ȳ���SIDС��; �����߶�/��ȣ� ����LID��
% % function order = getITEMorder(Item,whichSortItemOrder)
% % 
% % %��SID����: SID������˳������,���С����ǰ��
% % szRow = cellfun(@(x)size(x,1), Item.SID);
% % if (max(szRow)~=min(szRow)),  error('ͬһITEM��Ӧ���ж��SID');  end %ͬһItemӦ��ֻ��һ��SID
% % SIDorder = cell2mat(Item.SID);   %ֱ��cell2matת��; %ITEM��SID 1-n��˳�򷵻� 
% % 
% % %��LID����: LID��ָ��˳��, ����SID����ȫ��һ��,�ٰ�LID��С��������,��ʵû������(��SID/LID����ͬһITEM),��󿴸߶� 
% % szRow = cellfun(@(x)size(x,1), Item.LID);
% % if (max(szRow)~=min(szRow)),  error('ͬһITEM��Ӧ���ж��SID');  end %ͬһItemӦ��ֻ��һ��LID
% % LIDorder = cell2mat(Item.LID);   %ֱ��cell2matת��; %ITEM��SID 1-n��˳�򷵻� 
% %     
% % %��ITEM����/�����/���߶� ���������� ������ͬIDLU���ֿ�
% % % ����LUID 2: ȷ����ʹ������ȫ��ͬ ��LUID��ͬ�� Ҳ�����һ��
% % 
% % % V2: ********** ����isNonMixed
% % % �����ж�Item�ȳ���,��߶�����. % �����nbcol2����:
% % global ISisNonMixed
% % ISisNonMixed
% % % Ŀǰ˳�� : 1: SID ; 2: isNonMixed; һ�����濪ʼ: 3: Longth/Height; 4:Width; 5: LID; (3,4,5,����һ��) 6: Height
% % tmpItem = [SIDorder; Item.isNonMixed; Item.LWH(2,:); Item.LWH(1,:); LIDorder; Item.LWH(3,:); ];
% % if ISisNonMixed==1    
% %     [x,order] = sortrows(tmpItem',[1, 2, 3, 4, 5, 6 ],{'ascend','descend','descend','descend','descend','ascend'});
% % else
% %     [~,order] = sortrows(tmpItem',[1, 3, 4, 5, 6 ],{'ascend','descend','descend','descend','descend'});
% % end
% % 
% % 
% % % tmpItem = [SIDorder; LIDorder; Item.LWH; Item.isNonMixed]; % tmpItem = [SIDorder; LIDorder; Item.LWH; Item.isNonMixed; Item.isHeightFull];
% % % [~,order] = sortrows(tmpItem',[1,6, 4, 3, 2, 5 ],{'ascend','descend','descend','descend','descend','ascend'});
% % % [~,order] = sortrows(tmpItem',[1,6, 2, 4, 3, 5 ],{'ascend','descend','descend','descend','descend','descend'});  
% % 
% % 
% %         % V1: *********** ������isNonMixed
% %         % % tmpItem = [SIDorder; LIDorder; Item.LWH; ];  % tmpItem = [Item.SID; Item.LID; Item.LWH];  % tmpItem = [ Item.LWH];
% %         % % % [~,order] = sortrows(tmpItem',[1, 4, 2, 5],{'ascend','descend','descend','descend'}); 
% %         % % [~,order] = sortrows(tmpItem',[1, 4, 3, 2, 5],{'ascend','descend','descend','descend','descend'});  
% % 
% % 
% % if ~isrow(order), order=order'; end
% % end

%% ����2: getThisLevel
function [thisLevel,iLevel,sItem, LU] = getThisLevel( iItem, iLevel, sItem, Strip, LU,p)
% ��ͬwhichStripH��,��ù�ͬ��thisLevel
if p.whichStripH == 1 % 1 bestfit 2 firstfit 3 nextfit
    % ����Rotation���ӱ���
    if sItem.isRota(iItem) == 1 %��Item������ת
        %             if ParaArray.whichRotation == 1
        % �ҵ�����rotation�µ�level:��һ�ڷŷ���ɷ����iItem��level
        flag = find(Strip.LW(1, 1 : iLevel) >= sItem.LWH(1,iItem) |  ...
            Strip.LW(1, 1 : iLevel) >= sItem.LWH(2,iItem));
    else %��Item��������ת
        % ���������µ�ѡ��find����㹻�Ķ��level,����������Сʣ��ˮƽ��ȵ�
        flag = find(Strip.LW(1, 1 : iLevel) >= sItem.LWH(1,iItem));
    end
    if isempty(flag)
        iLevel = iLevel + 1;% �����Ȳ����㣬��level����
        [thisLevel,iLevel,sItem, LU] = getThisLevel( iItem, iLevel, sItem,  Strip, LU, p);
    else
        % ��ȡthisLevel: Ψһ��FF������⵽thisLevel�ļ��㣨ѡ��������������С�ģ�
        tmpLevels = Strip.LW(1,1:iLevel);   %��ȡ�����Ѱ��Ż��°��ŵ�level��ʣ��ˮƽƽ�������tmpLevels
        tepMinLeftWdith = min(tmpLevels(flag));                      %�ҳ�tepAvailableLevelArray�п����ɱ�iITem��ʣ��ˮƽ��ȵ���Сֵ������ΪtepMinLeftWdith
        thisLevel = find(tmpLevels==tepMinLeftWdith);            %�ҳ�����Сֵ��Ӧ���Ǹ�/Щlevel
        if ~all(ismember(thisLevel,flag)),     error('Not all thisLevel belongs to flag ');          end
        if length(thisLevel)>1
            thisLevel = thisLevel(1);
        end
    end
elseif p.whichStripH == 2 % firstfit  % firstfit���ܲ���ֱ������bestfit�Ĵ���?
    % ����Rotation���ӱ���
    if sItem.isRota(iItem) == 1 %��Item������ת
        % �ҵ�����rotation�µ�level:��һ�ڷŷ���ɷ����iItem��level
        flag = find(Strip.LW(1, 1 : iLevel) >= sItem.LWH(1,iItem) |  ...
            Strip.LW(1, 1 : iLevel) >= sItem.LWH(2,iItem));
    else
        % ���������µ�ѡ��find����㹻�Ķ��level,�������ڵ�һ�������� Ψһ������thisLevel�Ļ�ȡ
        flag = find(Strip.LW(1, 1 : iLevel) >= sItem.LWH(1,iItem));
    end
    if isempty(flag)
        iLevel = iLevel + 1;% �����Ȳ����㣬��level����
        [thisLevel,iLevel,sItem, LU] = getThisLevel( iItem, iLevel, sItem, Strip, LU, p);
    else
        thisLevel = flag(1);
        if ~all(ismember(thisLevel,flag)),     error('Not all thisLevel belongs to flag ');       end
    end
elseif p.whichStripH == 3 % nextfit
    % ����Rotation���ӱ���
    sItem.itemorder
    if sItem.isRota(iItem) == 1 %��Item������ת % nextfit�²���ֱ������bestfit�Ĵ���
        % V2: iItem�Ѿ�����ת���,��Ѱڷ�λ��,�����ж���ת����������
        % �ж���ǰlevel�Ƿ��������һ�ڷŷ���ɷ����iItem flaged: �����ݱ�ʾ���ԣ����򲻿���        
        flaged = find(Strip.LW(1,iLevel) >= sItem.LWH(1,iItem) );
        
% %         global ISsItemAdjust
% %         if ISsItemAdjust==1
% %         % �����iLevel�ǿ�,�ҷŵ��±���iITem, ����
% %         if  ~isempty(flaged) && Strip.LW(1,iLevel) < 2400 
% % %         Strip.LW(1,iLevel)%         sItem.LWH(1,iItem)%         sItem.isHeightFull(iItem)%         sItem.HLayer(iItem)%         sItem.LWH(3,iItem) 
% %         tmpLIDmat = cell2mat(sItem.LID);
% %         
% %         fIarray = sItem.Item_Strip(1,:) == iLevel; %�Ѱ�����iLvel�ڵ�Item
% %         fIscalar =  fIarray & sItem.Item_Strip(2,:) == max(sItem.Item_Strip(2,fIarray)) %�Ѱ�����iLvel�ڿ������(������)��һ��        
% %         
% %         fLIDfIscalar = tmpLIDmat == sItem.LID{fIscalar} % ��fIscalarͬLID��fֵ
% %         fLIDiItem = tmpLIDmat == sItem.LID{iItem}        % ��iItemͬLID��fֵ
% %         
% %         if sum(fIscalar)~=1, error('fff'); end
% %         % ������iLevel�������ŵ�Item���жԱ�, ���Ǳ������ڲ�ͬ��LIDֵ; ����ͬIdֵ, ֱ���˳�
% %         if sItem.LID{fIscalar} ~= sItem.LID{iItem} 
% %             % ��Ȼת��LID, �����fIscalar����û����ͬ������ͬID��
% %             if ismember(sItem.LID{fIscalar},tmpLIDmat(~fLIDfIscalar))    error('��Ȼת��LID, �����fIscalar����û����ͬ������ͬID��');   end 
% %             
% %                             %     iItemsLID = ismember(tmpLIDmat(~fIarray),sItem.LID{iItem}); % ����iItem�Ƿ�����ͬID, ���û��, ��ɶҲ��˵��, ����iItem����level
% %             if sum(fLIDiItem) > 0,  % �����iItem��ӦLID����1��, ɶҲ��˵��, ֱ���˳�
% %                     %     sItem.LWH(3,fIscalar) % sItem.LWH(3,iItem)                            
% %                 ItemHeightwithSameLIDofiItem = sItem.LWH(3,fLIDiItem)
% %                 
% %                 if ~issorted(ItemHeightwithSameLIDofiItem,'descend') || isempty(ItemHeightwithSameLIDofiItem)
% %                     error('ͬһLID�µ�Item�߶ȷǵݼ������Ϊ��ֵ, ��Ԥ�ڴ���')
% %                 end
% %                 
% %                 %���������Сֵ����, ���ͷ�ڷ�
% %                 if abs(sItem.LWH(3,fIscalar) - ItemHeightwithSameLIDofiItem(1)) > abs(sItem.LWH(3,fIscalar) - ItemHeightwithSameLIDofiItem(end)) 
% %                    ord = 1:length(sItem.Weight);
% %                    ord(:,fLIDiItem)= fliplr(ord(:,fLIDiItem)); % ������iItem��ͬLID��Items��˳��
% %                    ord(:,fLIDiItem)
% %                    sItem.LWH(3,iItem) 
% %                    % ����ItemҪ��, ��Ӧ��LUҲҪ��
% %                    zzz = sItem.itemorder
% %                    sItem = structfun(@(x) x(:,ord),sItem,'UniformOutput',false);
% %                    sItem.itemorder  = zzz               
% %                    
% %                    fLIDiItemIdx = find(fLIDiItem);
% %                    FLIPfLIDiItemIdx = fliplr(fLIDiItemIdx);
% %                    tmpLU_Item = LU.LU_Item(1,:);
% %                    for i=1:length(fLIDiItemIdx)
% %                        fLIDLU = tmpLU_Item == fLIDiItemIdx(i); %SAME : fLIDLU = ismember(LU.LU_Item(1,:) ,fLIDiItemIdx(i))                       
% %                        LU.LU_Item(1, fLIDLU )  = FLIPfLIDiItemIdx(i);
% %                    end
% %                    % ͬʱ����LU.DOC
% %                    LU.DOC(end-1:end,:) = LU.LU_Item;
% %                         LU
% %                         sItem
% %                         sItem.Item_Strip
% %                         sItem.CoordItemStrip
% %                    % �������
% %                    flaged = find(Strip.LW(1,iLevel) >= sItem.LWH(1,iItem) );
% %                    if isempty(flaged), error('ͬLID����Ȳ�ͬ, ��Ԥ�ڴ���'); end
% %                    
% %                 end
% %                 
% %             end
% %         end
% %         %     sItem.LWH(3,fIscalar)  %�Ѱ�����iLvel�� (�����һ����Item�߶�) %         sItem.LWH(1,iItem) %         sItem.isHeightFull(fI)         sItem.isHeightFull(iItem)
% %         end
% %         end % END OF ISsItemAdjust
        
        
                % V1 : �ж���ת���
                %                 flaged = find(Strip.LW(1,iLevel) >= sItem.LWH(1,iItem) |  ...
                %                                       Strip.LW(1,iLevel) >= sItem.LWH(2,iItem));
    else
        % ��ͬ�����µ�ѡ�������ǰitem�Ŀ�<=��ǰstrip�ĵ�ǰlevel�Ŀ�
        flaged = find(Strip.LW(1,iLevel) >= sItem.LWH(1,iItem) );
    end
    
        
    if  isempty(flaged)  %ע����֮ǰ~flag������
        iLevel = iLevel + 1;% �����Ȳ����㣬��level����
        [thisLevel,iLevel,sItem, LU] = getThisLevel( iItem, iLevel, sItem, Strip, LU, p) ;
    else
        if  isempty(flaged) ,   error(' �����ܵĴ��� ');      end
        thisLevel = iLevel; % ��ǰlevelһ���ŵ���
    end
    %             end
    
end
end

