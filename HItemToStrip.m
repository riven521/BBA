% function [d] = HItemToStrip(d,p)
function [Item,Strip]= HItemToStrip(LU,Item,Veh,p)
% ��Ҫ����:Item����Strip�� %  ����:�����(row);  ����:��������(coloum);
% Input ---  ITEM:  ID LWH Weight
% Output --- ITEM: itemorder Item_Strip itemRotaFlag CoordItemStrip
% Output --- Strip: LW Weight
% Item (1 LWH (��֪)
% Item (2 Item_Strip  (dim1:���item��ĳ��strip dim2:item����˳��(��->��) 
% Item (3 CoordItemStrip Item��strip������) 
% Item (4 itemo rder ������Item����˳��)
% Strip (1 LW )

%% ��ʼ��
% nDim Itemά��(2) nItem Item���� nStrip Strip���� 
% widthStrip Strip�����
nDim = size(Item.LWH,1); if nDim ==3, nDim = nDim-1;end
sz = size(Item.LWH);
nItem = sz(2);
nStrip = nItem;
wStrip = Veh.LWH(1,1);;

%% ���ж�ITEM����Horizontal/Vertical ��ʽ�ڷţ������Ƿ���ת�������жϽ����㷨��˳��
% �����Ƿ�������ת, ֻ���Ƿ���Ҫ��Horizontal/Vertical��ʽ�ڷ�

%  [ItemLWRota, ItemRotaed] = placeItemHori(Item,1);  %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��
% ���ԭ�ⲻ���ķ���ֵ:����ʼֵ
[Item.Rotaed] = placeItemHori(Item.LWH,Item.isRota,2);  %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��
Item.LWH = getRotaedLWH(Item.LWH, Item.Rotaed, LU.buff); 

    %% ITEM���� 555
    % getITEMorder - ��ȡItemLWRota��˳��(�ص��Ǹ߶ȵݼ�����) % ITEM��������ʽ �߶�/��̱�
    Item.itemorder = getITEMorder(Item,p.whichSortItemOrder );
    % getSortedITEM - ��ȡ��order������ITEM:sortedItemArray
    sItem = structfun(@(x) x(:,Item.itemorder),Item,'UniformOutput',false);
                                    % printstruct(d) ;printstruct(sItem)

    % 1��2���ڵ�, sortedItemArray��Ӧ��LWHRota��Rotaed������->���󷵻ص�ԭ����ItemArry��
    if p.whichRotationHori == 1 % �����ĸ�level,������horizontally��ʽ�ڷ�
        [ sItem.Rotaed] = placeItemHori(sItem.LWH,sItem.isRota,1);  %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��
    end
    if p.whichRotationHori == 2 % �����ĸ�level,������vertical��ʽ�ڷ�
        [ sItem.Rotaed] = placeItemHori(sItem.LWH,sItem.isRota,0);  %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��
    end
    sItem.LWH = getRotaedLWH(sItem.LWH, sItem.Rotaed, LU.buff); 
     
%% 55 LU->Item->Stripת�� 
% ��ȡ���� �˴�ֻʹ��LWHRota��Rotaed; ��ʹ��LWH
% ItemLWRotaSort = sItem.LWH(1:2,:); %ItemLWSortHori
% ItemRotaedSort = sItem.Rotaed; % ItemRotaSortHori % itemRotaSort = zeros(1,size(Item.LWH,2));
% ItemisRotaSort = sItem.isRota;
% ItemWeightSort = sItem.Weight;

% Itemorder = Item.itemorder;

%%
% ��ȡitemBeStripMatrixSort: ÿ�������Item���ĸ�Strip��  �Լ�˳��
% ��ȡLWStrip:  �����ɵ�Strip�ĳ���
% ��ȡCoordItemStripSort  Item��strip������ֵ
Strip.LW = zeros(2,nItem);   %strip���� dim2-����(����ߵļ���) (�߶Ƚ����ο�,current �߶�)
Strip.LW(1,:) = wStrip;   %dim1-���ʣ�� 
Strip.Weight = zeros(1,nStrip); % ��ʼ��ֵ
 
Strip_Item = zeros(1,nStrip);  % ÿ��Strip�ڵ�Item���� ���ڲ���

sItem_Strip = zeros(2,nItem); %dim1:���ڵڼ���level dim2:���ڸ�level�ڼ����ŷ� 555
sCoordItemStrip = zeros(2,nItem); %Item��strip������ֵ

% 55 ��ȡthisLevel - ��ǰitemҪ�����level���
% ѭ����strip�а���item,���̶�item,�仯ѡ��ͬlevel(thisLevel)
% ע�ͣ���ȡ FLAG        �ɷ��µ�ǰitem(iItem)������һ��level�ļ��ϣ� 
% ע�ͣ���ȡ thisLevel   ��FLAG���ҵ���������Ǹ�thisLevel, ��ִ�� insert����

iLevel = 1; iItem = 1; %iStrip����itemʵ��
while 1
    if iItem > nItem, break; end
    
    % ���ݲ�ͬ�����ҵ��ܷ��뵱ǰitem��strips/levels�е�һ��    
    [thisLevel,iLevel] = getThisLevel(iItem,iLevel,sItem, Strip, p);     %iLevel���ڴκ����ڲ��ϵ�������Զָʾ��ǰ���µ�level
    
    insertItemToStrip(thisLevel);
    
%     plot2DStrip(); %������ͼ    
    iItem = iItem + 1;
end

% plot2DStrip(); %һ���Ի�ͼ

% ���� ����ֵ��d
%Matalb code gerator use:
%         Item_Strip=sItem_Strip; CoordItemStrip=sCoordItemStrip;

% Item��أ����µİ�˳�򷵻أ��޸��µĲ��践�أ�
% sItem_Strip
% sCoordItemStrip
% sItem.Rotaed
% sItem.LWH(1:2,:)
% ��ȡitemBeStripMatrix : ÿ��Item���ĸ�Strip��  �Լ�˳��
% ��ȡCoordItemStrip : ÿ��Item��Strip������
% ��ȡRotaed : ÿ��item�Ƿ�Rotation�ı�־
% ��ȡLWHRota��ÿ��item���Rotaed��־���õ�LWH

    Item.Item_Strip(:,Item.itemorder) = sItem_Strip;  
    Item.CoordItemStrip(:,Item.itemorder) = sCoordItemStrip;    
    
    % ItemArray��ת���
    Item.Rotaed(:,Item.itemorder) = sItem.Rotaed;
    ItemLWRota(:,Item.itemorder) = sItem.LWH(1:2,:);
    Item.LWH = [ItemLWRota; Item.LWH(3,:)];             % ����ԭʼ˳�����ת���ItemArray

    % LUArray��ת���,��ʱ����    
    nbItem=length(Item.Rotaed);
    % ѭ��ÿ��item
    for idxItem=1:nbItem
        flagThisItem = (LU.LU_Item(1,:)==idxItem );
        % ��Ӧλ��LU.Rotaed����
        if Item.Rotaed(idxItem)
            LU.Rotaed(flagThisItem) = ~LU.Rotaed(flagThisItem);
            % ��Ӧλ��LU.LWH����
            LU.LWH(1, flagThisItem) = Item.LWH(1, idxItem);
            LU.LWH(2, flagThisItem) = Item.LWH(2, idxItem);
        end
    end

    
% Strip���: ��˳�����
% Strip.LW
% ��ȡLWStrip:  �����ɵ�strip�ĳ���
% ��ȡStripWeight:  �����ɵ�strip������

    Strip.LW = Strip.LW(:,Strip.LW(2,:)>0); % û��˳�� + ȥ��δʹ�õ�Strip    
    Strip.Weight = Strip.Weight(Strip.Weight(:)>0); % û��˳�� + ȥ��δʹ�õ�Strip    
    
    %% ����script
    % �����Ҫ���:���ÿ��level������ 
    printscript();
%     printstruct(d);
    
    %% Ƕ�׺���        

    function insertItemToStrip(thisLevel)
%         %Ϊ��matlab��coder ����
%         sCoordItemStrip=sCoordItemStrip;Strip.LW=Strip.LW;
%         ItemRotaSort=ItemRotaSort;ItemLWSort=ItemLWSort;
%         Strip_Item=Strip_Item;sItem_Strip=sItem_Strip;
        
        % 1 ����Item���Sort����
        %  1.1 ����CoordItemStripSort
        sCoordItemStrip(1,iItem) = wStrip - Strip.LW(1,thisLevel);  %����x����
        sCoordItemStrip(2,iItem) = sum(Strip.LW(2,1:thisLevel-1));      %����y���� %���iLevel=1,�����ߣ�����Ϊ0������Ϊ���
        
        % 2 ����Strip������ݣ�������
        %  2.1 ����LWStrip        
        if sItem.isRota(iItem) == 1 %��Item������ת
            % �ж����: 2��
            isflagCurr = Strip.LW(1,thisLevel) >=  sItem.LWH(1,iItem); %�ж��Ƿ�current's stripʣ���� >= ��ǰ�߶ȣ�����ת��
            isNewLevel = Strip.LW(1,thisLevel) == wStrip; % �ж��Ƿ� new Level            
            % ����strip��Ϣ
            if isNewLevel %�������,�����Է���,���ۺ��ְڷ�,���:ֱ�Ӹ��£�Item�Ѱ�Hori/Vert�ڷŹ���
                    updateLWStrip();
            else % �������level ����԰ڷ�,����; ����,����������ת�������
                    if isflagCurr
                        updateLWStrip();
                    else
                        rotateItem();
                        updateLWStrip();
                    end
             end            
        elseif sItem.isRota(iItem) == 0 %��Item��������ת
            % ����strip��Ϣ
            updateLWStrip();
        end
        %  2.2 ����stripBeItemArray
        Strip_Item(thisLevel) = Strip_Item(thisLevel) + 1; %ֻҪ��level����һ��item,����������1
        
        %  2.3 ���±�level��Ӧ��StripWeight: 
        Strip.Weight(thisLevel) =  Strip.Weight(thisLevel) + sItem.Weight(iItem);
        
        %  1.3 ����item����strip��ϢitemBeStripMatrixSort
        sItem_Strip(1,iItem) = thisLevel;    %�ڼ���level
        sItem_Strip(2,iItem) = Strip_Item(thisLevel); %��level�µڼ��ΰ���
        
        % 4 ����Ƕ�׺���
        function rotateItem()
            %  �������Rotaed�仯 ��Ҫ����Ʒ������rotate(��)��ȥ
            sItem.Rotaed(iItem) = ~sItem.Rotaed(iItem);
            tep = sItem.LWH(1,iItem);
            sItem.LWH(1,iItem) = sItem.LWH(2,iItem);
            sItem.LWH(2,iItem) = tep;
        end
        
        function updateLWStrip()
            Strip.LW(1,thisLevel) = Strip.LW(1,thisLevel) - sItem.LWH(1,iItem); %����wleft (�ڷŷ���ǰ��һ��)
            Strip.LW(2,thisLevel) = max(Strip.LW(2,thisLevel), sItem.LWH(2,iItem)); %����strip�߶�lleft(ȡ���ֵ)
        end
        
    end
    
    function printscript()
        % ���Դ���
        % % LWHStrip
        % % sItem_Strip
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
        hStrip = sum(Strip.LW(2,sItem_Strip(2,:)>0));        
        nstrip = sum(sItem_Strip(2,:)>0);

        nIDType = unique(sItem.ID);
        nColors = hsv(length(nIDType)); %��ͬ����LU���費ͬ��ɫ
        
        %% ��ͼ
        % 1 ��ͼ��������Strip
        DrawRectangle([wStrip/2 hStrip/2 wStrip hStrip 0],'--', [0.5 0.5 0.5]);
        hold on;
        % 2 ��ͼ�����strip/item ��ͼ
        for istrip = 1:nstrip
            % �ҳ���ǰistrip����Ʒ����
            idxDrawItem = find(sItem_Strip(1,:)==istrip);
            % ��ȡ�������µı���
            drawItemCoordMatrix = sCoordItemStrip(:,idxDrawItem);
            drawItemLWH = sItem.LWH(:,idxDrawItem);
            drawItemId = sItem.ID(:,idxDrawItem);

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

function order = getITEMorder(Item,whichSortItemOrder)        
        tmpLWHItem = Item.LWH(1:2,:);
        tmpIDItem =     Item.ID(1,:);
        if whichSortItemOrder == 1 %Descend of ��(��)
            tmpLWH = [tmpIDItem; tmpLWHItem]; %��������ITEM��ID����һ���γ���ʱ����
            [~,order] = sortrows(tmpLWH',[3 1 ],{'descend','descend'}); %���߶�,ID(��ͬ�߶�ʱ)�ݼ�����            
        end
        if whichSortItemOrder == 2  %Descend of shortest��̱�  -> ����Rotation���ӱ���
%             tmpLWH = [tmpLWHItem; min(tmpLWHItem(1:nDim,:))]; %����������̱ߵ��������γ���ʱ����tmpLWH
%             [~,order] = sortrows(tmpLWH',[3 2 1],{'descend','descend','descend'}); %way1 ����̱�,�߶�,��ȵݼ�����
             %BACKUP  [~,itemorder] = sort(tmpLWH(nDim+1,:),'descend'); %��ȡ��ʱ����������˳�� way2
        end
        if whichSortItemOrder == 3  %Descend of total area �ܱ����  ->
            %             printstruct(d);
       end        
        if ~isrow(order), order=order'; end
end

    function [thisLevel,iLevel] = getThisLevel( iItem, iLevel, sItem,  Strip, p)  
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
                [thisLevel,iLevel] = getThisLevel( iItem, iLevel, sItem,  Strip, p);
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
                 [thisLevel,iLevel] = getThisLevel( iItem, iLevel, sItem, Strip, p); 
            else
                thisLevel = flag(1);
                if ~all(ismember(thisLevel,flag)),     error('Not all thisLevel belongs to flag ');       end
            end
        elseif p.whichStripH == 3 % nextfit                 
            % ����Rotation���ӱ���
            if sItem.isRota(iItem) == 1 %��Item������ת % nextfit�²���ֱ������bestfit�Ĵ���
                % �ж���ǰlevel�Ƿ��������һ�ڷŷ���ɷ����iItem flaged: �����ݱ�ʾ���ԣ����򲻿���
                flaged = find(Strip.LW(1,iLevel) >= sItem.LWH(1,iItem) |  ...
                                      Strip.LW(1,iLevel) >= sItem.LWH(2,iItem));
            else
                % ��ͬ�����µ�ѡ�������ǰitem�Ŀ�<=��ǰstrip�ĵ�ǰlevel�Ŀ�
                flaged = find(Strip.LW(1,iLevel) >= sItem.LWH(1,iItem) );
            end
            if  isempty(flaged)  %ע����֮ǰ~flag������
                iLevel = iLevel + 1;% �����Ȳ����㣬��level����
                [thisLevel,iLevel] = getThisLevel( iItem, iLevel, sItem, Strip, p) ;
            else
                if  isempty(flaged) ,   error(' �����ܵĴ��� ');      end
                thisLevel = iLevel; % ��ǰlevelһ���ŵ���
            end
%             end
        
        end
    end


%% DEL ����

% %     function insertItemToStrip(thisLevel)
% % %         CoordItemStripSort=CoordItemStripSort;LWStrip=LWStrip;
% % %         %Ϊ��matlab��coder ����
% % %         ItemRotaSort=ItemRotaSort;ItemLWSort=ItemLWSort;
% % %         stripBeItemArray=stripBeItemArray;sItem_Strip=sItem_Strip;
% %         
% %         % 1 ����CoordItemStripSort
% %         CoordItemStripSort(1,iStrip) = widthStrip - LWStrip(1,thisLevel);  %����x����
% %         CoordItemStripSort(2,iStrip) = sum(LWStrip(2,1:thisLevel-1));      %����y���� %���iLevel=1,�����ߣ�����Ϊ0������Ϊ���
% %         
% %         % 2 ����LWStrip
% %         % 2 ����stripBeItemArray
% %         if ParaArray.whichRotation == 1
% %             % �ж����: 5��
% %             isflagHori = LWStrip(1,thisLevel) >=  ItemLWSort(1,iStrip); %�ж��Ƿ� wleft >= longest
% %             isflagVert = LWStrip(1,thisLevel) >= ItemLWSort(2,iStrip);  %�ж��Ƿ� wleft >= shortest
% %             isNewLevel = LWStrip(1,thisLevel) == widthStrip; % �ж��Ƿ� new Level
% %             
% %             % ����strip��Ϣ
% %             if isNewLevel
% %                 if ParaArray.whichRotationHori == 2 % ��level, �ض�����vertical��ʽ(����)
% %                     if ~isflagVert,  error('1'); end
% %                     updateStripVertical();
% %                 else % 0-1 ��% ��level, �ض�����horizontally��ʽ(��Ĭ�Ϸ�ʽ����)
% %                     if ~isflagHori,  error('1'); end
% %                     updateStripHorizontal();
% %                 end
% %             else % �������level
% %                 if ParaArray.whichRotationHori == 1 % ����level, ���Ȱ���horizontall��ʽ(����)
% %                     if isflagHori
% %                         updateStripHorizontal();
% %                     elseif isflagVert
% %                         updateStripVertical();
% %                     else
% %                         error('levelѡ�����,���ۺ������Ų���');
% %                     end
% %                 else % 0/2 : ����level, �ض�����vertical��ʽ(����) vһ�����ԣ�h��һ������ (ע��v���ԣ���h�����ԵĿ����ԣ����ʼ��δhotizntal orientation�Ͳ�����)
% %                     if ~isflagVert,  error('1');  end
% %                     updateStripVertical();
% %                 end
% %             end
% %         elseif ParaArray.whichRotation == 0
% %             % ����strip��Ϣ
% %             LWStrip(1,thisLevel) = LWStrip(1,thisLevel) - ItemLWSort(1,iStrip); %����wleft
% %             LWStrip(2,thisLevel) = max(LWStrip(2,thisLevel),ItemLWSort(2,iStrip)); % �����Ʒ�߶�>��strip�߶�->����strip�߶�lleft
% %         end
% %         
% %         % 3 ����item����strip��ϢitemBeStripMatrixSort + ����stripBeItemArray
% %         stripBeItemArray(thisLevel) = stripBeItemArray(thisLevel) + 1; %ֻҪ��level����һ��item,����������1
% %         itemBeStripMatrixSort(1,iStrip) = thisLevel;    %�ڼ���level
% %         itemBeStripMatrixSort(2,iStrip) = stripBeItemArray(thisLevel); %��level�µڼ��ΰ���
% %         
% %         % 4 ����Ƕ�׺���
% %         function updateStripVertical()
% %             LWStrip(1,thisLevel) = LWStrip(1,thisLevel) - ItemLWSort(2,iStrip); %����wleft by Vertically
% %             LWStrip(2,thisLevel) = max(LWStrip(2,thisLevel), ItemLWSort(1,iStrip)); %����strip�߶�lleft(ȡ���ֵ)
% %             % ��isflagVert==1 ������� ��Ҫ����Ʒ������rotation(��)��ȥ
% %             ItemRotaSort(iStrip) = ~ItemRotaSort(iStrip);
% %             tep = ItemLWSort(1,iStrip);
% %             ItemLWSort(1,iStrip) = ItemLWSort(2,iStrip);
% %             ItemLWSort(2,iStrip) = tep;
% %         end
% %         function updateStripHorizontal()
% %             LWStrip(1,thisLevel) = LWStrip(1,thisLevel) - ItemLWSort(1,iStrip); %����wleft by Horizontally
% %             LWStrip(2,thisLevel) = max(LWStrip(2,thisLevel), ItemLWSort(2,iStrip)); %����strip�߶�lleft(ȡ���ֵ)
% %         end
% %             
% %     end
