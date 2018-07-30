function [d] = HItemToStrip(d,ParaArray)
% ��Ҫ����:Item����Strip�� %  ����:�����(row);  ����:��������(coloum);
% Input ---  ITEM:  ID LWH Weight
% Output --- ITEM: itemorder itemBeStripMatrix itemRotaFlag CoordItemStrip
% Output --- StripArray: LW Weight
% d.ItemArray (1 LWH (��֪)
% d.ItemArray (2 itemBeStripMatrix  (dim1:���item��ĳ��strip dim2:item����˳��(��->��) 
% d.ItemArray (3 CoordItemStrip Item��strip������) 
% d.ItemArray (4 itemorder ������Item����˳��)
% d.StripArray (1 LW )

%% ��ʼ��
% nDim Itemά��(2) nItem Item���� nStrip Strip���� 
% widthStrip Strip�����
nDim = size(d.ItemArray.LWH,1);  if nDim ==3, nDim = nDim-1;end
nItem = size(d.ItemArray.LWH,2);
nStrip = nItem;

tmpUniqueBin = unique(d.BinArray.LWH(1:nDim,:)','rows')';
widthStrip = tmpUniqueBin(1);
clear tmpUniqueBin;

%% ���ж�ITEM����Horizontal/Vertical ��ʽ�ڷţ������Ƿ���ת�������жϽ����㷨��˳��
% �����Ƿ�������ת, ֻ���Ƿ���Ҫ��Horizontal/Vertical��ʽ�ڷ�

%  [ItemLWRota, ItemRotaed] = placeItemHori(d.ItemArray,1);  %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��
% ���ԭ�ⲻ���ķ���ֵ:����ʼֵ
[d.ItemArray.Rotaed] = placeItemHori(d.ItemArray.LWH,d.ItemArray.isRota,2);  %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��
d.ItemArray.LWH = getRotaedLWH(d.ItemArray.LWH, d.ItemArray.Rotaed, d.LUArray.BUFF); 

    %% ITEM���� 555
    % getITEMorder - ��ȡItemLWRota��˳��(�ص��Ǹ߶ȵݼ�����) % ITEM��������ʽ �߶�/��̱�
    d.ItemArray.itemorder = getITEMorder();
    % getSortedITEM - ��ȡ��order������ITEM:sortedItemArray
    sortedItemArray = getSortedITEM(d.ItemArray.itemorder);
                                    % printstruct(d) ;printstruct(sortedItemArray)

    % 1��2���ڵ�, sortedItemArray��Ӧ��LWHRota��Rotaed������->���󷵻ص�ԭ����ItemArry��
    if ParaArray.whichRotationHori == 1 % �����ĸ�level,������horizontally��ʽ�ڷ�
        [ sortedItemArray.Rotaed] = placeItemHori(sortedItemArray.LWH,sortedItemArray.isRota,1);  %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��
    end
    if ParaArray.whichRotationHori == 2 % �����ĸ�level,������vertical��ʽ�ڷ�
        [ sortedItemArray.Rotaed] = placeItemHori(sortedItemArray.LWH,sortedItemArray.isRota,0);  %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��
    end
     sortedItemArray.LWH = getRotaedLWH(sortedItemArray.LWH, sortedItemArray.Rotaed, d.LUArray.BUFF); 
     
%% 55 LU->Item->Stripת�� 
% ��ȡ���� �˴�ֻʹ��LWHRota��Rotaed; ��ʹ��LWH
ItemLWRotaSort = sortedItemArray.LWH(1:nDim,:); %ItemLWSortHori
ItemRotaedSort = sortedItemArray.Rotaed; % ItemRotaSortHori % itemRotaSort = zeros(1,size(d.ItemArray.LWH,2));
ItemisRotaSort = sortedItemArray.isRota;
Itemorder = d.ItemArray.itemorder;
                    % ��ȡ1 ItemLWSort ������ITEM�ĳ���� 2 itemRotaSort ������FLAG ȫ��Ϊ0 ����Ҫ
                    % if ParaArray.whichRotation == 0
                    %     ItemLWSort = sortedItemArray.LWH(1:nDim,:); 
                    %     itemRotaSort = zeros(1,size(d.ItemArray.LWH,2));
                    % end %itemRotaSort = sortedItemArray.itemRotaFlag; 
                    % ��ȡ1 ItemLWSort Hori orienting��LWH ������ƽ����bottom��
                    % 2 ItemRotaSort ԭʼitem�Ƿ�rotation  0 ����horizontal orientation 1 ���� vertical orientation
                    % if ParaArray.whichRotation == 1,
                    %     [ItemLWSort, ItemRotaSort] = horiOrient(ItemLWSort); 
                    % end
%%
% ��ȡitemBeStripMatrixSort: ÿ�������Item���ĸ�Strip��  �Լ�˳��
% ��ȡLWStrip:  �����ɵ�Strip�ĳ���
% ��ȡCoordItemStripSort  Item��strip������ֵ
LWStrip = zeros(nDim,nItem);   %strip���� dim2-����(����ߵļ���) (�߶Ƚ����ο�,current �߶�)
LWStrip(1,:) = widthStrip;   %dim1-���ʣ�� 
stripBeItemArray = zeros(1,nStrip);  % ÿ��Strip�ڵ�Item���� ���ڲ���
itemBeStripMatrixSort = zeros(2,nItem); %dim1:���ڵڼ���level dim2:���ڸ�level�ڼ����ŷ� 555
CoordItemStripSort = zeros(2,nItem); %Item��strip������ֵ

% 55 ��ȡthisLevel - ��ǰitemҪ�����level���
% ѭ����strip�а���item,���̶�item,�仯ѡ��ͬlevel(thisLevel)
% ע�ͣ���ȡ FLAG        �ɷ��µ�ǰitem(iItem)������һ��level�ļ��ϣ� 
% ע�ͣ���ȡ thisLevel   ��FLAG���ҵ���������Ǹ�thisLevel, ��ִ�� insert����

iLevel = 1; iItem = 1; %iStrip����itemʵ��
while 1
    if iItem > nItem, break; end
    
    % ���ݲ�ͬ�����ҵ��ܷ��뵱ǰitem��strips/levels�е�һ��    

    thisLevel = getThisLevel();     %iLevel���ڴκ����ڲ��ϵ�������Զָʾ��ǰ���µ�level
    insertItemToStrip(thisLevel);
    
%      plot2DStrip(); %������ͼ
    
    iItem = iItem + 1;
end

% plot2DStrip(); %һ���Ի�ͼ

% ���� ����ֵ��d
%Matalb code gerator use:
%         itemBeStripMatrix=itemBeStripMatrixSort; CoordItemStrip=CoordItemStripSort;

% Item��أ����µİ�˳�򷵻أ��޸��µĲ��践�أ�
% itemBeStripMatrixSort
% CoordItemStripSort
% ItemRotaedSort
% ItemLWRotaSort
% ��ȡitemBeStripMatrix : ÿ��Item���ĸ�Strip��  �Լ�˳��
% ��ȡCoordItemStrip : ÿ��Item��Strip������
% ��ȡRotaed : ÿ��item�Ƿ�Rotation�ı�־
% ��ȡLWHRota��ÿ��item���Rotaed��־���õ�LWH

    itemBeStripMatrix(:,Itemorder) = itemBeStripMatrixSort;
    d.ItemArray.itemBeStripMatrix = itemBeStripMatrix;
    
    CoordItemStrip(:,Itemorder) = CoordItemStripSort;    
    d.ItemArray.CoordItemStrip = CoordItemStrip;
    
    % ItemArray��ת���
    itemRotaed(:,Itemorder) = ItemRotaedSort;
    d.ItemArray.Rotaed = itemRotaed;
    ItemLWRota(:,Itemorder) = ItemLWRotaSort;
    d.ItemArray.LWH = [ItemLWRota; d.ItemArray.LWH(3,:)];  % ����ԭʼ˳�����ת���ItemArray

    % LUArray��ת���,��ʱ����    
    nbItem=length(d.ItemArray.Rotaed);
    % ѭ��ÿ��item
    for idxItem=1:nbItem
        flagThisItem = (d.LUArray.LUBeItemArray(1,:)==idxItem );
        % ��Ӧλ��LU.Rotaed����
        if d.ItemArray.Rotaed(idxItem)
            d.LUArray.Rotaed(flagThisItem) = ~d.LUArray.Rotaed(flagThisItem);
            % ��Ӧλ��LU.LWH����
            d.LUArray.LWH(1, flagThisItem) = d.ItemArray.LWH(1, idxItem);
            d.LUArray.LWH(2, flagThisItem) = d.ItemArray.LWH(2, idxItem);
        end
    end

    
% Strip���: ��˳�����
% LWStrip
% ��ȡLWStrip:  �����ɵ�strip�ĳ���

    d.StripArray.LW = LWStrip(:,LWStrip(2,:)>0); % û��˳�� + ȥ��δʹ�õ�Strip
    
                % % % ��ʹû��RotationFlage Ҳ�и����� �ж��Ƿ�Rotation TO BE DELE
                %     if ParaArray.whichRotation == 0,   itemRotaFlag(:,Itemorder) = itemRotaSort;  end
                %     if ParaArray.whichRotation == 1,   itemRotaFlag(:,Itemorder) = ItemRotaedSort;   end
                %     d.ItemArray.itemRotaFlag = itemRotaFlag;

    
    %% ����script
    % �����Ҫ���:���ÿ��level������ 
    printscript();
%     dap = rmfield(d, {'BinArray', 'LUArray'});  printstruct(dap);


    %% Ƕ�׺���        
    function order = getITEMorder()        
        tmpLWHItem = d.ItemArray.LWH(1:nDim,:);
        tmpIDItem =     d.ItemArray.ID(1,:);
        if ParaArray.whichSortItemOrder == 1 %Descend of ��(��)
            tmpLWH = [tmpIDItem; tmpLWHItem]; %��������ITEM��ID����һ���γ���ʱ����
            [~,order] = sortrows(tmpLWH',[3 1 ],{'descend','descend'}); %���߶�,ID(��ͬ�߶�ʱ)�ݼ�����            
        end
        if ParaArray.whichSortItemOrder == 2  %Descend of shortest��̱�  -> ����Rotation���ӱ���
%             tmpLWH = [tmpLWHItem; min(tmpLWHItem(1:nDim,:))]; %����������̱ߵ��������γ���ʱ����tmpLWH
%             [~,order] = sortrows(tmpLWH',[3 2 1],{'descend','descend','descend'}); %way1 ����̱�,�߶�,��ȵݼ�����
             %BACKUP  [~,itemorder] = sort(tmpLWH(nDim+1,:),'descend'); %��ȡ��ʱ����������˳�� way2
        end
        if ParaArray.whichSortItemOrder == 3  %Descend of total area �ܱ����  ->
            %             printstruct(d);
       end        
        if ~isrow(order), order=order'; end
    end

% %     function order = getITEMorder()
% %         ItemLWRota
% %         tmpLWHItem = d.ItemArray.LWH(1:nDim,:);
% %         tmpIDItem =     d.ItemArray.ID(1,:);
% %         if ParaArray.whichSortItemOrder == 1 %Descend of ��(��)
% %             tmpLWH = [tmpIDItem; tmpLWHItem]; %��������ITEM��ID����һ���γ���ʱ����
% %             [~,order] = sortrows(tmpLWH',[3 1 ],{'descend','descend'}); %���߶�,ID(��ͬ�߶�ʱ)�ݼ�����            
% %         end
% %         if ParaArray.whichSortItemOrder == 2  %Descend of shortest��̱�  -> ����Rotation���ӱ���
% %             tmpLWH = [tmpLWHItem; min(tmpLWHItem(1:nDim,:))]; %����������̱ߵ��������γ���ʱ����tmpLWH
% %             [~,order] = sortrows(tmpLWH',[3 2 1],{'descend','descend','descend'}); %way1 ����̱�,�߶�,��ȵݼ�����
% %              %BACKUP  [~,itemorder] = sort(tmpLWH(nDim+1,:),'descend'); %��ȡ��ʱ����������˳�� way2
% %         end
% %         if ParaArray.whichSortItemOrder == 3  %Descend of total area �ܱ����  ->
% % %             printstruct(d);
% %        end        
% %         if ~isrow(order), order=order'; end
% %     end

    function item = getSortedITEM(order)
        item = structfun(@(x) x(:,order),d.ItemArray,'UniformOutput',false);
    end

    function thisLevel = getThisLevel()
        % ��ͬwhichStripH��,��ù�ͬ��thisLevel
        if ParaArray.whichStripH == 1 % 1 bestfit 2 firstfit 3 nextfit
            % ����Rotation���ӱ���
            if ItemisRotaSort(iItem) == 1 %��Item������ת
                                        %             if ParaArray.whichRotation == 1
                                        % �ҵ�����rotation�µ�level:��һ�ڷŷ���ɷ����iItem��level
                flag = find(LWStrip(1, 1 : iLevel) >= ItemLWRotaSort(1,iItem) |  ...
                                  LWStrip(1, 1 : iLevel) >= ItemLWRotaSort(2,iItem));
            else %��Item��������ת
                % ���������µ�ѡ��find����㹻�Ķ��level,����������Сʣ��ˮƽ��ȵ�
                flag = find(LWStrip(1, 1 : iLevel) >= ItemLWRotaSort(1,iItem));
            end
            if isempty(flag)
                iLevel = iLevel + 1;% �����Ȳ����㣬��level����
                thisLevel = getThisLevel();
            else
                % ��ȡthisLevel: Ψһ��FF������⵽thisLevel�ļ��㣨ѡ��������������С�ģ�
                tmpLevels = LWStrip(1,1:iLevel);   %��ȡ�����Ѱ��Ż��°��ŵ�level��ʣ��ˮƽƽ�������tepAvailableLevelArray
                tepMinLeftWdith = min(tmpLevels(flag));                      %�ҳ�tepAvailableLevelArray�п����ɱ�iITem��ʣ��ˮƽ��ȵ���Сֵ������ΪtepMinLeftWdith
                thisLevel = find(tmpLevels==tepMinLeftWdith);            %�ҳ�����Сֵ��Ӧ���Ǹ�/Щlevel
                if ~all(ismember(thisLevel,flag)),     error('Not all thisLevel belongs to flag ');          end
                if length(thisLevel)>1
                    thisLevel = thisLevel(1);
                end
            end
        elseif ParaArray.whichStripH == 2 % firstfit  % firstfit���ܲ���ֱ������bestfit�Ĵ���?
            % ����Rotation���ӱ���
            if ItemisRotaSort(iItem) == 1 %��Item������ת
                % �ҵ�����rotation�µ�level:��һ�ڷŷ���ɷ����iItem��level
                flag = find(LWStrip(1, 1 : iLevel) >= ItemLWRotaSort(1,iItem) |  ...
                                  LWStrip(1, 1 : iLevel) >= ItemLWRotaSort(2,iItem));
            else
                % ���������µ�ѡ��find����㹻�Ķ��level,�������ڵ�һ�������� Ψһ������thisLevel�Ļ�ȡ
                flag = find(LWStrip(1, 1 : iLevel) >= ItemLWRotaSort(1,iItem));
            end
            if isempty(flag)
                iLevel = iLevel + 1;% �����Ȳ����㣬��level����
                 thisLevel = getThisLevel();
            else
                thisLevel = flag(1);
                if ~all(ismember(thisLevel,flag)),     error('Not all thisLevel belongs to flag ');       end
            end
        elseif ParaArray.whichStripH == 3 % nextfit 
            % ����Rotation���ӱ���
            if ItemisRotaSort(iItem) == 1 %��Item������ת % nextfit�²���ֱ������bestfit�Ĵ���
                % �ж���ǰlevel�Ƿ��������һ�ڷŷ���ɷ����iItem flaged: �����ݱ�ʾ���ԣ����򲻿���
                flaged = find(LWStrip(1,iLevel) >= ItemLWRotaSort(1,iItem) |  ...
                                      LWStrip(1,iLevel) >= ItemLWRotaSort(2,iItem));
            else
                % ��ͬ�����µ�ѡ�������ǰitem�Ŀ�<=��ǰstrip�ĵ�ǰlevel�Ŀ�
                flaged = find(LWStrip(1,iLevel) >= ItemLWRotaSort(1,iItem) );
            end
            if  isempty(flaged)  %ע����֮ǰ~flag������
                iLevel = iLevel + 1;% �����Ȳ����㣬��level����
                thisLevel = getThisLevel();
            else
                if  isempty(flaged) ,   error(' �����ܵĴ��� ');      end
                thisLevel = iLevel; % ��ǰlevelһ���ŵ���
            end
        end
    end

    function insertItemToStrip(thisLevel)
%         CoordItemStripSort=CoordItemStripSort;LWStrip=LWStrip;
%         %Ϊ��matlab��coder ����
%         ItemRotaSort=ItemRotaSort;ItemLWSort=ItemLWSort;
%         stripBeItemArray=stripBeItemArray;itemBeStripMatrixSort=itemBeStripMatrixSort;
        
        % 1 ����CoordItemStripSort
        CoordItemStripSort(1,iItem) = widthStrip - LWStrip(1,thisLevel);  %����x����
        CoordItemStripSort(2,iItem) = sum(LWStrip(2,1:thisLevel-1));      %����y���� %���iLevel=1,�����ߣ�����Ϊ0������Ϊ���
        
        % 2 ����LWStrip
        % 2 ����stripBeItemArray
        if ItemisRotaSort(iItem) == 1 %��Item������ת
            % �ж����: 2��
                                    %             isflagHori = LWStrip(1,thisLevel) >=  ItemLWSort(1,iItem); %�ж��Ƿ� wleft >= longest
                                    %             isflagVert = LWStrip(1,thisLevel) >= ItemLWSort(2,iItem);  %�ж��Ƿ� wleft >= shortest
            isflagCurr = LWStrip(1,thisLevel) >=  ItemLWRotaSort(1,iItem); %�ж��Ƿ�current's stripʣ���� >= ��ǰ�߶�
            isNewLevel = LWStrip(1,thisLevel) == widthStrip; % �ж��Ƿ� new Level            
            % ����strip��Ϣ
            if isNewLevel %�������,�����Է���,���ۺ��ְڷ�,���:ֱ�Ӹ��£�Item�Ѱ�Hori/Vert�ڷŹ���
%                 if ParaArray.whichRotationHori == 1 % ��level, �ض�����horizontally��ʽ�ڷ�
%                     updateLWStrip();               
%                 elseif ParaArray.whichRotationHori == 2 % ��level, �ض�����vertical��ʽ�ڷ�
                    updateLWStrip();
            else % �������level ����԰ڷ�,����; ����,������С����
%                 if ParaArray.whichRotationHori == 1 % ����level, ���Ȱ���horizontall��ʽ(����)
                    if isflagCurr
                        updateLWStrip();
                    else
                        rotateItem();
                        updateLWStrip();
                    end
    %                 elseif ParaArray.whichRotationHori == 2 % 2 : ����level, ���Ȱ���vertical��ʽ(����)
                        %vһ�����ԣ�h��һ������ (ע��v���ԣ���h�����ԵĿ����ԣ����ʼ��δhotizntal orientation�Ͳ�����)
                end            
        elseif ItemisRotaSort(iItem) == 0 %��Item��������ת
            % ����strip��Ϣ
            updateLWStrip();
        end
        
        % 3 ����item����strip��ϢitemBeStripMatrixSort + ����stripBeItemArray
        stripBeItemArray(thisLevel) = stripBeItemArray(thisLevel) + 1; %ֻҪ��level����һ��item,����������1
        itemBeStripMatrixSort(1,iItem) = thisLevel;    %�ڼ���level
        itemBeStripMatrixSort(2,iItem) = stripBeItemArray(thisLevel); %��level�µڼ��ΰ���
        
        % 4 ����Ƕ�׺���
        function rotateItem()
            %  �������Rotaed�仯 ��Ҫ����Ʒ������rotate(��)��ȥ
            ItemRotaedSort(iItem) = ~ItemRotaedSort(iItem);
            tep = ItemLWRotaSort(1,iItem);
            ItemLWRotaSort(1,iItem) = ItemLWRotaSort(2,iItem);
            ItemLWRotaSort(2,iItem) = tep;
        end
        
        function updateLWStrip()
            LWStrip(1,thisLevel) = LWStrip(1,thisLevel) - ItemLWRotaSort(1,iItem); %����wleft (�ڷŷ���ǰ��һ��)
            LWStrip(2,thisLevel) = max(LWStrip(2,thisLevel), ItemLWRotaSort(2,iItem)); %����strip�߶�lleft(ȡ���ֵ)
        end
        
    end
    
    function printscript()
        % ���Դ���
        % % LWHStrip
        % % itemBeStripMatrixSort
        % % ItemLWSort
        % % itemCoordMatrixSort
        % % itemBeStripMatrix
        % % LWHItem
        % % itemCoordMatrix
        %  printstruct(d);
        
        % �����Ҫ���:��ô�1��ʼÿ��strip����������
        for iStrip = 1:max(d.ItemArray.itemBeStripMatrix(1,:))
            [~,idx] = find(d.ItemArray.itemBeStripMatrix(1,:)==iStrip);
            fprintf('strip %d ��ʣ���+���Ϊ:  ',iStrip);
            fprintf('( %d ) ',d.StripArray.LW(:,iStrip));
            fprintf('\n');
            fprintf('strip %d ���� original Item ������(����)[��ת��־]{����}Ϊ  \n  ',iStrip);
            fprintf('%d ',idx);
            fprintf('( %d ) ', d.ItemArray.LWH(1:nDim,idx));fprintf('\n');
            fprintf('[ %d ] ', d.ItemArray.Rotaed(:,idx));fprintf('\n');  %ItemRotaedSort
            fprintf('{ %d } ', d.ItemArray.CoordItemStrip(:,idx));fprintf('\n');
            fprintf('\n');
        end
    end


    function plot2DStrip()
        %% ��ʼ��
        wStrip = widthStrip;        
        hStrip = sum(LWStrip(2,itemBeStripMatrixSort(2,:)>0));        
        nstrip = sum(itemBeStripMatrixSort(2,:)>0);

        nIDType = unique(sortedItemArray.ID);
        nColors = hsv(length(nIDType)); %��ͬ����LU���費ͬ��ɫ
        
        %% ��ͼ
        % 1 ��ͼ��������Strip
        DrawRectangle([wStrip/2 hStrip/2 wStrip hStrip 0],'--', [0.5 0.5 0.5]);
        hold on;
        % 2 ��ͼ�����strip/item ��ͼ
        for istrip = 1:nstrip
            % �ҳ���ǰistrip����Ʒ����
            idxDrawItem = find(itemBeStripMatrixSort(1,:)==istrip);
            % ��ȡ�������µı���
            drawItemCoordMatrix = CoordItemStripSort(:,idxDrawItem);
            drawItemLWH = ItemLWRotaSort(:,idxDrawItem);
            drawItemId = sortedItemArray.ID(:,idxDrawItem);

            % ��ͼ�����item
            nThisItem = size(drawItemLWH,2);
            for iplotItem = 1:nThisItem
                % ��ͼ��������iItem
                itemWidth = drawItemLWH(1,iplotItem);
                itemLength = drawItemLWH(2,iplotItem);
                itemCenter = [drawItemCoordMatrix(1,iplotItem)+itemWidth/2 ...
                    drawItemCoordMatrix(2,iplotItem)+itemLength/2 ];
                
                % ���Ӷ�rotation���ж� �ڲ�Ƕ�׺�������Ҫ
                %             if ParaArray.whichRotation == 1
                %                 drawItemRotaMatrix = ItemRotaSort(:,idxDrawItem); %����rotation����
                %             end
                
                % %         if ParaArray.whichRotation == 1 && drawItemRotaMatrix(iItem)
                % %             itemWidth = drawItemLWH(2,iItem);
                % %             itemLength = drawItemLWH(1,iItem);
                % %             itemCenter = [drawItemCoordMatrix(1,iItem)+itemWidth/2 ...
                % %                         drawItemCoordMatrix(2,iItem)+itemLength/2 ];
                % %         end
                
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

%% DEL ����

% %     function insertItemToStrip(thisLevel)
% % %         CoordItemStripSort=CoordItemStripSort;LWStrip=LWStrip;
% % %         %Ϊ��matlab��coder ����
% % %         ItemRotaSort=ItemRotaSort;ItemLWSort=ItemLWSort;
% % %         stripBeItemArray=stripBeItemArray;itemBeStripMatrixSort=itemBeStripMatrixSort;
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
