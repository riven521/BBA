function [da] = HItemToStrip(da,ParaArray)
% ��Ҫ����:Item����Strip�� %  ����:�����(row);  ����:��������(coloum);
% Input ---  ITEM:  ID LWH Weight
% Output --- ITEM: itemorder itemBeStripMatrix itemRotaFlag CoordItemStrip
% Output --- StripArray: LW Weight
% da.ItemArray (1 LWH (��֪)
% da.ItemArray (2 itemBeStripMatrix  (dim1:���item��ĳ��strip dim2:item����˳��(��->��) 
% da.ItemArray (3 CoordItemStrip Item��strip������) 
% da.ItemArray (4 itemorder ������Item����˳��)
% da.StripArray (1 LW )

%% ��ʼ��
% nDim Itemά��(2) nItem Item���� nStrip Strip���� 
% widthStrip Strip�����
nDim = size(da.ItemArray.LWH,1);  if nDim ==3, nDim = nDim-1;end
nItem = size(da.ItemArray.LWH,2);
nStrip = nItem;

tmpUniqueBin = unique(da.BinArray.LWH(1:nDim,:)','rows')';
widthStrip = tmpUniqueBin(1);
clear tmpUniqueBin;


    %% ITEM���� 555 
    % getITEMorder - ��ȡITEM��˳��(�ص��Ǹ߶ȵݼ�����) % ITEM��������ʽ �߶�/��̱�
    da.ItemArray.itemorder = getITEMorder(); 
    % getSortedITEM - ��ȡ��order������ITEM:sortedItemArray
    sortedItemArray = getSortedITEM(da.ItemArray.itemorder);
%     printstruct(da) ;printstruct(sortedItemArray)    

%% ����Rotation���ӱ��� 
itemorder = da.ItemArray.itemorder;
% ��ȡ1 LWHItemSort ������ITEM�ĳ���� 2 itemRotaSort ������FLAG ȫ��Ϊ0 ����Ҫ
if ParaArray.whichRotation == 0, LWHItemSort = sortedItemArray.LWH(1:nDim,:);; itemRotaSort = sortedItemArray.itemRotaFlag; end
% ��ȡ1 LWHItemSortHori Hori orienting��LWH 2 itemRotaSortHori ԭʼitem�Ƿ�rotation  0 ����horizontal orientation 1 ���� vertical orientation
if ParaArray.whichRotation == 1,     [LWHItemSortHori,itemRotaSortHori] = horiOrient(sortedItemArray.LWH(1:nDim,:)); end


%% 55 LU->Item->Stripת�� 
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
% ע�ͣ���ȡ FLAG        �ɷ��µ�ǰitem(iStrip)������һ��level�ļ��ϣ� 
% ע�ͣ���ȡ thisLevel   ��FLAG���ҵ���������Ǹ�thisLevel, ��ִ�� insert����

iLevel = 1; iStrip = 1;
while 1
    if iStrip > nItem, break; end
    
    thisLevel = getThisLevel();    
    insertItemToStrip(thisLevel);    
    
    iStrip = iStrip + 1;
end

% ���� ����ֵ��da
        %Matalb code gerator use:
        itemBeStripMatrix=itemBeStripMatrixSort;CoordItemStrip=CoordItemStripSort;

% ��ȡitemBeStripMatrix : ÿ��Item���ĸ�Strip��  �Լ�˳��
% ��ȡCoordItemStrip : ÿ��Item��Strip������
% ��ȡitemRotaFlag : ÿ��item�Ƿ�Rotation�ı�־
% ��ȡLWStrip:  �����ɵ�strip�ĳ���
                
    itemBeStripMatrix(:,itemorder) = itemBeStripMatrixSort;
    CoordItemStrip(:,itemorder) = CoordItemStripSort;
    da.ItemArray.itemBeStripMatrix = itemBeStripMatrix;
    da.ItemArray.CoordItemStrip = CoordItemStrip;

    % % % ��ʹû��RotationFlage Ҳ�и����� �ж��Ƿ�Rotation
    if ParaArray.whichRotation == 0,   itemRotaFlag(:,itemorder) = itemRotaSort;  end
    if ParaArray.whichRotation == 1,   itemRotaFlag(:,itemorder) = itemRotaSortHori;   end
    da.ItemArray.itemRotaFlag = itemRotaFlag;

    LWStrip = LWStrip(:,LWStrip(2,:)>0);
    da.StripArray.LW = LWStrip; % ȥ��δʹ�õ�Strip
    
    %% ����script
    % �����Ҫ���:���ÿ��level������ 
    printscript();
    dap = rmfield(da, {'BinArray', 'LUArray'});  printstruct(dap);


    %% Ƕ�׺���        
    function order = getITEMorder()
        tmpLWHItem = da.ItemArray.LWH(1:nDim,:);
        tmpIDItem =     da.ItemArray.ID(1,:);
        if ParaArray.whichSortItemOrder == 1 %Descend of ��(��)
            tmpLWH = [tmpIDItem; tmpLWHItem]; %��������ITEM��ID����һ���γ���ʱ����
            [~,order] = sortrows(tmpLWH',[3 1 ],{'descend','descend'}); %���߶�,ID(��ͬ�߶�ʱ)�ݼ�����            
        end        
        if ParaArray.whichSortItemOrder == 2  %Descend of shortest��̱�  -> ����Rotation���ӱ���
            tmpLWH = [tmpLWHItem; min(tmpLWHItem(1:nDim,:))]; %����������̱ߵ��������γ���ʱ����tmpLWH
            [~,order] = sortrows(tmpLWH',[3 2 1],{'descend','descend','descend'}); %way1 ����̱�,�߶�,��ȵݼ�����
             %BACKUP  [~,itemorder] = sort(tmpLWH(nDim+1,:),'descend'); %��ȡ��ʱ����������˳�� way2
        end
        if ~isrow(order), order=order'; end
    end

    function item = getSortedITEM(order)
        item = structfun(@(x) x(:,order),da.ItemArray,'UniformOutput',false);
    end

    function thisLevel = getThisLevel()
        % ��ͬwhichStripH��,��ù�ͬ��thisLevel
        if ParaArray.whichStripH == 1 % 1 bestfit 2 firstfit 3 nextfit
            % ����Rotation���ӱ���
            if ParaArray.whichRotation == 1
                % �ҵ�����rotation�µ�level:��һ�ڷŷ���ɷ����iItem��level
                flag = find(LWStrip(1,1:iLevel) >= LWHItemSortHori(1,iStrip) |  ...
                    LWStrip(1,1:iLevel) >= LWHItemSortHori(2,iStrip));
            else
                % ���������µ�ѡ��find����㹻�Ķ��level,����������Сʣ���ȵ�
                flag = find(LWStrip(1,1:iLevel) >= LWHItemSort(1,iStrip));
            end
            if isempty(flag)
                iLevel = iLevel + 1;% �����Ȳ����㣬��level����
                thisLevel = getThisLevel();
            else
                % ��ȡthisLevel: Ψһ��FF������⵽thisLevel�ļ��㣨ѡ��������������С�ģ�
                tepAvailableLevelArray = LWStrip(1,1:iLevel);   %��ȡ�����Ѱ��Ż��°��ŵ�level��ʣ��ˮƽƽ�������tepAvailableLevelArray
                tepMinLeftWdith = min(tepAvailableLevelArray(flag));                      %�ҳ�tepAvailableLevelArray�п����ɱ�iITem��ʣ��ˮƽ��ȵ���Сֵ������ΪtepMinLeftWdith
                thisLevel = find(tepAvailableLevelArray==tepMinLeftWdith);            %�ҳ�����Сֵ��Ӧ���Ǹ�(�飩level
                if length(thisLevel)>1
                    thisLevel = thisLevel(1);
                end
            end
        elseif ParaArray.whichStripH == 2 % firstfit  % firstfit���ܲ���ֱ������bestfit�Ĵ���?
            % ����Rotation���ӱ���
            if ParaArray.whichRotation == 1
                % �ҵ�����rotation�µ�level:��һ�ڷŷ���ɷ����iItem��level
                flag = find(LWStrip(1,1:iLevel) >= LWHItemSortHori(1,iStrip) |  ...
                    LWStrip(1,1:iLevel) >= LWHItemSortHori(2,iStrip));
            else
                % ���������µ�ѡ��find����㹻�Ķ��level,�������ڵ�һ�������� Ψһ������thisLevel�Ļ�ȡ
                flag = find(LWStrip(1,1:iLevel) >= LWHItemSort(1,iStrip));
            end
            if isempty(flag)
                iLevel = iLevel + 1;% �����Ȳ����㣬��level����
                 thisLevel = getThisLevel();
            else
                thisLevel = flag(1);
            end
        elseif ParaArray.whichStripH == 3 % nextfit
            % ����Rotation���ӱ���
            if ParaArray.whichRotation == 1 % nextfit�²���ֱ������bestfit�Ĵ��� FIXME
                % �ҵ�����rotation�µ�level:��һ�ڷŷ���ɷ����iItem��level
                flag = find(LWStrip(1,iLevel) >= LWHItemSortHori(1,iStrip) |  ...
                    LWStrip(1,iLevel) >= LWHItemSortHori(2,iStrip));
            else
                % ��ͬ�����µ�ѡ�������ǰitem�Ŀ�<=��ǰstrip�ĵ�ǰlevel�Ŀ�
                flag = find(LWStrip(1,iLevel) >= LWHItemSort(1,iStrip) );
            end
            if  isempty(flag) %ע����֮ǰ~flag������ FIXME ~flag
                iLevel = iLevel + 1;% �����Ȳ����㣬��level����
                 thisLevel = getThisLevel();
            else
                thisLevel = iLevel;
            end
        end
    end

    function insertItemToStrip(thisLevel)
        CoordItemStripSort=CoordItemStripSort;LWStrip=LWStrip;
        itemRotaSortHori=itemRotaSortHori;LWHItemSortHori=LWHItemSortHori;
        stripBeItemArray=stripBeItemArray;itemBeStripMatrixSort=itemBeStripMatrixSort;
        
        % 1 ����CoordItemStripSort
        CoordItemStripSort(1,iStrip) = widthStrip - LWStrip(1,thisLevel);  %����x����
        CoordItemStripSort(2,iStrip) = sum(LWStrip(2,1:thisLevel-1));      %����y���� %���iLevel=1,�����ߣ�����Ϊ0������Ϊ���
        
        % 2 ����stripBeItemArray
        % 2 ����LWStrip
        if ParaArray.whichRotation == 1 && ParaArray.whichRotationHori == 0
            % �ж����
            flagHori = LWStrip(1,thisLevel) >=  LWHItemSortHori(1,iStrip); %�ж��Ƿ� wleft >= longest
            flagVert = LWStrip(1,thisLevel) >= LWHItemSortHori(2,iStrip);  %�ж��Ƿ� wleft >= shortest
            if ~flagVert & ~flagVert
                LWStrip(1,thisLevel)
                LWHItemSortHori(1,iStrip)
                LWHItemSortHori(2,iStrip)
                1
            end
            isNewLevel = LWStrip(1,thisLevel) == widthStrip; % �ж��Ƿ� new Level
            flagLargerLengthHori = LWStrip(2,thisLevel) <  LWHItemSortHori(2,iStrip); %�ж��Ƿ� Ҫ����level��/�߶�
            flagLargerLengthVert = LWStrip(2,thisLevel) <  LWHItemSortHori(1,iStrip); %�ж��Ƿ� Ҫ����level��/�߶�
            
            % ����strip��Ϣ
            % �����level,���Ȱ���horizontally��ʽ(��Ĭ�Ϸ�ʽ����)
            if isNewLevel
                if flagHori
                    LWStrip(1,thisLevel) = LWStrip(1,thisLevel) - LWHItemSortHori(1,iStrip); %����wleft by Horizontally
                    if flagLargerLengthHori
                        LWStrip(2,thisLevel) = LWHItemSortHori(2,iStrip); %����strip�߶�lleft
                    end
                elseif flagVert
                    LWStrip(1,thisLevel) = LWStrip(1,thisLevel) ...
                        - LWHItemSortHori(2,iStrip); %����wleft by Vertically
                    if flagLargerLengthVert
                        LWStrip(2,thisLevel) = LWHItemSortHori(1,iStrip); %����strip�߶�lleft
                    end
                    % ��flagVert==1 ������� ��Ҫ����Ʒ������rotation(��)��ȥ
                    itemRotaSortHori(iStrip) = ~itemRotaSortHori(iStrip);
                    tep = LWHItemSortHori(1,iStrip);
                    LWHItemSortHori(1,iStrip) = LWHItemSortHori(2,iStrip);
                    LWHItemSortHori(2,iStrip) = tep;
                else      %���wleft <shorter edge
                    error('���ݴ���,���ۺ������Ų���');
                end
            else  % �����level,���Ȱ���vertically��ʽ��
                if flagVert
                    LWStrip(1,thisLevel) = LWStrip(1,thisLevel) ...
                        - LWHItemSortHori(2,iStrip); %����wleft by Vertically
                    if flagLargerLengthVert
                        LWStrip(2,thisLevel) = LWHItemSortHori(1,iStrip); %����strip�߶�lleft
                    end
                    % ��flagVert==1 ������� ��Ҫ����Ʒ������rotation(��)��ȥ
                    itemRotaSortHori(iStrip) = ~itemRotaSortHori(iStrip);
                    tep = LWHItemSortHori(1,iStrip);
                    LWHItemSortHori(1,iStrip) = LWHItemSortHori(2,iStrip);
                    LWHItemSortHori(2,iStrip) = tep;
                elseif flagHori
                    LWStrip(1,thisLevel) = LWStrip(1,thisLevel) ...
                        - LWHItemSortHori(1,iStrip); %����wleft by Horizontally
                    if flagLargerLengthHori
                        LWStrip(2,thisLevel) = LWHItemSortHori(2,iStrip); %����strip�߶�lleft
                    end
                else
                    error('levelѡ�����,���ۺ������Ų���');
                end               
            end
            stripBeItemArray(thisLevel) = stripBeItemArray(thisLevel) + 1; %ֻҪ��level����һ��item,����������1
        end
        
        if ParaArray.whichRotation == 1 && ParaArray.whichRotationHori == 1
           % �ж����
            flagHori = LWStrip(1,thisLevel) >=  LWHItemSortHori(1,iStrip); %�ж��Ƿ� wleft >= longest
            flagVert = LWStrip(1,thisLevel) >= LWHItemSortHori(2,iStrip);  %�ж��Ƿ� wleft >= shortest
            isNewLevel = LWStrip(1,thisLevel) == widthStrip; % �ж��Ƿ� new Level
            flagLargerLengthHori = LWStrip(2,thisLevel) <  LWHItemSortHori(2,iStrip); %�ж��Ƿ� Ҫ����level��/�߶�
            flagLargerLengthVert = LWStrip(2,thisLevel) <  LWHItemSortHori(1,iStrip); %�ж��Ƿ� Ҫ����level��/�߶�
            
            % ����strip��Ϣ
            % �����level,���Ȱ���horizontally��ʽ(��Ĭ�Ϸ�ʽ����)
            if isNewLevel
                if flagHori
                    LWStrip(1,thisLevel) = LWStrip(1,thisLevel) - LWHItemSortHori(1,iStrip); %����wleft by Horizontally
                    if flagLargerLengthHori
                        LWStrip(2,thisLevel) = LWHItemSortHori(2,iStrip); %����strip�߶�lleft
                    end
                elseif flagVert
                    LWStrip(1,thisLevel) = LWStrip(1,thisLevel) ...
                        - LWHItemSortHori(2,iStrip); %����wleft by Vertically
                    if flagLargerLengthVert
                        LWStrip(2,thisLevel) = LWHItemSortHori(1,iStrip); %����strip�߶�lleft
                    end
                    % ��flagVert==1 ������� ��Ҫ����Ʒ������rotation(��)��ȥ
                    itemRotaSortHori(iStrip) = ~itemRotaSortHori(iStrip);
                    tep = LWHItemSortHori(1,iStrip);
                    LWHItemSortHori(1,iStrip) = LWHItemSortHori(2,iStrip);
                    LWHItemSortHori(2,iStrip) = tep;
                else      %���wleft <shorter edge
                    error('���ݴ���,���ۺ������Ų���');
                end
            else  % �����level,���Ȱ���vertically��ʽ��
                if flagHori
                    LWStrip(1,thisLevel) = LWStrip(1,thisLevel) ...
                        - LWHItemSortHori(1,iStrip); %����wleft by Horizontally
                    if flagLargerLengthHori
                        LWStrip(2,thisLevel) = LWHItemSortHori(2,iStrip); %����strip�߶�lleft
                    end                    
                elseif flagVert
                    LWStrip(1,thisLevel) = LWStrip(1,thisLevel) ...
                        - LWHItemSortHori(2,iStrip); %����wleft by Vertically
                    if flagLargerLengthVert
                        LWStrip(2,thisLevel) = LWHItemSortHori(1,iStrip); %����strip�߶�lleft
                    end
                    % ��flagVert==1 ������� ��Ҫ����Ʒ������rotation(��)��ȥ
                    itemRotaSortHori(iStrip) = ~itemRotaSortHori(iStrip);
                    tep = LWHItemSortHori(1,iStrip);
                    LWHItemSortHori(1,iStrip) = LWHItemSortHori(2,iStrip);
                    LWHItemSortHori(2,iStrip) = tep;
                else
                    error('levelѡ�����,���ۺ������Ų���');
                end
            end
            stripBeItemArray(thisLevel) = stripBeItemArray(thisLevel) + 1; %ֻҪ��level����һ��item,����������1            
        end
        
        if ParaArray.whichRotation == 1 && ParaArray.whichRotationHori == 2
           % �ж����
            flagHori = LWStrip(1,thisLevel) >=  LWHItemSortHori(1,iStrip); %�ж��Ƿ� wleft >= longest
            flagVert = LWStrip(1,thisLevel) >= LWHItemSortHori(2,iStrip);  %�ж��Ƿ� wleft >= shortest
            isNewLevel = LWStrip(1,thisLevel) == widthStrip; % �ж��Ƿ� new Level
            flagLargerLengthHori = LWStrip(2,thisLevel) <  LWHItemSortHori(2,iStrip); %�ж��Ƿ� Ҫ����level��/�߶�
            flagLargerLengthVert = LWStrip(2,thisLevel) <  LWHItemSortHori(1,iStrip); %�ж��Ƿ� Ҫ����level��/�߶�
            
            % ����strip��Ϣ
            % �����level,���Ȱ���horizontally��ʽ(��Ĭ�Ϸ�ʽ����)
            if isNewLevel
                if flagVert
                    LWStrip(1,thisLevel) = LWStrip(1,thisLevel) ...
                        - LWHItemSortHori(2,iStrip); %����wleft by Vertically
                    if flagLargerLengthVert
                        LWStrip(2,thisLevel) = LWHItemSortHori(1,iStrip); %����strip�߶�lleft
                    end
                    % ��flagVert==1 ������� ��Ҫ����Ʒ������rotation(��)��ȥ
                    itemRotaSortHori(iStrip) = ~itemRotaSortHori(iStrip);
                    tep = LWHItemSortHori(1,iStrip);
                    LWHItemSortHori(1,iStrip) = LWHItemSortHori(2,iStrip);
                    LWHItemSortHori(2,iStrip) = tep;
                elseif flagHori
                    LWStrip(1,thisLevel) = LWStrip(1,thisLevel) - LWHItemSortHori(1,iStrip); %����wleft by Horizontally
                    if flagLargerLengthHori
                        LWStrip(2,thisLevel) = LWHItemSortHori(2,iStrip); %����strip�߶�lleft
                    end
                else      %���wleft <shorter edge
                    error('���ݴ���,���ۺ������Ų���');
                end
            else  % �����level,���Ȱ���vertically��ʽ��
                if flagVert
                    LWStrip(1,thisLevel) = LWStrip(1,thisLevel) ...
                        - LWHItemSortHori(2,iStrip); %����wleft by Vertically
                    if flagLargerLengthVert
                        LWStrip(2,thisLevel) = LWHItemSortHori(1,iStrip); %����strip�߶�lleft
                    end
                    % ��flagVert==1 ������� ��Ҫ����Ʒ������rotation(��)��ȥ
                    itemRotaSortHori(iStrip) = ~itemRotaSortHori(iStrip);
                    tep = LWHItemSortHori(1,iStrip);
                    LWHItemSortHori(1,iStrip) = LWHItemSortHori(2,iStrip);
                    LWHItemSortHori(2,iStrip) = tep;
                elseif flagHori
                    LWStrip(1,thisLevel) = LWStrip(1,thisLevel) ...
                        - LWHItemSortHori(1,iStrip); %����wleft by Horizontally
                    if flagLargerLengthHori
                        LWStrip(2,thisLevel) = LWHItemSortHori(2,iStrip); %����strip�߶�lleft
                    end                    
                else
                    error('levelѡ�����,���ۺ������Ų���');
                end
            end
            stripBeItemArray(thisLevel) = stripBeItemArray(thisLevel) + 1; %ֻҪ��level����һ��item,����������1            
        end
        
        if ParaArray.whichRotation == 0
            % 555����strip��Ϣ
            LWStrip(1,thisLevel) = LWStrip(1,thisLevel) - LWHItemSort(1,iStrip); %����wleft
            if LWHItemSort(2,iStrip) > LWStrip(2,thisLevel)         % �����Ʒ�߶�>��strip�߶�
                LWStrip(2,thisLevel) = LWHItemSort(2,iStrip); %����strip�߶�lleft
            end
            stripBeItemArray(thisLevel) = stripBeItemArray(thisLevel) + 1; %ֻҪ��level����һ��item,����������1
        end

        % 3 ����item����strip��Ϣ
        itemBeStripMatrixSort(1,iStrip) = thisLevel;    %�ڼ���level
        itemBeStripMatrixSort(2,iStrip) = stripBeItemArray(thisLevel); %��level�µڼ��ΰ���
    end
    
    function printscript()
        % ���Դ���
        % % LWHStrip
        % % itemBeStripMatrixSort
        % % LWHItemSort
        % % itemCoordMatrixSort
        % % itemBeStripMatrix
        % % LWHItem
        % % itemCoordMatrix
        %  printstruct(da);
        
        % �����Ҫ���:��ô�1��ʼÿ��strip����������
        for iStrip = 1:max(da.ItemArray.itemBeStripMatrix(1,:))
            [~,idx] = find(da.ItemArray.itemBeStripMatrix(1,:)==iStrip);
            fprintf('strip %d ��ʣ���+���Ϊ:  ',iStrip);
            fprintf('( %d ) ',da.StripArray.LW(:,iStrip));
            fprintf('\n');
            fprintf('strip %d ���� original Item ������(����)[��ת��־]{����}Ϊ  \n  ',iStrip);
            fprintf('%d ',idx);
            fprintf('( %d ) ', da.ItemArray.LWH(1:nDim,idx));fprintf('\n');
            fprintf('[ %d ] ', da.ItemArray.itemRotaFlag(:,idx));fprintf('\n');
            fprintf('{ %d } ', da.ItemArray.CoordItemStrip(:,idx));fprintf('\n');
            fprintf('\n');
        end
    end
end
