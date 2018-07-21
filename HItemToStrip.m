function [da] = HItemToStrip(da,ParaArray)
% ��Ҫ����:Item����Strip�� %  ����:�����(row);  ����:��������(coloum);
% Input ---  ITEM:  LWH
% Output --- ����ص�
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
%% Item���򲢳�ʼ��
% ��ȡitemorder
% ��ȡLWHItemSort
 % sort ��ȡItem���� %
LWHItem = da.ItemArray.LWH(1:nDim,:);
if ParaArray.whichSortItemOrder == 1 %Descend of ��(��)
    [~,itemorder] = sort(LWHItem(nDim,:),'descend');   % ��Item�ģ�����) �ݼ�����
elseif ParaArray.whichSortItemOrder == 2  %Descend of shortest��̱�  -> ����Rotation���ӱ���
    tmpLWH = [LWHItem; min(LWHItem(1:nDim,:))]; %����������̱��γ���ʱ����tmpLWH
    
    [~,itemorder] = sortrows(tmpLWH',[3 2 1],{'descend','descend','descend'}); %way1
    [~,itemorder] = sort(tmpLWH(nDim+1,:),'descend'); %��ȡ��ʱ����������˳�� way2
    clear tmpLWH;
else
    error('���ò�������');
end
LWHItemSort = LWHItem(:,itemorder);

%% ����Rotation���ӱ��� 
% ��ȡLWHItemSortHori
% ��ȡitemRotaSortHori
% itemRotaSortHori(ԭʼitem�Ƿ�rotation 0-1����) LWHItemSortHori(Hori orienting��LWH)
    itemRotaSortHori = zeros(1,nItem); % 0 ����horizontal orientation 1 ���� vertical orientation
if ParaArray.whichRotation == 1
%     itemRotaSortHori = zeros(1,nItem); % 0 ����horizontal orientation 1 ���� vertical orientation
    LWHItemSortHori = zeros(nDim,nItem);
    [LWHItemSortHori,itemRotaSortHori] = horiOrient(LWHItemSort);
end

%% 55 LU->Item->Stripת�� 
% ��ȡitemBeStripMatrixSort: ÿ�������Item���ĸ�Strip��  �Լ�˳��
% ��ȡLWStrip:  �����ɵ�Strip�ĳ���
% ��ȡCoordItemStripSort  Item��strip������ֵ
LWStrip = zeros(nDim,nItem);   %strip���� dim2-����(����ߵļ���) (�߶Ƚ����ο�,current �߶�)
LWStrip(1,:) = widthStrip;   %dim1-���ʣ�� 
stripBeItemArray = zeros(1,nStrip);  % ÿ��Strip�ڵ�Item���� ���ڲ���
itemBeStripMatrixSort = zeros(2,nItem); %dim1:���ڵڼ���level dim2:���ڸ�level�ڼ����ŷ� 555
CoordItemStripSort = zeros(2,nItem); %Item��strip������ֵ

iLevel = 1; iStrip = 1;
while 1
    if iStrip > nItem, break; end
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
            continue;
        else
            % ��ȡthisLevel: Ψһ��FF������⵽thisLevel�ļ��㣨ѡ��������������С�ģ�
            tepAvailableLevelArray = LWStrip(1,1:iLevel);   %��ȡ�����Ѱ��Ż��°��ŵ�level��ʣ��ˮƽƽ�������tepAvailableLevelArray
            tepMinLeftWdith = min(tepAvailableLevelArray(flag));                      %�ҳ�tepAvailableLevelArray�п����ɱ�iITem��ʣ��ˮƽ��ȵ���Сֵ������ΪtepMinLeftWdith
            thisLevel = find(tepAvailableLevelArray==tepMinLeftWdith);            %�ҳ�����Сֵ��Ӧ���Ǹ�(�飩level
            if length(thisLevel)>1
                thisLevel = thisLevel(1);
            end
        end
    elseif ParaArray.whichStripH == 2
        % ��ͬ�����µ�ѡ�����find����㹻�Ķ��level,�������ڵ�һ�������� Ψһ������thisLevel�Ļ�ȡ
        flag = find(LWStrip(1,1:iLevel) >= LWHItemSort(1,iStrip));
        if isempty(flag)
            iLevel = iLevel + 1;% �����Ȳ����㣬��level����
            continue;
        else
            thisLevel = flag(1);
        end
    elseif ParaArray.whichStripH == 3
        % ��ͬ�����µ�ѡ�������ǰitem�Ŀ�<=��ǰstrip�ĵ�ǰlevel�Ŀ�
        flag = LWHItemSort(1,iStrip) <= LWStrip(1,iLevel);
        if ~flag  %ע����ǰ��isempty(flag)������
            iLevel = iLevel + 1;% �����Ȳ����㣬��level����
            continue;
        else
            thisLevel = iLevel;
        end
    end    
    
    insertItemToStrip();
    iStrip = iStrip + 1;
end

%% ���� ����ֵ��da
% ��ȡitemBeStripMatrix : ÿ��Item���ĸ�Strip��  �Լ�˳��
% ��ȡCoordItemStrip : ÿ��Item��Strip������
% ��ȡitemRotaFlag : ÿ��item�Ƿ�Rotation�ı�־
% ��ȡLWStrip:  �����ɵ�strip�ĳ���
% ��ȡitemorder: item������
        %Matalb code gerator use:
        itemBeStripMatrix=itemBeStripMatrixSort;CoordItemStrip=CoordItemStripSort;itemRotaFlag=itemRotaSortHori;
        
    itemBeStripMatrix(:,itemorder) = itemBeStripMatrixSort;
    CoordItemStrip(:,itemorder) = CoordItemStripSort;
    da.ItemArray.itemBeStripMatrix = itemBeStripMatrix;
    da.ItemArray.CoordItemStrip = CoordItemStrip;
    
    % ��ʹû��RotationFlage Ҳ�и����� �ж��Ƿ�Rotation
        itemRotaFlag(:,itemorder) = itemRotaSortHori;
        da.ItemArray.itemRotaFlag = itemRotaFlag;

LWStrip = LWStrip(:,LWStrip(2,:)>0);
da.StripArray.LW = LWStrip; % ȥ��δʹ�õ�Strip
da.ItemArray.itemorder = itemorder;

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



    function insertItemToStrip()
        CoordItemStripSort=CoordItemStripSort;LWStrip=LWStrip;itemRotaSortHori=itemRotaSortHori;LWHItemSortHori=LWHItemSortHori;stripBeItemArray=stripBeItemArray;itemBeStripMatrixSort=itemBeStripMatrixSort;
        
        % 1 ����CoordItemStripSort
        CoordItemStripSort(1,iStrip) = widthStrip - LWStrip(1,thisLevel);  %����x����
        CoordItemStripSort(2,iStrip) = sum(LWStrip(2,1:thisLevel-1));      %����y���� %���iLevel=1,�����ߣ�����Ϊ0������Ϊ���
        
        % 2 ����stripBeItemArray
        % 2 ����LWStrip
        if ParaArray.whichRotation == 1 && ParaArray.whichRotationHori == 0
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
end
