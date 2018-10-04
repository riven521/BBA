%% GET ITEM �������
% 1: Item.isHeightFull (�Խ���>=�����϶, ��Ϊ����; Item�����LU�߶� >= �����϶, ��Ϊ����)
% 2: Item.Layer ����ÿ��Item�ڶѶ�Ĳ���
% 3: Item.isWeightFine ÿ��Item���Ƿ�������������
% 4: Item.isNonMixed  ����ÿ��Item�Ƿ�Ϊ����Ҫ��ƴ�Ŀ���

%% ����
function   [Item,LU] = cpuItem(Item,LU,Veh)
    %% ��ʼ��
    sz = size(Item.isRota);
    hVeh  = Veh.LWH(3,1);  

    Item.isHeightFull = ones(sz)*-1;    %Item���Ƿ�߶�����(��ʼΪ-1)
    Item.isWeightFine = ones(sz)*-1;    %Item���Ƿ���������(��ʼΪ-1)
    Item.isNonMixed = ones(sz)*-1;    %Item�Ƿ����Ҫ����ж�,��ż������Item��ǰ����Strip����

    %% SECTION 0 ����ITEM��PID,LID,SID
    % �ɻ�ϵ�LU.DOC����ITEM�ڰ�����PID,LID,SID������ 1808���� ����Item.PID,LID,SID��ʹ��
    LU.DOC=[LU.PID;LU.ID;LU.SID;zeros(size(LU.ID));zeros(size(LU.ID));...
        LU.LU_Item;];
    nItem = size(Item.LWH,2);
    for iItem=1:nItem
        tmp = LU.DOC([1,2,3], LU.DOC(6,:) == iItem);
        Item.PID(:,iItem) = num2cell(unique(tmp(1,:))',1); %unique(tmp(1,:))';
        Item.LID(:,iItem) = num2cell(unique(tmp(2,:))',1);
        Item.SID(:,iItem) =num2cell(unique(tmp(3,:))',1);
    end



    %% SECTION 1 ����ITEM��isHeightFull
    % (�Խ���>=�����϶, ��Ϊ����; Item�����LU�߶� >= �����϶, ��Ϊ����)
    for iItem=1:length(Item.isHeightFull)
        % diagItem: ITEM�ĶԽ��߳��� 
        diagItem = sqrt(Item.LWH(1,iItem)^2 + Item.LWH(2,iItem)^2);
        % maxHeightinLUofThisItem: Item��Lu����ߵĸ߶�
        flagLU = LU.LU_Item(1,:) == iItem;
        maxHeightinLUofThisItem = max(LU.LWH(3,flagLU));
        % hMargin: ITEM���복���ļ�϶
        hMargin = hVeh - Item.LWH(3,iItem);
        
        if diagItem >= hMargin,  Item.isHeightFull(iItem) = 1;  else Item.isHeightFull(iItem) = 0; end
        if maxHeightinLUofThisItem >= hMargin,   Item.isHeightFull(iItem) = 1;  else  Item.isHeightFull(iItem) = 0; end    
    end
       
    %% SECTION 2 ����ITEM��isWeightFine�������޸�
    % ****************** �������ص��ж�+�޸� ************ ����
    % isWeightUpDown: ITEM�����ж��Ƿ��������ص��ж�Item.isWeightFine
    Item = isWeightUpDown(Item,LU);
    % repairItemWeight: ��������������ص�case, �����޸�
    if ~all(Item.isWeightFine)
        [~,b] = find(Item.isWeightFine == 0);
        for iItem=1:length(b)
            LU = repairItemWeight(LU,b(iItem));
        end
    end
    Item = isWeightUpDown(Item,LU);
    if ~all(Item.isWeightFine),   error('������������casse, ����'); end

    %% SECTION 3 ����ITEM��isNonMixed�Ƿ�Ϊ����Ҫ��ƴ����
    % ****************** Iten���Ƿ�Ϊ����Ҫ��ƴ���� ************ ����
    % GET Item.isNonMixed: ����ÿ��Item�Ƿ�Ϊ����Ҫ��ƴ�Ŀ���
    ItemLID = cellfun(@(x) x(1), Item.LID); % arrayAllLID: ����ITEM��Ӧ��LIDֵ ������ʽ
    for iItem=1:length(unique(ItemLID))
        % Item i ���ڵ�LU flag���
        flagItem = ItemLID(:) == iItem;
        ItemWidth = unique(Item.LWH(1,flagItem));     % Item���    
        VehWidth = Veh.LWH(1,1);  % �������    
        maxWidthLayer= floor(VehWidth/ItemWidth); %Item�ɷſ�Ȳ���
        nb = sum(flagItem);
        if mod(nb,maxWidthLayer) == 0 %modΪ0���� ����Ҫ��� ����ϵ���ǰ��order����ǰ
            Item.isNonMixed(flagItem) = 1;
        else
            Item.isNonMixed(flagItem)= 0;
        end
    end
end

%% �ֲ����� %%

%% ����1 : isWeightFine: �ж�LU�Ƿ��������ع���
function Item = isWeightUpDown(Item,LU)
for iItem = 1:max(LU.LU_Item(1,:)) %��ITEM����ѭ��
    [~,idx] = find(LU.LU_Item(1,:)==iItem);
    nbLUinItem = length(idx);
    % ��ITME�ں�2������LU�Ľ����ж�
    if length(idx) > 1 %Item������ֻһ��Item,��Ҫ�ж��Ƿ������صı仯
        currLUWeight = zeros(1,nbLUinItem);
        for iIdx = 1:nbLUinItem
            currIdx = idx(LU.LU_Item(2,idx) == iIdx);
            currLUWeight(iIdx) = LU.Weight(:,currIdx);          %                 currLUHight(iIdx) = LU.LWH(3,currIdx);
        end
        if diff(currLUWeight) > 0 % ������������
            % �޸�LU.LU_Item��ֵ             1 5: idx  Ϊ 1 2: LU.LU_Item(2,idx) == iIdx��Ϊ 2 1
            Item.isWeightFine(1,iItem) = 0;
        else
            Item.isWeightFine(1,iItem) = 1;
        end
    else  %ITEM��ֻ��1��LU, �ض���������
        Item.isWeightFine(1,iItem) = 1; 
    end
end
end

%% ����2 : repairItemWeight: LU�粻����������,�����޸�
function LU = repairItemWeight(LU,itemIdx)
    [~,LUidx] = find(LU.LU_Item(1,:)==itemIdx); %�ҳ���item��Ӧ��lu��index
    nbLUinItem = length(LUidx);
    currLUWeight = zeros(1,nbLUinItem);
    for iIdx = 1:nbLUinItem
        currIdx = LUidx(LU.LU_Item(2,LUidx) == iIdx);
        currLUWeight(iIdx) = LU.Weight(:,currIdx);
    end
    % ��Item�ڵ�LU����������b; ������˳��LU_Item��˳��2�����޸�
    [~,b] = sort(currLUWeight,'descend');
    tt = LU.LU_Item(2,LUidx);
    LU.LU_Item(2,LUidx) = tt(b);
end


% ****************** �Խ��߿����޸� ************ �ر�
% % % V1: ����insert2Item�и����ʼֵ, �˴��ٽ����޸�
% % % repairItemFull: ������ڷ������case, ����΢��:Ϊ0�ĸ�Ϊ1,��ʣ��߶�С��ITEM�Խ��ߵĸ߶�, ��ΪFULL ����
% % if ~all(Item.isHeightFull)
% %     [~,b] = find(Item.isHeightFull == 0);
% %     for i=1:length(b)
% %            Item = repairItemFull(Item,hVeh,b(i)); %DONE 
% %     end
% % end