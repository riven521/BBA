%% GET ITEM �������
% 1: Item.isHeightFull (�Խ���>=�����϶, ��Ϊ����; Item�����LU�߶� >= �����϶, ��Ϊ����)
% 2: Item.Layer ����ÿ��Item�ڶѶ�Ĳ���
% 3: Item.isWeightFine ÿ��Item���Ƿ�������������
% 4: Item.isNonMixed  ����ÿ��Item�Ƿ�Ϊ����Ҫ��ƴ�Ŀ���
% 5: Item.nbItem ����ÿ��Item����ͬLU/ID������
%% ����
function   [Item,LU] = cpuItem(Item,LU,Veh)
    global ISdiagItem
    %% ��ʼ��
    sz = size(Item.isRota);
    hVeh  = Veh.LWH(3,1);  

    Item.isHeightFull = zeros(sz);    %Item���Ƿ�߶�����(��ʼΪ0)
%     Item.isWeightFine = ones(sz)*-1;    %Item���Ƿ���������(��ʼΪ-1) �޸�Ϊֱ���޸�ÿ��Item
    Item.isNonMixed = ones(sz)*-1;    %Item�Ƿ����Ҫ����ж�,��ż������Item��ǰ����Strip����
    Item.isMixedTile = zeros(sz);    %Item���,�ҳ����������Item��β�и�ֵΪ1

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
        % hMargin: ITEM���복���ļ�϶ ************* 1  �˷����������pingpuallʱ,��������Ҫ��˦β
        %         hMargin = hVeh - Item.LWH(3,iItem);
        % hMargin: ITEM��������Item�����ֵ�ļ�϶ ************* 2
        hMargin = max(Item.LWH(3,:)) - Item.LWH(3,iItem);
        %                 if abs(maxHeightinLUofThisItem - hMargin ) <=60
        %                 hMargin
        %                 maxHeightinLUofThisItem
        %                 end
            % V1: �໥��ͻ
                %         if ISdiagItem==1 && diagItem >= hMargin,  
                %             Item.isHeightFull(iItem) = 1;  else Item.isHeightFull(iItem) = 0; end
                %         if maxHeightinLUofThisItem >= hMargin,   Item.isHeightFull(iItem) = 1;  else  Item.isHeightFull(iItem) = 0; end    
            % V2: ��һ����, ��Ϊ����
        if ISdiagItem==1
            if maxHeightinLUofThisItem >= hMargin || diagItem >= hMargin, Item.isHeightFull(iItem) = 1;  else  Item.isHeightFull(iItem) = 0; end
        else
            if maxHeightinLUofThisItem >= hMargin,   Item.isHeightFull(iItem) = 1;  else  Item.isHeightFull(iItem) = 0; end
        end
    end

    %% SECTION 2 ����ITEM��isWeightFine�������޸�
    % ****************** �������ص��ж�+�޸� ************ ����
    % isWeightUpDown: ITEM�����ж��Ƿ��������ص��ж�Item.isWeightFine

    % V2�� �޸��������ء���check���̾��޸���
    LU  = repairItems(LU);

% V1�� �޸��������ء���
% %     Item = isWeightUpDown(Item,LU);
% %     Item.isWeightFine
% %     % repairItemWeight: ��������������ص�case, �����޸�
% %     if ~all(Item.isWeightFine)
% %         [~,order] = find(Item.isWeightFine == 0);
% %         for iItem=1:length(order)
% %             LU = repairItemWeight(LU,order(iItem));
% %         end
% %     end
% %     
% %     Item = isWeightUpDown(Item,LU);
% %     Item.isWeightFine
% %     checktLU(LU)
% %     if ~all(Item.isWeightFine),   error('������������casse, ����'); end

    %% SECTION 3 ����ITEM��isNonMixed/isMixedTile�Ƿ�Ϊ����Ҫ��ƴ/��ƴ������˦β����
    % ****************** Iten���Ƿ�Ϊ����Ҫ��ƴ���� ************ ����
    % GET Item.isNonMixed: ����ÿ��Item�Ƿ�Ϊ����Ҫ��ƴ�Ŀ���
    ItemLID = cellfun(@(x) x(1), Item.LID); % arrayAllLID: ����ITEM��Ӧ��LIDֵ ������ʽ
    uniItemLID = unique(ItemLID);
    % ѭ��: LID����
    for iItem=1:length(uniItemLID)
        % Item i ���ڵ�LU flag���
        flagItem = ItemLID(:) == iItem;  % flagItem = ItemLID(:) == uniItemLID(iItem);  idExchange(LU.ID)������,�����Ϊ���ߣ���Ӱ��޴� 
        
        ItemWidth = unique(Item.LWH(1,flagItem));     % Item���    
        VehWidth = Veh.LWH(1,1);  % �������    
        maxWidthLayer= floor(VehWidth/ItemWidth); %Item�ɷſ�Ȳ���
        nb = sum(flagItem);
        nbmod = mod(nb,maxWidthLayer);
        if nb ==0 || nbmod>nb, error('cpuItem�ּ���isNonMixed����'); end
        if nbmod == 0 %modΪ0���� ����Ҫ��� ����ϵ���ǰ��order����ǰ
            Item.isNonMixed(flagItem) = 1;
        else
            Item.isNonMixed(flagItem)= 0;
            % ����Item��isMixedTile
            tmpSort=[Item.isHeightFull;Item.HLayer;Item.LWH(3,:)];
            [~, order]=sortrows(tmpSort(:,flagItem)', [1,2,3],{'ascend','ascend','ascend'});
            flagItemIdx = find(flagItem);
            flagmodIdx = flagItemIdx(order(1:nbmod));
            Item.isMixedTile(flagmodIdx)=1;
        end
    end
    
    %% SECTION 4 ����ITEM��nbItem
    tmpItemLID = cell2mat(Item.LID);
    % V1
    %     for i=1:length(Item.Weight)
    %         Item.nbItem(i) = sum(tmpItemLID == tmpItemLID(i));
    %     end
    % V2
    Item.nbItem = sum(tmpItemLID == tmpItemLID');

end

%% �ֲ����� %%

%% ����1 : isWeightFine: �ж�LU�Ƿ��������ع���
%% v2 repairItems
function LU = repairItems(t)  % t�ض����������͵�table��struct
% 0 �ṹ��ת��Ϊtable Ԥ����
if isstruct(t)
    t = struct2table(structfun(@(x) x', t,'UniformOutput',false));  end

% 1.1 �����˲�LU_Item �粻���ڻ�ȡ��ITEMID 
if any(strcmp('LU_Item', t.Properties.VariableNames))
        t.ITEMID = t.LU_Item(:,1);
        t.ITEMSEQ = t.LU_Item(:,2); end

if any(strcmp('CoordLUBin', t.Properties.VariableNames))
    t.X = t.CoordLUBin(:,1);    t.Y = t.CoordLUBin(:,2);   t.Z = t.CoordLUBin(:,3);   end
    
%% 2 ��ͬITEMID�µ�CHEK
uniItemID = unique(t.ITEMID(:));
% 2.1 �������CoordLUBin����,������ϵ�ж���ITEMID�Ĳ���(1:XYֵ(�����Ƿ���ͬ); 2:������Zֵ)
for iItem = 1:length(uniItemID) %��ITEM����ѭ��
    flagIdx = t.ITEMID==uniItemID(iItem); %�Ե�һItemIdȥ�߼�ֵ;
    if any(strcmp('Y', t.Properties.VariableNames))
        vX = t{flagIdx,'X'};
        vY = t{flagIdx,'Y'};
        if any(vX ~= vX(1)) || any(vY ~= vY(1))
            error('��ͬITEM,��X��Y�����λ'); end
        
        v = t(flagIdx,{'Z','Weight','ITEMSEQ'});
    else
        v = t(flagIdx,{'Weight','ITEMSEQ'});
    end
    
    v = sortrows(v,'ITEMSEQ');
    
    if any(strcmp('Y', t.Properties.VariableNames))
        % ����������ǵݼ���Z�߶Ȳ��ǵ���
        if ~issorted(v.Z,'ascend') || ~issorted(v.Weight,'descend')
            issorted(v.Z,'ascend');  issorted(v.Weight,'descend');
            error('��ͬITEM,���������ǵݼ���Z�߶Ȳ��ǵ���'); end
        else
        if ~issorted(v.Weight,'descend')
                                                            %             uniItemID(iItem)
                                                            %             v
            t = repairItemWeight(t,uniItemID(iItem)); % ����t��LU_Item
                                                            %             v = t(flagIdx,{'Weight','LU_Item'});
                                                            %             v = sortrows(v,'LU_Item');
        end
    end
end
LU = t;
if istable(LU)
    LU = table2struct(LU,'ToScalar',true);
    LU = (structfun(@(x) x',LU,'UniformOutput',false));
end

end

%% V1�� isWeightUpDown �޷���ȫ�����������ص�case
% % function Item = isWeightUpDown(Item,LU)
% % for iItem = 1:max(LU.LU_Item(1,:)) %��ITEM����ѭ��
% %     [~,idx] = find(LU.LU_Item(1,:)==iItem);
% %     if iItem==21
% %         1
% %     end
% %     if isempty(idx), 
% %         1
% %     end
% %     nbLUinItem = length(idx);    
% %     % ��ITME�ں�2������LU�Ľ����ж�
% %     if nbLUinItem > 1 %Item������ֻһ��Item,��Ҫ�ж��Ƿ������صı仯
% %         currLUWeight = zeros(1,nbLUinItem);
% %         for iIdx = 1:nbLUinItem
% %             currIdx = idx(LU.LU_Item(2,idx) == iIdx);
% %             currLUWeight(iIdx) = LU.Weight(:,currIdx);          %                 currLUHight(iIdx) = LU.LWH(3,currIdx);
% %         end
% %         diff(currLUWeight)
% %         all(diff(currLUWeight) > 0) 
% %         if all(diff(currLUWeight) > 0) % ������������
% %             % �޸�LU.LU_Item��ֵ             1 5: idx  Ϊ 1 2: LU.LU_Item(2,idx) == iIdx��Ϊ 2 1
% %             Item.isWeightFine(1,iItem) = 0;
% %         else
% %             Item.isWeightFine(1,iItem) = 1;
% %         end
% %     else  %ITEM��ֻ��1��LU, �ض���������
% %         Item.isWeightFine(1,iItem) = 1; 
% %     end
% % end
% % end

%% ����2 : repairItemWeight: LU�粻����������,�����޸�
% V2: �����޸�LU.LU_Item�ĵڶ��е�ֵ LU.LU_Item(2,flagLU)
function LU = repairItemWeight(oLU,itemIdx)
if istable(oLU)
    LU = table2struct(oLU,'ToScalar',true);
    LU = (structfun(@(x) x',LU,'UniformOutput',false));
else
    LU = oLU;
end

flagLU = (LU.LU_Item(1,:)==itemIdx);   %�ҳ���item��Ӧ��lu��index
tmpWeight=LU.Weight(:,flagLU);
[~,b] = sort(tmpWeight,'descend');
[~,b] = sort(b);

LU.LU_Item(2,flagLU) = b;

if istable(oLU) %���ص�ҲҪ��table
    LU = struct2table(structfun(@(x) x',LU,'UniformOutput',false));
end
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