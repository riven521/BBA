%% GET ITEM �������
% 1: Item.isHeightFull ( Item�����LU�߶� * 0.95 >= �����϶, ��Ϊ����)  computeisHeightFullItem
% 2: Item.Layer ����ÿ��Item�ڶѶ�Ĳ���
% 3: Item.isWeightFine ÿ��Item���Ƿ�������������
% 4: Item.isNonMixed  ����ÿ��Item�Ƿ�Ϊ����Ҫ��ƴ�Ŀ���
% 5: Item.nbItem ����ÿ��Item����ͬLU/ID������
%% ����
function   [Item,LU] = cpuItem(Item,LU,Veh)

    %% ��ʼ��
    sz = size(Item.Weight);
    % hVeh  = Veh.LWH(3,1);  % ��ʱΪ��

    Item.isHeightFull = zeros(sz);    %Item���Ƿ�߶�����(��ʼΪ0)  %  Item.isWeightFine = ones(sz)*-1;    %Item���Ƿ���������(��ʼΪ-1) �޸�Ϊֱ���޸�ÿ��Item
    Item.isNonMixed = ones(sz)*-1;    %1:ITEM���������ŵ�����strip,���������ITEM���; 0: �б�Ҫ���; ��ż������Item��ǰ����Strip����
    Item.isMixedTile = zeros(sz);          % 0:��isNonMixed=1,ȫ��Ϊ0;��isNonMixed=0ʱ,�������±�����β��ITEM��ֵֵΪ1. 
    Item.nbItem = zeros(sz);  
    
    %% V2 SECTION 0 ����ITEM��PID,LID,SID,��TABLE����,����֪��ʲô��ʲô,����1,2,3�����滻
    % �ɻ�ϵ�LU.DOC����ITEM�ڰ�����PID,LID,SID������ 1808���� ����Item.PID,LID,SID��ʹ��
    t = struct2table(structfun(@(x) x',LU,'UniformOutput',false));
    
    nItem = size(Item.LWH,2);
    for iItem=1:nItem
        f = t.LU_Item(:,1) == iItem;                          %  f = logical(ones(height(t),1))
        Item.LID(:,iItem) = {unique(t.ID(f))};             % NOTE: ITEM���LID��LU��ID
        %         Item.LID(:,iItem) = {unique(t.LID(f))};
        Item.SID(:,iItem) = {unique(t.SID(f))};
        Item.EID(:,iItem) = {unique(t.EID(f))};
        Item.PID(:,iItem) = {unique(t.PID(f))};
    end    
    %  t2 = struct2table(structfun(@(x) x',Item,'UniformOutput',false));

    %% V1 : ��LU.DOC����
% %     %% SECTION 0 ����ITEM��PID,LID,SID
% %     % �ɻ�ϵ�LU.DOC����ITEM�ڰ�����PID,LID,SID������ 1808���� ����Item.PID,LID,SID��ʹ��
% %     LU.DOC=[LU.PID;LU.ID;LU.SID;zeros(size(LU.ID));zeros(size(LU.ID));...
% %         LU.LU_Item;];
% %     nItem = size(Item.LWH,2);
% %     for iItem=1:nItem
% %         tmp = LU.DOC([1,2,3], LU.DOC(6,:) == iItem);
% %         Item.PID(:,iItem) = num2cell(unique(tmp(1,:))',1); %unique(tmp(1,:))';
% %         Item.LID(:,iItem) = num2cell(unique(tmp(2,:))',1);
% %         Item.SID(:,iItem) =num2cell(unique(tmp(3,:))',1);
% %     end
    
    %% SECTION 1 ����ITEM��isHeightFull
    Item.isHeightFull = computeisHeightFullItem(Item,LU);

    %% SECTION 2 ����ITEM��isWeightFine�������޸�
    % V2�� �޸��������ء���check���̾��޸���
    LU.LU_Item(2,:)= repairItems(LU);    
    % �Ƿ��Ծ���������
    chktLU(LU)
    
    %% 	SECTION 3 Item.isNonMixed: and  Item.isMixedTile
    % ****************** Iten���Ƿ�Ϊ����Ҫ��ƴ���� ************ ����
    % GET Item.isNonMixed: ����ÿ��Item�Ƿ�Ϊ����Ҫ��ƴ�Ŀ���
    ItemLID = cellfun(@(x) x(1), Item.LID); % arrayAllLID: ����ITEM��Ӧ��LIDֵ ������ʽ
    uniItemLID = unique(ItemLID);
    % ѭ��: LID����
    for iItem=1:length(uniItemLID)
        % Item i ���ڵ�LU flag���
        flagItem = ItemLID(:) == iItem;  % flagItem = ItemLID(:) == uniItemLID(iItem);  idExchange(LU.ID)������,�����Ϊ���ߣ���Ӱ��޴� 
        
        ItemWidth = unique(Item.LWH(1,flagItem));     % Item����
        VehWidth = Veh.LWH(1,1);  % ��������
        
        nbmaxItem= floor(VehWidth/ItemWidth);     % Item�ɷ������Ȳ���
        nbItem = sum(flagItem);
        
        % GapWidth ����������������ʣ��Item���Gap����
        GapWidth = VehWidth - nbmaxItem*ItemWidth;
        
        % nbmod: ������������ �� �������� ȡ�� 
        nbmod = mod(nbItem,nbmaxItem);        if nbItem ==0 || nbmod>nbItem, error('cpuItem�ּ���isNonMixed����'); end
        
        % ���㣺flagNonMixedItem��1 ��������Item��������ϵģ���������0��������������������϶Ѷ�����Դ�
        % ������Ϊ0����һ�����������пգ����������Item��ϣ�
        if nbmod~=0 
            flagNonMixedItem = 0;
        else
            % ������Ϊ0����ֻ��һ��Strip��ʣ��Gap�����ֽϴ���ָ��Ϊ����������Item��ϣ�
            if nbItem==nbmaxItem && GapWidth > ItemWidth*0.5
                flagNonMixedItem = 0;
            else % ����ָ��Ϊ������������Item���
                flagNonMixedItem = 1;
            end
        end
        
        % ��ֵ������flagNonMixedItem����ֵItem��isNonMixed/isMixedTile ��Item�������Ҫ���ݣ�
        if   flagNonMixedItem 
            Item.isNonMixed(flagItem) = 1;  % ���ItemΪ�����ϵ�, ����˦β�ı�Ǿ�Ϊ0
            Item.isMixedTile(flagItem) = 0;
        else  
            Item.isNonMixed(flagItem)= 0;
            
            % Item.isMixedTile  ����Item��isMixedTile �����ҳ���������ϵ�Items�������Ա��
            if nbmod == 0 % �� nbItem==nbmaxItem && GapWidth > ItemWidth*0.5 ��ʵ������mix����Gap̫����ֻ��һ��
                Item.isMixedTile(flagItem)=1;  %NOTE: flagItem���е�Item�����
            else            
                tmpSort=[Item.isHeightFull;Item.HLayer;Item.LWH(3,:)];
                [~, order]=sortrows(tmpSort(:,flagItem)', [1,2,3],{'ascend','ascend','ascend'});
                flagItemIdx = find(flagItem);
                flagmodIdx = flagItemIdx(order(1:nbmod));
                Item.isMixedTile(flagmodIdx)=1;  %NOTE: flagmodIdx ����Item���
            end
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
function TF = computeisHeightFullItem(Item,LU)
    global ISdiagItem
    
    n = length(Item.Weight);
    TF = deal(ones(1,n)*-1);    
    
    % (�Խ���>=�����϶, ��Ϊ����; Item�����LU�߶� >= �����϶, ��Ϊ����)
    for iItem=1:n
        % ItemDiag : ITEM�ĶԽ��߳��� 
        ItemDiag = sqrt(Item.LWH(1,iItem)^2 + Item.LWH(2,iItem)^2);
        
        % maxHeightinLUofThisItem: Item��Lu����ߵĸ߶�
        flagLU = LU.LU_Item(1,:) == iItem;
        maxHeightinLUofThisItem = max(LU.LWH(3,flagLU));
        
        % hMargin: ITEM���복���ļ�϶ ************* V1  �˷����������pingpuallʱ,��������Ҫ��˦β
        %         hMargin = hVeh - Item.LWH(3,iItem);
        
        % hMargin: ��ITEM�������г�������Item�Ѷ�����ֵ�ļ�϶(���) ************* V2
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
        % ����1: ���Ѷ������LU�߶� > ���Ѷ�����߶Ѷ�ĸ߶Ȳ� -> ���� || �Ѷ�Խ��� >  ���Ѷ�����߶Ѷ�ĸ߶Ȳ�
        % ����2: ���Ѷ������LU�߶�*0.95 > ���Ѷ�����߶Ѷ�ĸ߶Ȳ� %����Ϊ2
        if ISdiagItem==1
            if maxHeightinLUofThisItem >= hMargin || ItemDiag >= hMargin, TF(iItem) = 1;  else  TF(iItem) = 0; end
        else
            if 0.95*maxHeightinLUofThisItem >= hMargin,   TF(iItem) = 1;  else  TF(iItem) = 0; end
        end
        
    end
    
end
    
%% ����1 : isWeightFine: �ж�LU�Ƿ��������ع���
% v3 repairItems
function ITEMSEQ = repairItems(LU)  % t�ض����������͵�table��struct
    
T = getTableLU(LU);

uniItemID = unique(T.ITEMID(:));
for iItem = 1:length(uniItemID)         %��ITEM����ѭ��
    
    flagLUIdx = T.ITEMID==uniItemID(iItem); %�Ե�һItemIdȥ�߼�ֵ;
    
    subT = T(flagLUIdx,{'Weight','ITEMSEQ'});
    
    if isUpDownWeight(subT) 
        T.ITEMSEQ(flagLUIdx) = repairItemWeight(T,flagLUIdx); % ����t��LU_Item
        ITEMSEQ = T.ITEMSEQ';
    end
    
end

% LU = getSturctT(T); %���ò����ؽṹ��LU
end


%% v2
% % function LU = repairItems(t)  % t�ض����������͵�table��struct
% %     chktLU(t)
% % % % 0 �ṹ��ת��Ϊtable Ԥ����
% % T = getTableLU(t)
% % T.Properties.VariableNames
% %    t=T
% % 
% % %% 2 ��ͬITEMID�µ�CHEK
% % 
% % uniItemID = unique(t.ITEMID(:));
% % % 2.1 �������CoordLUBin����,������ϵ�ж���ITEMID�Ĳ���(1:XYֵ(�����Ƿ���ͬ); 2:������Zֵ)
% % for iItem = 1:length(uniItemID) %��ITEM����ѭ��
% %     flagIdx = t.ITEMID==uniItemID(iItem); %�Ե�һItemIdȥ�߼�ֵ;
% %     if any(strcmp('Y', t.Properties.VariableNames))
% %         vX = t{flagIdx,'X'};
% %         vY = t{flagIdx,'Y'};
% %         if any(vX ~= vX(1)) || any(vY ~= vY(1))
% %             error('��ͬITEM,��X��Y�����λ'); end
% %         
% %         v = t(flagIdx,{'Z','Weight','ITEMSEQ'});
% %     else
% %         v = t(flagIdx,{'Weight','ITEMSEQ'});
% %     end
% %     
% %     v = sortrows(v,'ITEMSEQ');
% %     
% %     if any(strcmp('Y', t.Properties.VariableNames))
% %         % ����������ǵݼ���Z�߶Ȳ��ǵ���
% %         if ~issorted(v.Z,'ascend') || ~issorted(v.Weight,'descend')
% %             issorted(v.Z,'ascend')  
% %             issorted(v.Weight,'descend')
% %             error('��ͬITEM,���������ǵݼ���Z�߶Ȳ��ǵ���'); end
% %         else
% %         if ~issorted(v.Weight,'descend')
% %                                                             %             uniItemID(iItem)
% %                                                             %             v
% %             t = repairItemWeight(t,uniItemID(iItem)); % ����t��LU_Item
% %                                                             %             v = t(flagIdx,{'Weight','LU_Item'});
% %                                                             %             v = sortrows(v,'LU_Item');
% %         end
% %     end
% % end
% % LU = t;
% % 
% % % ɾ����ʼ����ȡ����
% % if any(strcmp('ITEMID', LU.Properties.VariableNames))
% %         LU.ITEMID = []; end
% % if any(strcmp('ITEMSEQ', LU.Properties.VariableNames))
% %         LU.ITEMSEQ = []; end
% %     if any(strcmp('X', LU.Properties.VariableNames))
% %         LU.X = []; end
% %     if any(strcmp('Y', LU.Properties.VariableNames))
% %         LU.Y = []; end
% %     if any(strcmp('Z', LU.Properties.VariableNames))
% %         LU.Z = []; end
% %     
% % if istable(LU)
% %     LU = table2struct(LU,'ToScalar',true);
% %     LU = (structfun(@(x) x',LU,'UniformOutput',false));
% % end
% % 
% % end
%% ����2 : repairItemWeight: LU�粻����������,�����޸�
% V3: ����struct������table�޸� �����޸�LU.LU_Item�ĵڶ��е�ֵ LU.LU_Item(2,flagLU)
function b = repairItemWeight(T,flagLUIdx)
    tmpWeight=T.Weight(flagLUIdx);
    [~,b] = sort(tmpWeight,'descend');
    [~,b] = sort(b);
    T.ITEMSEQ(flagLUIdx) = b;
end

% % %% ����2 : repairItemWeight: LU�粻����������,�����޸�
% % % V2: �����޸�LU.LU_Item�ĵڶ��е�ֵ LU.LU_Item(2,flagLU)
% % function LU = repairItemWeight(oLU,itemIdx)
% %     if istable(oLU)
% %         LU = table2struct(oLU,'ToScalar',true);
% %         LU = (structfun(@(x) x',LU,'UniformOutput',false));
% %     else
% %         LU = oLU;
% %     end
% % 
% %     flagLU = (LU.LU_Item(1,:)==itemIdx);   %�ҳ���item��Ӧ��lu��index
% %     tmpWeight=LU.Weight(:,flagLU);
% %     [~,b] = sort(tmpWeight,'descend');
% %     [~,b] = sort(b);
% % 
% %     LU.LU_Item(2,flagLU) = b;
% % 
% %     if istable(oLU) %���ص�ҲҪ��table
% %         LU = struct2table(structfun(@(x) x',LU,'UniformOutput',false));
% %     end
% % end


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

%% V1�� repairItemWeight ����insert2Item�и����ʼֵ, �˴��ٽ����޸�
% ****************** �Խ��߿����޸� ************ �ر�
% % % V1: ����insert2Item�и����ʼֵ, �˴��ٽ����޸�
% % % repairItemFull: ������ڷ������case, ����΢��:Ϊ0�ĸ�Ϊ1,��ʣ��߶�С��ITEM�Խ��ߵĸ߶�, ��ΪFULL ����
% % if ~all(Item.isHeightFull)
% %     [~,b] = find(Item.isHeightFull == 0);
% %     for i=1:length(b)
% %            Item = repairItemFull(Item,hVeh,b(i)); %DONE 
% %     end
% % end