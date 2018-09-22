function [LU,Item,ItemID] = HLUtoItem(LU,Veh)
% ��Ҫ����:LU�Ѷ���γ�Item %  ����:�����(row);  ����:��������(coloum);
% Input ---  LU: ID LWH Weight ��LU: ����ԭ��˳��
% Output --- LU: order LU_Item ��LU: ����ԭ��˳��(ORDER�ǽ���Item�㷨��LU˳��)
% Output --- Item: ID LWH Weight ...��ITEM:û��˳���㷨������˳��

% LU.LU_Item           (2,n) : ��1: LU�ڵڼ���item ��2:LU�����item��˳��(��-��))
% LU.order                (1,n):  LU����˳��)
% Item.LID                ��1,n): ITEM������ ��ͬ�ڲ�LU����)
% Item.LWH              (3,n): ITEM�ĳ����-������LU��ͬ,�߶�Ϊ�Ѷ��߶�)
% Item.Weight          (1,n): ITEM������
% tmpItem_LU         (1,n): ��1 ITEM��LU����

%% LU����
% ��ȡLU��˳��(�ص��Ǹ߶ȵݼ�����)

[LU.order]  = getLUorder(LU); %��ȡ LU����(��ID����,��߶ȵݼ�)
% printstruct(LU)
% ��ȡ��order������LU:sLU
if isSameCol(LU)
    sLU = structfun(@(x) x(:,LU.order),LU,'UniformOutput',false);
else
    error('����ʹ��structfun');
end
% LU.order(:,LU.order)
% os = sLU.order
% printstruct(sLU)

%% 55 LU->Itemת��

% ������sLU, �����Ѷ��ȡ��Item,�Լ�sLU��Item�ڵ�˳��
sz = size(sLU.ID);
nLU = sz(2);
hVeh  = Veh.LWH(3,1);  % tmpUniqueBin = unique(Veh.LWH(1:3,:)','rows')'; % hVeh = tmpUniqueBin(3);

% �����ʼ����Ҫ������fields
%     Item.LID = zeros(sz);             %Item��ID����
%     Item.SID = zeros(sz);
%     Item.UID = zeros(sz);
%     Item.PID = zeros(numel(unique(LU.PID)),sz(2));
    
    Item.isRota = ones(sz)*-1;    %Item�Ŀ���ת����(��ʼΪ-1)
    Item.Rotaed = ones(sz)*-1;
    Item.isFull = ones(sz)*-1;    %Item���Ƿ�����(��ʼΪ-1)
    Item.isWeightFine = ones(sz)*-1;    %Item���Ƿ���������(��ʼΪ-1)
Item.LWH = zeros(3,nLU); % Item.LWH(1,:) = wStrip;   %dim1-���ʣ��  % Item.LWH(3,:) = hVeh; % 
Item.Weight = zeros(1,nLU); %Item������
% ��ʱʹ��
tmpItem_LU = zeros(1,nLU);  % ��1��ÿ��Strip�ڵ�Item���� �� ��2��ÿ��Strip�ڵĲ�ͬLUID����
% sLU����
sLU.LU_Item = zeros(2,sz(2));     %dim1:���ڵڼ���Item dim2:���ڸ�Item�ڼ����ŷ�

iItem = 1; iLU = 1; %iStrip����itemʵ��
% �̶�LU, ѡ��ITEM; 
while 1
    if iLU > nLU, break; end
    [thisItem,iItem] = getThisItem(iItem);
    insertLUToItem(thisItem,iLU);
    iLU = iLU + 1;
end


% Get ITEM ��ؿ��Է� NEXT FIT 
    function [thisItem,iItem] = getThisItem(iItem)
    % isflagHeight : �Ƿ�ITEM�߶�����
    % isNewItem2 ���Ƿ�ITEM������
    % isSameID2 �� �Ƿ�ITEM�ڵ�ID��ͬ
    
        % ͬ��SID/UID ͬ��LUID Item�߶����� δ����Weight��
        isflagHeight =hVeh - Item.LWH(3,iItem) >= sLU.LWH(3,iLU); %�ж��Ƿ�current's itemʣ���� >= ��ǰiLU�߶�
        % 1 �����ж��Ƿ�������
        if hVeh - Item.LWH(3,iItem) >= sLU.LWH(3,iLU)*2
            Item.isFull(1,iItem) = 0;
        else
            Item.isFull(1,iItem) = 1;
        end
            
        flagLUinItem = sLU.LU_Item(1,:) == iItem;
        if ~any(flagLUinItem) %�����iItem�ڲ���������LU,����Item
            isNewItem2 = 1;
        else
            isNewItem2 = 0;
            if ~isscalar(unique(sLU.ID(flagLUinItem))),    error('��Ԥ�ڴ���');     end            
            isSameID2 = unique(sLU.ID(flagLUinItem)) ==  sLU.ID(iLU);  %����V2�汾:�ж�iLU��Item��LU�Ƿ�����ͬһ��ID
        end
        
            % �ϰ汾V1
                %         isSameID = Item.LID(iItem) == sLU.ID(iLU); %�ж�Item�ڲ�ID�Ƿ�=��ǰiLU��ID
                %         isNewItem = Item.LWH(3,iItem) == 0; % �ж��Ƿ� new Item �߶�==0
        
       % �������TIEM, һ���ɷţ���������߶����� �� �뱾ITEM�ڵ�ID��ͬ��Ҳ�ɷ�;
        if isNewItem2
                thisItem = iItem;
        else
            if isflagHeight && isSameID2 %����߶����� ��LU ID��ͬ %TODO ��������Ҫ���������ڱ�ITEM������
                 thisItem = iItem;
            else
               % 2 �����ж��Ƿ�������
                if isflagHeight && ~isSameID2 %����߶�������LU ID��ͬ, ������ITEM�Ƿ�����
                    Item.isFull(1,iItem) = 0;  % �����Ƿ�������
                end
                if ~isflagHeight && isSameID2 %����߶Ȳ�������LU ID��ͬ, ������ITEM������
                    Item.isFull(1,iItem) = 1;  % �����Ƿ�������
                end
                iItem = iItem + 1;
                [thisItem,iItem] = getThisItem(iItem);
            end
            

        end
    end

% Put LU into thisItem
    function insertLUToItem(thisItem,iLU)
        %����Height
        Item.LWH(3,thisItem) = Item.LWH(3,thisItem)  + sLU.LWH(3,iLU);
        Item.LWH(1:2,thisItem) = sLU.LWH(1:2,iLU);  %����item����
        Item.Weight(1,thisItem) = Item.Weight(1,thisItem) + sLU.Weight(1,iLU); %����item����
        
        tmpItem_LU(1,thisItem) = tmpItem_LU(1,thisItem) + 1;
        
        sLU.LU_Item(1,iLU) = thisItem;
        sLU.LU_Item(2,iLU) = tmpItem_LU(1,thisItem);
        
                    %         tmpLUThisItem = sLU.LU_Item(1,:) == thisItem;
                    %         tmpItem_LU(2,iItem) = numel(unique(sLU.PID(1,tmpLUThisItem)));

        Item.isRota(1,thisItem) = sLU.isRota(1,iLU);  %����ID����ת����
        Item.Rotaed(1,thisItem) = sLU.Rotaed(1,iLU);  %����ID��ת���
        
%         Item.LID(1,thisItem) = sLU.ID(1,iLU); %����ID����        
%         Item.SID(1,thisItem) = sLU.SID(1,iLU);   % Item.UID(1,thisItem) = sLU.UID(1,iLU);
        
%         Item.PID(sLU.PID(1,iLU),thisItem) = Item.PID(sLU.PID(1,iLU),thisItem) + 1;    % 555 ���¶���PID - ��ֵΪ���ִ���
%         Item.PID(sLU.PID(1,iLU),thisItem) = 1;      % 555 ���¶���PID - ��ֵΪ�������     
        
    end

% LU�ڲ�����,sLU����order�仯����
if isSameCol(sLU)
    LU = reorderStruct(LU.order, sLU);
else
    error('����ʹ��structfun');
end




% Itemȥ��δʹ�� %     Item.Rotaed(:,Item.itemorder) = sLU.Rotaed;
% ���ITEM������ȫ����ͬ
if isSameCol(Item)
    Item = structfun(@(x) x( : , Item.LWH(1,:)>0 ), Item, 'UniformOutput', false);
else
    error('����ʹ��structfun');
end

% ITEM�����ж��Ƿ��������ص��ж�Item.isWeightFine
Item = isWeightUpDown(Item,LU);
% ��������������ص�case, �����޸�
if ~all(Item.isWeightFine)
    [~,b] = find(Item.isWeightFine == 0);
    for i=1:length(b)
        LU = repairItemWeight(LU,b(i));
    end
end
Item = isWeightUpDown(Item,LU);
if ~all(Item.isWeightFine),   error('������������casse, ����'); end




% �ɻ�ϵ�LU.DOC����ITEM�ڰ�����PID,LID,SID������ 1808���� ��ʱδ��
LU.DOC=[LU.PID;LU.ID;LU.SID;zeros(size(LU.ID));zeros(size(LU.ID));...
    LU.LU_Item;];
nItem = size(Item.LWH,2);
for iItem=1:nItem
    tmp = LU.DOC([1,2,3], LU.DOC(6,:) == iItem);
    Item.PID(:,iItem) = num2cell(unique(tmp(1,:))',1); %unique(tmp(1,:))';
    Item.LID(:,iItem) = num2cell(unique(tmp(2,:))',1);
    Item.SID(:,iItem) =num2cell(unique(tmp(3,:))',1);
end

% ������� ItemID
% ItemID = getITEMIDArray(Item);
ItemID = [];
%% ����script TO BE FIX
% �����Ҫ���:���ÿ��item������ ԭʼ LU���
printscript(LU,Item);

end

% �ж�LU�Ƿ��������ع���
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

% ��LU�������ع��ɽ����޸�
function LU = repairItemWeight(LU,itemIdx)
    [~,LUidx] = find(LU.LU_Item(1,:)==itemIdx); %�ҳ���item��Ӧ��lu��index
    LU.LU_Item(:,LUidx)
    nbLUinItem = length(LUidx);
    currLUWeight = zeros(1,nbLUinItem);
    for iIdx = 1:nbLUinItem
        currIdx = LUidx(LU.LU_Item(2,LUidx) == iIdx);
        currLUWeight(iIdx) = LU.Weight(:,currIdx);
    end
    % ��Item�ڵ�LU����������b; ������˳��LU_Item��˳��2��������
    [~,b] = sort(currLUWeight,'descend');
    tt = LU.LU_Item(2,LUidx);
    LU.LU_Item(2,LUidx) = tt(b);
end

function printscript(LU,Item)
    for iItem = 1:max(LU.LU_Item(1,:))
%         [~,idx] = find(LU.LU_Item(1,:)==iItem);
%         fprintf('item %d �ĳ����Ϊ:  ',iItem);
%         fprintf('( %d ) ',Item.LWH(:,iItem));
%         fprintf('\n');
%         fprintf('item %d ���� original LU ������(�����)Ϊ  \n  ',iItem);
%         fprintf('%d ',idx);
%         fprintf('( %d ) ', LU.LWH(:,idx));
%         fprintf('\n');
%         fprintf('item %d ���� original LU ������(��)Ϊ  \n  ',iItem);
%         fprintf('%d ',idx);
%         fprintf('( %d ) ', LU.LWH(3,idx)); 
%         fprintf('\n');
%         fprintf('item %d ���� original LU ����Ϊ  \n  ',iItem);
%         fprintf('%d ',idx);
%         fprintf('( %d ) ', LU.Weight(:,idx));
%         fprintf('\n');
%                fprintf('item %d ���� original LU ***Ϊ  \n  ',iItem);
%         fprintf('%d ',idx);
%         fprintf('( %d ) ', LU.LU_Item(2,idx)); 
%         fprintf('\n'); 
%         isWeightUpDown
%         if length(idx) > 1 %Item������ֻһ��Item,��Ҫ�ж��Ƿ������صı仯
%             currLUWeight = zeros(1,length(idx));
%             currLUHight = zeros(1,length(idx));
%             for iIdx = 1:length(idx)
%                currIdx = idx(LU.LU_Item(2,idx) == iIdx);
%                currLUWeight(iIdx) = LU.Weight(:,currIdx);
%                currLUHight(iIdx) = LU.LWH(3,currIdx);      
%             end
%             if diff(currLUWeight) > 0 % ������������
%                 currLUWeight
%             end
%             if diff(currLUHight) >0  
%                 currLUHight
%             end
%         end       

        

    end
end


function [tepLUorder] = getLUorder(LU)
tmpLUMatrix = [LU.SID; LU.ID; LU.PID; LU.LWH; LU.Weight];
% [~,tepLUorder] = sortrows(tmpLUMatrix',[1, 2, 3, 6],{'ascend','ascend','ascend','descend'}); 
% ��Ӧ��; ���ȣ� ID��PID���߶ȣ�������
[~,tepLUorder] = sortrows(tmpLUMatrix',[1, 5, 2, 3, 6, 7],{'ascend','descend','ascend','ascend','descend','descend'}); 
% [~,tepLUorder] = sortrows(tmpLUMatrix',[1, 5, 2, 3, 7, 6],{'ascend','descend','ascend','ascend','descend','descend'}); 

% tmpLUMatrix = [LU.ID; LU.LWH; LU.SID; LU.PID];
% [~,tepLUorder] = sortrows(tmpLUMatrix',[5, 1, 6, 4],{'ascend','ascend','ascend','descend'}); %5:SID; 1:ID 4:Hight
%         tepLUorder = 1:length(LU.ID)'; %ֱ�Ӹ�ֵ1:n % tepLUorder = [2 3 4 1 5]';
if ~isrow(tepLUorder)
    tepLUorder = tepLUorder';
end
end


% ��LUת��ΪItem����Ҫ����
% % function [Item,LU_Item] = getItem(sLU,Veh)
% % %% ��ʼ��
% % % nDim LUά�� nLU LU���� nItem Item���� nLUid LU����
% % % heightBin Bin���߶�
% % 
% % sz = size(sLU.ID);
% % 
% % nLU = sz(2);
% % nLUid = size(unique(sLU.ID),2);
% % 
% % hVeh  = Veh.LWH(3,1);  % tmpUniqueBin = unique(Veh.LWH(1:3,:)','rows')'; % hVeh = tmpUniqueBin(3);
% % 
% % % �����ʼ����Ҫ������fields
% %     % Item.LID = zeros(sz);             %Item��ID����
% %     % Item.SID = zeros(sz);           
% %     % Item.UID = zeros(sz);           
% %     % Item.isRota = ones(sz)*2;    %Item�Ŀ���ת����(��ʼΪ2)
% % Item.Weight = zeros(sz);     %Item������
% % Item.LWH = zeros(3,sz(2));  %Item�Ŀ���
% % 
% % LU_Item = zeros(2,sz(2));     %dim1:���ڵڼ���Item dim2:���ڸ�Item�ڼ����ŷ�
% % 
% % iItem = 1;
% % Item_LU = zeros(sz);  % ÿ��Item�ڶѶ��LU���� ���ڲ���
% % for iLUid=1:nLUid
% %     hLeft = hVeh;
% %     for iLU=1:nLU
% %         if sLU.ID(iLU) == iLUid %���Ե�ǰLU��ӦLuid�ڸ�iLUid�ڵĽ��в���
% %             if hLeft < sLU.LWH(3,iLU) %�統ǰLU�߶Ȳ�����(�߶��ڵ�nDim��)
% %                 iItem =  iItem + 1;
% %                 hLeft = hVeh; 
% %             end
% %             
% %             hLeft = hLeft - sLU.LWH(3,iLU);                 %����ʣ��߶�
% %             Item.LWH(1:2,iItem) = sLU.LWH(1:2,iLU);  %����item����
% %             Item.LWH(3,iItem) = Item.LWH(3,iItem) + sLU.LWH(3,iLU); %����item�߶�
% %             Item.Weight(1,iItem) = Item.Weight(1,iItem) + sLU.Weight(1,iLU); %����item����
% %             
% %             Item_LU(iItem) = Item_LU(iItem) + 1;
% %             LU_Item(1,iLU) = iItem;
% %             LU_Item(2,iLU) = Item_LU(iItem);
% %             
% %             Item.LID(1,iItem) = sLU.ID(1,iLU);               %����ID����
% %             Item.SID(1,iItem) = sLU.SID(1,iLU);               
% %             Item.UID(1,iItem) = sLU.UID(1,iLU);               
% %             Item.isRota(1,iItem) = sLU.isRota(1,iLU);  %����ID����ת����
% %             Item.Roated(1,iItem) = sLU.Rotaed(1,iLU);  %����ID��ת���
% %         end
% %     end
% %     iItem =  iItem + 1;
% % end
% % 
% % end

%%  ��ȡITEMID�����������(ͬ����ID������������������item�Ƿ����ת)
% % function ItemID = getITEMIDArray(Item)
% % 
% % ItemID.ID = unique(Item.LID);
% % nItemID = numel(ItemID.ID);
% % 
% % for iID = 1:nItemID
% %     ItemID.Weight(iID) = sum(Item.Weight .* (Item.LID == ItemID.ID(iID)) );
% %     ItemID.Volume(iID) = sum(prod(Item.LWH) .* (Item.LID == ItemID.ID(iID)) );
% %     ItemID.Area(iID) = sum(prod(Item.LWH(1:2,:)) .* (Item.LID == ItemID.ID(iID)) );
% %     ItemID.isRota(iID) =  unique(Item.isRota(Item.LID == ItemID.ID(iID))); if ~isscalar(ItemID.isRota(iID)), error('��������'); end
% % end
% % 
% % end

