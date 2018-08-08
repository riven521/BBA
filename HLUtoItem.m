function [LU,Item,ItemID] = HLUtoItem(LU,Veh)
% ��Ҫ����:LU�Ѷ���γ�Item %  ����:�����(row);  ����:��������(coloum);
% Input ---  LU: ID LWH Weight ��LU: ����ԭ��˳��
% Output --- LU: order LU_Item ��LU: ����ԭ��˳��(ORDER�ǽ���Item�㷨��LU˳��)
% Output --- Item: ID LWH Weight ...��ITEM:û��˳���㷨������˳��

% LU.LU_Item           (2,n) : ��1: LU�ڵڼ���item ��2:LU�����item��˳��(��-��))
% LU.order                (1,n):  LU����˳��)
% Item.ID                ��1,n): ITEM������ ��ͬ�ڲ�LU����)
% Item.LWH              (3,n): ITEM�ĳ����-������LU��ͬ,�߶�Ϊ�Ѷ��߶�)
% Item.Weight          (1,n): ITEM������

%% LU����
% ��ȡLU��˳��(�ص��Ǹ߶ȵݼ�����)
LU.order = getLUorder(LU); %��ȡ LU����(��ID����,��߶ȵݼ�)
% printstruct(LU)
% ��ȡ��order������LU:sLU
sLU = structfun(@(x) x(:,LU.order),LU,'UniformOutput',false);
% printstruct(sLU)

%% 55 LU->Itemת��

% ������sLU, �����Ѷ��ȡ��Item,�Լ�sLU��Item�ڵ�˳��
sz = size(sLU.ID);
nLU = sz(2);
hVeh  = Veh.LWH(3,1);  % tmpUniqueBin = unique(Veh.LWH(1:3,:)','rows')'; % hVeh = tmpUniqueBin(3);

Item.ID = zeros(sz);             %Item��ID����
Item.SID = zeros(sz);           
% Item.UID = zeros(sz);           
Item.isRota = ones(sz)*2;    %Item�Ŀ���ת����(��ʼΪ2)

Item.LWH = zeros(3,nLU); % Item.LWH(1,:) = wStrip;   %dim1-���ʣ��  % Item.LWH(3,:) = hVeh; % 
Item.Weight = zeros(1,nLU); %Item������
Item.Item_LU = zeros(2,nLU);  % ��1��ÿ��Strip�ڵ�Item���� �� ��2��ÿ��Strip�ڵĲ�ͬLUID����
sLU_Item = zeros(2,sz(2));     %dim1:���ڵڼ���Item dim2:���ڸ�Item�ڼ����ŷ�

iItem = 1; iLU = 1; %iStrip����itemʵ��
while 1
    if iLU > nLU, break; end
    [thisItem,iItem] = getThisItem(iItem);
    insertLUToItem(thisItem,iLU);       
    iLU = iLU + 1;
end

% Get ITEM ��ؿ��Է�
    function [thisItem,iItem] = getThisItem(iItem)
        % ͬ��SID/UID ͬ��LUID Item�߶����� δ����Weight��
        isflagCurr =hVeh - Item.LWH(3,iItem) >= sLU.LWH(3,iLU); %�ж��Ƿ�current's itemʣ���� >= ��ǰiLU�߶�
        isSameID = Item.ID(iItem) == sLU.ID(iLU); %�ж�Item�ڲ�ID�Ƿ�=��ǰiLU��ID
        isNewItem = Item.LWH(3,iItem) == 0; % �ж��Ƿ� new Item �߶�==0            
        
        if isNewItem
                thisItem = iItem;
        else
            if isflagCurr && isSameID %����߶����� �� ID��ͬ
                 thisItem = iItem;
            else
                iItem = iItem + 1;
                [thisItem,iItem] = getThisItem(iItem);
            end
        end
    end

% Put LU into thisItem
    function insertLUToItem(thisItem,iLU)
        %����Height
        Item.LWH(3,iItem) = Item.LWH(3,iItem)  + sLU.LWH(3,iLU);
        Item.LWH(1:2,iItem) = sLU.LWH(1:2,iLU);  %����item����
        Item.Weight(1,iItem) = Item.Weight(1,iItem) + sLU.Weight(1,iLU); %����item����
        
        Item.Item_LU(1,iItem) = Item.Item_LU(1,iItem) + 1;
        
        sLU_Item(1,iLU) = iItem;
        sLU_Item(2,iLU) = Item.Item_LU(1,iItem);
        
        tmpLUThisItem = sLU_Item(1,:) == thisItem;
        Item.Item_LU(2,iItem) = numel(unique(sLU.PID(1,tmpLUThisItem)));
        
        Item.ID(1,iItem) = sLU.ID(1,iLU);               %����ID����
        Item.SID(1,iItem) = sLU.SID(1,iLU);
        Item.UID(1,iItem) = sLU.UID(1,iLU);
        Item.isRota(1,iItem) = sLU.isRota(1,iLU);  %����ID����ת����
        Item.Rotaed(1,iItem) = sLU.Rotaed(1,iLU);  %����ID��ת���
    end

% LU�ڲ�����
    LU.LU_Item(:,LU.order) = sLU_Item;
    
% Itemȥ��δʹ�� %     Item.Rotaed(:,Item.itemorder) = sLU.Rotaed;
    Item = structfun(@(x) x( : , Item.ID(1,:)>0 ), Item, 'UniformOutput', false);

% ������� ItemID
% ItemID = getITEMIDArray(Item);
 ItemID = [];
%% ����script TO BE FIX
% �����Ҫ���:���ÿ��item������ ԭʼ LU���
printscript(LU,Item);

end

function printscript(LU,Item)
    for iItem = 1:max(LU.LU_Item(1,:))
        [~,idx] = find(LU.LU_Item(1,:)==iItem);
        fprintf('item %d �ĳ����Ϊ:  ',iItem);
        fprintf('( %d ) ',Item.LWH(:,iItem));
        fprintf('\n');
        fprintf('item %d ���� original LU ������(�����)Ϊ  \n  ',iItem);
        fprintf('%d ',idx);
        fprintf('( %d ) ', LU.LWH(:,idx));
        fprintf('\n');
    end
end

function order = getLUorder(LU)
tmpLUMatrix = [LU.SID; LU.ID; LU.PID; LU.LWH];
% [~,tepLUorder] = sortrows(tmpLUMatrix',[1, 2, 3, 6],{'ascend','ascend','ascend','descend'}); 
[~,tepLUorder] = sortrows(tmpLUMatrix',[1, 5, 2, 3, 6],{'ascend','descend','ascend','ascend','descend'}); 

% tmpLUMatrix = [LU.ID; LU.LWH; LU.SID; LU.PID];
% [~,tepLUorder] = sortrows(tmpLUMatrix',[5, 1, 6, 4],{'ascend','ascend','ascend','descend'}); %5:SID; 1:ID 4:Hight
%         tepLUorder = 1:length(LU.ID)'; %ֱ�Ӹ�ֵ1:n % tepLUorder = [2 3 4 1 5]';
if ~isrow(tepLUorder)
    order = tepLUorder';
end
end


% ��LUת��ΪItem����Ҫ����
function [Item,LU_Item] = getItem(sLU,Veh)
%% ��ʼ��
% nDim LUά�� nLU LU���� nItem Item���� nLUid LU����
% heightBin Bin���߶�

sz = size(sLU.ID);

nLU = sz(2);
nLUid = size(unique(sLU.ID),2);

hVeh  = Veh.LWH(3,1);  % tmpUniqueBin = unique(Veh.LWH(1:3,:)','rows')'; % hVeh = tmpUniqueBin(3);

Item.ID = zeros(sz);             %Item��ID����
Item.SID = zeros(sz);           
Item.UID = zeros(sz);           
Item.isRota = ones(sz)*2;    %Item�Ŀ���ת����(��ʼΪ2)
Item.Weight = zeros(sz);     %Item������
Item.LWH = zeros(3,sz(2));  %Item�Ŀ���

LU_Item = zeros(2,sz(2));     %dim1:���ڵڼ���Item dim2:���ڸ�Item�ڼ����ŷ�

iItem = 1;
Item_LU = zeros(sz);  % ÿ��Item�ڶѶ��LU���� ���ڲ���
for iLUid=1:nLUid
    hLeft = hVeh;
    for iLU=1:nLU
        if sLU.ID(iLU) == iLUid %���Ե�ǰLU��ӦLuid�ڸ�iLUid�ڵĽ��в���
            if hLeft < sLU.LWH(3,iLU) %�統ǰLU�߶Ȳ�����(�߶��ڵ�nDim��)
                iItem =  iItem + 1;
                hLeft = hVeh; 
            end
            
            hLeft = hLeft - sLU.LWH(3,iLU);                 %����ʣ��߶�
            Item.LWH(1:2,iItem) = sLU.LWH(1:2,iLU);  %����item����
            Item.LWH(3,iItem) = Item.LWH(3,iItem) + sLU.LWH(3,iLU); %����item�߶�
            Item.Weight(1,iItem) = Item.Weight(1,iItem) + sLU.Weight(1,iLU); %����item����
            
            Item_LU(iItem) = Item_LU(iItem) + 1;
            LU_Item(1,iLU) = iItem;
            LU_Item(2,iLU) = Item_LU(iItem);
            
            Item.ID(1,iItem) = sLU.ID(1,iLU);               %����ID����
            Item.SID(1,iItem) = sLU.SID(1,iLU);               
            Item.UID(1,iItem) = sLU.UID(1,iLU);               
            Item.isRota(1,iItem) = sLU.isRota(1,iLU);  %����ID����ת����
        end
    end
    iItem =  iItem + 1;
end

end

%%  ��ȡITEMID�����������(ͬ����ID������������������item�Ƿ����ת)
function ItemID = getITEMIDArray(Item)

ItemID.ID = unique(Item.ID);
nItemID = numel(ItemID.ID);

for iID = 1:nItemID
    ItemID.Weight(iID) = sum(Item.Weight .* (Item.ID == ItemID.ID(iID)) );
    ItemID.Volume(iID) = sum(prod(Item.LWH) .* (Item.ID == ItemID.ID(iID)) );
    ItemID.Area(iID) = sum(prod(Item.LWH(1:2,:)) .* (Item.ID == ItemID.ID(iID)) );
    ItemID.isRota(iID) =  unique(Item.isRota(Item.ID == ItemID.ID(iID))); if ~isscalar(ItemID.isRota(iID)), error('��������'); end
end

end

