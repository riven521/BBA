%% GET STRIP �������
% 1 Strip.isMixed % 1����ϲ㣻 0��������
% 2 Strip.isHeightFull % 1��ȫ���������㣻 0������������ ����Item.isHeightFull �ж� 
% 3 Strip.isWidthFull % % 1��ȫ�����ǿ�����㣻 0����ȷ����� ���� Strip���ʣ�� �Ƿ� ������ ��С��Item����ж� 
% 4 Strip.maxHeight % Strip������Item�����߶�

% 5 Strip.nbItem % ����������ֵ, ����ITEM�ĶѶ���� ��ͷ�ڷ����� -1�����strip
% 6 Strip.isAllPured % ��1�������Ҹ�ITEMû�л���ͣ� 0������������strip�л�ϵģ� -1�����strip
% 7 Strip.isSingleItem %: 1:�����ҵ����� 0�������Ҷ����-1�����strip

% 8 Strip.loadingrate   % ÿ��strip��װ�ر���
% 9 Strip.loadingrateLimit   % ÿ��strip������װ�ر���

% Strip.seqSW
% Strip.LID

%% ����
function   [Strip,LU] = cpuStrip(Strip,Item,LU,Veh)
%% ��ʼ��
    Strip.isMixed = ones(size(Strip.Weight))*-1;   %�Ƿ�Ϊ�����,�������LID 
    Strip.isHeightFull = ones(size(Strip.Weight))*-1;   %�Ƿ������Full��Item
    Strip.nbItem = ones(size(Strip.Weight))*-1;   %��STRIP�ڲ�ITEM���͸���, �����Ĭ��Ϊ-1
    Strip.isAllPured = ones(size(Strip.Weight))*-1;   %��STRIP��ӦLID�Ƿ�������STRIP, ���������Ĭ��Ϊ-1
    Strip.isSingleItem = ones(size(Strip.Weight))*-1;   %��Strip�ڶ�Ӧֻ��1��ITEM
    Strip.isWidthFull = ones(size(Strip.Weight))*-1;     %�Ƿ�Ϊ��ȷ�Full��Item
    Strip.maxHeight = ones(size(Strip.Weight))*-1;     %Strip����߸߶�.
%     Strip.seqSW = ones(size(Strip.Weight))*-1;     %Strip��??? ��δ��

    Strip.Stripvolume = ones(size(Strip.Weight))*-1;  %ÿ��strip�Ŀ������ = �߶�*���(�����Ŀ��)
    Strip.StripvolumeLimit = ones(size(Strip.Weight))*-1; % %ÿ��strip�����޿������ = �߶�*���(stripʹ�ÿ��=�������-stripʣ����)
    Strip.Itemvolume = ones(size(Strip.Weight))*-1; %ÿ��strip������Itemװ�����
    Strip.loadingrate = ones(size(Strip.Weight))*-1; % ÿ��strip��װ�ر���
    Strip.loadingrateLimit = ones(size(Strip.Weight))*-1;     % ÿ��strip������װ�ر���
    
    

%% 0.0 ����LU_Strip, Strip�ڵ�PID,LID,SID
% �ɻ�ϵ�LU.DOC����LU_STRIP, ����STRIP�ڰ�����PID,LID,SID������ 1808����
nbLU = size(LU.LWH,2);
LU.LU_Strip = [zeros(1,nbLU);zeros(1,nbLU)];
for iLU=1:nbLU
    theItem = LU.LU_Item(1,iLU);   %iLU���ڵڼ���Item
    LU.LU_Strip(1,iLU)= Item.Item_Strip(1,theItem);
end

LU.DOC=[LU.DOC; LU.LU_Strip];
nStrip = size(Strip.LW,2);
for iStrip=1:nStrip
    tmp = LU.DOC([1,2,3], LU.DOC(8,:) == iStrip);
    Strip.PID(:,iStrip) = num2cell(unique(tmp(1,:))',1);
    Strip.LID(:,iStrip) = num2cell(unique(tmp(2,:))',1);
    Strip.SID(:,iStrip) = num2cell(unique(tmp(3,:))',1);
end
    
%% 0: ����stripװ����
Strip = computeLoadingRateStrip(Strip,Item,Veh); 
    
%% 1: STRIP.isMixed: STRIP�����ж��Ƿ񵥴���/������ж�
Strip = isMixedStrip(Strip);

%% 2: STRIP.isHeightFull: STRIP�����ж��Ƿ������Full��Item. 
Strip = isFullStrip(Strip,Item);

%% 3: STRIP.isWidthFull : STRIP�����ж��Ƿ�������width��Full��Item. 
Strip = isWidthFullStrip(Strip,Item);

%% 4: STRIP.maxHeight : ����Strip�����߶�
for i=1:length(Strip.maxHeight)
    % �������ֵ
    % Item.Item_Strip(1,:) == i) : Strip i �ڲ���Item flag
    Strip.maxHeight(i) = max(Item.LWH(3, Item.Item_Strip(1,:) == i));
end

%% 5,6,7
%Strip.isAllPured�����:-1; ����: 1 (���strip��û�и�ID) ; 0 (���strip�ں��иõ���strip��ID)
%Strip.nbItem: ���:-1; ����: ��ӦStrip�ڲ���Item��nbLID���͸���,��ֵԽ��,����LU����Խ��
%Strip.isSingleItem: ���: -1; ����: Strip�ڽ���һ��Item,�ض��ǵ�����.
LIDinItemsArray = cellfun(@(x) x(1), Item.LID); % arrayAllLID: ����ITEM��Ӧ��LIDֵ ������ʽ
mixedStrip = find(Strip.isMixed(1,:) == 1);
mixedLID = [];
for m=1:length(mixedStrip)
     f = Item.Item_Strip(1,:) == mixedStrip(m);
     mixedLID = [mixedLID, LIDinItemsArray(f)];
end
mixedLID = unique(mixedLID);

uniItem = unique(Item.Item_Strip(1,:));
for i=1:length(Strip.isAllPured)
    if ~Strip.isMixed(1,i) %���ǵ�����        
        cellLID = Item.LID(Item.Item_Strip(1,:) == uniItem(i)); % cellLID: ��Strip�ڵ�ITEM��Ӧ��LIDֵ
            %         cellLID = Item.LID(Item.Item_Strip(1,:) == i); % cellLID: ��Strip�ڵ�ITEM��Ӧ��LIDֵ
        LIDinThisItemArray = cellfun(@(x)x(1), cellLID);
        if isscalar(unique(LIDinThisItemArray)) 
            
            if isscalar(LIDinThisItemArray)
                Strip.isSingleItem(1,i) = 1;
            else
                Strip.isSingleItem(1,i) = 0;
            end
            
            if ismember(unique(LIDinThisItemArray),mixedLID)
                Strip.isAllPured(1,i) = 0;
            else
                Strip.isAllPured(1,i) = 1;
            end
            
            Strip.nbItem(1,i) = sum(LIDinItemsArray == unique(LIDinThisItemArray));
            
        else
             error('������STRIP�ڵ�ITEM�����Ͳ�ͬ'); %arrayLID            
        end
    end
end

end

%% �ֲ����� %%

%% ����1: �ж�STRIP�Ƿ������HeightFull��Item
function Strip = isFullStrip(Strip,Item)
    % ѭ���ж�Strip�Ƿ�full
    uniItem = unique(Item.Item_Strip(1,:));
    for i=1:length(Strip.isHeightFull)
%          if all(Item.isHeightFull(Item.Item_Strip(1,:) == i)) %�����STRIP��ӦITEM��isFull��Ϊ1,��STRIPҲΪfull
         if all(Item.isHeightFull(Item.Item_Strip(1,:) == uniItem(i))) %�����STRIP��ӦITEM��isFull��Ϊ1,��STRIPҲΪfull
             Strip.isHeightFull(i) = 1;
         else
             Strip.isHeightFull(i) = 0;
         end
    end
end

%% ����2: �ж�STRIP�Ƿ������WidthFull��Item
% ****************** Strip���Ƿ����Width��Full��Item���� ************ ����
function Strip = isWidthFullStrip(Strip,Item) 
    % ѭ���ж�Strip�Ƿ�Widthfull, ��ȼ�϶ >= ��Strip������Item�Ŀ����Сֵ
    uniItem = unique(Item.Item_Strip(1,:));
    for i=1:length(Strip.isWidthFull)
        flagItem = Item.Item_Strip(1,:) == uniItem(i); 
        max(Item.LWH(1, flagItem))
        if Strip.LW(1, i) >= min(Item.LWH(1, flagItem))
            Strip.isWidthFull(i) = 0;
        else
            Strip.isWidthFull(i) = 1;
        end        
    end
end

%% ����3: �ж�STRIP�Ƿ�����
function Strip = isMixedStrip(Strip)
    % ѭ���ж�Strip�Ƿ�Ϊ�����
    for i=1:length(Strip.isMixed)
         if numel(Strip.LID{i}) > 1
             Strip.isMixed(i) = 1;
         else
             Strip.isMixed(i) = 0;
         end
    end
end

function Strip = computeLoadingRateStrip(Strip,Item,Veh)
    % ��ʼ��
    nStrip = size(Strip.LW,2);

    % ����ÿ��strip��װ����
    %ÿ��strip�Ŀ������ = �߶�*���(�����Ŀ��)
    Strip.Stripvolume = Strip.LW(2,:)*Veh.LWH(1,1);
    %ÿ��strip�����޿������ = �߶�*���(stripʹ�ÿ��=�������-stripʣ����)
    Strip.StripvolumeLimit = Strip.LW(2,:) .* (Veh.LWH(1,1) - Strip.LW(1,:));
    a = Item.LWH;
    b = Item.Item_Strip;
    uniItem = unique(Item.Item_Strip(1,:));
    for iStrip =1:nStrip
        %ÿ��strip������Itemװ�����
        Strip.Itemvolume(iStrip)= sum(a(1, (b(1,:)==uniItem(iStrip))) .* a(2, (b(1,:)==uniItem(iStrip))));
    end
    %ÿ��strip��װ�ر���
    Strip.loadingrate =  Strip.Itemvolume ./ Strip.Stripvolume;
    %ÿ��strip������װ�ر���
    Strip.loadingrateLimit =  Strip.Itemvolume ./ Strip.StripvolumeLimit;
end
