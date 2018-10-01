%% GET STRIP �������
% 1 Strip.isMixed % 1����ϲ㣻 0��������
% 2 Strip.isHeightFull % 1��ȫ���������㣻 0������������ ����Item.isHeightFull �ж� 
% 3 Strip.isWidthFull % % 1��ȫ�����ǿ�����㣻 0����ȷ����� ���� Strip���ʣ�� �Ƿ� ������ ��С��Item����ж� 
% 4 Strip.maxHeight % Strip������Item�����߶�

% 5 Strip.nbItem % ����������ֵ, ����ITEM�ĶѶ���� ��ͷ�ڷ����� -1�����strip
% 6 Strip.isAllPured % ��1�������Ҹ�ITEMû�л���ͣ� 0��������Ҳ�ɻ�ϵģ� -1�����strip
% 7 Strip.isSingleItem %: 1:�����ҵ����� 0�������Ҷ����-1�����strip

% 8 Strip.loadingrate   % ÿ��strip��װ�ر���
% 9 Strip.loadingrateLimit   % ÿ��strip������װ�ر���

% Strip.seqSW
% Strip.LID

%% ����
function   [Strip] = cpuStrip(Strip,Item,Veh)
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
%Strip.isAllPured����STRIP��ӦLID�Ƿ�������STRIP, ���������Ĭ��Ϊ-1
%Strip.nbItem: STRIP�����ڲ���ITEM��nbLID���͸���,��ֵԽ��,����LU����Խ��
%Strip.isSingleItem: Strip
arrayAllLID = cellfun(@(x) x(1), Item.LID); % arrayAllLID: ����ITEM��Ӧ��LIDֵ ������ʽ
mixedStrip = find(Strip.isMixed(1,:) == 1);
mixedLID = [];
for m=1:length(mixedStrip)
     flagitemIdx = Item.Item_Strip(1,:) == mixedStrip(m);
     mixedLID = [mixedLID, arrayAllLID(flagitemIdx)];
end
mixedLID = unique(mixedLID);
for i=1:length(Strip.isAllPured)
    if ~Strip.isMixed(1,i) %���ǵ�����        
        cellLID = Item.LID(Item.Item_Strip(1,:) == i); % cellLID: ��Strip�ڵ�ITEM��Ӧ��LIDֵ
        arrayLID = cellfun(@(x)x(1), cellLID);
        if isscalar(unique(arrayLID)) 
            if isscalar(arrayLID)
                Strip.isSingleItem(1,i) = 1;
            else
                Strip.isSingleItem(1,i) = 0;
            end
            if ismember(unique(arrayLID),mixedLID)
                Strip.isAllPured(1,i) = 0;
            else
                Strip.isAllPured(1,i) = 1;
            end
            Strip.nbItem(1,i) = sum(arrayAllLID == unique(arrayLID));
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
    for i=1:length(Strip.isHeightFull)
         if all(Item.isHeightFull(Item.Item_Strip(1,:) == i)) %�����STRIP��ӦITEM��isFull��Ϊ1,��STRIPҲΪfull
             Strip.isHeightFull(i) = 1;
         else
             Strip.isHeightFull(i) = 0;
         end
    end    
end

%% ����2: �ж�STRIP�Ƿ������WightFull��Item
% ****************** Strip���Ƿ����Width��Full��Item���� ************ ����
function Strip = isWidthFullStrip(Strip,Item) 
    % ѭ���ж�Strip�Ƿ�Widthfull, ��ȼ�϶ >= ��Strip������Item�Ŀ����Сֵ
    for i=1:length(Strip.isWidthFull)
        flagItem = Item.Item_Strip(1,:) == i;        
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
    for iStrip =1:nStrip
        %ÿ��strip������Itemװ�����
        Strip.Itemvolume(iStrip)= sum(a(1, (b(1,:)==iStrip)) .* a(2, (b(1,:)==iStrip)));
    end
    %ÿ��strip��װ�ر���
    Strip.loadingrate =  Strip.Itemvolume ./ Strip.Stripvolume;
    %ÿ��strip������װ�ر���
    Strip.loadingrateLimit =  Strip.Itemvolume ./ Strip.StripvolumeLimit;
end
