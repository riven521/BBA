%% GET STRIP �������
% LU.LU_Strip, LU.CoordLUStrip
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
    Strip.isMixedSID = ones(size(Strip.Weight))*-1;   %�Ƿ�Ϊ�����,�������SID 
    Strip.isMixedEID = ones(size(Strip.Weight))*-1;   %�Ƿ�Ϊ�����,�������EID 
    
    Strip.nbLID = ones(size(Strip.Weight))*-1;   %STRIP��LID����
    Strip.nbSID= ones(size(Strip.Weight))*-1;   %SID����
    Strip.nbEID= ones(size(Strip.Weight))*-1;   %EID����
    
    
    Strip.GapValue = ones(size(Strip.Weight))*-1;         %���STRIP��isMixed==1������Ȳ���STRIP����ʱ��С����Ϊ0���ڵ�GAPֵ
    Strip.isGapBalance = ones(size(Strip.Weight))*-1;   %GAPֵ�����ITEM��1/3��Ϊ��balance��
    
    Strip.isHeightFull = ones(size(Strip.Weight))*-1;   %�Ƿ������Full��Item
    Strip.isHeightBalance = ones(size(Strip.Weight))*-1;   %�Ƿ�Item�ĸ߶Ȳ��첻��
    
    Strip.nbItem = ones(size(Strip.Weight))*-1;            %��STRIP�ڲ�ITEM���͸���, �����Ĭ��Ϊ-1
    Strip.isAllPured = ones(size(Strip.Weight))*-1;       %��STRIP��ӦLID�Ƿ�������STRIP, ���������Ĭ��Ϊ-1
    Strip.isSingleItem = ones(size(Strip.Weight))*-1;   %��Strip�ڶ�Ӧֻ��1��ITEM
    Strip.isWidthFull = ones(size(Strip.Weight))*-1;     %�Ƿ�Ϊ��ȷ�Full��Item  �� ��ȼ�϶ >= ��Strip������Item�Ŀ����Сֵ�� ���Ȳ�����ȡֵ0
    
    Strip.maxHeight = ones(size(Strip.Weight))*-1;     %Strip����߸߶�.
    Strip.lowestHeight = ones(size(Strip.Weight))*-1;     %Strip����͸߶�.
    Strip.meanHeight = ones(size(Strip.Weight))*-1;     %Strip����͸߶�.
    
    Strip.maxLULength= ones(size(Strip.Weight))*-1;     %Strip����߸߶�.
    Strip.lowestLULength = ones(size(Strip.Weight))*-1;     %Strip����͸߶�.
%     Strip.meanLUWidth = ones(size(Strip.Weight))*-1;     %Strip����͸߶�.
%     Strip.seqSW = ones(size(Strip.Weight))*-1;     %Strip��??? ��δ��

    Strip.Stripvolume = ones(size(Strip.Weight))*-1;  %ÿ��strip�Ŀ������ = �߶�*���(�����Ŀ��)
    Strip.StripvolumeLimit = ones(size(Strip.Weight))*-1; % %ÿ��strip�����޿������ = �߶�*���(stripʹ�ÿ��=�������-stripʣ����)
    Strip.Itemvolume = ones(size(Strip.Weight))*-1; %ÿ��strip������Itemװ�����
    Strip.loadingrate = ones(size(Strip.Weight))*-1; % ÿ��strip��װ�ر���
    Strip.loadingrateLimit = ones(size(Strip.Weight))*-1;     % ÿ��strip������װ�ر���

%% 0.0 ����LU.LU_Strip, LU.CoordLUStrip
% �ɻ�ϵ�LU.DOC����LU_STRIP, ����STRIP�ڰ�����PID,LID,SID������ 1808����
nbLU = size(LU.LWH,2);
LU.LU_Strip = zeros(2,nbLU);
LU.CoordLUStrip = zeros(3,nbLU);

% ����LU_Strip
for iLU=1:nbLU
    % ����LU_Strip��һ��
    iItem = LU.LU_Item(1,iLU);   %iLU���ڵڼ���Item, Item���ڵڼ���Strip,��Lu���ڵڼ���Strip
    LU.LU_Strip(1,iLU)= Item.Item_Strip(1,iItem);
    % ����LU_Strip�ڶ���
    fiItem = find(Item.Item_Strip(1,:) == Item.Item_Strip(1,iItem) & Item.Item_Strip(2,:) < Item.Item_Strip(2,iItem));
    nbLUfiItem = sum(ismember(LU.LU_Item(1,:),fiItem));
    LU.LU_Strip(2,iLU) = nbLUfiItem+LU.LU_Item(2,iLU); % ����Strip˳��: ͬһStrip����ǰ�������nbLUfiItem + ��iLU��Item��˳��
    % ����LU.CoordLUStrip
    LU.CoordLUStrip(1,iLU) = Item.CoordItemStrip(1,iItem);
    LU.CoordLUStrip(2,iLU) = Item.CoordItemStrip(2,iItem);
        % fLU: ��iLUͬ��iItem �� ˳�����ڱ�iLU; ����Ϊ��, ��Ӱ��.
    fLU = LU.LU_Item(1,:) == iItem & LU.LU_Item(2,:) < LU.LU_Item(2,iLU);
    LU.CoordLUStrip(3,iLU) = sum(LU.LWH(3,fLU));
end

%% V2
    %% SECTION 0 ����STRIP��PID,LID,SID,��TABLE����,����֪��ʲô��ʲô,����1,2,3�����滻
    t = struct2table(structfun(@(x) x',LU,'UniformOutput',false));
    
    nStrip = size(Strip.LW,2);
    for iStrip=1:nStrip
        f = t.LU_Strip(:,1) == iStrip;   
        Strip.LID(:,iStrip) = {unique(t.ID(f))};           % NOTE: STRIP���LID��LU��ID
        %         Item.LID(:,iItem) = {unique(t.LID(f))};
        Strip.SID(:,iStrip) = {unique(t.SID(f))};
        Strip.EID(:,iStrip) = {unique(t.EID(f))};
        Strip.PID(:,iStrip) = {unique(t.PID(f))};
    end    
    %  t2 = struct2table(structfun(@(x) x',Strip,'UniformOutput',false));
    
%% V1
% % LU.DOC=[LU.DOC; LU.LU_Strip];
% % nStrip = size(Strip.LW,2);
% % for iStrip=1:nStrip
% %     tmp = LU.DOC([1,2,3], LU.DOC(8,:) == iStrip);
% %     Strip.PID(:,iStrip) = num2cell(unique(tmp(1,:))',1);
% %     Strip.LID(:,iStrip) = num2cell(unique(tmp(2,:))',1);
% %     Strip.SID(:,iStrip) = num2cell(unique(tmp(3,:))',1);
% % end
    
%% 0: ����stripװ����
Strip = computeLoadingRateStrip(Strip,Item,Veh); 
    
%% 1: STRIP.isMixed: STRIP�����ж��Ƿ񵥴���/������ж�
Strip = isMixedStrip(Strip);

%% 2: STRIP.isHeightFull: STRIP�����ж��Ƿ������Full��Item. 
Strip = isHeightBalanceStrip(Strip,Item);
Strip = isFullStrip(Strip,Item);

%% 3: STRIP.isWidthFull : STRIP�����ж��Ƿ�������width��Full��Item. 
Strip = isWidthFullStrip(Strip,Item);

%% 4: STRIP.maxHeight/maxLULength : ����Strip�ĸ߶ȺͿ��
for i=1:length(Strip.maxHeight)
    % �����ڲ��߶�
    % Item.Item_Strip(1,:) == i) : Strip i �ڲ���Item flag
    Strip.maxHeight(i) = max(Item.LWH(3, Item.Item_Strip(1,:) == i));
    Strip.lowestHeight(i) = min(Item.LWH(3, Item.Item_Strip(1,:) == i));
    Strip.meanHeight(i) = Strip.maxHeight(i) - Strip.lowestHeight(i); 
    Strip.meanHeight(i) = (Strip.maxHeight(i) + Strip.lowestHeight(i))/2; 
    if any(Strip.meanHeight(i)<0), error('��ֵ������;'); end
    % �����ڲ�LU���
    Strip.maxLULength(i) = max(Item.LWH(2, Item.Item_Strip(1,:) == i));
    Strip.lowestLULength(i) = min(Item.LWH(2, Item.Item_Strip(1,:) == i));
end

%% 4.1 STRIP.GapValue �� isGAPBalance�ļ���

fmix = Strip.isMixed==1;
fwid = Strip.isWidthFull==0;

% StripCheck �ɶ�ε���
chkStrip(Strip);

% Calc GapValue
Strip.GapValue(~fmix) =0;
Strip.GapValue(fmix) = Strip.maxLULength(fmix) - Strip.lowestLULength(fmix);
Strip.GapValue(fwid) = Strip.maxLULength(fwid) - 0;

% Calc isGapBalance
Strip.isGapBalance = Strip.GapValue < 0.33*Strip.maxLULength;

% fgap=Strip.isGapBalance ==0;
% Strip.LID
% printstruct(Strip)

%% 5,6
%Strip.isAllPured�����:-1; ����: 1 (���strip��û�и�ID) ; 0 (���strip�ں��иõ���strip��ID)
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
            
            % Strip.nbItem(1,i) = sum(LIDinItemsArray ==
            % unique(LIDinThisItemArray)); %�������뺯������
            
        else
             error('������STRIP�ڵ�ITEM�����Ͳ�ͬ'); %arrayLID            
        end
    else % strip���ǻ����     
        Strip.isSingleItem(1,i) = 0;
        Strip.isAllPured(1,i) = 0;
    end
end

%% 7
%Strip.nbItem: ���:-; ����: ��ӦStrip�ڲ���Item��nbLID���͸���,��ֵԽ��,����LU����Խ��
[Strip.nbItem, Strip.nbLU, Strip.nbLULID] = cpuStripnbItem(Strip,Item,LU);

end

%% �ֲ����� %%

%% ����1: �ж�STRIP�Ƿ������HeightFull��Item
function Strip = isHeightBalanceStrip(Strip,Item) 
global parBalance
    % 1 ѭ���ж�Strip�Ƿ����ItemΪfull��,�����,��StripΪfull
    uniStrip = unique(Item.Item_Strip(1,:));

    % 2 ѭ���ж�Strip�ڲ�Item֮�������ֵ, �Ƿ�<= ��С�ĶԽ��߻�һ������ֵ, ����,��ΪFull; 
    for i=1:length(Strip.isHeightFull)
        fItem = Item.Item_Strip(1,:) == uniStrip(i);        
        % �߶ȼ�϶
        maxHeightDiff = max(abs(diff(Item.LWH(3,fItem))));
        % �Ա�2: ��ߵ�Item��1/3
        oneThirdsHeightItem = max(Item.LWH(3,fItem))*parBalance;
        
        % �Ա�2: ����ֵ, �Ƿ�<= 1/3���Item (Item�߶�ƽ��, ��ʹ�ܵ�,  Ҳ��Ϊ������)
        if ~isempty(maxHeightDiff)
            if maxHeightDiff <= oneThirdsHeightItem                 
                     Strip.isHeightBalance(i) = 1;
            else
                     Strip.isHeightBalance(i) = 0;       
            end
        else %��Ϊ��ֵ, ��ֻ��Stripһ���Ѷ�, ���Ǿ����
            if sum(fItem)~=1, error('��strip��������һ���Ѷ�!'); end
            Strip.isHeightBalance(i) = 1;
        end
    end
    
    % �������
    if any(Strip.isHeightBalance==-1), error('����Strip.isHeightBalanceδ����!'); end
end

%% V1 ����������isHeightBalance�߶Ⱦ���, ��ȫ
% function Strip = isFullStrip(Strip,Item)
%     % 1 ѭ���ж�Strip�Ƿ����ItemΪfull��,�����,��StripΪfull
%     uniStrip = unique(Item.Item_Strip(1,:));
% %     for i=1:length(Strip.isHeightFull)
% %         %          if all(Item.isHeightFull(Item.Item_Strip(1,:) == i)) %�����STRIP��ӦITEM��isFull��Ϊ1,��STRIPҲΪfull
% %          if all(Item.isHeightFull(Item.Item_Strip(1,:) == uniStrip(i))) %�����STRIP��ӦITEM��isFull��Ϊ1,��STRIPҲΪfull
% %              Strip.isHeightFull(i) = 1;
% %          else
% %              Strip.isHeightFull(i) = 0;
% %          end
% %     end
%     
% 
%     % 2 ѭ���ж�Strip�ڲ�Item֮�������ֵ, �Ƿ�<= ��С�ĶԽ��߻�һ������ֵ, ����,��ΪFull; 
%     for i=1:length(Strip.isHeightFull)
%         fItem = Item.Item_Strip(1,:) == uniStrip(i);        
%         diagItem = sqrt(Item.LWH(1,fItem).^2 + Item.LWH(2,fItem).^2);
%         % �߶ȼ�϶
%         maxHeightDiff = max(abs(diff(Item.LWH(3,fItem))));
%         
%         % �Ա�1: ��С��Item�Խ���
%         minDiagItem = min(diagItem);
%         % �Ա�2: ��ߵ�Item��1/3
%         oneThirdsHeightItem = max(Item.LWH(3,fItem))*1/3
%         % �Ա�3: ����ֵ
%         absHeight = 300;
%         
%         % �Ա�2: ����ֵ, �Ƿ�<= 1/3���Item (Item�߶�ƽ��, ��ʹ�ܵ�,  Ҳ��Ϊ������)
%         if ~isempty(maxHeightDiff)
%             if maxHeightDiff <= oneThirdsHeightItem                 
%                      Strip.isHeightFull(i) = 1;
%                      % 3 ���Ӽ�ʹmaxHeightDiff��С, ������߶ȵ�, Ҳ��Ϊ��HeightFull
%                      if ~all(Item.isHeightFull(Item.Item_Strip(1,:) == uniStrip(i))) %�����STRIP��ӦITEM��isFull���Ǿ�Ϊ1,��STRIP��full
%                          Strip.isHeightFull(i) = 0;        
%                      end
%             else
%                      Strip.isHeightFull(i) = 0;       
%                      % 3 ���Ӽ�ʹmaxHeightDiff�ܴ�, ������߶ȸ�, Ҳ��ΪHeightFull
%                      if all(Item.isHeightFull(Item.Item_Strip(1,:) == uniStrip(i))) %�����STRIP��ӦITEM��isFull��Ϊ1,��STRIPҲΪfull
%                          Strip.isHeightFull(i) = 1;
%                      end
%             end
%         else
%             if all(Item.isHeightFull(Item.Item_Strip(1,:) == uniStrip(i))) %�����STRIP��ӦITEM��isFull��Ϊ1,��STRIPҲΪfull
%                 Strip.isHeightFull(i) = 1;
%             else
%                 Strip.isHeightFull(i) = 0;
%             end
%         end
%     end
%     
%     % �������
%     if any(Strip.isHeightFull==-1), error('����Strip.isHeightFullδ����!'); end
%      
% end

%% V2 ����isHeightBalance�߶Ⱦ���
function Strip = isFullStrip(Strip,Item)
    % 1 ѭ���ж�Strip�Ƿ����ItemΪfull��,�����,��StripΪfull
    uniStrip = unique(Item.Item_Strip(1,:));

    % 2 ѭ���ж�Strip�ڲ�Item֮�������ֵ, �Ƿ�<= ��С�ĶԽ��߻�һ������ֵ, ����,��ΪFull; 
    for i=1:length(Strip.isHeightFull)
        if Strip.isHeightBalance(i) == 0 %����߶Ȳ�����,һ���Ǹ߶Ȳ��� (�����Ǹ߶ȶ��ܸ�,������������)
            Strip.isHeightFull(i) = 0;
        else %����߶Ⱦ���,��Item�Ƿ�������,�粻��,Strip�߶�Ҳ����
            Strip.isHeightFull(i) = 1; %Strip��һItemҲ�Ǹ߶Ⱦ���, ������Item�ǲ���, ��StripҲ�ǲ���
            if ~all(Item.isHeightFull(Item.Item_Strip(1,:) == uniStrip(i))) %�����STRIP��ӦITEM��isFull���Ǿ�Ϊ1,��STRIP��full
                Strip.isHeightFull(i) = 0;
            end
        end
    end

    % �������
    if any(Strip.isHeightFull==-1), error('����Strip.isHeightFullδ����!'); end

end

%% ����2: �ж�STRIP�Ƿ������WidthFull��Item
% ****************** Strip���Ƿ����Width��Full��Item���� ************ ����
function Strip = isWidthFullStrip(Strip,Item) 
    % ѭ���ж�Strip�Ƿ�Widthfull, ��ȼ�϶ >= ��Strip������Item�Ŀ����Сֵ
    uniItem = unique(Item.Item_Strip(1,:));
    for i=1:length(Strip.isWidthFull)
        flagItem = Item.Item_Strip(1,:) == uniItem(i); 
                    max(Item.LWH(1, flagItem));
        if Strip.LW(1, i) >= min(Item.LWH(1, flagItem))
            Strip.isWidthFull(i) = 0;
        else
            Strip.isWidthFull(i) = 1;
        end
    end
end

%% ����3: �ж�STRIP�Ƿ�����
function Strip = isMixedStrip(Strip)
    % ѭ���ж�Strip�Ƿ�Ϊ��ͬLU.ID�Ļ����
    for i=1:length(Strip.isMixed)
        
         if numel(Strip.LID{i}) > 1
             Strip.isMixed(i) = 1;             
         else
             Strip.isMixed(i) = 0;
         end
         Strip.nbLID(i) = numel(Strip.LID{i});

         if numel(Strip.SID{i}) > 1
             Strip.isMixedSID(i) = 1;
         else
             Strip.isMixedSID(i) = 0;
         end
         Strip.nbSID(i) = numel(Strip.LID{i});
         
         if numel(Strip.EID{i}) > 1
             Strip.isMixedEID(i) = 1;
         else
             Strip.isMixedEID(i) = 0;
         end
         Strip.nbEID(i) = numel(Strip.LID{i});
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
