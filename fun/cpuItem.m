%% GET ITEM �������
% 1: Item.isHeightFull ( Item�����LU�߶� * 0.95 >= �����϶, ��Ϊ����)  computeisHeightFullItem
% 4: Item.isNonMixed  ����ÿ��Item�Ƿ�Ϊ����Ҫ��ƴ�Ŀ���
% 5: Item.nbItem ����ÿ��Item����ͬLU/ID������
%% ����
function   [Item] = cpuItem(Item,LU,Veh)

    %% ��ʼ��
    sz = size(Item.Weight);
    % hVeh  = Veh.LWH(3,1);  % ��ʱΪ��

    Item.isHeightFull = zeros(sz);    %Item���Ƿ�߶�����(��ʼΪ0)  %  Item.isWeightFine = ones(sz)*-1;    %Item���Ƿ���������(��ʼΪ-1) �޸�Ϊֱ���޸�ÿ��Item
    Item.MixOrder = ones(sz)*-1;    %1:ITEM���������ŵ�����strip,���������ITEM���; 0: �б�Ҫ���; ��ż������Item��ǰ����Strip����
        %Item.isMixedTile = zeros(sz);          % 0:��isNonMixed=1,ȫ��Ϊ0;��isNonMixed=0ʱ,�������±�����β��ITEM��ֵֵΪ1. 
    Item.nbItem = zeros(sz);  
    
    %% V2 SECTION 0 ����ITEM��PID,LID,SID,��TABLE����,����֪��ʲô��ʲô,����1,2,3�����滻
    % �ɻ�ϵ�LU.DOC����ITEM�ڰ�����PID,LID,SID������ 1808���� ����Item.PID,LID,SID��ʹ��
    T = struct2table(structfun(@(x) x',LU,'UniformOutput',false));
    
    nItem = size(Item.LWH,2);
    for iLuID=1:nItem
        f = T.LU_Item(:,1) == iLuID;                          %  f = logical(ones(height(t),1))
        Item.LID(:,iLuID) = {unique(T.ID(f))};             % NOTE: ITEM���LID��LU��ID
        %         Item.LID(:,iItem) = {unique(t.LID(f))};
        Item.SID(:,iLuID) = {unique(T.SID(f))};
        Item.EID(:,iLuID) = {unique(T.EID(f))};
        Item.PID(:,iLuID) = {unique(T.PID(f))};
    end    
    %  t2 = struct2table(structfun(@(x) x',Item,'UniformOutput',false));


    %% SECTION 1 ����ITEM��isHeightFull
    Item.isHeightFull = computeisHeightFullItem(Item,LU);
    
    %% 	SECTION 2 Item.isNonMixed: Item��������
    [Item.MixOrder] = computeMixandTile(Item,Veh);  
    
    %% SECTION 3 ����ITEM��nbItem
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


% v3��Item����ָ���Ϊ1��, nm����Ҫָ�ꡣ
function [nm] = computeMixandTile(Item,Veh) 

    LuIDs = cellfun(@(x) x(1), Item.LID);   % arrayAllLID: ����ITEM��Ӧ��LIDֵ ������ʽ
    uniLuLID = unique(LuIDs);
    
    n = length(LuIDs);
    [nm] = deal(zeros(1,n));
    
    % ѭ��: ���̵�ID/���͸���
    for iLuID=1:length(uniLuLID)
       
        % 1 ����nbmod���Ƿ���Ҫ��ϵ���Ҫָ�꣩
        % ����1 ���ڵ�LU flag���
        flagItem = LuIDs ==  uniLuLID(iLuID);  % flagItem = ItemLID(:) == uniItemLID(iItem);  idExchange(LU.ID)������,�����Ϊ���ߣ���Ӱ��޴� 
        
        ItemWidth = unique(Item.LWH(1,flagItem));     % Item���
        VehWidth = Veh.LWH(1,1);  % �������        
        nbmaxItem= floor(VehWidth/ItemWidth);       % Item�ɷ�����Ȳ���
        
        nbItem = sum(flagItem);                                     % Item ����
        
        % nbmod: ���������� �� ������� ȡ�� ��Ϊ0�������������������ͶѶ�Ŀ����ԣ�1 ���Ȼ ����ֻ��strip��
        nbmod = mod(nbItem,nbmaxItem);        if nbItem ==0 || nbmod>nbItem, error('cpuItem�ּ���isNonMixed����'); end
        
        % GapWidth = VehWidth - nbmaxItem*ItemWidth; % GapWidth ����������������ʣ��Item���Gap��� ����ʱ���ã�
        
        % 2 ������Ϊ0�����պ�����������
        if nbmod==0
            % ������Ϊ0����ֻ��һ��Strip�����������ȼ�����Ϊ�������Ȳ����㣬������������ϣ�Ҳ����ͨ��gap�����ı䣩
            if nbItem==nbmaxItem                                        % && GapWidth > ItemWidth*0.5
                nm(flagItem)  = 2;
            else                                                                         % ����ָ��Ϊ�����������Item���
                nm(flagItem)  = 1;
            end            
        elseif nbmod~=0            
            nm(flagItem)  = 1;  %��Ĭ��Ϊ1���̶��Ե����ı����ϵı��Ϊ3����������û��2
            
            tmpSort=[Item.isHeightFull;Item.HLayer;Item.LWH(3,:)];
            [~, order]=sortrows(tmpSort(:,flagItem)', [1,2,3],{'ascend','ascend','ascend'});
            flagItemIdx = find(flagItem);
            flagmodIdx = flagItemIdx(order(1:nbmod));
            nm(flagmodIdx)=3;  %NOTE: flagmodIdx ����Item���
        end
        
    end
end

%% v2�� v3����2������Ϊ1��
% % function [nm, mt] = computeMixandTile(Item,Veh) 
% % 
% %     
% %     LuIDs = cellfun(@(x) x(1), Item.LID);   % arrayAllLID: ����ITEM��Ӧ��LIDֵ ������ʽ
% %     uniLuLID = unique(LuIDs);
% %     
% %     n = length(LuIDs);
% %     [nm, mt] = deal(zeros(1,n));
% %     
% %     % ѭ��: ���̵�ID/���͸���
% %     for iLuID=1:length(uniLuLID)
% %        
% %         % ����1 ���ڵ�LU flag���
% %         flagItem = LuIDs ==  uniLuLID(iLuID);  % flagItem = ItemLID(:) == uniItemLID(iItem);  idExchange(LU.ID)������,�����Ϊ���ߣ���Ӱ��޴� 
% %         
% %         ItemWidth = unique(Item.LWH(1,flagItem));     % Item���
% %         VehWidth = Veh.LWH(1,1);  % �������        
% %         nbmaxItem= floor(VehWidth/ItemWidth);       % Item�ɷ�����Ȳ���
% %         
% %         nbItem = sum(flagItem);                                     % Item ����
% %         
% %         % nbmod: ���������� �� ������� ȡ�� ��Ϊ0�������������������ͶѶ�Ŀ����ԣ�1 ���Ȼ ����ֻ��strip��
% %         nbmod = mod(nbItem,nbmaxItem);        if nbItem ==0 || nbmod>nbItem, error('cpuItem�ּ���isNonMixed����'); end
% %         
% %         % GapWidth ����������������ʣ��Item���Gap��� ����ʱ���ã�
% %         GapWidth = VehWidth - nbmaxItem*ItemWidth;
% %         
% %         % ���㣺flagNonMixedItem��1 ��������Item�������ϵģ���������0��������������������϶Ѷ�����Դ�
% %         % ������Ϊ0����һ�����������пգ����������Item��ϣ�
% %         if nbmod~=0 
% %             flagNonMixedItem = 0;
% %         else
% %             % ������Ϊ0����ֻ��һ��Strip��ʣ��Gap����ֽϴ���ָ��Ϊ���������Item��ϣ�
% %             if nbItem==nbmaxItem && GapWidth > ItemWidth*0.5
% %                 flagNonMixedItem = 0; % flagNonMixedItem = 0;
% %             else % ����ָ��Ϊ�����������Item���
% %                 flagNonMixedItem = 1;
% %             end
% %         end
% %         
% %         % ��ֵ������flagNonMixedItem����ֵItem��isNonMixed/isMixedTile ��Item�������Ҫ���ݣ�
% %         if   flagNonMixedItem 
% %             nm(flagItem) = 1;  % ���ItemΪ�����ϵ�, ����˦β�ı�Ǿ�Ϊ0
% %             mt(flagItem) = 0;
% %         else  
% %             nm(flagItem)= 0;
% %             
% %             % mt  ����Item��isMixedTile �����ҳ��������ϵ�Items�������Ա��
% %             if nbmod == 0 % �� nbItem==nbmaxItem && GapWidth > ItemWidth*0.5 ��ʵ������mix����Gap̫����ֻ��һ��
% %                 mt(flagItem)=1;  %NOTE: flagItem���е�Item�����
% %             else            
% %                 tmpSort=[Item.isHeightFull;Item.HLayer;Item.LWH(3,:)];
% %                 [~, order]=sortrows(tmpSort(:,flagItem)', [1,2,3],{'ascend','ascend','ascend'});
% %                 flagItemIdx = find(flagItem);
% %                 flagmodIdx = flagItemIdx(order(1:nbmod));
% %                 mt(flagmodIdx)=2;  %NOTE: flagmodIdx ����Item���
% %             end
% %         end
% %     end
% % end


    %% v1: computeMixandTile
% % % %     %% 	SECTION 2 Item.isNonMixed: and  Item.isMixedTile
% % % %     %     [Item.isNonMixed, Item.isMixedTile] = computeMixandTile(Item) todo �ĵ�һ��������
% % % %     % ****************** Iten���Ƿ�Ϊ����Ҫ��ƴ���� ************ ����
% % % %     % GET Item.isNonMixed: ����ÿ��Item�Ƿ�Ϊ����Ҫ��ƴ�Ŀ���
% % % %     ItemLID = cellfun(@(x) x(1), Item.LID); % arrayAllLID: ����ITEM��Ӧ��LIDֵ ������ʽ
% % % %     uniItemLID = unique(ItemLID);
% % % %     % ѭ��: LID����
% % % %     for iItem=1:length(uniItemLID)
% % % %         % Item i ���ڵ�LU flag���
% % % %         flagItem = ItemLID(:) == iItem;  % flagItem = ItemLID(:) == uniItemLID(iItem);  idExchange(LU.ID)������,�����Ϊ���ߣ���Ӱ��޴� 
% % % %         
% % % %         ItemWidth = unique(Item.LWH(1,flagItem));     % Item���
% % % %         VehWidth = Veh.LWH(1,1);  % �������
% % % %         
% % % %         nbmaxItem= floor(VehWidth/ItemWidth);     % Item�ɷ�����Ȳ���
% % % %         nbItem = sum(flagItem);
% % % %         
% % % %         % GapWidth ����������������ʣ��Item���Gap���
% % % %         GapWidth = VehWidth - nbmaxItem*ItemWidth;
% % % %         
% % % %         % nbmod: ���������� �� ������� ȡ�� 
% % % %         nbmod = mod(nbItem,nbmaxItem);        if nbItem ==0 || nbmod>nbItem, error('cpuItem�ּ���isNonMixed����'); end
% % % %         
% % % %         % ���㣺flagNonMixedItem��1 ��������Item�������ϵģ���������0��������������������϶Ѷ�����Դ�
% % % %         % ������Ϊ0����һ�����������пգ����������Item��ϣ�
% % % %         if nbmod~=0 
% % % %             flagNonMixedItem = 0;
% % % %         else
% % % %             % ������Ϊ0����ֻ��һ��Strip��ʣ��Gap����ֽϴ���ָ��Ϊ���������Item��ϣ�
% % % %             if nbItem==nbmaxItem && GapWidth > ItemWidth*0.5
% % % %                 flagNonMixedItem = 0;
% % % %             else % ����ָ��Ϊ�����������Item���
% % % %                 flagNonMixedItem = 1;
% % % %             end
% % % %         end
% % % %         
% % % %         % ��ֵ������flagNonMixedItem����ֵItem��isNonMixed/isMixedTile ��Item�������Ҫ���ݣ�
% % % %         if   flagNonMixedItem 
% % % %             Item.isNonMixed(flagItem) = 1;  % ���ItemΪ�����ϵ�, ����˦β�ı�Ǿ�Ϊ0
% % % %             Item.isMixedTile(flagItem) = 0;
% % % %         else  
% % % %             Item.isNonMixed(flagItem)= 0;
% % % %             
% % % %             % Item.isMixedTile  ����Item��isMixedTile �����ҳ��������ϵ�Items�������Ա��
% % % %             if nbmod == 0 % �� nbItem==nbmaxItem && GapWidth > ItemWidth*0.5 ��ʵ������mix����Gap̫����ֻ��һ��
% % % %                 Item.isMixedTile(flagItem)=1;  %NOTE: flagItem���е�Item�����
% % % %             else            
% % % %                 tmpSort=[Item.isHeightFull;Item.HLayer;Item.LWH(3,:)];
% % % %                 [~, order]=sortrows(tmpSort(:,flagItem)', [1,2,3],{'ascend','ascend','ascend'});
% % % %                 flagItemIdx = find(flagItem);
% % % %                 flagmodIdx = flagItemIdx(order(1:nbmod));
% % % %                 Item.isMixedTile(flagmodIdx)=1;  %NOTE: flagmodIdx ����Item���
% % % %             end
% % % %         end
% % % %     end


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
    