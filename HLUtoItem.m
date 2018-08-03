function [d] = HLUtoItem(d,ParaArray)
% ��Ҫ����:LU�Ѷ���γ�Item %  ����:�����(row);  ����:��������(coloum);
% Input ---  LU: ID LWH Weight ��LU: ����ԭ��˳��
% Output --- LU: order LU_Item ��LU: ����ԭ��˳��(ORDER�ǽ���Item�㷨��LU˳��)
% Output --- Item: ID LWH Weight��ITEM:û��˳���㷨������˳��

% d.LU.LU_Item (2��*n��: ��1: LU�ڵڼ���item ��2:LU�����item��˳��(��-��)) 
% d.LU.order (1��*n��: LU����˳��)
% d.Item.ID ��1��*n��: ITEM������ ��ͬ�ڲ�LU����)
% d.Item.LWH (3��*n��: ITEM�ĳ����-������LU��ͬ,�߶�Ϊ�Ѷ��߶�)
% d.Item.Weight (1��*n��: ITEM������)
% Ƕ�׺���(������ ��ת��)
% getLUorder
% getSortedLU
% getItem
% 

    %% ��ʼ��
% nDim LUά�� nLU LU���� nItem Item���� nLUid LU���� 
% heightBin Bin���߶�
nDim = size(d.LU.LWH,1);
nLU = size(d.LU.LWH,2);
nItem = nLU;
nLUid = size(unique(d.LU.ID),2);

tmpUniqueBin = unique(d.Veh.LWH(1:nDim,:)','rows')';
heightBin = tmpUniqueBin(3);
clear tmpUniqueBin;

    %% LU����
    % getLUorder - ��ȡLU��˳��(�ص��Ǹ߶ȵݼ�����)
    % getSortedLU - ��ȡ��order������LU:sortedLUArray
    d.LU.order = getLUorder(); %��ȡ LU����(��ID����,��߶ȵݼ�)
    sortedLUArray = getSortedLU(d.LU.order);

    %% 55 LU->Itemת��
    ItemID = zeros(1,nItem);   %Item��ID����
    ItemisRota = ones(1,nItem)*2;   %Item�Ŀ���ת����(��ʼΪ2)
    
    ItemWeight = zeros(1,nItem);   %Item������
    ItemLWH = zeros(nDim,nItem);   %Item�Ŀ���
    LUBeItemSort = zeros(2,nLU);     %dim1:���ڵڼ���Item dim2:���ڸ�Item�ڼ����ŷ� 555
                % RotaFlagItem = ones(1,nItem)*5;   %Item����תflag��0��1,����Ϊδʹ�ã� TOBE DELE    
    
     getItem();     %555 LU��Item��ת��
    
    d.Item.ID = ItemID(:,ItemID(1,:)>0);                                    % ȥ��δʹ�õ�  
    d.Item.isRota = ItemisRota(:,ItemisRota(1,:)<=1);               % ȥ��δʹ�õ�     
    d.Item.Weight = ItemWeight(:,ItemWeight(1,:)>0);            % ȥ��δʹ�õ� 
    d.Item.LWH = ItemLWH(:,ItemLWH(1,:)>0);                        % ȥ��δʹ�õ� 
    d.LU.LU_Item(:,d.LU.order) = LUBeItemSort; % d.LU.LU_Item : ÿ�������LU���ĸ�Item��  �Լ�˳��
                %  d.Item.itemRotaFlag = RotaFlagItem(:,RotaFlagItem(1,:)<=1); %TOBE DELE


    %% �������+����script
    % �������
    getITEMIDArray();
    % �����Ҫ���:���ÿ��item������ ԭʼ LU���
     printscript();
    
    %% Ƕ�׺���
    function order = getLUorder()
        tmpLUMatrix = [d.LU.ID; d.LU.LWH(1:nDim,:)];
        [~,tepLUorder] = sortrows(tmpLUMatrix',[1 4],{'ascend','descend'}); %1:ID 4:Hight
%         tepLUorder = 1:length(d.LU.ID)'; %ֱ�Ӹ�ֵ1:n
        % tepLUorder = [2 3 4 1 5]';
        order = tepLUorder';
    end

    function lu = getSortedLU(LUorder)
        lu = structfun(@(x) x(:,LUorder),d.LU,'UniformOutput',false);
    end

    % ��LUת��ΪItem����Ҫ����
    function getItem()
        iItem = 1;
        tmpitemBeLU = zeros(1,nItem);  % ÿ��Item�ڶѶ��LU���� ���ڲ���
        for iLUid=1:nLUid
            heightLeft = heightBin;
            for iLU=1:nLU
                if sortedLUArray.ID(iLU) == iLUid %���Ե�ǰLU��ӦLuid�ڸ�iLUid�ڵĽ��в���
                    if heightLeft < sortedLUArray.LWH(nDim,iLU) %�統ǰLU�߶�����(�߶��ڵ�nDim��)
                        iItem =  iItem + 1;
                        heightLeft = heightBin; %
                    end
                    
                    heightLeft = heightLeft - sortedLUArray.LWH(nDim,iLU); %����ʣ��߶�
                    ItemLWH(1:2,iItem) = sortedLUArray.LWH(1:2,iLU); %����item����
                    ItemLWH(3,iItem) = ItemLWH(3,iItem) + sortedLUArray.LWH(nDim,iLU); %����item�߶�
                    ItemWeight(1,iItem) = ItemWeight(1,iItem) + sortedLUArray.Weight(1,iLU); %����item����
                    
                    tmpitemBeLU(iItem) = tmpitemBeLU(iItem) + 1;
                    LUBeItemSort(1,iLU) = iItem;
                    LUBeItemSort(2,iLU) = tmpitemBeLU(iItem);
                    
                    ItemID(1,iItem) = sortedLUArray.ID(1,iLU); %����ID����
                    ItemisRota(1,iItem) = sortedLUArray.isRota(1,iLU);  %����ID����ת����
                                        %  RotaFlagItem(1,iItem) = sortedLUArray.RotaFlag(1,iLU); %����RotaFlag����                    
                end
            end
            iItem =  iItem + 1;
        end
    end

    %%  ��ȡITEMID�����������(ͬ����ID������������������item�Ƿ����ת)
    function getITEMIDArray()
                %  printstruct(d);  
        d.ItemID.ID = unique(d.Item.ID);
        nItemID = numel(d.ItemID.ID);        

        for iID = 1:nItemID
        d.ItemID.Weight(iID) = sum(d.Item.Weight .* (d.Item.ID == d.ItemID.ID(iID)) );
        d.ItemID.Volume(iID) = sum(prod(d.Item.LWH) .* (d.Item.ID == d.ItemID.ID(iID)) );
        d.ItemID.Area(iID) = sum(prod(d.Item.LWH(1:2,:)) .* (d.Item.ID == d.ItemID.ID(iID)) );
        d.ItemID.isRota(iID) =  unique(d.Item.isRota(d.Item.ID == d.ItemID.ID(iID))); if ~isscalar(d.ItemID.isRota(iID)), error('��������'); end
        end
        printstruct(d);
    end

    function printscript()
        for iItem = 1:max(d.LU.LU_Item(1,:))
            [~,idx] = find(d.LU.LU_Item(1,:)==iItem);
            fprintf('item %d �ĳ����Ϊ:  ',iItem);
            fprintf('( %d ) ',d.Item.LWH(:,iItem));
            fprintf('\n');
            fprintf('item %d ���� original LU ������(�����)Ϊ  \n  ',iItem);
            fprintf('%d ',idx);
            fprintf('( %d ) ', d.LU.LWH(:,idx));
            fprintf('\n');
        end
    end
end




%% ���� ����ֵ��da
% ��ȡLUBeItemArray : ÿ�������LU���ĸ�Item��  �Լ�˳��
% ��ȡLWHItem:  �����ɵ�Item�ĳ����
% ��ȡLUorder: LU������
% d.Item.LWH = LWHItem(:,LWHItem(1,:)>0); % ȥ��δʹ�õ�
% LU_Item=LUBeItemArraySort;
% LU_Item(:,LUorder) = LUBeItemArraySort;
% d.LU.LUorder = LUorder;
% d.LU.LU_Item = LU_Item;

% % function outStructArray = SortArrayofStruct( structArray, fieldName )
% %     %UNTITLED2 Summary of this function goes here
% %     %   Detailed explanation goes here
% %     if ( ~isempty(structArray) &&  length (structArray) > 0)
% %       [~,I] = sort(arrayfun (@(x) x.(fieldName), structArray)) ;
% %       outStructArray = structArray(I) ;        
% %     else 
% %         disp ('Array of struct is empty');
% %     end      
% % end
