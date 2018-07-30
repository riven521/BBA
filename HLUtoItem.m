function [d] = HLUtoItem(d,ParaArray)
% ��Ҫ����:LU�Ѷ���γ�Item %  ����:�����(row);  ����:��������(coloum);
% Input ---  LUArray: ID LWH Weight ��LU: ����ԭ��˳��
% Output --- LUArray: order LUBeItemArray ��LU: ����ԭ��˳��(ORDER�ǽ���Item�㷨��LU˳��)
% Output --- ItemArray: ID LWH Weight��ITEM:û��˳���㷨������˳��

% d.LUArray.LUBeItemArray (2��*n��: ��1: LU�ڵڼ���item ��2:LU�����item��˳��(��-��)) 
% d.LUArray.order (1��*n��: LU����˳��)
% d.ItemArray.ID ��1��*n��: ITEM������ ��ͬ�ڲ�LU����)
% d.ItemArray.LWH (3��*n��: ITEM�ĳ����-������LU��ͬ,�߶�Ϊ�Ѷ��߶�)
% d.ItemArray.Weight (1��*n��: ITEM������)
% Ƕ�׺���(������ ��ת��)
% getLUorder
% getSortedLU
% getItem
% 

    %% ��ʼ��
% nDim LUά�� nLU LU���� nItem Item���� nLUid LU���� 
% heightBin Bin���߶�
nDim = size(d.LUArray.LWH,1);
nLU = size(d.LUArray.LWH,2);
nItem = nLU;
nLUid = size(unique(d.LUArray.ID),2);

tmpUniqueBin = unique(d.BinArray.LWH(1:nDim,:)','rows')';
heightBin = tmpUniqueBin(3);
clear tmpUniqueBin;

    %% LU����
    % getLUorder - ��ȡLU��˳��(�ص��Ǹ߶ȵݼ�����)
    % getSortedLU - ��ȡ��order������LU:sortedLUArray
    d.LUArray.order = getLUorder(); %��ȡ LU����(��ID����,��߶ȵݼ�)
    sortedLUArray = getSortedLU(d.LUArray.order);

    %% 55 LU->Itemת��
    ItemID = zeros(1,nItem);   %Item��ID����
    ItemisRota = ones(1,nItem)*2;   %Item�Ŀ���ת����(��ʼΪ2)
    
    ItemWeight = zeros(1,nItem);   %Item������
    ItemLWH = zeros(nDim,nItem);   %Item�Ŀ���
    LUBeItemSort = zeros(2,nLU);     %dim1:���ڵڼ���Item dim2:���ڸ�Item�ڼ����ŷ� 555
                % RotaFlagItem = ones(1,nItem)*5;   %Item����תflag��0��1,����Ϊδʹ�ã� TOBE DELE    
    
     getItem();     %555 LU��Item��ת��
    
    d.ItemArray.ID = ItemID(:,ItemID(1,:)>0);                                    % ȥ��δʹ�õ�  
    d.ItemArray.isRota = ItemisRota(:,ItemisRota(1,:)<=1);               % ȥ��δʹ�õ�     
    d.ItemArray.Weight = ItemWeight(:,ItemWeight(1,:)>0);            % ȥ��δʹ�õ� 
    d.ItemArray.LWH = ItemLWH(:,ItemLWH(1,:)>0);                        % ȥ��δʹ�õ� 
    d.LUArray.LUBeItemArray(:,d.LUArray.order) = LUBeItemSort; % d.LUArray.LUBeItemArray : ÿ�������LU���ĸ�Item��  �Լ�˳��
                %  d.ItemArray.itemRotaFlag = RotaFlagItem(:,RotaFlagItem(1,:)<=1); %TOBE DELE


    %% �������+����script
    % �������
    getITEMIDArray();
    % �����Ҫ���:���ÿ��item������ ԭʼ LU���
     printscript();
    
    %% Ƕ�׺���
    function order = getLUorder()
        tmpLUMatrix = [d.LUArray.ID; d.LUArray.LWH(1:nDim,:)];
        [~,tepLUorder] = sortrows(tmpLUMatrix',[1 4],{'ascend','descend'}); %1:ID 4:Hight
%         tepLUorder = 1:length(d.LUArray.ID)'; %ֱ�Ӹ�ֵ1:n
        % tepLUorder = [2 3 4 1 5]';
        order = tepLUorder';
    end

    function lu = getSortedLU(LUorder)
        lu = structfun(@(x) x(:,LUorder),d.LUArray,'UniformOutput',false);
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
        d.ItemIDArray.ID = unique(d.ItemArray.ID);
        nItemID = numel(d.ItemIDArray.ID);        

        for iID = 1:nItemID
        d.ItemIDArray.Weight(iID) = sum(d.ItemArray.Weight .* (d.ItemArray.ID == d.ItemIDArray.ID(iID)) );
        d.ItemIDArray.Volume(iID) = sum(prod(d.ItemArray.LWH) .* (d.ItemArray.ID == d.ItemIDArray.ID(iID)) );
        d.ItemIDArray.Area(iID) = sum(prod(d.ItemArray.LWH(1:2,:)) .* (d.ItemArray.ID == d.ItemIDArray.ID(iID)) );
        d.ItemIDArray.isRota(iID) =  unique(d.ItemArray.isRota(d.ItemArray.ID == d.ItemIDArray.ID(iID))); if ~isscalar(d.ItemIDArray.isRota(iID)), error('��������'); end
        end
        printstruct(d);
    end

    function printscript()
        for iItem = 1:max(d.LUArray.LUBeItemArray(1,:))
            [~,idx] = find(d.LUArray.LUBeItemArray(1,:)==iItem);
            fprintf('item %d �ĳ����Ϊ:  ',iItem);
            fprintf('( %d ) ',d.ItemArray.LWH(:,iItem));
            fprintf('\n');
            fprintf('item %d ���� original LU ������(�����)Ϊ  \n  ',iItem);
            fprintf('%d ',idx);
            fprintf('( %d ) ', d.LUArray.LWH(:,idx));
            fprintf('\n');
        end
    end
end




%% ���� ����ֵ��da
% ��ȡLUBeItemArray : ÿ�������LU���ĸ�Item��  �Լ�˳��
% ��ȡLWHItem:  �����ɵ�Item�ĳ����
% ��ȡLUorder: LU������
% d.ItemArray.LWH = LWHItem(:,LWHItem(1,:)>0); % ȥ��δʹ�õ�
% LUBeItemArray=LUBeItemArraySort;
% LUBeItemArray(:,LUorder) = LUBeItemArraySort;
% d.LUArray.LUorder = LUorder;
% d.LUArray.LUBeItemArray = LUBeItemArray;

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
