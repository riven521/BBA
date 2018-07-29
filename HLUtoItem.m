function [da] = HLUtoItem(da,ParaArray)
% ��Ҫ����:LU�Ѷ���γ�Item %  ����:�����(row);  ����:��������(coloum);
% Input ---  LUArray: ID LWH Weight
% Output --- LUArray: order LUBeItemArray
% Output --- ItemArray: ID LWH Weight
% da.LUArray.LUBeItemArray (2��*n��: ��1: LU�ڵڼ���item ��2:LU�����item��˳��(��-��)) 
% da.LUArray.order (1��*n��: LU����˳��)
% da.ItemArray.ID ��1��*n��: ITEM������ ��ͬ�ڲ�LU����)
% da.ItemArray.LWH (3��*n��: ITEM�ĳ����-������LU��ͬ,�߶�Ϊ�Ѷ��߶�)
% da.ItemArray.Weight (1��*n��: ITEM������)
% Ƕ�׺���(������ ��ת��)
% getLUorder
% getSortedLU
% getItem

    %% ��ʼ��
% nDim LUά�� nLU LU���� nItem Item���� nLUid LU���� 
% heightBin Bin���߶�
nDim = size(da.LUArray.LWH,1);
nLU = size(da.LUArray.LWH,2);
nItem = nLU;
nLUid = size(unique(da.LUArray.ID),2);

tmpUniqueBin = unique(da.BinArray.LWH(1:nDim,:)','rows')';
heightBin = tmpUniqueBin(3);
clear tmpUniqueBin;

    %% LU����
    % getLUorder - ��ȡLU��˳��(�ص��Ǹ߶ȵݼ�����)
    % getSortedLU - ��ȡ��order������LU:sortedLUArray
    da.LUArray.order = getLUorder(); %��ȡ LU����(��ID����,��߶ȵݼ�)
    sortedLUArray = getSortedLU(da.LUArray.order);

    %% 55 LU->Itemת��
    IDItem = zeros(1,nItem);   %Item��ID����
    WeightItem = zeros(1,nItem);   %Item������
    LWHItem = zeros(nDim,nItem);   %Item�Ŀ���
% %     RotaFlagItem = ones(1,nItem)*5;   %Item����תflag��0��1,����Ϊδʹ�ã� TOBE
% DELE
    LUBeItemArraySort = zeros(2,nLU); %dim1:���ڵڼ���Item dim2:���ڸ�Item�ڼ����ŷ� 555
    
    getItem(); %555 ת��
    
    da.ItemArray.ID = IDItem(:,IDItem(1,:)>0); % ȥ��δʹ�õ�     
    da.ItemArray.Weight = WeightItem(:,WeightItem(1,:)>0); % ȥ��δʹ�õ� 
% %     da.ItemArray.itemRotaFlag = RotaFlagItem(:,RotaFlagItem(1,:)<=1); %
% TO BE DELE

    da.ItemArray.LWH = LWHItem(:,LWHItem(1,:)>0); % ȥ��δʹ�õ� 
    da.LUArray.LUBeItemArray(:,da.LUArray.order) = LUBeItemArraySort; % da.LUArray.LUBeItemArray : ÿ�������LU���ĸ�Item��  �Լ�˳��

    %% �������+����script
    % �������
    getITEMIDArray();
    % �����Ҫ���:���ÿ��item������ ԭʼ LU���
    printscript();
    
    %% Ƕ�׺���
    function order = getLUorder()
        tmpLUMatrix = [da.LUArray.ID; da.LUArray.LWH(1:nDim,:)];
        [~,tepLUorder] = sortrows(tmpLUMatrix',[1 4],{'ascend','descend'}); %1:ID 4:Hight
        % tepLUorder = 1:length(da.LUArray.ID)'; %ֱ�Ӹ�ֵ1:n
        % tepLUorder = [2 3 4 1 5]';
        order = tepLUorder';
    end

    function lu = getSortedLU(LUorder)
        lu = structfun(@(x) x(:,LUorder),da.LUArray,'UniformOutput',false);
    end

    function getItem()
        iItem = 1;
        itemBeLUArray = zeros(1,nItem);  % ÿ��Item�ڶѶ��LU���� ���ڲ���
        for iLUid=1:nLUid
            heightLeft = heightBin;
            for iLU=1:nLU
                if sortedLUArray.ID(iLU) == iLUid %���Ե�ǰLU��ӦLuid�ڸ�iLUid�ڵĽ��в���
                    if heightLeft < sortedLUArray.LWH(nDim,iLU) %�統ǰLU�߶�����(�߶��ڵ�nDim��)
                        iItem =  iItem + 1;
                        heightLeft = heightBin;%
                    end
                    heightLeft = heightLeft - sortedLUArray.LWH(nDim,iLU); %����ʣ��߶�
                    LWHItem(1:2,iItem) = sortedLUArray.LWH(1:2,iLU); %����item����
                    LWHItem(3,iItem) = LWHItem(3,iItem) + sortedLUArray.LWH(nDim,iLU); %����item�߶�
                    WeightItem(1,iItem) = WeightItem(1,iItem) + sortedLUArray.Weight(1,iLU); %����item����
                    IDItem(1,iItem) = sortedLUArray.ID(1,iLU); %����ID����
%                     RotaFlagItem(1,iItem) = sortedLUArray.RotaFlag(1,iLU); %����RotaFlag����                    
                    itemBeLUArray(iItem) = itemBeLUArray(iItem) + 1;
                    LUBeItemArraySort(1,iLU) = iItem;
                    LUBeItemArraySort(2,iLU) = itemBeLUArray(iItem);
                end
            end
            iItem =  iItem + 1;
        end
    end

    %%  ��ȡITEMID�����������(ͬ����ID����������������)
    function getITEMIDArray()
%          printstruct(da);
        da.ItemIDArray.ID = unique(da.ItemArray.ID);
        nItemID = numel(da.ItemIDArray.ID);        

        for iID = 1:nItemID
        da.ItemIDArray.Weight(iID) = sum(da.ItemArray.Weight .* (da.ItemArray.ID == da.ItemIDArray.ID(iID)) );
        da.ItemIDArray.Volume(iID) = sum(prod(da.ItemArray.LWH) .* (da.ItemArray.ID == da.ItemIDArray.ID(iID)) );
        da.ItemIDArray.Area(iID) = sum(prod(da.ItemArray.LWH(1:2,:)) .* (da.ItemArray.ID == da.ItemIDArray.ID(iID)) );
        end
        
    end

    function printscript()
        for iItem = 1:max(da.LUArray.LUBeItemArray(1,:))
            [~,idx] = find(da.LUArray.LUBeItemArray(1,:)==iItem);
            fprintf('item %d �ĳ����Ϊ:  ',iItem);
            fprintf('( %d ) ',da.ItemArray.LWH(:,iItem));
            fprintf('\n');
            fprintf('item %d ���� original LU ������(�����)Ϊ  \n  ',iItem);
            fprintf('%d ',idx);
            fprintf('( %d ) ', da.LUArray.LWH(:,idx));
            fprintf('\n');
        end
    end
end




%% ���� ����ֵ��da
% ��ȡLUBeItemArray : ÿ�������LU���ĸ�Item��  �Լ�˳��
% ��ȡLWHItem:  �����ɵ�Item�ĳ����
% ��ȡLUorder: LU������
% da.ItemArray.LWH = LWHItem(:,LWHItem(1,:)>0); % ȥ��δʹ�õ�
% LUBeItemArray=LUBeItemArraySort;
% LUBeItemArray(:,LUorder) = LUBeItemArraySort;
% da.LUArray.LUorder = LUorder;
% da.LUArray.LUBeItemArray = LUBeItemArray;

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
