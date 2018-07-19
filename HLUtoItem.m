function [da] = HLUtoItem(da,ParaArray)
% ����: LU: ID LWH
% ���: 
% da.LUArray (1 LWH (��֪)
% da.LUArray (2 LUBeItemArray (dim1:���LU��ĳ��item dim2:lu����˳��(�ײ�-�߲�)) 
% da.LUArray (3 LUorder ������LU����˳��)
% da.ItemArray (1 LWH ) 

%% ��ʼ��
% nDim LUά�� nLU LU���� nItem Item���� nLUid LU���� tmpheightBin Bin���߶�
nDim = size(da.LUArray.LWH,1);
nLU = size(da.LUArray.LWH,2);
nItem = nLU;
nLUid = size(unique(da.LUArray.ID),2);

tmpUniqueBin = unique(da.BinArray.LWH(1:nDim,:)','rows')';
heightBin = tmpUniqueBin(3);
clear tmpUniqueBin;

%% LU���򲢳�ʼ��
% ��ȡLUorder
% ��ȡLWHLuSort
% ��ȡIDLuSort
% sort ��ȡ��id��߶ȵ����� %

tmpLUMatrix = [da.LUArray.ID; da.LUArray.LWH(1:nDim,:)]; %1:ID 4:Hight
[~,tepLUorder] = sortrows(tmpLUMatrix',[1 4],{'ascend','descend'});
% tepLUorder = 1:length(da.LUArray.ID)';
% tepLUorder = [2 3 4 1 5]';
LUorder = tepLUorder'; %��ȡ LU����(��ID��߶�)���˳��
LWHLuSort = da.LUArray.LWH(:,LUorder);
IDLuSort = da.LUArray.ID(:,LUorder);
clear tmpLUMatrix tepLUorder;

%% 55 LU->Itemת�� 
% ��ȡLUBeItemArraySort : ÿ�������LU���ĸ�Item��  �Լ�˳��
% ��ȡLWHItem:  �����ɵ�Item�ĳ����
LWHItem = zeros(nDim,nItem);   %Item�Ŀ���
LUBeItemArraySort = zeros(2,nLU); %dim1:���ڵڼ���Item dim2:���ڸ�Item�ڼ����ŷ� 555
itemBeLUArray = zeros(1,nItem);  % ÿ��Item�ڶѶ��LU���� ���ڲ���
iItem = 1;

for iLUid=1:nLUid
    heightLeft = heightBin;
    for iLU=1:nLU
        if IDLuSort(iLU) == iLUid %���Ե�ǰLU��ӦLuid�ڸ�iLUid�ڵĽ��в���
            if heightLeft < LWHLuSort(nDim,iLU) %�統ǰLU�߶�����(�߶��ڵ�nDim��)
                iItem =  iItem + 1;
                heightLeft = heightBin;%  
            end              
            heightLeft = heightLeft - LWHLuSort(nDim,iLU); %����ʣ��߶�
            LWHItem(1:2,iItem) = LWHLuSort(1:2,iLU); %����item����
            LWHItem(3,iItem) = LWHItem(3,iItem) + LWHLuSort(nDim,iLU); %����item�߶�
            itemBeLUArray(iItem) = itemBeLUArray(iItem) + 1;
            LUBeItemArraySort(1,iLU) = iItem;
            LUBeItemArraySort(2,iLU) = itemBeLUArray(iItem);
        end
    end
    iItem =  iItem + 1;
end

%% ���� ����ֵ��da
% ��ȡLUBeItemArray : ÿ�������LU���ĸ�Item��  �Լ�˳��
% ��ȡLWHItem:  �����ɵ�Item�ĳ����
% ��ȡLUorder: LU������
da.ItemArray.LWH = LWHItem(:,LWHItem(1,:)>0); % ȥ��δʹ�õ�
LUBeItemArray=LUBeItemArraySort;
LUBeItemArray(:,LUorder) = LUBeItemArraySort;
da.LUArray.LUorder = LUorder;
da.LUArray.LUBeItemArray = LUBeItemArray;

% ���Դ���
% da.ItemArray.LWH 
% LUBeItemArray
% da.LUArray.LWH

% �����Ҫ���:���ÿ��item������ ԭʼ LU���
% % for iItem = 1:max(da.LUArray.LUBeItemArray(1,:))
% %     [~,idx] = find(da.LUArray.LUBeItemArray(1,:)==iItem);
% %     fprintf('item %d �ĳ����Ϊ:  ',iItem);
% %     fprintf('( %d ) ',da.ItemArray.LWH(:,iItem));
% %     fprintf('\n');
% %     fprintf('item %d ���� original LU ������(�����)Ϊ  \n  ',iItem);
% %     fprintf('%d ',idx);
% %     fprintf('( %d ) ', da.LUArray.LWH(:,idx));
% %     fprintf('\n');
% % end

end