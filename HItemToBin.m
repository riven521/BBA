function [da] = HItemToBin(da)
% ��Ҫ����:Item����Bin��
%% ��ȡitem��bin�Ĺ�ϵ itemBeBinMatrix  + ���� CoordItemBin
%% ��ȡLU��bin�Ĺ�ϵ LUBeBinMatrix  + ���� CoordLUBin
%% ��ȡLU��item�Ĺ�ϵ LURotaFlag

%% ��ʼ��
nDim = size(da.ItemArray.LWH,1);  
nItem = size(da.ItemArray.LWH,2); %����ʹ�õ�Item������
nLU= size(da.LUArray.LWH,2); %����ʹ�õ�LU������
uniBinDataMatrix = unique((da.BinArray.LWH(1:nDim,:))','rows')';
%% �ṹ����ȡ
LUBeItemArray = da.LUArray.LUBeItemArray;
stripBeBinMatrix = da.StripArray.stripBeBinMatrix;
itemBeStripMatrix = da.ItemArray.itemBeStripMatrix;

itemRotaed = da.ItemArray.Rotaed;
% % itemRotaFlag = da.ItemArray.itemRotaFlag;
CoordItemStrip = da.ItemArray.CoordItemStrip;

LWStrip = da.StripArray.LW; % LWStripSort = da.StripArray.LW(:,da.StripArray.striporder);
LWHLU = da.LUArray.LWH;
    
itemBeBinMatrix=zeros(2,nItem);
CoordItemBin=zeros(2,nItem);

LUBeBinMatrix=zeros(2,nLU);
CoordLUBin=zeros(2,nLU);
LURotaed=zeros(1,nLU);
    
% ѭ��ÿ��bin
for iBin=1:max(stripBeBinMatrix(1,:))
    tmpItemSeq=1;  %ÿ��bin��˳���1��ʼ
    tmpLUSeq=1;
    tmpStrip = []; %����CoordItemBin��(��)��ʹ�� %       tmpLWStrip = LWStrip(:,stripBeBinMatrix(1,:)==iBin); %��bin��ʣ���Ⱥͳ�(��)��
    nbStrip=numel(find(stripBeBinMatrix(1,:)==iBin));    %nbStrip����iBin��strip����    
    % ѭ��bin��ÿ��strip
    for iStrip=1:nbStrip %�Ӱ��ŵĵ�һ��Strip��ʼ�𲽷���
        [~,thisStrip] = find(stripBeBinMatrix(1,:)==iBin & stripBeBinMatrix(2,:)==iStrip ); %��bin�ڵ�iStrip��strip      
        if ~isscalar(thisStrip), error('�������');  end
        nbItem=numel(find(itemBeStripMatrix(1,:)==thisStrip));
        % ѭ��stripÿ��item
        for iItem=1:nbItem %ͬһstrip�� ѭ��ȡֵiItem
            [~,thisItem] = find(itemBeStripMatrix(1,:)==thisStrip & itemBeStripMatrix(2,:)==iItem );%��Strip�ڵ�iItem��item
            if ~isscalar(thisItem), error('�������');  end
            
            % ����itemBeBinMatrix 555
            itemBeBinMatrix(1,thisItem)=iBin;  % ��thisItem����iBin�����
            itemBeBinMatrix(2,thisItem)=tmpItemSeq;  % ��thisItem�ڱ�iBin�ڵ�˳�� ��1��ʼ
            tmpItemSeq=tmpItemSeq+1;     
            % ����CoordItemBin(��strip���¶���Item) 555 HAVE checked?
            CoordItemBin(1,thisItem) = CoordItemStrip(1,thisItem); % 55 ��thisItem��x���겻��
                % iStrip = 1�� �׸�strip�߶�һ��Ϊ0��iStrip=2���߶�һ����Ϊ0
            CoordItemBin(2,thisItem) = sum(LWStrip(2,tmpStrip)); % 555 ����ԭ��2��thisStrip˳��ȷ�� % 555 LWStripSort����ԭ��:iStripÿ��bin����ȡ��;LWStripSort�����е�Strip %             CoordItemBin(2,thisItem) = sum(tmpLWStrip(2,1:iStrip-1)); 
            CoordItemBin(3,thisItem) = 0; %Item����0��ʼ
            
            % ���Ӷ�LU�ĸ���
            tmpLU = []; %����CoordLUBin�߶�ʹ�� % tmpLWLU = LWHLU(:,LUBeItemArray(1,:)==thisItem); %����
            nbLU=numel(find(LUBeItemArray(1,:)==thisItem));            
            % ѭ��itemÿ��LU
            for iLU=1:nbLU
                [~,thisLU] = find(LUBeItemArray(1,:)==thisItem & LUBeItemArray(2,:)==iLU);
                if ~isscalar(thisLU), error('�������');  end
                
                % ����LURotaed 555
                LURotaed(1,thisLU)=itemRotaed(1,thisItem);
                % ����LUBeBinMatrix 555
                LUBeBinMatrix(1,thisLU)=iBin;
                LUBeBinMatrix(2,thisLU)=tmpLUSeq;
                tmpLUSeq=tmpLUSeq+1;               
                % ����CoordLUBin 555 HAVE checked
                CoordLUBin(1,thisLU) = CoordItemBin(1,thisItem);
                CoordLUBin(2,thisLU) = CoordItemBin(2,thisItem);                
                CoordLUBin(3,thisLU) = sum(LWHLU(3,tmpLU));  %555 ����ԭ��ͬ��% CoordLUBin(3,thisLU) = sum(tmpLWLU(3,1:iLU-1));
                tmpLU = [tmpLU thisLU];
            end
        end
        tmpStrip = [tmpStrip thisStrip];
    end
end

%% �ṹ�帳ֵ
da.ItemArray.itemBeBinMatrix=itemBeBinMatrix;
da.ItemArray.CoordItemBin=CoordItemBin;

da.LUArray.LUBeBinMatrix= LUBeBinMatrix;
da.LUArray.CoordLUBin = CoordLUBin;

%NOTE: �˴���LU��ITEM�׶ε���ת��ֵ��LU��   
%GcheekInput����ĶԱȵ�����ԭʼֵ���˴�Ҫ�滻Ϊ��һ�����ԭʼֵ.
%% % ���治��Ҫ��, ��ΪLUArray.Rotaed�Ѿ���HItemToStrip���¹���
% da.LUArray.Rotaed = LURotaed;


% printstruct(da);
% �����Ҫ���:��ô�1��ʼÿ��bin������item����
% CoordItemBin itemBeBinMatrix
for iBin = 1:max(da.ItemArray.itemBeBinMatrix(1,:))
    [~,idx] = find(da.ItemArray.itemBeBinMatrix(1,:)==iBin); %��iBin�µ�item������
    idxSeq = da.ItemArray.itemBeBinMatrix(2,idx); %��iBin��item����˳��Seq
    fprintf('bin �Ŀ�+��+��Ϊ: ' );
    fprintf(' %d  ',uniBinDataMatrix);
    fprintf('\n');
    fprintf('bin %d ��ʣ���+ʣ�೤Ϊ:  ',iBin);fprintf('\n');
    fprintf('( %d ) ',da.BinSArray.LW(:,iBin));fprintf('\n');
    fprintf('\n');
    
    fprintf('bin %d ���� original item ������{˳��}(����)[��ת��־]{����}Ϊ  \n  ',iBin);
    fprintf('%d ',idx);fprintf('\n');
    fprintf('{%d} ',idxSeq);fprintf('\n');
    fprintf(' (%d %d %d) ', da.ItemArray.LWH(1:nDim,idx));fprintf('\n');
    fprintf(' [%d]     ', da.ItemArray.Rotaed(:,idx));fprintf('\n');
    fprintf(' {%d %d %d} ', da.ItemArray.CoordItemBin(:,idx));fprintf('\n');
    fprintf('\n');
    
    [~,idxLU] = find(da.LUArray.LUBeBinMatrix(1,:)==iBin); %��iBin�µ�item������
    fprintf('bin %d ���� original LU ������{˳��}[item���](����)[��ת��־]{����}Ϊ  \n  ',iBin);
    idxLUSeq = da.LUArray.LUBeBinMatrix(2,idxLU); %��iBin��item����˳��Seq
    idxLUItem = da.LUArray.LUBeItemArray(1,idxLU);
    fprintf('%d ',idxLU);fprintf('\n');
    fprintf('{%d} ',idxLUSeq);fprintf('\n');
    fprintf('[%d] ',idxLUItem);fprintf('\n');
    fprintf(' (%d %d %d) ', da.LUArray.LWH(1:nDim,idxLU));fprintf('\n');
    fprintf(' [%d]     ', da.LUArray.Rotaed(:,idxLU));fprintf('\n');
    fprintf(' {%d %d %d} ', da.LUArray.CoordLUBin(:,idxLU));fprintf('\n');
    fprintf('\n');

    % ������˳��չʾ
    % %     [~,x]=sort(da.LUArray.LUBeBinMatrix(2,idxLU));
    % %     idxLUSeq = idxLUSeq(x); %��iBin��item����˳��Seq
    % %     idxLUItem = idxLUItem(x);
    % %     fprintf('%d ',idxLU);fprintf('\n');
    % %     fprintf('{%d} ',idxLUSeq);fprintf('\n');
    % %     fprintf('[%d] ',idxLUItem);fprintf('\n');
    % %     fprintf(' (%d %d %d) ', da.LUArray.LWH(1:nDim,idxLU(x)));fprintf('\n');
    % %     fprintf(' [%d]     ', da.LUArray.LURotaFlag(:,idxLU(x)));fprintf('\n');
    % %     fprintf(' {%d %d %d} ', da.LUArray.CoordLUBin(:,idxLU(x)));fprintf('\n');
    % %     fprintf('\n');
end
end
