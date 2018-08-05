function [d] = HItemToBin(d)
% ��Ҫ����:Item����Bin��
%% ��ȡitem��bin�Ĺ�ϵ Item_Bin  + ���� CoordItemBin
%% ��ȡLU��bin�Ĺ�ϵ LU_Bin  + ���� CoordLUBin
%% ��ȡLU��item�Ĺ�ϵ LURotaFlag

%% ��ʼ��
nDim = size(d.Item.LWH,1);  
nItem = size(d.Item.LWH,2); %����ʹ�õ�Item������
nLU= size(d.LU.LWH,2); %����ʹ�õ�LU������

%% �ṹ����ȡ
LU_Item = d.LU.LU_Item;
stripBeBinMatrix = d.Strip.Strip_Bin;
Item_Strip = d.Item.Item_Strip;

itemRotaed = d.Item.Rotaed;
% % itemRotaFlag = d.Item.itemRotaFlag;
CoordItemStrip = d.Item.CoordItemStrip;

LWStrip = d.Strip.LW; % LWStripSort = d.Strip.LW(:,d.Strip.striporder);
LWHLU = d.LU.LWH;
    
Item_Bin=zeros(2,nItem);
CoordItemBin=zeros(2,nItem);

LU_Bin=zeros(2,nLU);
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
        nbItem=numel(find(Item_Strip(1,:)==thisStrip));
        % ѭ��stripÿ��item
        for iItem=1:nbItem %ͬһstrip�� ѭ��ȡֵiItem
            [~,thisItem] = find(Item_Strip(1,:)==thisStrip & Item_Strip(2,:)==iItem );%��Strip�ڵ�iItem��item
            if ~isscalar(thisItem), error('�������');  end
            
            % ����itemBeBinMatrix 555
            Item_Bin(1,thisItem)=iBin;  % ��thisItem����iBin�����
            Item_Bin(2,thisItem)=tmpItemSeq;  % ��thisItem�ڱ�iBin�ڵ�˳�� ��1��ʼ
            tmpItemSeq=tmpItemSeq+1;     
            % ����CoordItemBin(��strip���¶���Item) 555 HAVE checked?
            CoordItemBin(1,thisItem) = CoordItemStrip(1,thisItem); % 55 ��thisItem��x���겻��
                % iStrip = 1�� �׸�strip�߶�һ��Ϊ0��iStrip=2���߶�һ����Ϊ0
            CoordItemBin(2,thisItem) = sum(LWStrip(2,tmpStrip)); % 555 ����ԭ��2��thisStrip˳��ȷ�� % 555 LWStripSort����ԭ��:iStripÿ��bin����ȡ��;LWStripSort�����е�Strip %             CoordItemBin(2,thisItem) = sum(tmpLWStrip(2,1:iStrip-1)); 
            CoordItemBin(3,thisItem) = 0; %Item����0��ʼ
            
            % ���Ӷ�LU�ĸ���
            tmpLU = []; %����CoordLUBin�߶�ʹ�� % tmpLWLU = LWHLU(:,LU_Item(1,:)==thisItem); %����
            nbLU=numel(find(LU_Item(1,:)==thisItem));            
            % ѭ��itemÿ��LU
            for iLU=1:nbLU
                [~,thisLU] = find(LU_Item(1,:)==thisItem & LU_Item(2,:)==iLU);
                if ~isscalar(thisLU), error('�������');  end
                
                % ����LURotaed 555
                LURotaed(1,thisLU)=itemRotaed(1,thisItem);
                % ����LUBeBinMatrix 555
                LU_Bin(1,thisLU)=iBin;
                LU_Bin(2,thisLU)=tmpLUSeq;
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
d.Item.Item_Bin=Item_Bin;
d.Item.CoordItemBin=CoordItemBin;

d.LU.LU_Bin= LU_Bin;
d.LU.CoordLUBin = CoordLUBin;

%NOTE: �˴���LU��ITEM�׶ε���ת��ֵ��LU��   
%GcheekInput����ĶԱȵ�����ԭʼֵ���˴�Ҫ�滻Ϊ��һ�����ԭʼֵ.
%% % ���治��Ҫ��, ��ΪLUArray.Rotaed�Ѿ���HItemToStrip���¹���
% d.LU.Rotaed = LURotaed;


% printstruct(d);
% �����Ҫ���:��ô�1��ʼÿ��bin������item����
% CoordItemBin Item_Bin
for iBin = 1:max(d.Item.Item_Bin(1,:))
    [~,idx] = find(d.Item.Item_Bin(1,:)==iBin); %��iBin�µ�item������
    idxSeq = d.Item.Item_Bin(2,idx); %��iBin��item����˳��Seq
    fprintf('bin �Ŀ�+��+��Ϊ: ' );
    fprintf(' %d  ',d.Veh.LWH);
    fprintf('\n');
    fprintf('bin %d ��ʣ���+ʣ�೤Ϊ:  ',iBin);fprintf('\n');
    fprintf('( %d ) ',d.Bin.LW(:,iBin));fprintf('\n');
    fprintf('\n');
    
    fprintf('bin %d ���� original item ������{˳��}(����)[��ת��־]{����}Ϊ  \n  ',iBin);
    fprintf('%d ',idx);fprintf('\n');
    fprintf('{%d} ',idxSeq);fprintf('\n');
    fprintf(' (%d %d %d) ', d.Item.LWH(1:nDim,idx));fprintf('\n');
    fprintf(' [%d]     ', d.Item.Rotaed(:,idx));fprintf('\n');
    fprintf(' {%d %d %d} ', d.Item.CoordItemBin(:,idx));fprintf('\n');
    fprintf('\n');
    
    [~,idxLU] = find(d.LU.LU_Bin(1,:)==iBin); %��iBin�µ�item������
    fprintf('bin %d ���� original LU ������{˳��}[item���](����)[��ת��־]{����}Ϊ  \n  ',iBin);
    idxLUSeq = d.LU.LU_Bin(2,idxLU); %��iBin��item����˳��Seq
    idxLUItem = d.LU.LU_Item(1,idxLU);
    fprintf('%d ',idxLU);fprintf('\n');
    fprintf('{%d} ',idxLUSeq);fprintf('\n');
    fprintf('[%d] ',idxLUItem);fprintf('\n');
    fprintf(' (%d %d %d) ', d.LU.LWH(1:nDim,idxLU));fprintf('\n');
    fprintf(' [%d]     ', d.LU.Rotaed(:,idxLU));fprintf('\n');
    fprintf(' {%d %d %d} ', d.LU.CoordLUBin(:,idxLU));fprintf('\n');
    fprintf('\n');

    % ������˳��չʾ
    % %     [~,x]=sort(d.LU.LU_Bin(2,idxLU));
    % %     idxLUSeq = idxLUSeq(x); %��iBin��item����˳��Seq
    % %     idxLUItem = idxLUItem(x);
    % %     fprintf('%d ',idxLU);fprintf('\n');
    % %     fprintf('{%d} ',idxLUSeq);fprintf('\n');
    % %     fprintf('[%d] ',idxLUItem);fprintf('\n');
    % %     fprintf(' (%d %d %d) ', d.LU.LWH(1:nDim,idxLU(x)));fprintf('\n');
    % %     fprintf(' [%d]     ', d.LU.LURotaFlag(:,idxLU(x)));fprintf('\n');
    % %     fprintf(' {%d %d %d} ', d.LU.CoordLUBin(:,idxLU(x)));fprintf('\n');
    % %     fprintf('\n');
end
end
