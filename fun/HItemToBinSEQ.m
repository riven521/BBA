%% V2
function [LU,Item,Strip] = HItemToBinSEQ(LU,Item,Strip)

nItem = size(Item.LWH,2); %����ʹ�õ�Item������
nLU= size(LU.LWH,2); %����ʹ�õ�LU������

% ����������ʼ��
Item.Item_Bin=zeros(2,nItem);
LU.LU_Bin=zeros(2,nLU);

[LU.leftA,LU.bottomA,LU.rightA,LU.topA] = deal(zeros(1,nLU)); %LU���ĸ�����

iiStrip = 0;
% ѭ��ÿ��bin   
for iBin=1:max(Strip.Strip_Bin(1,:))
    tmpItemSeq=1;  %ÿ��bin��˳���1��ʼ
    tmpLUSeq=1;
    tmpStrip = []; %����CoordItemBin��(��)��ʹ�� %       tmpLWStrip = Strip.LW(:,Strip.Strip_Bin(1,:)==iBin); %��bin��ʣ���Ⱥͳ�(��)��
    
    nbStrip=numel(find(Strip.Strip_Bin(1,:)==iBin));    %nbStrip����iBin��strip����   
    
    % ѭ��bin��ÿ��strip
    for iStrip=1:nbStrip %�Ӱ��ŵĵ�һ��Strip��ʼ�𲽷���
        iiStrip = iiStrip+1;
        tmpLUSeqinStrip=1;
        [~,thisStrip] = find(Strip.Strip_Bin(1,:)==iBin & Strip.Strip_Bin(2,:)==iStrip ); %��bin�ڵ�iStrip��strip        

 
        if ~isscalar(thisStrip),     error('�������');     end
        
        nbItem=numel(find(Item.Item_Strip(1,:)==thisStrip));
        % ѭ��stripÿ��item,  �ӵ�һ��˳��ʼ�𲽻�ȡITEM��BIN�ڵ�����, ����ITEM_STRIP��˳��
        for iItem=1:nbItem %ͬһstrip�� ѭ��ȡֵiItem
            
            [~,thisItem] = find(Item.Item_Strip(1,:)==thisStrip & Item.Item_Strip(2,:)==iItem );    %��Strip�ڵ�iItem��item
            if ~isscalar(thisItem), error('�������');  end
            
            % ����itemBeBinMatrix 555
            Item.Item_Bin(1,thisItem)=iBin;  % ��thisItem����iBin�����
            Item.Item_Bin(2,thisItem)=tmpItemSeq;  % ��thisItem�ڱ�iBin�ڵ�˳�� ��1��ʼ
            tmpItemSeq=tmpItemSeq+1;    
                    
            % ���Ӷ�LU�ĸ���
            tmpLU = [];             %����CoordLUBin�߶�ʹ�� % tmpLWLU = LU.LWH(:,LU.LU_Item(1,:)==thisItem); 
            nbLU=numel(find(LU.LU_Item(1,:)==thisItem));     
            
            % ѭ��itemÿ��LU, �ӵ�һ��˳��ʼ�𲽻�ȡLU��BIN�ڵ�����
            for iLU=1:nbLU                
                [~,thisLU] = find(LU.LU_Item(1,:)==thisItem & LU.LU_Item(2,:)==iLU);
                if ~isscalar(thisLU), error('�������');  end
                
                                            % ����LU_Strip % NOTE: ��plot3DBPPʱ ����, ���ִ˴���ֵ����, so ɾ�������
                                            %                 LU.LU_Strip(1,thisLU)=iiStrip;
                                            %                 LU.LU_Strip(2,thisLU)=tmpLUSeqinStrip;

                tmpLUSeqinStrip=tmpLUSeqinStrip+1;
                                            % ����LURotaed 555  % LURotaed(1,thisLU)=Item.Rotaed(1,thisItem);
                % ����LU_Bin 555
                LU.LU_Bin(1,thisLU)=iBin;
                LU.LU_Bin(2,thisLU)=tmpLUSeq;
                tmpLUSeq=tmpLUSeq+1;
                
          tmpLU = [tmpLU thisLU];
                
           
            end
        end  % EOF ITEM
        tmpStrip = [tmpStrip thisStrip];
    end  %EOF STRIP
end  %EOF BIN




end % EOF function
