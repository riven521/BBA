function [LU,Item] = HItemToBin(LU,Item,Strip)
% ��Ҫ����:Item����Bin��
%% ��ȡitem��bin�Ĺ�ϵ Item_Bin  + ���� CoordItemBin
%% ��ȡLU��bin�Ĺ�ϵ LU_Bin  + ���� CoordLUBin

%% ��ʼ��
%nDim = 3;
nItem = size(Item.LWH,2); %����ʹ�õ�Item������
nLU= size(LU.LWH,2); %����ʹ�õ�LU������
    
% ����������ʼ��
Item.Item_Bin=zeros(2,nItem);
Item.CoordItemBin=zeros(2,nItem);
LU.LU_Bin=zeros(2,nLU);
LU.CoordLUBin=zeros(2,nLU); % LU.LU_Strip=zeros(2,nLU); %% NOTE: ��plot3DBPPʱ ����, ���ִ˴���ֵ����, so ɾ�������
[LU.leftA,LU.bottomA,LU.rightA,LU.topA] = deal(zeros(1,nLU));

iiStrip = 0;
% ѭ��ÿ��bin 5555555555 �ǳ���Ҫ�ĺ��� 55555555555555
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

        %%
        % ƽ��Ҫ�� 1:����strip�����ж�� 2: �����ǻ�ϵ� 3: Item�Ĳ���Ҫ����,����>=2
%         if Strip.isSingleItem(thisStrip) ~= 1 && Strip.isMixed(thisStrip) ~= 1
% %             flagItem = Item.Item_Strip(1,:) == thisStrip;
% %             idxItem = find(flagItem)
% %             Item.Layer(flagItem)
% %             difLayer = diff( Item.Layer(flagItem))
% %             di = find(abs(difLayer)>=2);
% %              if ~isempty(di) && ~isscalar(di), error('EEE'); end
% %              flagItem(di)
% %              flagItem(di+1)
% %             if abs (diff( Item.Layer(flagItem)) ) > 1 % ���Item�Ĳ��� >= 2 ����Ҫƽ��
% %                 LU.LU_Item(1, : )
% %                 LU.LU_Item(2, : )
% %                 Item.Layer(flagItem)
% %                 for i=1:sum(flagItem)
% %                     
% %                 end
%             end
%             [~,y] = find(Item.Item_Strip(1,:) == thisStrip); %y��Item������,��ǰstrip��
%             for i=1:length(y) % ÿ��Item��ѭ��
%                 
%                 sum(LU.LU_Item(1,:) == y(i))
%                 Item.Item_Strip(2,flagItem)
%                 Item.LWH(:,flagItem)
% 
%             end
%             Strip.isFull(thisStrip)
%             Strip.isMixed(thisStrip)
%             Strip.isSingleItem(thisStrip)
%             Strip.loadingrateLimit(thisStrip)
%         end
        
        % ƽ��Ҫ�� 1:����stripֻ��1��, ����������, �����ܱ�ƽ��
        
        
        % ƽ��Ҫ�� 1:����stripֻ�ж��, ����Ŀ�����ǰƽ��
        
 %%       
        if ~isscalar(thisStrip), 
            error('�������'); 
        end
        nbItem=numel(find(Item.Item_Strip(1,:)==thisStrip));
        % ѭ��stripÿ��item,  �ӵ�һ��˳��ʼ�𲽻�ȡITEM��BIN�ڵ�����, ����ITEM_STRIP��˳��
        for iItem=1:nbItem %ͬһstrip�� ѭ��ȡֵiItem
            [~,thisItem] = find(Item.Item_Strip(1,:)==thisStrip & Item.Item_Strip(2,:)==iItem );    %��Strip�ڵ�iItem��item
            if ~isscalar(thisItem), error('�������');  end
            
            % ����itemBeBinMatrix 555
            Item.Item_Bin(1,thisItem)=iBin;  % ��thisItem����iBin�����
            Item.Item_Bin(2,thisItem)=tmpItemSeq;  % ��thisItem�ڱ�iBin�ڵ�˳�� ��1��ʼ
            tmpItemSeq=tmpItemSeq+1;     
            % ����CoordItemBin(��strip���¶���Item) 555 HAVE checked?
            Item.CoordItemBin(1,thisItem) = Item.CoordItemStrip(1,thisItem); % 55 ��thisItem��x���겻��
                % iStrip = 1�� �׸�strip�߶�һ��Ϊ0��iStrip=2���߶�һ����Ϊ0
            Item.CoordItemBin(2,thisItem) = sum(Strip.LW(2,tmpStrip)); % 555 ����ԭ��2��thisStrip˳��ȷ�� % 555 LWStripSort����ԭ��:iStripÿ��bin����ȡ��;LWStripSort�����е�Strip %             Item.CoordItemBin(2,thisItem) = sum(tmpLWStrip(2,1:iStrip-1)); 
            Item.CoordItemBin(3,thisItem) = 0; %Item����0��ʼ
            
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
                
                % ����CoordLUBin 555 HAVE checked
                LU.CoordLUBin(1,thisLU) = Item.CoordItemBin(1,thisItem);
                LU.CoordLUBin(2,thisLU) = Item.CoordItemBin(2,thisItem);                
                LU.CoordLUBin(3,thisLU) = sum(LU.LWH(3,tmpLU));  %555 ����ԭ��ͬ��% LU.CoordLUBin(3,thisLU) = sum(tmpLWLU(3,1:iLU-1));
                tmpLU = [tmpLU thisLU];
                
                % ������LU���ĸ�����
                LU.leftA(thisLU) = LU.CoordLUBin(1,thisLU);
                LU.bottomA(thisLU) = LU.CoordLUBin(2,thisLU);
                LU.rightA(thisLU) = LU.leftA(thisLU) + LU.LWH(1,thisLU);
                LU.topA(thisLU) = LU.bottomA(thisLU) + LU.LWH(2,thisLU);        
                
            end
        end
        tmpStrip = [tmpStrip thisStrip];
    end
end

            %NOTE: �˴���LU��ITEM�׶ε���ת��ֵ��LU��   
            %GcheekInput����ĶԱȵ�����ԭʼֵ���˴�Ҫ�滻Ϊ��һ�����ԭʼֵ.
            % LU.Rotaed = LURotaed; % ���治��Ҫ��, ��ΪLUArray.Rotaed�Ѿ���HItemToStrip���¹���


% printstruct(d);
%     printscript();
    
% �����Ҫ���:��ô�1��ʼÿ��bin������item����
% CoordItemBin Item_Bin
% %     function printscript()
% %         nBin = max(Item.Item_Bin(1,:));
% %         for iBin = 1: nBin
% %             [~,idx] = find(Item.Item_Bin(1,:)==iBin); %��iBin�µ�item������
% %             idxSeq = Item.Item_Bin(2,idx); %��iBin��item����˳��Seq
% %             fprintf('bin �Ŀ�+��+��Ϊ: ' );
% %             fprintf(' %d  ',Veh.LWH);
% %             fprintf('\n');
% %             fprintf('bin %d ��ʣ���+ʣ�೤Ϊ:  ',iBin);fprintf('\n');
% %             fprintf('( %d ) ',Bin.LW(:,iBin));fprintf('\n');
% %             fprintf('\n');
% %             
% %             fprintf('bin %d ���� original item ������{˳��}(����)[��ת��־]{����}Ϊ  \n  ',iBin);
% %             fprintf('%d ',idx);fprintf('\n');
% %             fprintf('{%d} ',idxSeq);fprintf('\n');
% %             fprintf(' (%d %d %d) ', Item.LWH(1:nDim,idx));fprintf('\n');
% %             fprintf(' [%d]     ', Item.Rotaed(:,idx));fprintf('\n');
% %             fprintf(' {%d %d %d} ', Item.CoordItemBin(:,idx));fprintf('\n');
% %             fprintf('\n');
% %             
% %             [~,idxLU] = find(LU.LU_Bin(1,:)==iBin); %��iBin�µ�item������
% %             fprintf('bin %d ���� original LU ������{˳��}[item���](����)[��ת��־]{����}Ϊ  \n  ',iBin);
% %             idxLUSeq = LU.LU_Bin(2,idxLU); %��iBin��item����˳��Seq
% %             idxLUItem = LU.LU_Item(1,idxLU);
% %             fprintf('%d ',idxLU);fprintf('\n');
% %             fprintf('{%d} ',idxLUSeq);fprintf('\n');
% %             fprintf('[%d] ',idxLUItem);fprintf('\n');
% %             fprintf(' (%d %d %d) ', LU.LWH(1:nDim,idxLU));fprintf('\n');
% %             fprintf(' [%d]     ', LU.Rotaed(:,idxLU));fprintf('\n');
% %             fprintf(' {%d %d %d} ', LU.CoordLUBin(:,idxLU));fprintf('\n');
% %             fprintf('\n');
% %             
% %             % ������˳��չʾ
% %             % %     [~,x]=sort(LU.LU_Bin(2,idxLU));
% %             % %     idxLUSeq = idxLUSeq(x); %��iBin��item����˳��Seq
% %             % %     idxLUItem = idxLUItem(x);
% %             % %     fprintf('%d ',idxLU);fprintf('\n');
% %             % %     fprintf('{%d} ',idxLUSeq);fprintf('\n');
% %             % %     fprintf('[%d] ',idxLUItem);fprintf('\n');
% %             % %     fprintf(' (%d %d %d) ', LU.LWH(1:nDim,idxLU(x)));fprintf('\n');
% %             % %     fprintf(' [%d]     ', LU.LURotaFlag(:,idxLU(x)));fprintf('\n');
% %             % %     fprintf(' {%d %d %d} ', LU.CoordLUBin(:,idxLU(x)));fprintf('\n');
% %             % %     fprintf('\n');
% %         end
% %     end
end
