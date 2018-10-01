function [d] = RunAlgorithm(d,p)
        
        %% Ԥ����:����Input��������
        d = GcheckInput(d); %���Բ��� 
        
                %         pgon = getPolyshape(d.LU.LWH);    maxX = sum(d.LU.LWH(1,:))+10;    maxY = max(max(d.LU.LWH'))+10;  maxV = max(maxX,maxY);
                %                    plot(pgon);        axis equal;   axis ([0 maxX 0 maxY]);   
                %                    plot3Dshape(d.LU.LWH);

        % ����Ԥ�����ص㣺��ȡLU.Rotaed,�����Ƿ�����
        [d.LU, d.Veh] = Gpreproc(d.LU, d.Veh,p.whichSortItemOrder); %������
        
        %% ����ʽ: LU��Item���㷨    
        [d.LU,d.Item,d.ItemID] = HLUtoItem(d.LU,d.Veh); %Item����ID������򣨵���һ�������仯˳��
        
        [d.Item] = cpuItem(d.Item,d.LU,d.Veh);

                        %         printstruct(d.LU);         printstruct(d.Item);
                        %         pgon = getPolyshape(d.Item.LWH);    %          figure; plot(pgon);  axis equal;  axis ([0 maxX 0 maxY]);
        %% �����½�
%         lb = computerLB(d.Item,d.Veh);   fprintf('LB = %d \n', lb); %��ĳ��bin����Ϊ׼
        %% ����ʽItem��Strip���㷨
        [d.LU,d.Item,d.Strip] = HItemToStrip(d.LU,d.Item,d.Veh,p);                %         printstruct(d);   %  printstruct(d.Item);
        
        [d.Strip] = cpuStrip(d.Strip,d.Item,d.Veh);

                %% ��Strip�н���һ���Ҹ�>���Item����ѡ�񲢸�����Ӧ����
%         d = modifyStripWithOneItem(d);
%         function d = modifyStripWithOneItem(d)
%             stripheight = d.Strip.LW(2,:);
%             binwidth = d.Veh.LWH(1,1);
%             stripleftwidth = d.Strip.LW(1,:);
%             stripwidth = ( binwidth - stripleftwidth );
%             [tmpset] = find(stripheight > stripwidth);
%             if ~isempty(tmpset)
%                 if isscalar(tmpset) %�Ը�strip�����ڲ�����1��Item����,��������漰CoordItemStrip
%                     d.Strip.LW(:,tmpset) = [binwidth-stripheight(tmpset),stripwidth(tmpset)];    %strip�ĳ������
%                     %�ڲ�Item��itemRotaFlag���� 
%                     idxItem = find(d.Item.Item_Strip(1,:)==tmpset );
%                     if isscalar(idxItem)
%                         d.Item.itemRotaFlag(idxItem) = ~d.Item.Rotaed(idxItem);
%                     end                    
%                     %�ڲ�LU��LURotaFlag ���{ ߀δ�� %�ڲ�Item��CoordItemStrip���{                    
%                 end
%             end
%         end

        %% ����ʽ��Strip��Bin���㷨
        [d.Strip,d.Bin]= HStripToBin(d.Strip,d.Veh,d.LU,p);
        
        %% Item��bin����Ϣ��ȡ:
        [d.LU,d.Item] = HItemToBin(d.LU,d.Item,d.Strip);          printstruct(d.Item);

        [d.Bin] = cpuBin(d.Bin,d.Item,d.Veh);
        
     
%     printOut(d.Bin,d.Strip,d.Item,d.LU,d.Veh); %����,��ʱע��
        function printOut(Bin,Strip,Item,LU,Veh)
            nBin = size(Bin.LW,2);
            for iBin = 1: nBin
                [~,idx] = find(Item.Item_Bin(1,:)==iBin); %��iBin�µ�item������
                idxSeq = Item.Item_Bin(2,idx); %��iBin��item����˳��Seq
                fprintf('bin �Ŀ�+��+��Ϊ: ' );
                fprintf(' %d  ',Veh.LWH);
                fprintf('\n');
                fprintf('bin %d ��ʣ���+ʣ�೤Ϊ:  ',iBin);fprintf('\n');
                fprintf('( %d ) ',Bin.LW(:,iBin));fprintf('\n');
                fprintf('\n');

                fprintf('bin %d ���� original item ������{˳��}(����)[��ת��־]{����}Ϊ  \n  ',iBin);
                fprintf('%d ',idx);fprintf('\n');
                fprintf('{%d} ',idxSeq);fprintf('\n');
                fprintf(' (%d %d %d) ', Item.LWH(1:3,idx));fprintf('\n');
                fprintf(' [%d]     ', Item.Rotaed(:,idx));fprintf('\n');
                fprintf(' {%d %d %d} ', Item.CoordItemBin(:,idx));fprintf('\n');
                fprintf('\n');

                [~,idxLU] = find(LU.LU_Bin(1,:)==iBin); %��iBin�µ�item������
                fprintf('bin %d ���� original LU ������{˳��}[item���](����)[��ת��־]{����}Ϊ  \n  ',iBin);
                idxLUSeq = LU.LU_Bin(2,idxLU); %��iBin��item����˳��Seq
                idxLUItem = LU.LU_Item(1,idxLU);
                fprintf('%d ',idxLU);fprintf('\n');
                fprintf('{%d} ',idxLUSeq);fprintf('\n');
                fprintf('[%d] ',idxLUItem);fprintf('\n');
                fprintf(' (%d %d %d) ', LU.LWH(1:3,idxLU));fprintf('\n');
                fprintf(' [%d]     ', LU.Rotaed(:,idxLU));fprintf('\n');
                fprintf(' {%d %d %d} ', LU.CoordLUBin(:,idxLU));fprintf('\n');
                fprintf('\n');

                % ������˳��չʾ
                % %     [~,x]=sort(LU.LU_Bin(2,idxLU));
                % %     idxLUSeq = idxLUSeq(x); %��iBin��item����˳��Seq
                % %     idxLUItem = idxLUItem(x);
                % %     fprintf('%d ',idxLU);fprintf('\n');
                % %     fprintf('{%d} ',idxLUSeq);fprintf('\n');
                % %     fprintf('[%d] ',idxLUItem);fprintf('\n');
                % %     fprintf(' (%d %d %d) ', LU.LWH(1:nDim,idxLU(x)));fprintf('\n');
                % %     fprintf(' [%d]     ', LU.LURotaFlag(:,idxLU(x)));fprintf('\n');
                % %     fprintf(' {%d %d %d} ', LU.CoordLUBin(:,idxLU(x)));fprintf('\n');
                % %     fprintf('\n');
            end
        end
 
end

%% �����½�
% % % % %     function [ lb ] = computerLB(Item,Veh)
% % % % %         sum1 = sum(prod(Item.LWH,1));        
% % % % %         % todo �����ж��Ƿ����е�BinArray�����е�bin����ͬ�� ����� �����ִ��
% % % % %         sum2 = prod(Veh.LWH(:,1));
% % % % %         lb = ceil(sum1/sum2);
% % % % %         if lb <=0, error('EEE');end
% % % % %     end
    
    
%%     d.LU = doc(d.LU,d.Item,d.Strip,d.Bin);
% % % %     function LU = doc(LU,Item,Strip,Bin)
% % % %         
% % % %         LU.DOC=[LU.PID;LU.ID;LU.SID;zeros(size(LU.ID));zeros(size(LU.ID));...
% % % %             LU.LU_Item;LU.LU_Strip;LU.LU_Bin];
% % % %         
% % % %         nItem = size(Item.LWH,2);
% % % %         for iItem=1:nItem
% % % %             tmp = LU.DOC([1,2,3], LU.DOC(6,:) == iItem);
% % % %             Item.PID1(:,iItem) = num2cell(unique(tmp(1,:))',1); %unique(tmp(1,:))';
% % % %             Item.LID1(:,iItem) = num2cell(unique(tmp(2,:))',1);
% % % %             Item.SID1(:,iItem) =num2cell(unique(tmp(3,:))',1);
% % % %         end
% % % %         
% % % %         nStrip = size(Strip.LW,2);
% % % %         for iStrip=1:nStrip
% % % %             tmp = LU.DOC([1,2,3], LU.DOC(8,:) == iStrip);
% % % %             Strip.PID1(:,iStrip) = num2cell(unique(tmp(1,:))',1);
% % % %             Strip.LID1(:,iStrip) = num2cell(unique(tmp(2,:))',1);
% % % %             Strip.SID1(:,iStrip) = num2cell(unique(tmp(3,:))',1);
% % % %         end
% % % %         
% % % %         nBin = size(Bin.LW,2);
% % % %         for iBin=1:nBin
% % % %             tmp = LU.DOC([1,2,3], LU.DOC(10,:) == iBin);
% % % %             Bin.PID1(:,iBin) = num2cell(unique(tmp(1,:))',1);
% % % %             Bin.LID1(:,iBin) = num2cell(unique(tmp(2,:))',1);
% % % %             Bin.SID1(:,iBin) = num2cell(unique(tmp(3,:))',1);
% % % %         end
% % % %         
% % % %         LU.DOC
% % % %         LU.DOC([1,2,3,8],:)
% % % %     end
