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
        
        [d.Item,d.LU] = cpuItem(d.Item,d.LU,d.Veh);

                        %         pgon = getPolyshape(d.Item.LWH);    %          figure; plot(pgon);  axis equal;  axis ([0 maxX 0 maxY]);
        %% �����½�
%         lb = computerLB(d.Item,d.Veh);   fprintf('LB = %d \n', lb); %��ĳ��bin����Ϊ׼
        %% ����ʽItem��Strip���㷨
        [d.LU,d.Item,d.Strip] = HItemToStrip(d.LU,d.Item,d.Veh,p);                %         printstruct(d);   %  printstruct(d.Item);
        
        [d.Strip,d.LU] = cpuStrip(d.Strip,d.Item,d.LU,d.Veh);

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
        
        % ����ͷ����2: ÿ��ʣ��stripȫ���ڱȽ��� better than ����1
        [d.Strip,d.Bin]= HreStripToBin(d.Bin,d.Strip,d.Item,d.LU,d.Veh,p);
        
        % ����ͷ����1: ÿ��Bin��strip�Ƚ���
        %     [d.Strip,d.Bin]= HreStripToEachBin(d.Bin,d.Strip,d.Item,d.LU,d.Veh,p);
        
        %% Item��bin����Ϣ��ȡ:
        [d.LU,d.Item] = HItemToBin(d.LU,d.Item,d.Strip);          printstruct(d.Item);
        [d.Bin,d.LU] = cpuBin(d.Bin,d.Strip,d.Item,d.LU,d.Veh);

%     printOut(d.Bin,d.Strip,d.Item,d.LU,d.Veh); %����,��ʱע��
        function printOut(Bin,Strip,Item,LU,Veh)
            nBin = size(Bin.LW,2);
            for iBin = 1: nBin
                [~,ibin] = find(Item.Item_Bin(1,:)==iBin); %��iBin�µ�item������
                idxSeq = Item.Item_Bin(2,ibin); %��iBin��item����˳��Seq
                fprintf('bin �Ŀ�+��+��Ϊ: ' );
                fprintf(' %d  ',Veh.LWH);
                fprintf('\n');
                fprintf('bin %d ��ʣ���+ʣ�೤Ϊ:  ',iBin);fprintf('\n');
                fprintf('( %d ) ',Bin.LW(:,iBin));fprintf('\n');
                fprintf('\n');

                fprintf('bin %d ���� original item ������{˳��}(����)[��ת��־]{����}Ϊ  \n  ',iBin);
                fprintf('%d ',ibin);fprintf('\n');
                fprintf('{%d} ',idxSeq);fprintf('\n');
                fprintf(' (%d %d %d) ', Item.LWH(1:3,ibin));fprintf('\n');
                fprintf(' [%d]     ', Item.Rotaed(:,ibin));fprintf('\n');
                fprintf(' {%d %d %d} ', Item.CoordItemBin(:,ibin));fprintf('\n');
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

%% �ֲ����� %%

%% ����1: ����ͷ����2: 
function [Strip,Bin] = HreStripToBin(Bin,Strip,Item,LU,Veh,p)
        % DONE: ����ͷ����2: 
        % Ŀ��: �����LU���������Ϊ��С���԰ڷų�ͷ;
        % ����: ÿ��Bin����ʣ��Strip�������򲢷ֱ�ִ��HStripToBin�㷨.
        nbBin = max(Strip.Strip_Bin(1,:));
        if nbBin>1
            ibin=2;
            while 1
                % 1 ��ȡf and fidx: �ų��׸�bin���ʣ��Strip�߼�ֵ and ����ֵ
                Strip.Strip_Bin(1, Strip.Strip_Bin(1,:) >= ibin ) = -1; % ���з���bin��bin��Ÿ�ֵ-1
                % �ҳ�δ��ȷ������Strip
                f = Strip.Strip_Bin(1,:) == -1; 
                fidx = find(f);
                if ~any(f), break; end
                
                % 2 ��ȡs1 : ����f ��ȡ��ʣ��Strip�ṹ��
                s1 = rmfield(Strip,{'striporder','Strip_Bin'});
                s1 = structfun(@(x) x(:,f),s1,'UniformOutput',false);
                
                % 3 ��ȡi1 : ����ʣ��Strip��ȡ��ʣ��Item�ṹ��
                i1  = Item;
                f2 = zeros(1,length(i1.Weight));
                iIdxs= find(f);
                for j=1:length(iIdxs)
                    f2 = f2+(i1.Item_Strip(1,:) == iIdxs(j));
                end                
                i1 = structfun(@(x) x(:,logical(f2)), i1, 'UniformOutput',false);
                
                % 4 ��ȡ��s2��b2 : ����s1��i1, ���¼���Strip.nbItem and ����ִ������ʽS2B
                [s1] = cpuStripnbItem(s1,i1,Veh);
                
                [s2,b2]= HStripToBin(s1,Veh,LU,p);
                
                % 5 ��ȡfs: ʣ��Strip�ڷź���׸�bin�ڵ� so �Ҳ�==1
                % �滻ԭʼStrip_Bin��ֵ
                fs = s2.Strip_Bin(1,:)==1;
                Strip.Strip_Bin(1,fidx(fs)) = ibin;
                Strip.Strip_Bin(2,fidx(fs)) = s2.Strip_Bin(2,fs);
  
                % 6 ��ֵ���
                Bin.Weight(ibin) = b2.Weight(1);
                Bin.LW(:,ibin) = b2.LW(1);
                ibin = ibin+1;                
            end
        end
end
        
%% ����2: ����ͷ����1: 
        % Ŀ��: �����LU���������Ϊ��С���԰ڷų�ͷ;
        % ����: ÿ��Bin���ֱ𵫸�Bin��Stripִ��HStripToBin�㷨.
function [Strip,Bin] = HreStripToEachBin(Bin,Strip,Item,LU,Veh,p)
       % ��Ժ��������������򲢰��ŵ���Bin��
        nbBin = max(Strip.Strip_Bin(1,:));
        if nbBin>100
            for ibin=2:nbBin
                f = Strip.Strip_Bin(1,:) == ibin;
                
                s1 = rmfield(Strip,{'striporder','Strip_Bin'});
                s1 = structfun(@(x) x(:,f),s1,'UniformOutput',false);
                i1  = Item;  
                
                iIdxs= find(f);
                f2 = zeros(1,length(i1.Weight));
                for j=1:length(iIdxs)
                    f2 = f2+(i1.Item_Strip(1,:) == iIdxs(j));
                end
                
                i1 = structfun(@(x) x(:,logical(f2)), i1, 'UniformOutput',false);
                
                % 1 ���¼���Strip.nbItem and ����ִ������ʽS2B
                [s1] = cpuStripnbItem(s1,i1,Veh);
                
                [s2,b2]= HStripToBin(s1,Veh,LU,p);
                
                % 2 �滻��ԭʼStrip_Bin��ibin�Ľ���˳��
                Strip.Strip_Bin(2,f) = s2.Strip_Bin(2,:);
                
                % 3 �������
                if any(s2.Strip_Bin(1,:)~=1), error('���·����bin�Ų���'); end
                if size(b2.Weight,2) ~=1, error('��������ͬ'); end
                if size(b2.LW,2) ~=1, error('������ͬ'); end
                if Bin.Weight(ibin) ~= b2.Weight, error('��������ͬ'); end
                if Bin.LW(:,ibin) ~= b2.LW, error('������ͬ'); end                
                Bin.Weight(ibin) = b2.Weight;
                Bin.LW(:,ibin) = b2.LW;
                
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
