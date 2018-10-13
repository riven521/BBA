function [d] = RunAlgorithm(d,p)
        global ISshuaiwei ISplotStrip ISreStripToBin %����ִ��HStripToBin����
                
        %% Ԥ����:����Input��������
        d = GcheckInput(d);    %���Բ���   printstruct(d,'sortfields',1,'PRINTCONTENTS',0);

        % ����Ԥ�����ص㣺��ȡLU.Rotaed,�����Ƿ�����
        [d.LU, d.Veh] = Gpreproc(d.LU, d.Veh,p.whichSortItemOrder); %������ ��cpuLU

        %% ����ʽ: LU��Item���㷨    
        [d.LU,d.Item] = HLUtoItem(d.LU,d.Veh);   %Item����ID������򣨵���һ�������仯˳��

        % Item.isNonMixed Item.isMixedTile isHeightFull
        [d.Item,d.LU] = cpuItem(d.Item,d.LU,d.Veh);        % printstruct(d,'sortfields',1,'PRINTCONTENTS',0);

        %% ����ʽ: Item��Strip���㷨        
        %  *******  *******  *******
        % 1 Item����: %1: SID ; 2: isNonMixed; һ�����濪ʼ: 
        %                       3: Longth/Height; 4:Width; 5: LID; (3,4,5,����һ��) 6: Height
        % 2 ��Item����˳������»�ɵ�Strip
        [d.LU,d.Item,d.Strip] = HItemToStrip(d.LU,d.Item,d.Veh,p);     %   printstruct(d);   %  printstruct(d.Item);

        % ����LU.LU_Strip, LU.CoordLUStrip
        [d.Strip,d.LU] = cpuStrip(d.Strip,d.Item,d.LU,d.Veh);

        if ISplotStrip==1,      figure(222);     plot3DStrip(d.LU,d.Item,d.Veh,'LU');        end    % igure(111);   plot3DStrip(d.LU,d.Item,d.Veh,'Item');  % ����LU.CoordLUBin %         figure(222);     plot3DStrip(d.LU,d.Item,d.Veh,'LU');         %    igure(111);          plot3DStrip(d.LU,d.Item,d.Veh,'Item'); 

        %% ����ʽ��Strip��Bin���㷨        
        % ********* 1 Strip����: % 1: SID 2: priorityofLID
        [d.Strip,d.Bin] = HStripToBin(d.Strip,d.Veh,d.LU,p);
        
        [d.LU,d.Item] = HItemToBin(d.LU,d.Item,d.Strip); % ����LU��Bin������and˳��   %  Item.Item_Bin  Item.CoordItemBin LU.LU_Bin LU.CoordLUBin
        [d.Bin,d.LU] = cpuBin(d.Bin,d.Strip,d.Item,d.LU,d.Veh);  %����Bin��������� % ����isTileNeed
        if ISplotStrip==1,      plot3DBPP(d,p);      end    % igure(111);   plot3DStrip(d.LU,d.Item,d.Veh,'Item');  % ����LU.CoordLUBin 
        
        % ********* 2 ����ͷ����2: ÿ��ʣ��stripȫ���ڱȽ��� better than ����1
        if ISreStripToBin==1,   [d.Strip,d.Bin] = HreStripToBin(d.Bin,d.Strip,d.Item,d.LU,d.Veh,p); end     % ����ͷ����1: ÿ��Bin��strip�Ƚ���      % [d.Strip,d.Bin]= HreStripToEachBin(d.Bin,d.Strip,d.Item,d.LU,d.Veh,p);
        
        [d.LU,d.Item] = HItemToBin(d.LU,d.Item,d.Strip); % ����LU��Bin������and˳��   %  Item.Item_Bin  Item.CoordItemBin LU.LU_Bin LU.CoordLUBin
        [d.Bin,d.LU] = cpuBin(d.Bin,d.Strip,d.Item,d.LU,d.Veh);  %����Bin��������� % ����isTileNeed
        if ISplotStrip==1,      plot3DBPP(d,p);      end    % igure(111);   plot3DStrip(d.LU,d.Item,d.Veh,'Item');  % ����LU.CoordLUBin
        
        % ********* 3 ����Strip��˦β�Ż� *********** �޸� Strip.Strip_Bin
        if ISshuaiwei==1,      [d.Strip] = HStripSW(d.Strip);      end
        
        [d.LU,d.Item] = HItemToBin(d.LU,d.Item,d.Strip); % ����LU��Bin������and˳��   %  Item.Item_Bin  Item.CoordItemBin LU.LU_Bin LU.CoordLUBin
        [d.Bin,d.LU] = cpuBin(d.Bin,d.Strip,d.Item,d.LU,d.Veh);  %����Bin��������� % ����isTileNeed
        if ISplotStrip==1,      plot3DBPP(d,p);      end    % igure(111);   plot3DStrip(d.LU,d.Item,d.Veh,'Item');  % ����LU.CoordLUBin
        
                %         printstruct(d.Item,'sortfields',1,'PRINTCONTENTS',0)  
                %         printstruct(d.Bin,'sortfields',1,'PRINTCONTENTS',1)  
                %         printstruct(d,'sortfields',1,'PRINTCONTENTS',0)  

        %% ��ӡ���
%     printOut(d.Bin,d.Strip,d.Item,d.LU,d.Veh); %����,��ʱע��
% %         function printOut(Bin,Strip,Item,LU,Veh)
% %             nBin = size(Bin.LW,2);
% %             for iBin = 1: nBin
% %                 [~,ibin] = find(Item.Item_Bin(1,:)==iBin); %��iBin�µ�item������
% %                 idxSeq = Item.Item_Bin(2,ibin); %��iBin��item����˳��Seq
% %                 fprintf('bin �Ŀ�+��+��Ϊ: ' );
% %                 fprintf(' %d  ',Veh.LWH);
% %                 fprintf('\n');
% %                 fprintf('bin %d ��ʣ���+ʣ�೤Ϊ:  ',iBin);fprintf('\n');
% %                 fprintf('( %d ) ',Bin.LW(:,iBin));fprintf('\n');
% %                 fprintf('\n');
% % 
% %                 fprintf('bin %d ���� original item ������{˳��}(����)[��ת��־]{����}Ϊ  \n  ',iBin);
% %                 fprintf('%d ',ibin);fprintf('\n');
% %                 fprintf('{%d} ',idxSeq);fprintf('\n');
% %                 fprintf(' (%d %d %d) ', Item.LWH(1:3,ibin));fprintf('\n');
% %                 fprintf(' [%d]     ', Item.Rotaed(:,ibin));fprintf('\n');
% %                 fprintf(' {%d %d %d} ', Item.CoordItemBin(:,ibin));fprintf('\n');
% %                 fprintf('\n');
% % 
% %                 [~,idxLU] = find(LU.LU_Bin(1,:)==iBin); %��iBin�µ�item������
% %                 fprintf('bin %d ���� original LU ������{˳��}[item���](����)[��ת��־]{����}Ϊ  \n  ',iBin);
% %                 idxLUSeq = LU.LU_Bin(2,idxLU); %��iBin��item����˳��Seq
% %                 idxLUItem = LU.LU_Item(1,idxLU);
% %                 fprintf('%d ',idxLU);fprintf('\n');
% %                 fprintf('{%d} ',idxLUSeq);fprintf('\n');
% %                 fprintf('[%d] ',idxLUItem);fprintf('\n');
% %                 fprintf(' (%d %d %d) ', LU.LWH(1:3,idxLU));fprintf('\n');
% %                 fprintf(' [%d]     ', LU.Rotaed(:,idxLU));fprintf('\n');
% %                 fprintf(' {%d %d %d} ', LU.CoordLUBin(:,idxLU));fprintf('\n');
% %                 fprintf('\n');
% % 
% %                 % ������˳��չʾ
% %                 % %     [~,x]=sort(LU.LU_Bin(2,idxLU));
% %                 % %     idxLUSeq = idxLUSeq(x); %��iBin��item����˳��Seq
% %                 % %     idxLUItem = idxLUItem(x);
% %                 % %     fprintf('%d ',idxLU);fprintf('\n');
% %                 % %     fprintf('{%d} ',idxLUSeq);fprintf('\n');
% %                 % %     fprintf('[%d] ',idxLUItem);fprintf('\n');
% %                 % %     fprintf(' (%d %d %d) ', LU.LWH(1:nDim,idxLU(x)));fprintf('\n');
% %                 % %     fprintf(' [%d]     ', LU.LURotaFlag(:,idxLU(x)));fprintf('\n');
% %                 % %     fprintf(' {%d %d %d} ', LU.CoordLUBin(:,idxLU(x)));fprintf('\n');
% %                 % %     fprintf('\n');
% %             end
% %         end
 
end

%% ������Ҫע��
        % %         %  *******  ����д���,�ͱ���; 
        %         ti = d.Item;  tl = d.LU;
        %         [d.LU,d.Item] = HItemToBin(d.LU,d.Item,d.Strip);
        %         [match, er1, er2] = comp_struct(ti,d.Item,1);
        %         list_struct(er1)        
        %         if isequal(ti,d.Item)~=1, warning('eeeee'); end
        %         if isequal(tl,d.LU)~=1, warning('eeeee'); end
        
                        %         pgon = getPolyshape(d.LU.LWH);    maxX = sum(d.LU.LWH(1,:))+10;    maxY = max(max(d.LU.LWH'))+10;  maxV = max(maxX,maxY);
                %                    plot(pgon);        axis equal;   axis ([0 maxX 0 maxY]);   
                %                    plot3Dshape(d.LU.LWH);
                
                                        %  pgon = getPolyshape(d.Item.LWH);    % figure; plot(pgon);  axis equal;  axis ([0 maxX 0 maxY]);
                                        
                                        % �����½�
                                        %         lb = computerLB(d.Item,d.Veh);   fprintf('LB = %d \n', lb); %��ĳ��bin����Ϊ׼
                                        
                                                % ��Strip�н���һ���Ҹ�>���Item����ѡ�񲢸�����Ӧ����
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
                                                
                                                                % ********** ƽ���Ż� **********
                                                                % %         if any(d.Strip.isWidthFull)
                                                                % %             fw = d.Strip.isWidthFull == 0;
                                                                % %             fm = d.Strip.isMixed == 0;
                                                                % %             %             fi = find(fw);
                                                                % %             nbbin = numel(d.Bin.Weight);
                                                                % %             %    for i=1:sum(fi)
                                                                % %             for i=1:nbbin
                                                                % %                 f2 = d.Strip.Strip_Bin(1,:) == i;
                                                                % %                 f3 = f2 & fw & fm;   %����bin i��, ����isWidthFull��strip            
                                                                % % %                 f2 = d.Strip.Strip_Bin(1,:) == d.Strip.Strip_Bin(1,fi(i)); % �뱾fi(i)����ͬһbin�ڵ�isWidthFull�ĸ���
                                                                % % %                 f3 = f2 & f;                
                                                                % %                 if sum(f3) > 1                    
                                                                % %                     error('�뱾fi(i)����ͬһbin�ڵ�isWidthFull�ĸ���>1');
                                                                % %                 end
                                                                % %                 fi3 = find(f3);
                                                                % %                 
                                                                % %                 d.Strip
                                                                % %                 d.Strip.LW
                                                                % %                 
                                                                % %                 % ��Ҫƽ�̵�strip���ڵ�item��lu��logicalֵ
                                                                % %                 fitem =d.Item.Item_Strip(1,:) == fi3;
                                                                % %                 flu =d.LU.LU_Strip(1,:) == fi3;
                                                                % %                 d.LU
                                                                % %                 d.Item
                                                                % %                 1
                                                                % %                 % �޸�LU.CoordLUBin��ֵ
                                                                % %                 % �޸�LU.LU_Item�仯�Ķ��ֵ
                                                                % %                 % �޸�Item.LWH��ֵ
                                                                % %                 % �޸�Item.isFull1��ֵ,Weight,isHeightFull,isWeightFine,Layer,PID,LID,SID
                                                                % %                 % CoordItemBin, CoordItemStrip, Item_Strip(���µ�����仯)
                                                                % %                 
                                                                % %                 %             d.Strip
                                                                % %                 d.Strip.isMixed
                                                                % %                 d.Strip.isAllPured
                                                                % %                 d.Strip.maxHeight
                                                                % %             end    
                                                                % %         1
                                                                % %         end

%% �ֲ����� %%

%% ����1: ����ͷ����2: V1: s1 :ȫ��Ϊ����Strip�ų��Ѱ���Bin���ʣ��Strip;;cpuStripnbItemΪ����Strip����
function [Strip,Bin] = HreStripToBin(Bin,Strip,Item,LU,Veh,p)
        % DONE: ����ͷ����2: 
        % Ŀ��: �����LU���������Ϊ��С���԰ڷų�ͷ;
        % ����: ÿ��Bin�����ų�ǰ��Strip���ʣ��Strip�������򲢷ֱ�ִ��HStripToBin�㷨.
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
                [s1] = cpuStripnbItem(s1,i1,LU);
                %[s1,~] = cpuStrip(s1,i1,LU,Veh);  % �˺�������������Ӻ���, �Ƿ��б�Ҫ����ȫ��, ����ֻ���������Ӻ��� TODO
                
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

%% ����1: ����ͷ����2: V2: s1 :HStripToBin����ʱ:Ϊ����Strip�ų��Ѱ���Bin���ʣ��Strip;cpuStripnbItemΪȫ��Strip����
% % function [Strip,Bin] = HreStripToBin(Bin,Strip,Item,LU,Veh,p)
% %         % DONE: ����ͷ����2: 
% %         % Ŀ��: �����LU���������Ϊ��С���԰ڷų�ͷ;
% %         % ����: ÿ��Bin�����ų�ǰ��Strip���ʣ��Strip�������򲢷ֱ�ִ��HStripToBin�㷨.
% %         nbBin = max(Strip.Strip_Bin(1,:));
% %         if nbBin>1
% %             ibin=2;
% %             while 1
% %                 % 1 ��ȡf and fidx: �ų��׸�bin���ʣ��Strip�߼�ֵ and ����ֵ
% %                 Strip.Strip_Bin(1, Strip.Strip_Bin(1,:) == 1 ) = -1; 
% %                 Strip.Strip_Bin(1, Strip.Strip_Bin(1,:) >= ibin ) = 0; % ���з���bin��bin��Ÿ�ֵ-1
% %                 % �ҳ�δ��ȷ������Strip
% %                 f = Strip.Strip_Bin(1,:) == 0; 
% %                 fidx=find(f);
% %                 if ~any(f), break; end
% %                 
% %                 % 2 ��ȡs1 : �ų�����field,����f���½ṹ��
% %                 s1 = rmfield(Strip,{'striporder','Strip_Bin'});
% %                 s1.f = f;
% %                 
% %                 % 3 ��ȡi1 : ����ʣ��Strip��ȡ��ʣ��Item�ṹ��
% %                 i1  = Item;
% %                 % 4 ��ȡ��s2��b2 : ����s1��i1, ���¼���Strip.nbItem and ����ִ������ʽS2B
% %                 [s1] = cpuStripnbItem(s1,i1,LU);
% %                 %[s1,~] = cpuStrip(s1,i1,LU,Veh);  % �˺�������������Ӻ���, �Ƿ��б�Ҫ����ȫ��, ����ֻ���������Ӻ��� TODO
% %                 
% %                 s1 = structfun(@(x) x(:,f),s1,'UniformOutput',false); %����HStripToBinֻ�ܸ����ֵ� x(:,f)
% %                 [s2,b2]= HStripToBin(s1,Veh,LU,p);
% %                 
% %                 % 5 ��ȡfs: ʣ��Strip�ڷź���׸�bin�ڵ� so �Ҳ�==1
% %                 % �滻ԭʼStrip_Bin��ֵ
% %                 fs = s2.Strip_Bin(1,:)==1; %�ҳ���i1��bin���߼�ֵ                
% %                Strip.Strip_Bin(1,fidx(fs)) = ibin;
% %                Strip.Strip_Bin(2,fidx(fs)) = s2.Strip_Bin(2,fs);
% % 
% %                 % 6 ��ֵ���
% %                 Bin.Weight(ibin) = b2.Weight(1);
% %                 Bin.LW(:,ibin) = b2.LW(1);
% %                 ibin = ibin+1;                
% %             end
% %             Strip.Strip_Bin(1, Strip.Strip_Bin(1,:) == -1 ) = 1;             
% %             % sortrows(Strip.Strip_Bin',[1,2],{'ascend','ascend'})'
% %         end
% % end

        
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
