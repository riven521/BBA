function [d] = RunAlgorithm(d,p)
        
        % ����Input��������
%         d.LU.isRota = [1 0 0]
        printstruct(d.LU);
        d = GcheckInput(d); %���Բ��� 
        pgon = getPolyshape(d.LU.LWH);    maxX = sum(d.LU.LWH(1,:))+10;    maxY = max(max(d.LU.LWH'))+10;  maxV = max(maxX,maxY);
%                   plot(pgon);        axis equal;   axis ([0 maxX 0 maxY]);   
%                   plot3Dshape(d.LU.LWH);

        % ����Ԥ�����ص㣺��ȡLU.Rotaed,�����Ƿ�����
        [d.LU, d.Veh] = Gpreproc(d.LU, d.Veh,p.whichSortItemOrder); %�����Բ��� 
        
        %% ����ʽ: LU��Item���㷨    
        [d.LU,d.Item,d.ItemID] = HLUtoItem(d.LU,d.Veh); %Item����ID������򣨵���һ�������仯˳��
        printstruct(d.LU);
        printstruct(d.Item);
        pgon = getPolyshape(d.Item.LWH);
%          figure; plot(pgon);  axis equal;  axis ([0 maxX 0 maxY]);
        %% �����½�
%         lb = computerLB(d.Item,d.Veh);   fprintf('LB = %d \n', lb); %��ĳ��bin����Ϊ׼
        %% ����ʽ��Item��Strip���㷨
%         printstruct(d);
%         printstruct(d.Item);
        [d.LU,d.Item,d.Strip] = HItemToStrip(d.LU,d.Item,d.Veh,p);
        
        %% ����stripװ����
%         printstruct(d);
        d = computeLoadingRateStrip(d);
        function d = computeLoadingRateStrip(d)
            % ��ʼ��
            nStrip = size(d.Strip.LW,2);
            d.Strip.Stripvolume = zeros(1,nStrip);
            d.Strip.StripvolumeLimit = zeros(1,nStrip);
            d.Strip.Itemvolume = zeros(1,nStrip);
            d.Strip.loadingrate = zeros(1,nStrip);
            d.Strip.loadingrateLimit = zeros(1,nStrip);
            
            % ����ÿ��strip��װ����
            %ÿ��strip�Ŀ������ = �߶�*���(�����Ŀ��)
            d.Strip.Stripvolume = d.Strip.LW(2,:)*d.Veh.LWH(1,1);
            %ÿ��strip�����޿������ = �߶�*���(stripʹ�ÿ��=�������-stripʣ����)
            d.Strip.StripvolumeLimit = d.Strip.LW(2,:) .* (d.Veh.LWH(1,1) - d.Strip.LW(1,:));
            a = d.Item.LWH;
            b = d.Item.Item_Strip;
            for iStrip =1:nStrip
                %ÿ��strip��װ�����
                d.Strip.Itemvolume(iStrip)= sum(a(1, (b(1,:)==iStrip)) .* a(2, (b(1,:)==iStrip)));
            end
            %ÿ��strip��װ�ر���
            d.Strip.loadingrate =  d.Strip.Itemvolume ./ d.Strip.Stripvolume;
            %ÿ��strip������װ�ر���
            d.Strip.loadingrateLimit =  d.Strip.Itemvolume ./ d.Strip.StripvolumeLimit;
        end
        %% ��Strip�н���һ���Ҹ�>���Item����ѡ�񲢸�����Ӧ����
%         d = modifyStripWithOneItem(d);
        function d = modifyStripWithOneItem(d)
            stripheight = d.Strip.LW(2,:);
            binwidth = d.Veh.LWH(1,1);
            stripleftwidth = d.Strip.LW(1,:);
            stripwidth = ( binwidth - stripleftwidth );
            [tmpset] = find(stripheight > stripwidth);
            if ~isempty(tmpset)
                if isscalar(tmpset) %�Ը�strip�����ڲ�����1��Item����,��������漰CoordItemStrip
                    d.Strip.LW(:,tmpset) = [binwidth-stripheight(tmpset),stripwidth(tmpset)];    %strip�ĳ������
                    %�ڲ�Item��itemRotaFlag���� 
                    idxItem = find(d.Item.Item_Strip(1,:)==tmpset );
                    if isscalar(idxItem)
                        d.Item.itemRotaFlag(idxItem) = ~d.Item.Rotaed(idxItem);
                    end                    
                    %�ڲ�LU��LURotaFlag ���{ ߀δ�� %�ڲ�Item��CoordItemStrip���{                    
                end
            end
        end
        %% ����ʽ��Strip��Bin���㷨
%         printstruct(d);
        [d.Strip,d.Bin]= HStripToBin(d.Strip,d.Veh,d.LU,p);
printstruct(d.Strip)
        %% Item��bin����Ϣ��ȡ:
%         printstruct(d);
%         [d] = HItemToBin(d);
        [d.LU,d.Item] = HItemToBin(d.LU,d.Item,d.Strip);
        printstruct(d.Item);
         %% ����binװ����
         % ItemloadingrateLimit - ÿ��bin��Item�������/ÿ��binȥ��ʣ���ߺ�������
         % Itemloadingrate - ÿ��bin��Item�������/ÿ��bin���������
         d = computeLoadingRate2DBin(d);
        function d = computeLoadingRate2DBin(d)
            % ��ʼ��
            nBin = size(d.Bin.LW,2);
            d.Bin.Binarea = zeros(1,nBin);
            d.Bin.Itemarea = zeros(1,nBin);
            d.Bin.Itemloadingrate = zeros(1,nBin);
            d.Bin.ItemloadingrateLimit = zeros(1,nBin);
            % ����ÿ��Bin��װ����            
            BinWidth = d.Veh.LWH(1,1);
            BinHeight = d.Veh.LWH(2,1);
            BinArea = BinWidth .* BinHeight;
            %ÿ��Bin�Ŀ������ = �����߶�*�������
            d.Bin.Binarea = repmat(BinArea,1,nBin);            
            %ÿ��Bin �����޿������ = ���(binʹ�ÿ��=�������-binʣ����) *�߶�(binʹ�ø߶�=�����߶�-binʣ��߶�)
            d.Bin.BinareaLimit = (BinWidth - d.Bin.LW(1,:)) .* (BinHeight - d.Bin.LW(2,:));
            
            a = d.Item.LWH;
            b = d.Item.Item_Bin;
            for iBin =1:nBin
                %ÿ��Bin��װ�����
                d.Bin.Itemarea(iBin)= sum(a(1, (b(1,:)==iBin)) .* a(2, (b(1,:)==iBin)));
            end
            %ÿ��bin��װ�ر���
            d.Bin.loadingrate =  d.Bin.Itemarea ./ d.Bin.Binarea;
            %ÿ��bin������װ�ر���
            d.Bin.loadingrateLimit =  d.Bin.Itemarea ./ d.Bin.BinareaLimit;
        end
        printstruct(d);

        
% % %     d.LU = doc(d.LU,d.Item,d.Strip,d.Bin);
    function LU = doc(LU,Item,Strip,Bin)
        
        LU.DOC=[LU.PID;LU.ID;LU.SID;zeros(size(LU.ID));zeros(size(LU.ID));...
            LU.LU_Item;LU.LU_Strip;LU.LU_Bin];
        
        nItem = size(Item.LWH,2);
        for iItem=1:nItem
            tmp = LU.DOC([1,2,3], LU.DOC(6,:) == iItem);
            Item.PID1(:,iItem) = num2cell(unique(tmp(1,:))',1); %unique(tmp(1,:))';
            Item.LID1(:,iItem) = num2cell(unique(tmp(2,:))',1);
            Item.SID1(:,iItem) =num2cell(unique(tmp(3,:))',1);
        end
        
        nStrip = size(Strip.LW,2);
        for iStrip=1:nStrip
            tmp = LU.DOC([1,2,3], LU.DOC(8,:) == iStrip);
            Strip.PID1(:,iStrip) = num2cell(unique(tmp(1,:))',1);
            Strip.LID1(:,iStrip) = num2cell(unique(tmp(2,:))',1);
            Strip.SID1(:,iStrip) = num2cell(unique(tmp(3,:))',1);
        end
        
        nBin = size(Bin.LW,2);
        for iBin=1:nBin
            tmp = LU.DOC([1,2,3], LU.DOC(10,:) == iBin);
            Bin.PID1(:,iBin) = num2cell(unique(tmp(1,:))',1);
            Bin.LID1(:,iBin) = num2cell(unique(tmp(2,:))',1);
            Bin.SID1(:,iBin) = num2cell(unique(tmp(3,:))',1);
        end
        
        LU.DOC
        LU.DOC([1,2,3,8],:)
    end

%         isV(d.Bin,d.Strip,d.Item,d.LU,d.Veh);
        function isV(Bin,Strip,Item,LU,Veh) %���һ��Bin�ڵ�Item�Ƿ���Ա���С��Bin����
            nBin = size(Bin.LW,2);
            maxH = max(Item.LWH(3,Item.Item_Bin(1,:)==nBin)); %Item��λ�����һ��Bin����Ʒ�߶ȵ����ֵ
            flag1 = Veh.LWH(3,:) >= maxH
            Veh.area = Veh.volume./Veh.LWH(3,:);
            flag2 =Veh.area >= Bin.Itemarea(1,nBin)
            flag = flag1 & flag2;
            thisVeh = max(find(flag==1))
            
            %�ҳ�nBin��Ӧ��LU���
            tLU = find(LU.LU_Bin(1,:) == nBin);
            
            if isSameCol(LU)
                tLU = structfun(@(x) x( : , LU.LU_Bin(1,:) == nBin ), LU, 'UniformOutput', false);
                printstruct(tLU);
            end
            if isSameCol(Veh)
                tVeh = structfun(@(x) x( : , thisVeh ), Veh, 'UniformOutput', false);
                printstruct(tVeh);
            end
            
            [tLU,tItem,~] = HLUtoItem(tLU,tVeh); 
            [tItem,tStrip] = HItemToStrip(tLU,tItem,tVeh,p);
            printstruct(tLU);
            printstruct(tItem);
            
%             d.LU,d.Item,d.ItemID] = HLUtoItem(tLU,d.Veh); 
%             tItem = find(Item.Item_Bin(1,:) == nBin);
            1
            %�ҳ�nBin��Ӧ��Item���
            
            %�ҳ�nBin��Ӧ��Strip���
            1
        end
    
%         isR(d.Bin,d.Strip,d.Item,d.LU,d.Veh);
        function isR(Bin,Strip,Item,LU,Veh)
            
        end
        
         printOut(d.Bin,d.Strip,d.Item,d.LU,d.Veh); %����,��ʱע��
        %% �����һ�����滻Ϊ��С�ĳ�����С�����滻��
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

        
        
%         printstruct(d.Veh);
%         isReAssignVeh(d.Veh,d.Bin)
%         isRegetItem()
%         printstruct(d.Bin);


end

    function [ lb ] = computerLB(Item,Veh)
        sum1 = sum(prod(Item.LWH,1));        
        % todo �����ж��Ƿ����е�BinArray�����е�bin����ͬ�� ����� �����ִ��
        sum2 = prod(Veh.LWH(:,1));
        lb = ceil(sum1/sum2);
        if lb <=0, error('EEE');end
    end
