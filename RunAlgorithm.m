function [d] = RunAlgorithm(d,p)
        
        %% ����Input��������
        printstruct(d);
        d = GcheckInput(d);
        pgon = getPolyshape(d.LU.LWH);
                maxX = sum(d.LU.LWH(1,:))+10;
                maxY = max(max(d.LU.LWH'))+10;
                maxV = max(maxX,maxY);
%                 plot(pgon);        axis equal;   axis ([0 maxX 0 maxY]);
%         plot3Dshape(d.LU.LWH);

 

        %% ����Ԥ����
        d = Gpreproc(d);
        %% ����ʽ: LU��Item���㷨    
         printstruct(d);
        [d.LU,d.Item,d.ItemID] = HLUtoItem(d.LU,d.Veh); %Item����ID������򣨵���һ�������仯˳��
        printstruct(d);
                pgon = getPolyshape(d.Item.LWH);
%                 figure; plot(pgon);  axis equal;  axis ([0 maxX 0 maxY]);
        %% �����½�
        lb = computerLB(d.Item,d.Veh);   fprintf('LB = %d \n', lb); %��ĳ��bin����Ϊ׼
        %% ����ʽ��Item��Strip���㷨
%         printstruct(d);
%         printstruct(d.Item);
        [d] = HItemToStrip(d,p);
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
         d = modifyStripWithOneItem(d);
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
        [d.Strip,d.Bin]= HStripToBin(d.Strip,d.Veh,p);
        %% Item��bin����Ϣ��ȡ:
%         printstruct(d);
        [d] = HItemToBin(d);
         %% ����binװ����
         % ItemloadingrateLimit - ÿ��bin��Item�������/ÿ��binȥ��ʣ���ߺ�������
         % Itemloadingrate - ÿ��bin��Item�������/ÿ��bin���������
         d = computeLoadingRateBin(d);
        function d = computeLoadingRateBin(d)
            % ��ʼ��
            nBin = size(d.Bin.LW,2);
            d.Bin.Binvolume = zeros(1,nBin);
            d.Bin.Itemvolume = zeros(1,nBin);
            d.Bin.Itemloadingrate = zeros(1,nBin);
            d.Bin.ItemloadingrateLimit = zeros(1,nBin);
            % ����ÿ��Bin��װ����            
            BinWidth = d.Veh.LWH(1,:);
            BinHeight = d.Veh.LWH(2,:);
            BinVolume = BinWidth .* BinHeight;
            %ÿ��Bin�Ŀ������ = �����߶�*�������
            d.Bin.Binvolume = repmat(BinVolume,1,nBin);            
            %ÿ��Bin �����޿������ = ���(binʹ�ÿ��=�������-binʣ����) *�߶�(binʹ�ø߶�=�����߶�-binʣ��߶�)
            d.Bin.BinvolumeLimit = (BinWidth - d.Bin.LW(1,:)) .* (BinHeight - d.Bin.LW(2,:));
            
            a = d.Item.LWH;
            b = d.Item.Item_Bin;
            for iBin =1:nBin
                %ÿ��Bin��װ�����
                d.Bin.Itemvolume(iBin)= sum(a(1, (b(1,:)==iBin)) .* a(2, (b(1,:)==iBin)));
            end
            %ÿ��bin��װ�ر���
            d.Bin.loadingrate =  d.Bin.Itemvolume ./ d.Bin.Binvolume;
            %ÿ��bin������װ�ر���
            d.Bin.loadingrateLimit =  d.Bin.Itemvolume ./ d.Bin.BinvolumeLimit;
        end

end

    function [ lb ] = computerLB(Item,Veh)
        sum1 = sum(prod(Item.LWH,1));        
        % todo �����ж��Ƿ����е�BinArray�����е�bin����ͬ�� ����� �����ִ��
        sum2 = prod(Veh.LWH(:,1));
        lb = ceil(sum1/sum2);
        if lb <=0, error('EEE');end
    end
