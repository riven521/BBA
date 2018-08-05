function [d] = RunAlgorithm(d,p)
        
        %% 检验Input输入数据
        printstruct(d);
        d = GcheckInput(d);
        pgon = getPolyshape(d.LU.LWH);
                maxX = sum(d.LU.LWH(1,:))+10;
                maxY = max(max(d.LU.LWH'))+10;
                maxV = max(maxX,maxY);
%                 plot(pgon);        axis equal;   axis ([0 maxX 0 maxY]);
%         plot3Dshape(d.LU.LWH);

 

        %% 数据预处理
        d = Gpreproc(d);
        %% 启发式: LU到Item的算法    
         printstruct(d);
        [d.LU,d.Item,d.ItemID] = HLUtoItem(d.LU,d.Veh); %Item将按ID序号排序（但下一操作将变化顺序）
        printstruct(d);
                pgon = getPolyshape(d.Item.LWH);
%                 figure; plot(pgon);  axis equal;  axis ([0 maxX 0 maxY]);
        %% 计算下届
        lb = computerLB(d.Item,d.Veh);   fprintf('LB = %d \n', lb); %以某个bin类型为准
        %% 启发式：Item到Strip的算法
%         printstruct(d);
%         printstruct(d.Item);
        [d] = HItemToStrip(d,p);
        %% 计算strip装载率
%         printstruct(d);
        d = computeLoadingRateStrip(d);
        function d = computeLoadingRateStrip(d)
            % 初始化
            nStrip = size(d.Strip.LW,2);
            d.Strip.Stripvolume = zeros(1,nStrip);
            d.Strip.StripvolumeLimit = zeros(1,nStrip);
            d.Strip.Itemvolume = zeros(1,nStrip);
            d.Strip.loadingrate = zeros(1,nStrip);
            d.Strip.loadingrateLimit = zeros(1,nStrip);
            
            % 计算每个strip的装载率
            %每个strip的可用体积 = 高度*宽度(车辆的宽度)
            d.Strip.Stripvolume = d.Strip.LW(2,:)*d.Veh.LWH(1,1);
            %每个strip的有限可用体积 = 高度*宽度(strip使用宽度=车辆宽度-strip剩余宽度)
            d.Strip.StripvolumeLimit = d.Strip.LW(2,:) .* (d.Veh.LWH(1,1) - d.Strip.LW(1,:));
            a = d.Item.LWH;
            b = d.Item.Item_Strip;
            for iStrip =1:nStrip
                %每个strip的装载体积
                d.Strip.Itemvolume(iStrip)= sum(a(1, (b(1,:)==iStrip)) .* a(2, (b(1,:)==iStrip)));
            end
            %每个strip的装载比率
            d.Strip.loadingrate =  d.Strip.Itemvolume ./ d.Strip.Stripvolume;
            %每个strip的有限装载比率
            d.Strip.loadingrateLimit =  d.Strip.Itemvolume ./ d.Strip.StripvolumeLimit;
        end
        %% 对Strip中仅有一个且高>宽的Item进行选择并更新相应数据
         d = modifyStripWithOneItem(d);
        function d = modifyStripWithOneItem(d)
            stripheight = d.Strip.LW(2,:);
            binwidth = d.Veh.LWH(1,1);
            stripleftwidth = d.Strip.LW(1,:);
            stripwidth = ( binwidth - stripleftwidth );
            [tmpset] = find(stripheight > stripwidth);
            if ~isempty(tmpset)
                if isscalar(tmpset) %对该strip调换内部仅有1个Item方可,多个调整涉及CoordItemStrip
                    d.Strip.LW(:,tmpset) = [binwidth-stripheight(tmpset),stripwidth(tmpset)];    %strip的长宽调整
                    %内部Item的itemRotaFlag调整 
                    idxItem = find(d.Item.Item_Strip(1,:)==tmpset );
                    if isscalar(idxItem)
                        d.Item.itemRotaFlag(idxItem) = ~d.Item.Rotaed(idxItem);
                    end                    
                    %内部LU的LURotaFlag 不調 還未到 %内部Item的CoordItemStrip不調                    
                end
            end
        end
        %% 启发式：Strip到Bin的算法
%         printstruct(d);
        [d.Strip,d.Bin]= HStripToBin(d.Strip,d.Veh,p);
        %% Item到bin的信息获取:
%         printstruct(d);
        [d] = HItemToBin(d);
         %% 计算bin装载率
         % ItemloadingrateLimit - 每个bin内Item的体积和/每个bin去除剩余宽高后的总体积
         % Itemloadingrate - 每个bin内Item的体积和/每个bin可用总体积
         d = computeLoadingRateBin(d);
        function d = computeLoadingRateBin(d)
            % 初始化
            nBin = size(d.Bin.LW,2);
            d.Bin.Binvolume = zeros(1,nBin);
            d.Bin.Itemvolume = zeros(1,nBin);
            d.Bin.Itemloadingrate = zeros(1,nBin);
            d.Bin.ItemloadingrateLimit = zeros(1,nBin);
            % 计算每个Bin的装载率            
            BinWidth = d.Veh.LWH(1,:);
            BinHeight = d.Veh.LWH(2,:);
            BinVolume = BinWidth .* BinHeight;
            %每个Bin的可用体积 = 车辆高度*车辆宽度
            d.Bin.Binvolume = repmat(BinVolume,1,nBin);            
            %每个Bin 的有限可用体积 = 宽度(bin使用宽度=车辆宽度-bin剩余宽度) *高度(bin使用高度=车辆高度-bin剩余高度)
            d.Bin.BinvolumeLimit = (BinWidth - d.Bin.LW(1,:)) .* (BinHeight - d.Bin.LW(2,:));
            
            a = d.Item.LWH;
            b = d.Item.Item_Bin;
            for iBin =1:nBin
                %每个Bin的装载体积
                d.Bin.Itemvolume(iBin)= sum(a(1, (b(1,:)==iBin)) .* a(2, (b(1,:)==iBin)));
            end
            %每个bin的装载比率
            d.Bin.loadingrate =  d.Bin.Itemvolume ./ d.Bin.Binvolume;
            %每个bin的有限装载比率
            d.Bin.loadingrateLimit =  d.Bin.Itemvolume ./ d.Bin.BinvolumeLimit;
        end

end

    function [ lb ] = computerLB(Item,Veh)
        sum1 = sum(prod(Item.LWH,1));        
        % todo 增加判断是否所有的BinArray中所有的bin是相同的 如果是 则继续执行
        sum2 = prod(Veh.LWH(:,1));
        lb = ceil(sum1/sum2);
        if lb <=0, error('EEE');end
    end
