function [d] = RunAlgorithm(d,p)
        
        % 检验Input输入数据
%         d.LU.isRota = [1 0 0]
        printstruct(d.LU);
        d = GcheckInput(d); %可以不做 
        pgon = getPolyshape(d.LU.LWH);    maxX = sum(d.LU.LWH(1,:))+10;    maxY = max(max(d.LU.LWH'))+10;  maxV = max(maxX,maxY);
%                   plot(pgon);        axis equal;   axis ([0 maxX 0 maxY]);   
%                   plot3Dshape(d.LU.LWH);

        % 数据预处理：重点：获取LU.Rotaed,托盘是否排序
        [d.LU, d.Veh] = Gpreproc(d.LU, d.Veh,p.whichSortItemOrder); %不可以不做 
        
        %% 启发式: LU到Item的算法    
        [d.LU,d.Item,d.ItemID] = HLUtoItem(d.LU,d.Veh); %Item将按ID序号排序（但下一操作将变化顺序）
        printstruct(d.LU);
        printstruct(d.Item);
        pgon = getPolyshape(d.Item.LWH);
%          figure; plot(pgon);  axis equal;  axis ([0 maxX 0 maxY]);
        %% 计算下届
%         lb = computerLB(d.Item,d.Veh);   fprintf('LB = %d \n', lb); %以某个bin类型为准
        %% 启发式：Item到Strip的算法
%         printstruct(d);
%         printstruct(d.Item);
        [d.LU,d.Item,d.Strip] = HItemToStrip(d.LU,d.Item,d.Veh,p);
        
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
%         d = modifyStripWithOneItem(d);
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
        [d.Strip,d.Bin]= HStripToBin(d.Strip,d.Veh,d.LU,p);
printstruct(d.Strip)
        %% Item到bin的信息获取:
%         printstruct(d);
%         [d] = HItemToBin(d);
        [d.LU,d.Item] = HItemToBin(d.LU,d.Item,d.Strip);
        printstruct(d.Item);
         %% 计算bin装载率
         % ItemloadingrateLimit - 每个bin内Item的体积和/每个bin去除剩余宽高后的总体积
         % Itemloadingrate - 每个bin内Item的体积和/每个bin可用总体积
         d = computeLoadingRate2DBin(d);
        function d = computeLoadingRate2DBin(d)
            % 初始化
            nBin = size(d.Bin.LW,2);
            d.Bin.Binarea = zeros(1,nBin);
            d.Bin.Itemarea = zeros(1,nBin);
            d.Bin.Itemloadingrate = zeros(1,nBin);
            d.Bin.ItemloadingrateLimit = zeros(1,nBin);
            % 计算每个Bin的装载率            
            BinWidth = d.Veh.LWH(1,1);
            BinHeight = d.Veh.LWH(2,1);
            BinArea = BinWidth .* BinHeight;
            %每个Bin的可用体积 = 车辆高度*车辆宽度
            d.Bin.Binarea = repmat(BinArea,1,nBin);            
            %每个Bin 的有限可用体积 = 宽度(bin使用宽度=车辆宽度-bin剩余宽度) *高度(bin使用高度=车辆高度-bin剩余高度)
            d.Bin.BinareaLimit = (BinWidth - d.Bin.LW(1,:)) .* (BinHeight - d.Bin.LW(2,:));
            
            a = d.Item.LWH;
            b = d.Item.Item_Bin;
            for iBin =1:nBin
                %每个Bin的装载体积
                d.Bin.Itemarea(iBin)= sum(a(1, (b(1,:)==iBin)) .* a(2, (b(1,:)==iBin)));
            end
            %每个bin的装载比率
            d.Bin.loadingrate =  d.Bin.Itemarea ./ d.Bin.Binarea;
            %每个bin的有限装载比率
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
        function isV(Bin,Strip,Item,LU,Veh) %最后一个Bin内的Item是否可以被最小的Bin放下
            nBin = size(Bin.LW,2);
            maxH = max(Item.LWH(3,Item.Item_Bin(1,:)==nBin)); %Item的位于最后一个Bin的物品高度的最大值
            flag1 = Veh.LWH(3,:) >= maxH
            Veh.area = Veh.volume./Veh.LWH(3,:);
            flag2 =Veh.area >= Bin.Itemarea(1,nBin)
            flag = flag1 & flag2;
            thisVeh = max(find(flag==1))
            
            %找出nBin对应的LU序号
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
            %找出nBin对应的Item序号
            
            %找出nBin对应的Strip序号
            1
        end
    
%         isR(d.Bin,d.Strip,d.Item,d.LU,d.Veh);
        function isR(Bin,Strip,Item,LU,Veh)
            
        end
        
         printOut(d.Bin,d.Strip,d.Item,d.LU,d.Veh); %可用,暂时注释
        %% 将最后一个车替换为更小的车（从小到大替换）
        function printOut(Bin,Strip,Item,LU,Veh)
            nBin = size(Bin.LW,2);
            for iBin = 1: nBin
                [~,idx] = find(Item.Item_Bin(1,:)==iBin); %本iBin下的item索引号
                idxSeq = Item.Item_Bin(2,idx); %本iBin内item放入顺序Seq
                fprintf('bin 的宽+长+高为: ' );
                fprintf(' %d  ',Veh.LWH);
                fprintf('\n');
                fprintf('bin %d 的剩余宽+剩余长为:  ',iBin);fprintf('\n');
                fprintf('( %d ) ',Bin.LW(:,iBin));fprintf('\n');
                fprintf('\n');

                fprintf('bin %d 包含 original item 索引号{顺序}(长宽)[旋转标志]{坐标}为  \n  ',iBin);
                fprintf('%d ',idx);fprintf('\n');
                fprintf('{%d} ',idxSeq);fprintf('\n');
                fprintf(' (%d %d %d) ', Item.LWH(1:3,idx));fprintf('\n');
                fprintf(' [%d]     ', Item.Rotaed(:,idx));fprintf('\n');
                fprintf(' {%d %d %d} ', Item.CoordItemBin(:,idx));fprintf('\n');
                fprintf('\n');

                [~,idxLU] = find(LU.LU_Bin(1,:)==iBin); %本iBin下的item索引号
                fprintf('bin %d 包含 original LU 索引号{顺序}[item序号](长宽)[旋转标志]{坐标}为  \n  ',iBin);
                idxLUSeq = LU.LU_Bin(2,idxLU); %本iBin内item放入顺序Seq
                idxLUItem = LU.LU_Item(1,idxLU);
                fprintf('%d ',idxLU);fprintf('\n');
                fprintf('{%d} ',idxLUSeq);fprintf('\n');
                fprintf('[%d] ',idxLUItem);fprintf('\n');
                fprintf(' (%d %d %d) ', LU.LWH(1:3,idxLU));fprintf('\n');
                fprintf(' [%d]     ', LU.Rotaed(:,idxLU));fprintf('\n');
                fprintf(' {%d %d %d} ', LU.CoordLUBin(:,idxLU));fprintf('\n');
                fprintf('\n');

                % 按安放顺序展示
                % %     [~,x]=sort(LU.LU_Bin(2,idxLU));
                % %     idxLUSeq = idxLUSeq(x); %本iBin内item放入顺序Seq
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
        % todo 增加判断是否所有的BinArray中所有的bin是相同的 如果是 则继续执行
        sum2 = prod(Veh.LWH(:,1));
        lb = ceil(sum1/sum2);
        if lb <=0, error('EEE');end
    end
