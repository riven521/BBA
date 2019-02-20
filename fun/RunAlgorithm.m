function [d] = RunAlgorithm(d,p)
        global  ISshuaiwei ISreStripToBin  %重新执行HStripToBin方案
        global  ISplotshuaiwei  ISplotStripToBinAgain   ISplotRunAlgo  ISplotRunLIS

        %% 预处理:检验Input输入数据
        fprintf(1,'   Running RunAlgorithm  ...\n');
        
        %         d = GcheckInput(d);    %可以不做   因为main做一下即可
                
        % 数据预处理：重点：获取LU.Rotaed,托盘是否排序 类似cpuLU函数 增加了许多LU的属性  （进入之前LWH已经时包含margin的）        
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        [d.LU, d.Veh] = cpuVehLU(d.LU, d.Veh);            %必须做>??  in case 意外错误

        %% 启发式: cpuLU 
        % 计算LU的isNonMixed和LU的isMixedTile 属性 用途似乎不大, 目的为了生成Item时,满垛的一起,不满垛的后面
         
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % ADD: isNonmixed/ismixedtile
        %         [d.LU] = cpuLU(d.LU,d.Veh);     % REMOVE 因为isNonmixed/ismixedtile 不再必要
        
        %% 启发式: HLUtoItem LU到Item的算法 
        %fprintf(1,'     Running HLUtoItem...\n');
        
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % isNonmixed/ismixedtile
        % ADD: order/LU_Item
        % Item:
        % ADD: isRota/Rotaed/HLayer/LWH/Weight
        [d.LU,d.Item] = HLUtoItem(d.LU,d.Veh);          %Item将按ID序号排序（但下一操作将变化顺序）
        
        %% 启发式: repairItems 修复LU_Item内是否上轻下重(修改顺序而已,坐标此时还未计算)
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % isNonmixed/ismixedtile
        % order/LU_Item
        % REPAIR: LU_Item
        % Item:
        % isRota/Rotaed/HLayer/LWH/Weight
        d.LU.LU_Item(2,:)= repairItems(d.LU);
        
        % 核查是否还有上轻下重的情况
        chktLU(d.LU)
    
        %% 启发式: cpuItem - 计算Item的 MixOrder isHeightFull nbItem 三个参数
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % isNonmixed/ismixedtile
        % order/LU_Item
        % Item:
        % isRota/Rotaed/HLayer/LWH/Weight
        % ADD: /isHeightFull/MixOrder/nbItem/IDs
        [d.Item] = cpuItem(d.Item,d.LU,d.Veh);        % printstruct(d,'sortfields',1,'PRINTCONTENTS',0);
        %        printstruct(d.LU,'sortfields',1,'PRINTCONTENTS',1);

        %% 启发式: HItemToStrip ： Item到Strip的算法        
        % 1 Item排序: %1: SID ; 2: MixOrder; 一般正真开始: 
        %                       3: Longth/Height; 4:Width; 5: LID; (3,4,5,多数一样) 6: Height
        % 2 按Item给定顺序插入新或旧的Strip
        
        %fprintf(1,'     Running HItemToStrip...\n');     
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % isNonmixed/ismixedtile
        % order/LU_Item
        % Item:
        % isRota/Rotaed/HLayer/LWH/Weight
        % isHeightFull/MixOrder/nbItem/IDs
        % ADD: itemorder/Item_Strip/CoordItemStrip
        % Strip:
        % ADD: LW/Weight
        [d.Item,d.Strip] = HItemToStrip(d.LU,d.Item,d.Veh,p);     %   printstruct(d);   %  printstruct(d.Item);
      
        %% 启发式: cpuLU2Strip ： 计算LU_Strip, CoordLUStrip  
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % isNonmixed/ismixedtile
        % order/LU_Item
        % ADD: LU_Strip, CoordLUStrip
        % Item:
        % isRota/Rotaed/HLayer/LWH/Weight
        % isHeightFull/MixOrder/nbItem/IDs
        % itemorder/Item_Strip/CoordItemStrip
        % Strip:
        % LW/Weight        
        [d.LU] = cpuLU2Strip(d.LU,d.Item);     %   printstruct(d);   %  printstruct(d.Item);
        
        %% 启发式: cpuStrip 大量关键strip属性计算         
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % isNonmixed/ismixedtile
        % order/LU_Item
        % LU_Strip, CoordLUStrip
        % Item:
        % isRota/Rotaed/HLayer/LWH/Weight
        % isHeightFull/MixOrder/nbItem/IDs
        % itemorder/Item_Strip/CoordItemStrip
        % Strip:
        % LW/Weight
        % ADD: isMixedIDs/nbIDs/isHeightFull/isHeightBalance/isWidthFull/nbItem/nbLU/nbLULID
        % ADD: isAllPured/isSingleItem/mlmdHeight/IDs/~Limit
        % ADD: nLUID / nLULID
         [d.Strip] = cpuStrip(d.Strip,d.Item,d.LU,d.Veh); %          [d.Strip,d.LU] = cpuStrip(d.Strip,d.Item,d.LU,d.Veh);

        %% 启发式：HStripToBin - Strip到Bin的算法
        % ********* 1 Strip排序: % 1: SID 2: priorityofLID
        %fprintf(1,'     Running HStripToBin...\n');
        
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % isNonmixed/ismixedtile
        % order/LU_Item
        % LU_Strip, CoordLUStrip
        % Item:
        % isRota/Rotaed/HLayer/LWH/Weight
        % isHeightFull/MixOrder/nbItem/IDs
        % itemorder/Item_Strip/CoordItemStrip
        % Strip:
        % LW/Weight
        % isMixedIDs/nbIDs/isHeightFull/isHeightBalance/isWidthFull/nbItem/nbLU/nbLULID
        % isAllPured/isSingleItem/mlmdHeight/IDs/~Limit
        % ADD: striporder/strip_bin
        % Bin:
        % ADD: LW/Weight
        [d.Strip,d.Bin] = HStripToBin(d.Strip,d.Veh,p);
        
         %% 启发式：HItemToBin - % 计算LU在Bin内坐标and顺序等数据计算 
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % isNonmixed/ismixedtile
        % order/LU_Item
        % LU_Strip, LU.CoordLUStrip
        % ADD: LU_Bin/CoordLUBin/ lbrt_A / nIDBin / nLIDBin 
        % Item:
        % isRota/Rotaed/HLayer/LWH/Weight
        % isHeightFull/MixOrder/nbItem/IDs
        % itemorder/Item_Strip/CoordItemStrip
        % ADD: Item_bin/CoordItemBin
        % Strip:
        % LW/Weight
        % isMixedIDs/nbIDs/isHeightFull/isHeightBalance/isWidthFull/nbItem/nbLU/nbLULID
        % isAllPured/isSingleItem/mlmdHeight/IDs/~Limit
        % striporder/strip_bin
        % ADD: nLUID / nLULID（暂时无用)  / nLUIDBin / nLULIDBin  
        % Bin:
        % LW/Weight
        [d.LU,d.Item,d.Strip] = HItemToBin(d.LU,d.Item,d.Strip); % 计算LU在Bin内坐标and顺序   %  Item.Item_Bin  Item.CoordItemBin LU.LU_Bin LU.CoordLUBin

        %% 启发式：cpuBin - 计算isTileNed 是否需要甩尾平铺等bin指标
        % LU:
        % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
        % Roated/maxHLayer/nbID/nbLID
        % isNonmixed/ismixedtile
        % order/LU_Item
        % LU_Strip, LU.CoordLUStrip
        % LU_Bin/CoordLUBin/ lbrt_A
        % Item:
        % isRota/Rotaed/HLayer/LWH/Weight
        % isHeightFull/MixOrder/nbItem/IDs
        % itemorder/Item_Strip/CoordItemStrip
        % Item_bin/CoordItemBin
        % Strip:
        % LW/Weight
        % isMixedIDs/nbIDs/isHeightFull/isHeightBalance/isWidthFull/nbItem/nbLU/nbLULID
        % isAllPured/isSingleItem/mlmdHeight/IDs/~Limit
        % striporder/strip_bin
        % Bin:
        % LW/Weight        
        % ADD: isTileNed/IDs
        [d.Bin] = cpuBin(d.Bin,d.Strip,d.Item,d.LU,d.Veh);  %计算Bin内isTileNeed （高度，宽度不满）
        
        
        if ISplotRunAlgo && ISplotRunLIS
            %             plotSolutionT(d.LU,d.Veh,1,0,0,0,3,'原顺序LU'); %LU当前顺序
                        plotSolutionT(d.LU,d.Veh,3,0,0,0,3,'排序后LU');     %LU排序后，因为Lu.order此时已知
                        plotSolutionT(d.LU,d.Veh,0,1,0,0,3,'排序前Item',[],[],[]);   %Item排序前，因为Item->Strip前的排序LU_Item内没有
                        plotSolutionT(d.LU,d.Veh,0,2,0,0,3,'排序后Item',[],d.Item.itemorder,[]);
                        plotSolutionT(d.LU,d.Veh,0,0,1,0,3,'排序前Strip');     %Strip排序前，因为Strip->Bin前的排序LU_Strip内没有
            plotSolutionT(d.LU,d.Veh,0,0,2,0,1,'排序后Strip',[],[],d.Strip.striporder);     %Strip排序后
            plotSolutionT(d.LU,d.Veh,0,0,0,1,1,'ID作图Bin');                  %     Bin排序后
            %             plotSolutionT(d.LU,d.Veh,0,0,0,1,1,'LID作图Bin');            % Bin排序后 
        end
                                                    
        %% 启发式：量大车头 方案2: 每个剩余strip全体内比较量 better than 方案1 有故障? % 特别是对第三辆车及以后
        if ISreStripToBin==1
                fprintf(1,'     Running HStripToBinAgain...\n');                
                % LU:
                % IDs/OIDs/LWH/Weight/Index/margin/isRota/maxL
                % Roated/maxHLayer/nbID/nbLID
                % isNonmixed/ismixedtile
                % order/LU_Item
                % LU_Strip, LU.CoordLUStrip
                % LU_Bin/CoordLUBin/ lbrt_A
                % Item:
                % isRota/Rotaed/HLayer/LWH/Weight
                % isHeightFull/MixOrder/nbItem/IDs
                % itemorder/Item_Strip/CoordItemStrip
                % Item_bin/CoordItemBin
                % Strip:
                % LW/Weight
                % isMixedIDs/nbIDs/isHeightFull/isHeightBalance/isWidthFull/nbItem/nbLU/nbLULID
                % isAllPured/isSingleItem/mlmdHeight/IDs/~Limit
                % striporder/strip_bin
                % Bin:
                % LW/Weight
                % isTileNed/IDs
                [d.Strip,d.Bin,TFHStripToBinAgain] = HStripToBinAgain(d.Bin,d.Strip,d.Item,d.LU,d.Veh,p); 
                
                % sum(~d.Strip.isHeightBalance & d.Strip.isHeightFull)
                
                [d.LU,d.Item] = HItemToBin(d.LU,d.Item,d.Strip); % 计算LU在Bin内坐标and顺序   %  Item.Item_Bin  Item.CoordItemBin LU.LU_Bin LU.CoordLUBin

                [d.Bin] = cpuBin(d.Bin,d.Strip,d.Item,d.LU,d.Veh);  %计算Bin内相关属性 % 计算isTileNeed

                if TFHStripToBinAgain &&  ISplotStripToBinAgain && ISplotRunAlgo
                      plotSolutionT(d.LU,d.Veh,0,0,0,1,1,'量大车头后Bin'); % Bin排序后
                end
        end
                                                        
        %% 启发式：增加Strip的甩尾优化 *********** 修改 Strip.Strip_Bin 第二行
        if ISshuaiwei==1   
                fprintf(1,'     Running HStripSW...\n');

                % 返回strip的strip_bin顺序, strip是否甩尾,strip甩尾顺序 
                [d.Strip.Strip_Bin(2,:),  d.Strip.isShuaiWei, d.LU.isShuaiWei, TFHStripSW] = HStripSW(d.Strip,d.LU);

                % 计算LU/Item 在Bin内坐标and顺序LU_Bin  NOTE: Lu_Strip and Item_Strip 不受甩尾影响
                [d.LU,d.Item] = HItemToBin(d.LU,d.Item,d.Strip);    %  Item.Item_Bin  Item.CoordItemBin LU.LU_Bin LU.CoordLUBin
       
                [d.Bin] = cpuBin(d.Bin,d.Strip,d.Item,d.LU,d.Veh);  %计算Bin的isTileNeed

                if TFHStripSW && ISplotshuaiwei && ISplotRunAlgo
                        plotSolutionT(d.LU,d.Veh,0,0,0,1,1,'甩尾后Bin'); % Bin排序后 
%                     plotSolutionT(d.LU,d.Veh,0,0,0,1,8,'甩尾后Bin'); % Bin排序后 
                end                
        end

        %% 启发式：最后一部：赋值车型
       d.LU.LU_VehType = ones(size(d.LU.ID)) * d.Veh.order(1); % 针对车型选择,增加变量LU_VehType : 由于Veh内部按体积递减排序,获取order的第一个作为最大值
       fprintf(1,'   Running RunAlgorithm  DONE ...\n');                

end

%         printstruct(d.Item,'sortfields',1,'PRINTCONTENTS',0)
%         printstruct(d.Bin,'sortfields',1,'PRINTCONTENTS',1)
%         printstruct(d,'sortfields',1,'PRINTCONTENTS',0)

% cpuLU2Strip 计算LU.LU_Strip, LU.CoordLUStrip
function [LU] = cpuLU2Strip(LU,Item)

nbLU = size(LU.LWH,2);
LU.LU_Strip = zeros(2,nbLU);
LU.CoordLUStrip = zeros(3,nbLU);

% 更新LU_Strip
for iLU=1:nbLU
    % 更新LU_Strip第一行
    iItem = LU.LU_Item(1,iLU);   %iLU属于第几个Item, Item属于第几个Strip,则Lu属于第几个Strip
    LU.LU_Strip(1,iLU)= Item.Item_Strip(1,iItem);
    % 更新LU_Strip第二行
    fiItem = find(Item.Item_Strip(1,:) == Item.Item_Strip(1,iItem) & Item.Item_Strip(2,:) < Item.Item_Strip(2,iItem));
    nbLUfiItem = sum(ismember(LU.LU_Item(1,:),fiItem));
    LU.LU_Strip(2,iLU) = nbLUfiItem+LU.LU_Item(2,iLU); % 进入Strip顺序: 同一Strip内先前进入个数nbLUfiItem + 本iLU在Item的顺序
    % 更新LU.CoordLUStrip
    LU.CoordLUStrip(1,iLU) = Item.CoordItemStrip(1,iItem);
    LU.CoordLUStrip(2,iLU) = Item.CoordItemStrip(2,iItem);
    % fLU: 与iLU同属iItem 且 顺序晚于本iLU; 可能为空, 不影响.
    fLU = LU.LU_Item(1,:) == iItem & LU.LU_Item(2,:) < LU.LU_Item(2,iLU);
    LU.CoordLUStrip(3,iLU) = sum(LU.LWH(3,fLU));
end

end



%% 启发式：特殊：堆垛均衡设置 在cpuStrip之后
        % ********* 1 5555 堆垛均衡设置 *********** 修改 d.LU.maxHLayer(luidxPP) 
%          if ISstripbalance==0
%              [d.Strip,d.Item,d.LU,TFStripBalance] = HStripBalance(d.Strip,d.Item,d.LU,d.Veh,p);  end
%          if ISstripbalance && TFStripBalance && ISplotstripbalance && ISplotRunAlgo
%              plotSolutionT(d.LU,d.Veh,0,0,1,0,3,'均衡后Strip'); %Strip:HStripBalance后
%          end     % V1:  if ISplotStrip==1,      figure(222);      plot3DStrip(d.LU,d.Item,d.Veh,'LU');        end    %  plot3DStrip(d.LU,d.Item,d.Veh,'Item');  % 基于LU.CoordLUBin %   figure(222);     plot3DStrip(d.LU,d.Item,d.Veh,'LU');         %    igure(111);          plot3DStrip(d.LU,d.Item,d.Veh,'Item'); 

       %% 后处理: PID/SID等返回 改为返回顺序，此类不动的值从原始记录中获取
%        d.LU.SID = d.LU.OSID;
%        d.LU.PID = d.LU.OPID;

        %% 打印输出
%     printOut(d.Bin,d.Strip,d.Item,d.LU,d.Veh); %可用,暂时注释
% %         function printOut(Bin,Strip,Item,LU,Veh)
% %             nBin = size(Bin.LW,2);
% %             for iBin = 1: nBin
% %                 [~,ibin] = find(Item.Item_Bin(1,:)==iBin); %本iBin下的item索引号
% %                 idxSeq = Item.Item_Bin(2,ibin); %本iBin内item放入顺序Seq
% %                 fprintf('bin 的宽+长+高为: ' );
% %                 fprintf(' %d  ',Veh.LWH);
% %                 fprintf('\n');
% %                 fprintf('bin %d 的剩余宽+剩余长为:  ',iBin);fprintf('\n');
% %                 fprintf('( %d ) ',Bin.LW(:,iBin));fprintf('\n');
% %                 fprintf('\n');
% % 
% %                 fprintf('bin %d 包含 original item 索引号{顺序}(长宽)[旋转标志]{坐标}为  \n  ',iBin);
% %                 fprintf('%d ',ibin);fprintf('\n');
% %                 fprintf('{%d} ',idxSeq);fprintf('\n');
% %                 fprintf(' (%d %d %d) ', Item.LWH(1:3,ibin));fprintf('\n');
% %                 fprintf(' [%d]     ', Item.Rotaed(:,ibin));fprintf('\n');
% %                 fprintf(' {%d %d %d} ', Item.CoordItemBin(:,ibin));fprintf('\n');
% %                 fprintf('\n');
% % 
% %                 [~,idxLU] = find(LU.LU_Bin(1,:)==iBin); %本iBin下的item索引号
% %                 fprintf('bin %d 包含 original LU 索引号{顺序}[item序号](长宽)[旋转标志]{坐标}为  \n  ',iBin);
% %                 idxLUSeq = LU.LU_Bin(2,idxLU); %本iBin内item放入顺序Seq
% %                 idxLUItem = LU.LU_Item(1,idxLU);
% %                 fprintf('%d ',idxLU);fprintf('\n');
% %                 fprintf('{%d} ',idxLUSeq);fprintf('\n');
% %                 fprintf('[%d] ',idxLUItem);fprintf('\n');
% %                 fprintf(' (%d %d %d) ', LU.LWH(1:3,idxLU));fprintf('\n');
% %                 fprintf(' [%d]     ', LU.Rotaed(:,idxLU));fprintf('\n');
% %                 fprintf(' {%d %d %d} ', LU.CoordLUBin(:,idxLU));fprintf('\n');
% %                 fprintf('\n');
% % 
% %                 % 按安放顺序展示
% %                 % %     [~,x]=sort(LU.LU_Bin(2,idxLU));
% %                 % %     idxLUSeq = idxLUSeq(x); %本iBin内item放入顺序Seq
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

%% 部分重要注释
        % %         %  *******  如果有错误,就保留; 
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
                                        
                                        % 计算下届
                                        %         lb = computerLB(d.Item,d.Veh);   fprintf('LB = %d \n', lb); %以某个bin类型为准
                                        
                                                % 对Strip中仅有一个且高>宽的Item进行选择并更新相应数据
                                                %         d = modifyStripWithOneItem(d);
                                                %         function d = modifyStripWithOneItem(d)
                                                %             stripheight = d.Strip.LW(2,:);
                                                %             binwidth = d.Veh.LWH(1,1);
                                                %             stripleftwidth = d.Strip.LW(1,:);
                                                %             stripwidth = ( binwidth - stripleftwidth );
                                                %             [tmpset] = find(stripheight > stripwidth);
                                                %             if ~isempty(tmpset)
                                                %                 if isscalar(tmpset) %对该strip调换内部仅有1个Item方可,多个调整涉及CoordItemStrip
                                                %                     d.Strip.LW(:,tmpset) = [binwidth-stripheight(tmpset),stripwidth(tmpset)];    %strip的长宽调整
                                                %                     %内部Item的itemRotaFlag调整 
                                                %                     idxItem = find(d.Item.Item_Strip(1,:)==tmpset );
                                                %                     if isscalar(idxItem)
                                                %                         d.Item.itemRotaFlag(idxItem) = ~d.Item.Rotaed(idxItem);
                                                %                     end                    
                                                %                     %内部LU的LURotaFlag 不調 還未到 %内部Item的CoordItemStrip不調                    
                                                %                 end
                                                %             end
                                                %         end
                                                
                                                                % ********** 平铺优化 **********
                                                                % %         if any(d.Strip.isWidthFull)
                                                                % %             fw = d.Strip.isWidthFull == 0;
                                                                % %             fm = d.Strip.isMixed == 0;
                                                                % %             %             fi = find(fw);
                                                                % %             nbbin = numel(d.Bin.Weight);
                                                                % %             %    for i=1:sum(fi)
                                                                % %             for i=1:nbbin
                                                                % %                 f2 = d.Strip.Strip_Bin(1,:) == i;
                                                                % %                 f3 = f2 & fw & fm;   %即在bin i内, 又是isWidthFull的strip            
                                                                % % %                 f2 = d.Strip.Strip_Bin(1,:) == d.Strip.Strip_Bin(1,fi(i)); % 与本fi(i)处在同一bin内的isWidthFull的个数
                                                                % % %                 f3 = f2 & f;                
                                                                % %                 if sum(f3) > 1                    
                                                                % %                     error('与本fi(i)处在同一bin内的isWidthFull的个数>1');
                                                                % %                 end
                                                                % %                 fi3 = find(f3);
                                                                % %                 
                                                                % %                 d.Strip
                                                                % %                 d.Strip.LW
                                                                % %                 
                                                                % %                 % 需要平铺的strip对于的item和lu的logical值
                                                                % %                 fitem =d.Item.Item_Strip(1,:) == fi3;
                                                                % %                 flu =d.LU.LU_Strip(1,:) == fi3;
                                                                % %                 d.LU
                                                                % %                 d.Item
                                                                % %                 1
                                                                % %                 % 修改LU.CoordLUBin的值
                                                                % %                 % 修改LU.LU_Item变化的多个值
                                                                % %                 % 修改Item.LWH的值
                                                                % %                 % 修改Item.isFull1的值,Weight,isHeightFull,isWeightFine,Layer,PID,LID,SID
                                                                % %                 % CoordItemBin, CoordItemStrip, Item_Strip(上下调整会变化)
                                                                % %                 
                                                                % %                 %             d.Strip
                                                                % %                 d.Strip.isMixed
                                                                % %                 d.Strip.isAllPured
                                                                % %                 d.Strip.maxHeight
                                                                % %             end    
                                                                % %         1
                                                                % %         end

                                                                

%% 计算下届
% % % % %     function [ lb ] = computerLB(Item,Veh)
% % % % %         sum1 = sum(prod(Item.LWH,1));        
% % % % %         % todo 增加判断是否所有的BinArray中所有的bin是相同的 如果是 则继续执行
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
