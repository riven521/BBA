function [d] = RunAlgorithm(d,p)
        global ISshuaiwei ISplotStrip ISreStripToBin %重新执行HStripToBin方案
                
        %% 预处理:检验Input输入数据
        d = GcheckInput(d);    %可以不做   printstruct(d,'sortfields',1,'PRINTCONTENTS',0);

        % 数据预处理：重点：获取LU.Rotaed,托盘是否排序
        [d.LU, d.Veh] = Gpreproc(d.LU, d.Veh,p.whichSortItemOrder); %必须做 含cpuLU

        %% 启发式: LU到Item的算法    
        [d.LU,d.Item] = HLUtoItem(d.LU,d.Veh);   %Item将按ID序号排序（但下一操作将变化顺序）

        % Item.isNonMixed Item.isMixedTile isHeightFull
        [d.Item,d.LU] = cpuItem(d.Item,d.LU,d.Veh);        % printstruct(d,'sortfields',1,'PRINTCONTENTS',0);

        %% 启发式: Item到Strip的算法        
        %  *******  *******  *******
        % 1 Item排序: %1: SID ; 2: isNonMixed; 一般正真开始: 
        %                       3: Longth/Height; 4:Width; 5: LID; (3,4,5,多数一样) 6: Height
        % 2 按Item给定顺序插入新或旧的Strip
        [d.LU,d.Item,d.Strip] = HItemToStrip(d.LU,d.Item,d.Veh,p);     %   printstruct(d);   %  printstruct(d.Item);

        % 计算LU.LU_Strip, LU.CoordLUStrip
        [d.Strip,d.LU] = cpuStrip(d.Strip,d.Item,d.LU,d.Veh);

        if ISplotStrip==1,      figure(222);     plot3DStrip(d.LU,d.Item,d.Veh,'LU');        end    % igure(111);   plot3DStrip(d.LU,d.Item,d.Veh,'Item');  % 基于LU.CoordLUBin %         figure(222);     plot3DStrip(d.LU,d.Item,d.Veh,'LU');         %    igure(111);          plot3DStrip(d.LU,d.Item,d.Veh,'Item'); 

        %% 启发式：Strip到Bin的算法        
        % ********* 1 Strip排序: % 1: SID 2: priorityofLID
        [d.Strip,d.Bin] = HStripToBin(d.Strip,d.Veh,d.LU,p);
        
        [d.LU,d.Item] = HItemToBin(d.LU,d.Item,d.Strip); % 计算LU在Bin内坐标and顺序   %  Item.Item_Bin  Item.CoordItemBin LU.LU_Bin LU.CoordLUBin
        [d.Bin,d.LU] = cpuBin(d.Bin,d.Strip,d.Item,d.LU,d.Veh);  %计算Bin内相关属性 % 计算isTileNeed
        if ISplotStrip==1,      plot3DBPP(d,p);      end    % igure(111);   plot3DStrip(d.LU,d.Item,d.Veh,'Item');  % 基于LU.CoordLUBin 
        
        % ********* 2 量大车头方案2: 每个剩余strip全体内比较量 better than 方案1
        if ISreStripToBin==1,   [d.Strip,d.Bin] = HreStripToBin(d.Bin,d.Strip,d.Item,d.LU,d.Veh,p); end     % 量大车头方案1: 每个Bin内strip比较量      % [d.Strip,d.Bin]= HreStripToEachBin(d.Bin,d.Strip,d.Item,d.LU,d.Veh,p);
        
        [d.LU,d.Item] = HItemToBin(d.LU,d.Item,d.Strip); % 计算LU在Bin内坐标and顺序   %  Item.Item_Bin  Item.CoordItemBin LU.LU_Bin LU.CoordLUBin
        [d.Bin,d.LU] = cpuBin(d.Bin,d.Strip,d.Item,d.LU,d.Veh);  %计算Bin内相关属性 % 计算isTileNeed
        if ISplotStrip==1,      plot3DBPP(d,p);      end    % igure(111);   plot3DStrip(d.LU,d.Item,d.Veh,'Item');  % 基于LU.CoordLUBin
        
        % ********* 3 增加Strip的甩尾优化 *********** 修改 Strip.Strip_Bin
        if ISshuaiwei==1,      [d.Strip] = HStripSW(d.Strip);      end
        
        [d.LU,d.Item] = HItemToBin(d.LU,d.Item,d.Strip); % 计算LU在Bin内坐标and顺序   %  Item.Item_Bin  Item.CoordItemBin LU.LU_Bin LU.CoordLUBin
        [d.Bin,d.LU] = cpuBin(d.Bin,d.Strip,d.Item,d.LU,d.Veh);  %计算Bin内相关属性 % 计算isTileNeed
        if ISplotStrip==1,      plot3DBPP(d,p);      end    % igure(111);   plot3DStrip(d.LU,d.Item,d.Veh,'Item');  % 基于LU.CoordLUBin
        
                %         printstruct(d.Item,'sortfields',1,'PRINTCONTENTS',0)  
                %         printstruct(d.Bin,'sortfields',1,'PRINTCONTENTS',1)  
                %         printstruct(d,'sortfields',1,'PRINTCONTENTS',0)  

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
 
end

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

%% 局部函数 %%

%% 函数1: 量大车头方案2: V1: s1 :全部为基于Strip排除已安排Bin后的剩余Strip;;cpuStripnbItem为部分Strip进入
function [Strip,Bin] = HreStripToBin(Bin,Strip,Item,LU,Veh,p)
        % DONE: 量大车头方案2: 
        % 目的: 量大的LU被车辆拆分为量小但仍摆放车头;
        % 方法: 每个Bin都对排除前任Strip后的剩余Strip重新排序并分别执行HStripToBin算法.
        nbBin = max(Strip.Strip_Bin(1,:));
        if nbBin>1
            ibin=2;
            while 1
                % 1 获取f and fidx: 排除首个bin后的剩余Strip逻辑值 and 索引值
                Strip.Strip_Bin(1, Strip.Strip_Bin(1,:) >= ibin ) = -1; % 所有非首bin的bin序号赋值-1
                % 找出未正确排序后的Strip
                f = Strip.Strip_Bin(1,:) == -1; 
                fidx = find(f);
                if ~any(f), break; end
                
                % 2 获取s1 : 基于f 获取的剩余Strip结构体
                s1 = rmfield(Strip,{'striporder','Strip_Bin'});
                s1 = structfun(@(x) x(:,f),s1,'UniformOutput',false);
                
                % 3 获取i1 : 基于剩余Strip获取的剩余Item结构体
                i1  = Item;
                f2 = zeros(1,length(i1.Weight));
                iIdxs= find(f);
                for j=1:length(iIdxs)
                    f2 = f2+(i1.Item_Strip(1,:) == iIdxs(j));
                end                
                i1 = structfun(@(x) x(:,logical(f2)), i1, 'UniformOutput',false);
                
                % 4 获取新s2和b2 : 基于s1和i1, 重新计算Strip.nbItem and 重新执行启发式S2B
                [s1] = cpuStripnbItem(s1,i1,LU);
                %[s1,~] = cpuStrip(s1,i1,LU,Veh);  % 此函数包含上面的子函数, 是否有必要计算全部, 还是只计算上面子函数 TODO
                
                [s2,b2]= HStripToBin(s1,Veh,LU,p);
                
                % 5 获取fs: 剩余Strip摆放后的首个bin内的 so 右侧==1
                % 替换原始Strip_Bin的值
                fs = s2.Strip_Bin(1,:)==1;
                Strip.Strip_Bin(1,fidx(fs)) = ibin;
                Strip.Strip_Bin(2,fidx(fs)) = s2.Strip_Bin(2,fs);
  
                % 6 赋值语句
                Bin.Weight(ibin) = b2.Weight(1);
                Bin.LW(:,ibin) = b2.LW(1);
                ibin = ibin+1;                
            end
        end
end

%% 函数1: 量大车头方案2: V2: s1 :HStripToBin计算时:为基于Strip排除已安排Bin后的剩余Strip;cpuStripnbItem为全部Strip进入
% % function [Strip,Bin] = HreStripToBin(Bin,Strip,Item,LU,Veh,p)
% %         % DONE: 量大车头方案2: 
% %         % 目的: 量大的LU被车辆拆分为量小但仍摆放车头;
% %         % 方法: 每个Bin都对排除前任Strip后的剩余Strip重新排序并分别执行HStripToBin算法.
% %         nbBin = max(Strip.Strip_Bin(1,:));
% %         if nbBin>1
% %             ibin=2;
% %             while 1
% %                 % 1 获取f and fidx: 排除首个bin后的剩余Strip逻辑值 and 索引值
% %                 Strip.Strip_Bin(1, Strip.Strip_Bin(1,:) == 1 ) = -1; 
% %                 Strip.Strip_Bin(1, Strip.Strip_Bin(1,:) >= ibin ) = 0; % 所有非首bin的bin序号赋值-1
% %                 % 找出未正确排序后的Strip
% %                 f = Strip.Strip_Bin(1,:) == 0; 
% %                 fidx=find(f);
% %                 if ~any(f), break; end
% %                 
% %                 % 2 获取s1 : 排除部分field,增加f的新结构体
% %                 s1 = rmfield(Strip,{'striporder','Strip_Bin'});
% %                 s1.f = f;
% %                 
% %                 % 3 获取i1 : 基于剩余Strip获取的剩余Item结构体
% %                 i1  = Item;
% %                 % 4 获取新s2和b2 : 基于s1和i1, 重新计算Strip.nbItem and 重新执行启发式S2B
% %                 [s1] = cpuStripnbItem(s1,i1,LU);
% %                 %[s1,~] = cpuStrip(s1,i1,LU,Veh);  % 此函数包含上面的子函数, 是否有必要计算全部, 还是只计算上面子函数 TODO
% %                 
% %                 s1 = structfun(@(x) x(:,f),s1,'UniformOutput',false); %进入HStripToBin只能给部分的 x(:,f)
% %                 [s2,b2]= HStripToBin(s1,Veh,LU,p);
% %                 
% %                 % 5 获取fs: 剩余Strip摆放后的首个bin内的 so 右侧==1
% %                 % 替换原始Strip_Bin的值
% %                 fs = s2.Strip_Bin(1,:)==1; %找出第i1个bin的逻辑值                
% %                Strip.Strip_Bin(1,fidx(fs)) = ibin;
% %                Strip.Strip_Bin(2,fidx(fs)) = s2.Strip_Bin(2,fs);
% % 
% %                 % 6 赋值语句
% %                 Bin.Weight(ibin) = b2.Weight(1);
% %                 Bin.LW(:,ibin) = b2.LW(1);
% %                 ibin = ibin+1;                
% %             end
% %             Strip.Strip_Bin(1, Strip.Strip_Bin(1,:) == -1 ) = 1;             
% %             % sortrows(Strip.Strip_Bin',[1,2],{'ascend','ascend'})'
% %         end
% % end

        
%% 函数2: 量大车头方案1: 
        % 目的: 量大的LU被车辆拆分为量小但仍摆放车头;
        % 方法: 每个Bin都分别但各Bin内Strip执行HStripToBin算法.
function [Strip,Bin] = HreStripToEachBin(Bin,Strip,Item,LU,Veh,p)
       % 针对后续车辆重新排序并安排到新Bin内
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
                
                % 1 重新计算Strip.nbItem and 重新执行启发式S2B
                [s1] = cpuStripnbItem(s1,i1,Veh);
                
                [s2,b2]= HStripToBin(s1,Veh,LU,p);
                
                % 2 替换回原始Strip_Bin的ibin的进入顺序
                Strip.Strip_Bin(2,f) = s2.Strip_Bin(2,:);
                
                % 3 防错语句
                if any(s2.Strip_Bin(1,:)~=1), error('重新分配后bin放不下'); end
                if size(b2.Weight,2) ~=1, error('重量不相同'); end
                if size(b2.LW,2) ~=1, error('长宽不相同'); end
                if Bin.Weight(ibin) ~= b2.Weight, error('重量不相同'); end
                if Bin.LW(:,ibin) ~= b2.LW, error('长宽不相同'); end                
                Bin.Weight(ibin) = b2.Weight;
                Bin.LW(:,ibin) = b2.LW;
                
            end
        end
end
        
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
