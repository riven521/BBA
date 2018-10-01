function [d] = RunAlgorithm(d,p)
        
        %% 预处理:检验Input输入数据
        d = GcheckInput(d); %可以不做 
        
                %         pgon = getPolyshape(d.LU.LWH);    maxX = sum(d.LU.LWH(1,:))+10;    maxY = max(max(d.LU.LWH'))+10;  maxV = max(maxX,maxY);
                %                    plot(pgon);        axis equal;   axis ([0 maxX 0 maxY]);   
                %                    plot3Dshape(d.LU.LWH);

        % 数据预处理：重点：获取LU.Rotaed,托盘是否排序
        [d.LU, d.Veh] = Gpreproc(d.LU, d.Veh,p.whichSortItemOrder); %必须做
        
        %% 启发式: LU到Item的算法    
        [d.LU,d.Item,d.ItemID] = HLUtoItem(d.LU,d.Veh); %Item将按ID序号排序（但下一操作将变化顺序）
        
        [d.Item] = cpuItem(d.Item,d.LU,d.Veh);

                        %         printstruct(d.LU);         printstruct(d.Item);
                        %         pgon = getPolyshape(d.Item.LWH);    %          figure; plot(pgon);  axis equal;  axis ([0 maxX 0 maxY]);
        %% 计算下届
%         lb = computerLB(d.Item,d.Veh);   fprintf('LB = %d \n', lb); %以某个bin类型为准
        %% 启发式Item到Strip的算法
        [d.LU,d.Item,d.Strip] = HItemToStrip(d.LU,d.Item,d.Veh,p);                %         printstruct(d);   %  printstruct(d.Item);
        
        [d.Strip] = cpuStrip(d.Strip,d.Item,d.Veh);

                %% 对Strip中仅有一个且高>宽的Item进行选择并更新相应数据
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

        %% 启发式：Strip到Bin的算法
        [d.Strip,d.Bin]= HStripToBin(d.Strip,d.Veh,d.LU,p);
        
        %% Item到bin的信息获取:
        [d.LU,d.Item] = HItemToBin(d.LU,d.Item,d.Strip);          printstruct(d.Item);

        [d.Bin] = cpuBin(d.Bin,d.Item,d.Veh);
        
     
%     printOut(d.Bin,d.Strip,d.Item,d.LU,d.Veh); %可用,暂时注释
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
