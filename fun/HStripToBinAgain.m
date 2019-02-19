%% 函数1: 量大车头方案2: V3: 
    % s1 :HStripToBin计算时:为基于Strip排除已安排Bin后的剩余Strip; cpuStripnbItem为全部Strip进入
function [Strip,Bin,TF] = HStripToBinAgain(Bin,Strip,Item,LU,Veh,p)
        % 目的: 解决量大的LU被车辆拆分为量小但仍摆放车头;
        % 方法: ****** 每个Bin都对排除前任已安排Strip后的剩余Strip, 重新排序并分别执行HStripToBin算法. ******
        nbBin = max(Strip.Strip_Bin(1,:));
        TF = false;
        if nbBin>1
            fprintf(1,'       Exsiting 量大车头 in HStripToBinAgain (nBin>1)...\n');
            ibin=2;
            while 1
                % 1 获取f and fidx: 排除首个bin后的剩余Strip逻辑值 and 索引值
                Strip.Strip_Bin(1, Strip.Strip_Bin(1,:) >= ibin ) = -1; % 所有非首bin的bin序号赋值-1, 第一个bin永远不会动
                % 找出需要重新排序计算的strip,lu,item等的flag
                fStrip = Strip.Strip_Bin(1,:) == -1; 
                fidx=find(fStrip);    if ~any(fStrip), break; end
                
                % V1 错误原因是LU.LU_Bin会发生变化, 但此处考虑不到
                %                 fItem = Item.Item_Bin(1,:) >= ibin ; %
                %                 fLu = LU.LU_Bin(1,:) >= ibin ;    

                % V2 修改引入对LU.LU_Bin变化的优化
                Item.Item_Bin(1, Item.Item_Bin(1,:) >= ibin |  Item.Item_Bin(1,:) == 0) = -1;
                fItem = Item.Item_Bin(1,:) == -1;
                                
                LU.LU_Bin(1, LU.LU_Bin(1,:) >= ibin |  LU.LU_Bin(1,:) == 0) = -1; 
                fLu = LU.LU_Bin(1,:) == -1;

                % 2 获取s1,i1,l1; 分别增加位于剩余bin内的标记
                s1 = rmfield(Strip,{'striporder','Strip_Bin','nbLU','nbItem'});   s1.f = fStrip; %Strip
                i1  = Item;  i1.f = fItem;
                l1 = LU;      l1.f = fLu;

                % 3 进入cpuStripnbItem给了全部,增加了flag标记
                [s1.nbItem,s1.nbLU, s1.nbLULID] = cpuStripnbItem(s1,i1,l1);           
                
%                 [s1.isMixed,s1.isMixedSID,s1.isMixedEID] = isMixedStrip(s1);  % 增加了也没鸟用？
                
                % 4 获取新s2和b2 and 重新执行启发式S2B (进入HStripToBin只能给剩余部分的 x(:,f) )
                s1 = structfun(@(x) x(:,fStrip),s1,'UniformOutput',false);    % 排除已经确定入bin的剩余Strip
                [s2,b2]= HStripToBin(s1,Veh,p);
        
                % 5 获取fs: 剩余Strip摆放后的首个bin内的 so 右侧==1
                % 替换原始Strip_Bin的值
                fs = s2.Strip_Bin(1,:)==1; %找出第i1个bin的逻辑值                
                Strip.Strip_Bin(1,fidx(fs)) = ibin;
                Strip.Strip_Bin(2,fidx(fs)) = s2.Strip_Bin(2,fs);
                
                % 6 Bin更新赋值语句
                Bin.Weight(ibin) = b2.Weight(1);
                Bin.LW(:,ibin) = b2.LW(1);

                % 7 主要增加对LU_Bin的计算,
                [LU,Item] = HItemToBin(LU,Item,Strip); % 计算LU在Bin内坐标and顺序   %  Item.Item_Bin  Item.CoordItemBin LU.LU_Bin LU.CoordLUBin
                [Bin] = cpuBin(Bin,Strip,Item,LU,Veh);  %计算Bin内相关属性 % 计算isTileNeed
                
                % 8 ibin自增
                ibin = ibin+1;
                TF = true;
            end
        end
end

                                                                
%% 函数1: 量大车头方案2: V1: s1 :全部为基于Strip排除已安排Bin后的剩余Strip;;cpuStripnbItem为部分Strip进入
% % function [Strip,Bin] = HStripToBinAgain(Bin,Strip,Item,LU,Veh,p)
% %         % DONE: 量大车头方案2: 
% %         % 目的: 量大的LU被车辆拆分为量小但仍摆放车头;
% %         % 方法: 每个Bin都对排除前任Strip后的剩余Strip重新排序并分别执行HStripToBin算法.
% %         nbBin = max(Strip.Strip_Bin(1,:));
% %         if nbBin>1
% %             ibin=2;
% %             while 1
% %                 % 1 获取f and fidx: 排除首个bin后的剩余Strip逻辑值 and 索引值
% %                 Strip.Strip_Bin(1, Strip.Strip_Bin(1,:) >= ibin ) = -1; % 所有非首bin的bin序号赋值-1
% %                 % 找出未正确排序后的Strip
% %                 f = Strip.Strip_Bin(1,:) == -1; 
% %                 fidx = find(f);
% %                 if ~any(f), break; end
% %                 
% %                 % 2 获取s1 : 基于f 获取的剩余Strip结构体
% %                 s1 = rmfield(Strip,{'striporder','Strip_Bin'});
% %                 s1 = structfun(@(x) x(:,f),s1,'UniformOutput',false);
% %                 
% %                 % 3 获取i1 : 基于剩余Strip获取的剩余Item结构体
% %                 i1  = Item;
% %                 f2 = zeros(1,length(i1.Weight));
% %                 iIdxs= find(f);
% %                 for j=1:length(iIdxs)
% %                     f2 = f2+(i1.Item_Strip(1,:) == iIdxs(j));
% %                 end                
% %                 i1 = structfun(@(x) x(:,logical(f2)), i1, 'UniformOutput',false);
% %                 
% %                 % 4 获取新s2和b2 : 基于s1和i1, 重新计算Strip.nbItem and 重新执行启发式S2B
% %                 [s1] = cpuStripnbItem(s1,i1,LU);
% %                 %[s1,~] = cpuStrip(s1,i1,LU,Veh);  % 此函数包含上面的子函数, 是否有必要计算全部, 还是只计算上面子函数 TODO
% %                 
% %                 [s2,b2]= HStripToBin(s1,Veh,LU,p);
% %                 
% %                 % 5 获取fs: 剩余Strip摆放后的首个bin内的 so 右侧==1
% %                 % 替换原始Strip_Bin的值
% %                 fs = s2.Strip_Bin(1,:)==1;
% %                 Strip.Strip_Bin(1,fidx(fs)) = ibin;
% %                 Strip.Strip_Bin(2,fidx(fs)) = s2.Strip_Bin(2,fs);
% %   
% %                 % 6 赋值语句
% %                 Bin.Weight(ibin) = b2.Weight(1);
% %                 Bin.LW(:,ibin) = b2.LW(1);
% %                 ibin = ibin+1;                
% %             end
% %         end
% % end


%% %% 函数1: 量大车头方案2: V2: s1 :HStripToBin计算时:为基于Strip排除已安排Bin后的剩余Strip;cpuStripnbItem为全部Strip进入
% % function [Strip,Bin] = HStripToBinAgain(Bin,Strip,Item,LU,Veh,p)
% %         % DONE: 量大车头方案2: 
% %         % 目的: 解决量大的LU被车辆拆分为量小但仍摆放车头;
% %         % 方法: 每个Bin都对排除前任已安排Strip后的剩余Strip, 重新排序并分别执行HStripToBin算法.
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
% % function [Strip,Bin] = HreStripToEachBin(Bin,Strip,Item,LU,Veh,p)
% %        % 针对后续车辆重新排序并安排到新Bin内
% %         nbBin = max(Strip.Strip_Bin(1,:));
% %         if nbBin>100
% %             for ibin=2:nbBin
% %                 f = Strip.Strip_Bin(1,:) == ibin;
% %                 
% %                 s1 = rmfield(Strip,{'striporder','Strip_Bin'});
% %                 s1 = structfun(@(x) x(:,f),s1,'UniformOutput',false);
% %                 i1  = Item;  
% %                 
% %                 iIdxs= find(f);
% %                 f2 = zeros(1,length(i1.Weight));
% %                 for j=1:length(iIdxs)
% %                     f2 = f2+(i1.Item_Strip(1,:) == iIdxs(j));
% %                 end
% %                 
% %                 i1 = structfun(@(x) x(:,logical(f2)), i1, 'UniformOutput',false);
% %                 
% %                 % 1 重新计算Strip.nbItem and 重新执行启发式S2B
% %                 [s1] = cpuStripnbItem(s1,i1,Veh);
% %                 
% %                 [s2,b2]= HStripToBin(s1,Veh,LU,p);
% %                 
% %                 % 2 替换回原始Strip_Bin的ibin的进入顺序
% %                 Strip.Strip_Bin(2,f) = s2.Strip_Bin(2,:);
% %                 
% %                 % 3 防错语句
% %                 if any(s2.Strip_Bin(1,:)~=1), error('重新分配后bin放不下'); end
% %                 if size(b2.Weight,2) ~=1, error('重量不相同'); end
% %                 if size(b2.LW,2) ~=1, error('长宽不相同'); end
% %                 if Bin.Weight(ibin) ~= b2.Weight, error('重量不相同'); end
% %                 if Bin.LW(:,ibin) ~= b2.LW, error('长宽不相同'); end                
% %                 Bin.Weight(ibin) = b2.Weight;
% %                 Bin.LW(:,ibin) = b2.LW;
% %                 
% %             end
% %         end
% % end
        