%% 函数: Strip均衡
function   [Strip,Item,LU] = HStripBalance(Strip,Item,LU,Veh,p)
        flagBalanced = 0;
        while 1
            % 找到strip不混合的且不平衡的
            fstrip = Strip.isMixed==0 & Strip.isHeightBalance==0;
            if ~any(fstrip) || flagBalanced==1, break; end
            idxstr = find(fstrip);
            nbstrip = numel(idxstr);
            
            for idx=1:nbstrip
                luidxPP = ismember(LU.LU_Strip(1,:), idxstr(idx));
                
                % FINDOUT 本strip对应Item的LU个数,通过Item
                tmpItems = LU.LU_Item(1,luidxPP);
                [~,~,ic] = unique(tmpItems);
                nbLuArray= accumarray(ic,1);
                
                if length(unique(nbLuArray)) >2
                                %                     nbLuArray
                                %                     diff(nbLuArray)
                                                    %          plotSolutionT(LU,d.Veh);  %???
                                                    %          为何会有nbLuArray为标量 
                                                    %          为何会有Strip内只有1个Item的层进来? 需要CHECK
                    error('存在unblanced的堆垛小于2或唯一数超过2，超预期错误');        
                end
                if isscalar(nbLuArray) 
                    warning('nbLuArray仅有1个值,问题待查');
                end
                
               
                if  isscalar(nbLuArray)                         %此strip没有balance的必要 只有1个Item
                    flagBalanced = 1;
                    continue;
                elseif abs(diff(nbLuArray)) < 2          %此strip没有balance的必要 层数差不到2个
                    flagBalanced = 1;
                    continue;
                else
                    % 直接降到目前最大Item堆放的层数，后续再-1
                    LU.maxHLayer(luidxPP) = max(nbLuArray);
                end

                LU.maxHLayer(luidxPP) = min( max(LU.maxL(3,luidxPP)), max(LU.maxHLayer(luidxPP))) - 1;

                % 若降低层数到0 或 降低到原先LU的最小层数, 即可能错无论, 跳出本strip的高度均衡
                if any(LU.maxHLayer(luidxPP)<=0) || any(LU.maxHLayer(luidxPP) <= min(nbLuArray))
                    LU.maxHLayer(LU.maxHLayer(luidxPP)<=0) = 1;
                    flagBalanced = 1;
                    continue;
                end
                
                % 555 对降低层数的LU重新获取新的Strip
                %   [LU,Veh] = cpuLUVeh(LU,Veh);
                [LU] = cpuLU(LU,Veh);
                [LU,Item] = HLUtoItem(LU,Veh);          %Item将按ID序号排序（但下一操作将变化顺序）
                [Item,LU] = cpuItem(Item,LU,Veh);        % printstruct(d,'sortfields',1,'PRINTCONTENTS',0);
                [LU,Item,Strip] = HItemToStrip(LU,Item,Veh,p);
                [Strip,LU] = cpuStrip(Strip,Item,LU,Veh);
                % fstrip = Strip.isMixed==0 & Strip.isHeightBalance==0;
                flagBalanced = 0;
            end
        end           
end
