%% V2 
function [flagTiledArray,do2Array,do3Array] = HBinpingpu(maind,do,p)
% HBinpingpu ==> 从输入数据maind，主输出数据do，参数p 获取整车平铺do3Array/甩尾平铺do2Array的输出数据
%   函数内
%   d2inIbin=d3inIbin：ibin内的输入数据，来自maind
%   do2inIbin=do3inIbin：ibin内的输出数据，来自do  
        
        global ISpingpuShuaiWei ISpingpuAll ISshuaiwei
        global ISplotEachPingPuAll ISplotEachPingPuShuaiWei 
        
        % 3.1 初始化5个输出数据
        nbin = length(do.Bin.Weight);

        flagTiledArray = zeros(1,nbin);  %1代表整车平铺 2代表甩尾平铺
        d2Array(1:nbin) = maind;    % fixme 修改maind初始化 容易造成误解  甩尾平铺
        do2Array(1:nbin) = structfun(@(x) [], do, 'UniformOutput', false);
        d3Array(1:nbin) = maind;    %整车平铺
        do3Array(1:nbin) = structfun(@(x) [], do, 'UniformOutput', false); %do;

        %% 3.2 进入bin循环平铺
        bidx = 1:nbin; % NOTE: 修改只考虑甩尾平铺到全部Bin纳入考虑, 对非甩尾平铺的进行平铺判断 % bidx = find(do.Bin.isTileNeed);
        
        % 循环: 每个bin分别尝试平铺(整车平铺和甩尾平铺只能二选一，优先整车平铺)
        for i=1:numel(bidx)
            ibin = bidx(i);
            
            % 1 GET d2inIbin/do2inIbin(d3inIbin/do3inIbin)
            d2inIbin = getPartDinThisBin(do, ibin, maind);  % V2  % V1 %    luIdx = do.LU.LU_Bin(1,:) == ibin;   %      d2inIbin = getPartDinThisBin(maind,luIdx,do); %maind提取不含margin %修改成从maind提取IuIdx个输入,而非从运算后的do中提取 %         d2inIbin = getdinThisVeh(do,luIdx)
            doinIbin = getPartDinThisBin(do,ibin);             % v2  % v1  %    luIdx = do.LU.LU_Bin(1,:) == ibin; %                 stripidx = do.Strip.Strip_Bin(1,:) == ibin; %do.LU.LU_Strip(1,:) == istrip %                 itemidx = do.Item.Item_Bin(1,:) == ibin; %do.LU.LU_Strip(1,:) == istrip %                 do2inIbin.Veh = do.Veh; %                 do2inIbin.Bin = structfun(@(x) x(:,ibin),do.Bin,'UniformOutput',false); %                 do2inIbin.LU = structfun(@(x) x(:,luIdx),do.LU,'UniformOutput',false);                %                 do2inIbin.Strip = structfun(@(x) x(:,stripidx),do.Strip,'UniformOutput',false); %                 do2inIbin.Item = structfun(@(x) x(:,itemidx),do.Item,'UniformOutput',false);

            d3inIbin = d2inIbin;  
            
            %% 2 如果允许全部平铺-待核验(可能是非甩尾平铺), 观察本ibin内是否可以全部平铺,如可以,就取消甩尾平铺; 否则,进入甩尾平铺
            if ISpingpuAll==1
                
                % 2.1 确定minmaxLayer：
                % 1 特殊case：如果最后一车 且 一种完全一样的ID的方可1-2-3的递增平铺
                if ibin == numel(bidx) && range(d3inIbin.LU.ID) == 0 
                    % minmaxLayer = min(max(d3inIbin.LU.maxHLayer), max(d3inIbin.LU.maxL(3,:))); %maxHLayer：指定层数; maxL(3,:)：计算最大层数;
                    minmaxLayer = min(max(d3inIbin.LU.maxHLayer), max(d3inIbin.LU.maxL(3,:)));  % fixme 获取可以删除 maxL的判定, 可能是个雷
                else % 2一般case：整车平铺即仅铺一层
                    minmaxLayer = 1;
                end
                
                iLayer = 1;
                
                while 1
                    if iLayer > minmaxLayer, break; end
                    
                    % update d3inIbin -> 赋值到d3array
                    d3inIbin.LU.maxHLayer(:) = iLayer;   %d3内全部LU的层数设定为1  55555 全部平铺的重要条件                    
                    d3Array(ibin) = d3inIbin;
                    
                    iLayer=iLayer+1;
                    
                    % 2.2  reRunAlgorithm do3是d3运算后的结果                    
                    d3inIbin.LU = setLULWHwithbuff(d3inIbin.LU, d3inIbin.Veh);
                    do3inIbin = RunAlgorithm(d3inIbin,p);   % do3Array(ibin) = do3; (仅对平铺成功的放入do3Array，其它的放不放无所谓，用不到
                    [do3inIbin.LU,do3inIbin.Item] = setLCwithoutbuff(do3inIbin.LU,do3inIbin.Item);               %  margin替换回来
                    
                    % 2.3   当全部平铺没有问题,做后处理
                    if max(do3inIbin.LU.LU_Bin(1,:)) == 1
                        fprintf(1,'       Exsiting 整车平铺 in HBinpingpu (do3)...\n');
                        
                        flagTiledArray(ibin)=1;   %1代表整车平铺                       %do3.LU.LU_VehType = ones(size(d3inIbin.LU.ID)) * do3.Veh.order(1); % 针对车型选择,增加变量LU_VehType : 由于Veh内部按体积递减排序,获取order的第一个作为最大值                        
                        
                        do3Array(ibin) = do3inIbin;
                                   
                        if ISplotEachPingPuAll  % plot 整车平铺前和整车平铺后的bin图
%                             plotSolutionT(doinIbin.LU,doinIbin.Veh, 0, 0, 0 , 1 ,3,'整车平铺前 Bin'); 
                            plotSolutionT(do3inIbin.LU,do3inIbin.Veh,  0, 0, 0 , 1 ,3,'整车平铺后 Bin');   % plot整车平铺后的bin % do3 修改到d中？？？目前保留到do2Array中，未与d合并
                        end
                       
                        break;    % 整车平铺成功，跳出while后循环下一个bin
                    end
                    
                end % END WHILE
                
            end % END ISpingpuAll
            
                %% 判定：若整车平铺成功，迭代下一个bin；否则尝试甩尾平铺
            if flagTiledArray(ibin)==1
                continue;
            end
            
            %% 3 若整车平铺失败 进入甩尾平铺（前提：有甩尾操作） 更为复杂  
            if ISpingpuShuaiWei==1 && ISshuaiwei==1
                
                % 甩尾平铺: 步骤1:找出需要甩尾平铺的bin,即bin.isTileNeed==1的bin
                %                步骤2:从车尾逐步找出对应高度宽度不满Strip,降低其最高堆垛的堆垛层数,进行reLoading,直到空间不够或无需甩尾了
                
                do2inIbin = doinIbin;  % 每个条带要重新看,刚开始没有条带,或者重新允许Runalgorithm,没必要,就赋值之前maind的do结果就可以
                if ~isscalar(do2inIbin.Bin.isTileNeed), error('1'); end
                                               
                while 1   
                    
                    %do2inIbin循环判定,直到Run后没有需要甩尾平铺或甩尾平铺空间不够(需要多余1个bin)
                    if do2inIbin.Bin.isTileNeed ~= 1, break; end                    
                    if ~isscalar(unique(do2inIbin.LU.LU_Bin(1,:))), break; end  %if max(do2inIbin.LU.LU_Bin(1,:)) ~= 1, break; end
                    
                    nStrip = numel(do2inIbin.Strip.Weight);      if unique(do2inIbin.Strip.Strip_Bin(2, :)) ~= nStrip,    error('超预期错误');    end
                                        
                    % 不断更新:d2inIbin的LU.maxHLayer 以进行reLoading
                    freLoading = 0;                   
                
                    for istrip = nStrip : -1 : 1  % 从车尾条带开始

                        fstrip = find(do2inIbin.Strip.Strip_Bin(2,:) == istrip); %real strip 序号
                        i1 = ~do2inIbin.Strip.isWidthFull(fstrip); %宽度不满
                        i2 = ~do2inIbin.Strip.isHeightFull(fstrip) && ~do2inIbin.Strip.isHeightBalance(fstrip); %高度不满且非均衡
                        
                        i3 = do2inIbin.Strip.isHeightFull(fstrip) &&  ~do2inIbin.Strip.isHeightBalance(fstrip); %高度满层 , 但高度非均衡
                        if i3
                            error('impossible');
                        end
                        if ~(i1 || i2 ) % 如二者满足任意一个, 就可有进行降低层数算法
                            continue;
                        end

                        % 所有等于 该 istrip的托盘索引
                        luidxPP = getLuIdx(istrip,do2inIbin);                        
                        [curLuLayer,curLuItemHeight] = getLuLayer(luidxPP,do2inIbin);
                        
                        if all(curLuLayer(luidxPP) <=1) % 如层数均为1,不可
                            continue;
                        else
                            
                            % 获取要降低层数的idx: 1:单堆垛 2:层数>1 3:最高的
                            idx = (curLuLayer > 1 & luidxPP);  
                            idx = curLuItemHeight == max(curLuItemHeight(idx));  %若多个,找出最高层堆垛 , 降低指定层数
                            
                            if ~isscalar(unique(do2inIbin.LU.LU_Item(1,idx)))
                                uniItem = unique(do2inIbin.LU.LU_Item(1,idx));
                                idx = do2inIbin.LU.LU_Item(1,:) == uniItem(1);                % 若出现,选择第一个相同高度堆垛降低层数                
                                warning('本层堆垛高度一样,高度均衡,按道理只出现在宽度不满层,也可能出现在3个堆垛,2个一样高'); 
                            end
                            
                            d2inIbin.LU.maxHLayer(idx) = curLuLayer(idx) - 1;  if any(d2inIbin.LU.maxHLayer<1), error('1'); end % 相比目前降低一层,且不可能=0.
                            freLoading = 1;
                        end
           
                        if freLoading
                            break; 
                        end
                        
                    end

                    % 循环结束: 没有i降低层数,无需reloading
                    if ~freLoading
                        break;
                    end
                    
                    %% 运行主算法及后处理 $5 reRunAlgorithm
                    tmpd2 = d2inIbin;
                    [tmpd2.LU] = setLULWHwithbuff(tmpd2.LU, tmpd2.Veh);          
                    do2inIbin = RunAlgorithm(tmpd2,p);        %  do2 = RunAlgorithmPP(d2inIbin,p);  %do2.LU.LU_VehType = ones(size(d2inIbin.LU.ID)) * do2.Veh.order(1); % 针对车型选择,增加变量LU_VehType : 由于Veh内部按体积递减排序,获取order的第一个作为最大值                   
                    [do2inIbin.LU,do2inIbin.Item] = setLCwithoutbuff(do2inIbin.LU,do2inIbin.Item);  % do2Array(ibin) = do2; 必须注释，因为是个循环

                    % $6 后处理 甩尾平铺一次不够，需要循环再进行甩尾平铺 直到超出一车容量
                    if isscalar(unique(do2inIbin.LU.LU_Bin(1,:)))  % 若所有LU在一个车内,对的 % if max(do2inIbin.LU.LU_Bin(1,:)) == 1  %    do2.LU.LU_VehType = ones(size(do2.LU.ID))*do.Veh.order(1);
                        fprintf(1,'       Exsiting 甩尾平铺 in HBinpingpu (do2)...\n');
                        
                        flagTiledArray(ibin)=2; %2代表甩尾平铺
              
                        do2Array(ibin) = do2inIbin;                        % todo -> do2 数据不进入d 仅在return2bba中修改    % do2 数据进入d???? return2bba不修改？？？
                        
                        if ISplotEachPingPuShuaiWei  % plot 甩尾平铺前和甩尾平铺后的bin图
%                             plotSolutionT(doinIbin.LU,doinIbin.Veh, 0, 0, 0 , 1 ,3,'甩尾平铺前 Bin');
                            plotSolutionT(do2inIbin.LU,do2inIbin.Veh, 0, 0, 0 , 1 ,3,'甩尾平铺后 Bin');  % plot甩尾平铺后的bin
                        end
                    end                    
                
                end % END OF do2 WHILE 本车甩尾
                
            end % END OF ISpingpuShuaiwei
            
        end% END OF FOR
        
        %% chk
        if any(flagTiledArray)
            for ibin=1:nbin
                if flagTiledArray(ibin)==1 %整车平铺
                    %                 t3 = struct2table(structfun(@(x) x',d3Array(ibin).LU,'UniformOutput',false));
                    %                 to3 = struct2table(structfun(@(x) x',do3Array(ibin).LU,'UniformOutput',false));
                    chkLUnewold(d3Array(ibin).LU,do3Array(ibin).LU);
                    chktLU(d3Array(ibin).LU);
                    chktLU(do3Array(ibin).LU);
                    
                end
                if flagTiledArray(ibin)==2 %甩尾平铺
                    %                 t2 = struct2table(structfun(@(x) x',d2Array(ibin).LU,'UniformOutput',false));
                    %                 to2 = struct2table(structfun(@(x) x',do2Array(ibin).LU,'UniformOutput',false));
%                     chkLUnewold(d2Array(ibin).LU,do2Array(ibin).LU);
                    chktLU(d2Array(ibin).LU);
                    chktLU(do2Array(ibin).LU);
                end
                
%                 flagTiledArray;
%                 d2Array;
%                 do2Array;
%                 d3Array;
%                 do3Array;
            end
        end

end

%% 局部函数
 
                                                            
%% 计算istrip对应Lu的逻辑值
function LuIdx = getLuIdx(istrip,do)

if istrip == 0,  error('此bin不存在tileneed,超预期错误');   end

u=unique(do.LU.LU_Strip(1,:));                                % 获取Strip序号的唯一排序值

fi = find(do.Strip.Strip_Bin(2,:) == istrip );              % 同bin的strip序号 and 顺序>=istrip fi = find(do.Strip.Strip_Bin(2,:) >= istrip );          

LuIdx = ismember(do.LU.LU_Strip(1,:), u(fi));          %% fi->u(fi) 真正的序号 *********************

if ~any(LuIdx),  error('LuIdx全部为空, 不存在u(fi)对应的Lu逻辑判断'); end


end

%% 计算LuIdx对应Lu的层数/高度等
function   [curLuLayer,curLuItemHeight] = getLuLayer(LuIdx,do)

    [curLuLayer,curLuItemHeight] = deal(zeros(size(LuIdx)));
      
    fidx = find(LuIdx);    
    for i=1:length(fidx)
        curLuLayer(fidx(i)) = sum(do.LU.LU_Item(1,:) == do.LU.LU_Item(1, fidx(i))); % 本idx在本bin内对应堆垛的层数
        curLuItemHeight(fidx(i)) = sum(do.LU.LWH(3, (do.LU.LU_Item(1,:) == do.LU.LU_Item(1, fidx(i))))); % 本idx在本bin内对应堆垛的层数
    end 
end


%%  V1 HBinpingpu V2: 删除注释
% % function [flagTiledArray,do2Array,do3Array] = HBinpingpu(maind,do,p)
% % % HBinpingpu ==> 从输入数据maind，主输出数据do，参数p 获取整车平铺do3Array/甩尾平铺do2Array的输出数据
% % %   函数内
% % %   d2inIbin=d3inIbin：ibin内的输入数据，来自maind
% % %   do2inIbin=do3inIbin：ibin内的输出数据，来自do  
% %         
% %         global ISpingpuShuaiWei ISpingpuAll ISshuaiwei
% %         global ISplotEachPingPuAll ISplotEachPingPuShuaiWei 
% %         
% %         % 3.1 初始化5个输出数据
% %         nbin = length(do.Bin.Weight);
% % 
% %         flagTiledArray = zeros(1,nbin);  %1代表整车平铺 2代表甩尾平铺
% %         d2Array(1:nbin) = maind;    % fixme 修改maind初始化 容易造成误解  甩尾平铺
% %         do2Array(1:nbin) = structfun(@(x) [], do, 'UniformOutput', false);
% %         d3Array(1:nbin) = maind;    %整车平铺
% %         do3Array(1:nbin) = structfun(@(x) [], do, 'UniformOutput', false); %do;
% % 
% %         %% 3.2 进入bin循环平铺
% %         bidx = 1:nbin; % NOTE: 修改只考虑甩尾平铺到全部Bin纳入考虑, 对非甩尾平铺的进行平铺判断 % bidx = find(do.Bin.isTileNeed);
% %         
% %         % 循环: 每个bin分别尝试平铺(整车平铺和甩尾平铺只能二选一，优先整车平铺)
% %         for i=1:numel(bidx)
% %             ibin = bidx(i);
% %             
% %             % 1 GET d2inIbin/do2inIbin(d3inIbin/do3inIbin)
% %             d2inIbin = getPartDinThisBin(do, ibin, maind);  % V2  % V1 %    luIdx = do.LU.LU_Bin(1,:) == ibin;   %      d2inIbin = getPartDinThisBin(maind,luIdx,do); %maind提取不含margin %修改成从maind提取IuIdx个输入,而非从运算后的do中提取 %         d2inIbin = getdinThisVeh(do,luIdx)
% %             doinIbin = getPartDinThisBin(do,ibin);             % v2  % v1  %    luIdx = do.LU.LU_Bin(1,:) == ibin; %                 stripidx = do.Strip.Strip_Bin(1,:) == ibin; %do.LU.LU_Strip(1,:) == istrip %                 itemidx = do.Item.Item_Bin(1,:) == ibin; %do.LU.LU_Strip(1,:) == istrip %                 do2inIbin.Veh = do.Veh; %                 do2inIbin.Bin = structfun(@(x) x(:,ibin),do.Bin,'UniformOutput',false); %                 do2inIbin.LU = structfun(@(x) x(:,luIdx),do.LU,'UniformOutput',false);                %                 do2inIbin.Strip = structfun(@(x) x(:,stripidx),do.Strip,'UniformOutput',false); %                 do2inIbin.Item = structfun(@(x) x(:,itemidx),do.Item,'UniformOutput',false);
% % 
% %             d3inIbin = d2inIbin;  
% %             
% %             %% 2 如果允许全部平铺-待核验(可能是非甩尾平铺), 观察本ibin内是否可以全部平铺,如可以,就取消甩尾平铺; 否则,进入甩尾平铺
% %             if ISpingpuAll==1
% %                 
% %                 % 2.1 确定minmaxLayer：
% %                 % 1 特殊case：如果最后一车 且 一种完全一样的ID的方可1-2-3的递增平铺
% %                 if ibin == numel(bidx) && range(d3inIbin.LU.ID) == 0 
% %                     % minmaxLayer = min(max(d3inIbin.LU.maxHLayer), max(d3inIbin.LU.maxL(3,:))); %maxHLayer：指定层数; maxL(3,:)：计算最大层数;
% %                     minmaxLayer = min(max(d3inIbin.LU.maxHLayer), max(d3inIbin.LU.maxL(3,:)));  % fixme 获取可以删除 maxL的判定, 可能是个雷
% %                 else % 2一般case：整车平铺即仅铺一层
% %                     minmaxLayer = 1;
% %                 end
% %                 
% %                 iLayer = 1;
% %                 
% %                 while 1
% %                     if iLayer > minmaxLayer, break; end
% %                     
% %                     % update d3inIbin -> 赋值到d3array
% %                     d3inIbin.LU.maxHLayer(:) = iLayer;   %d3内全部LU的层数设定为1  55555 全部平铺的重要条件                    
% %                     d3Array(ibin) = d3inIbin;
% %                     
% %                     iLayer=iLayer+1;
% %                     
% %                     % 2.2  reRunAlgorithm do3是d3运算后的结果                    
% %                     d3inIbin.LU = setLULWHwithbuff(d3inIbin.LU, d3inIbin.Veh);
% %                     do3inIbin = RunAlgorithm(d3inIbin,p);   % do3Array(ibin) = do3; (仅对平铺成功的放入do3Array，其它的放不放无所谓，用不到
% %                     [do3inIbin.LU,do3inIbin.Item] = setLCwithoutbuff(do3inIbin.LU,do3inIbin.Item);               %  margin替换回来
% %                     
% %                     % 2.3   当全部平铺没有问题,做后处理
% %                     if max(do3inIbin.LU.LU_Bin(1,:)) == 1
% %                         fprintf(1,'       Exsiting 整车平铺 in HBinpingpu (do3)...\n');
% %                         
% %                         flagTiledArray(ibin)=1;   %1代表整车平铺                       %do3.LU.LU_VehType = ones(size(d3inIbin.LU.ID)) * do3.Veh.order(1); % 针对车型选择,增加变量LU_VehType : 由于Veh内部按体积递减排序,获取order的第一个作为最大值                        
% %                         
% %                         do3Array(ibin) = do3inIbin;
% %                                    
% %                         if ISplotEachPingPuAll  % plot 整车平铺前和整车平铺后的bin图
% %                             plotSolutionT(doinIbin.LU,doinIbin.Veh, 0, 0, 0 , 1 ,3,'整车平铺前 Bin'); 
% %                             plotSolutionT(do3inIbin.LU,do3inIbin.Veh,  0, 0, 0 , 1 ,3,'整车平铺后 Bin');   % plot整车平铺后的bin % do3 修改到d中？？？目前保留到do2Array中，未与d合并
% %                         end
% %                        
% %                         break;    % 整车平铺成功，跳出while后循环下一个bin
% %                     end
% %                     
% %                 end % END WHILE
% %                 
% %             end % END ISpingpuAll
% %             
% % %% V1 ISpingpuAll
% % % %          if ISpingpuAll==1
% % % %                 
% % % %                 % 2.1 确定minmaxLayer：
% % % %                 % 1 特殊case：如果最后一车 且 一种完全一样的ID的方可1-2-3的递增平铺
% % % %                 if ibin == numel(bidx) && range(d3inIbin.LU.ID) == 0 
% % % %                     minmaxLayer = min(max(d3inIbin.LU.maxHLayer), max(d3inIbin.LU.maxL(3,:))); %maxHLayer：指定层数; maxL(3,:)：计算最大层数;
% % % %                 else % 2一般case：整车平铺即仅铺一层
% % % %                     minmaxLayer = 1;
% % % %                 end
% % % %                 
% % % %                 iLayer = 1;
% % % %                 
% % % %                 while 1
% % % %                     if iLayer > minmaxLayer, break; end
% % % %                     
% % % %                     % update d3inIbin -> 赋值到d3array
% % % %                     d3inIbin.LU.maxHLayer(:) = iLayer;   %d3内全部LU的层数设定为1  55555 全部平铺的重要条件                    
% % % %                     d3Array(ibin) = d3inIbin;
% % % %                     
% % % %                     iLayer=iLayer+1;
% % % %                     
% % % %                     % 2.2  reRunAlgorithm do3是d3运算后的结果                    
% % % %                     d3inIbin.LU = setLULWHwithbuff(d3inIbin.LU, d3inIbin.Veh);
% % % %                     do3inIbin = RunAlgorithm(d3inIbin,p);   % do3Array(ibin) = do3; (仅对平铺成功的放入do3Array，其它的放不放无所谓，用不到
% % % %                     [do3inIbin.LU,do3inIbin.Item] = setLCwithoutbuff(do3inIbin.LU,do3inIbin.Item);               %  margin替换回来
% % % %                     
% % % %                     % 2.3   当全部平铺没有问题,做后处理
% % % %                     if max(do3inIbin.LU.LU_Bin(1,:)) == 1
% % % %                         fprintf(1,'       Exsiting 整车平铺 in HBinpingpu (do3)...\n');
% % % %                         
% % % %                         flagTiledArray(ibin)=1;   %1代表整车平铺                       %do3.LU.LU_VehType = ones(size(d3inIbin.LU.ID)) * do3.Veh.order(1); % 针对车型选择,增加变量LU_VehType : 由于Veh内部按体积递减排序,获取order的第一个作为最大值                        
% % % %                         
% % % %                         do3Array(ibin) = do3inIbin;
% % % %                                    
% % % %                         if ISplotEachPingPuAll  % plot 整车平铺前和整车平铺后的bin图
% % % %                             plotSolutionT(doinIbin.LU,doinIbin.Veh, 0, 0, 0 , 1 ,3,'整车平铺前 Bin'); 
% % % %                             plotSolutionT(do3inIbin.LU,do3inIbin.Veh,  0, 0, 0 , 1 ,3,'整车平铺后 Bin');   % plot整车平铺后的bin % do3 修改到d中？？？目前保留到do2Array中，未与d合并
% % % %                         end
% % % %                        
% % % %                         break;    % 整车平铺成功，跳出while后循环下一个bin
% % % %                     end
% % % %                     
% % % %                 end % END WHILE
% % % %                 
% % % %             end % END ISpingpuAll
% %             
% %                 %% 判定：若整车平铺成功，迭代下一个bin；否则尝试甩尾平铺
% %             if flagTiledArray(ibin)==1
% %                 continue;
% %             end
% %             
% %             %% 3 若整车平铺失败 进入甩尾平铺（前提：有甩尾操作） 更为复杂  
% %             if ISpingpuShuaiWei==1 && ISshuaiwei==1
% %                 
% %                 % 甩尾平铺: 步骤1:找出需要甩尾平铺的bin,即bin.isTileNeed==1的bin
% %                 %                步骤2:从车尾逐步找出对应高度宽度不满Strip,降低其最高堆垛的堆垛层数,进行reLoading,直到空间不够
% %                 
% %                 do2inIbin = doinIbin;  % 每个条带要重新看,刚开始没有条带,或者重新允许Runalgorithm,没必要,就赋值之前maind的do结果就可以
% %                 if ~isscalar(do2inIbin.Bin.isTileNeed), error('1'); end
% %                 
% %                 plotSolutionT(doinIbin.LU,doinIbin.Veh, 0, 0, 0 , 1 ,3,'甩尾平铺前 Bin');
% %                 plotSolutionT(doinIbin.LU,doinIbin.Veh, 0, 0, 0 , 1 ,8,'甩尾平铺前 Bin');
% %                                
% %                 while 1   
% %                     
% %                     %do2inIbin循环判定,直到Run后没有需要甩尾平铺或甩尾平铺空间不够(需要多余1个bin)
% %                     if do2inIbin.Bin.isTileNeed ~= 1, break; end                    
% %                     if ~isscalar(unique(do2inIbin.LU.LU_Bin(1,:))), break; end  %if max(do2inIbin.LU.LU_Bin(1,:)) ~= 1, break; end
% %                     
% %                     nStrip = numel(do2inIbin.Strip.Weight);      if unique(do2inIbin.Strip.Strip_Bin(2, :)) ~= nStrip,    error('超预期错误');    end
% %                                         
% %                     % 不断更新:d2inIbin的LU.maxHLayer 以进行reLoading
% %                     freLoading = 0;
% %                     
% % %                     uniStrip = unique(do2inIbin.LU.LU_Strip(1,:));
% %                     
% %                     for istrip = nStrip : -1 : 1  % 从车尾条带开始
% % 
% % %                         do2inIbin.LU.LU_Strip(1,:) == istrip
% % %                         do2inIbin.Strip.Strip_Bin
% %                         fstrip = find(do2inIbin.Strip.Strip_Bin(2,:) == istrip)
% %                         
% %                         
% %                         i1 = ~do2inIbin.Strip.isWidthFull(fstrip); %宽度不满
% %                         i2 = ~do2inIbin.Strip.isHeightFull(fstrip) && ~do2inIbin.Strip.isHeightBalance(fstrip); %高度不满且非均衡
% %                         
% %                         if ~(i1 || i2) % 如二者有一个不满足, 就继续需要甩尾平铺
% %                             continue;
% %                         end
% % 
% %                         % 所有等于 该 istrip的托盘索引
% %                         luidxPP = getLuIdx(istrip,do2inIbin);                        
% %                         [curLuLayer,curLuItemHeight] = getLuLayer(luidxPP,do2inIbin);
% %                         
% %                         if all(curLuLayer(luidxPP) <=1)
% %                             continue;
% %                         else                            
% %                             idx = (curLuLayer > 1 & luidxPP);  % 肯定有curLuLayer>=2
% %                             % 不属于stripbalance的
% %                             idx2 = curLuItemHeight == max(curLuItemHeight(idx));  %若多个,找出最高层,降低指定层数
% %                             
% %                             
% %                                     if ~isscalar(unique(do2inIbin.LU.LU_Item(1,idx2))), 
% %                                         uniItem = unique(do2inIbin.LU.LU_Item(1,idx2))
% %                                         idx2 = do2inIbin.LU.LU_Item(1,:) == uniItem(1)
% %                                         
% %                                         warning('本层堆垛高度一样,高度均衡,按道理只出现在宽度不满层'); end
% %                                     
% %                             d2inIbin.LU.maxHLayer(idx2) = curLuLayer(idx2) - 1; % 相比目前降低一层,且不可能=0.
% %                                                     if any(d2inIbin.LU.maxHLayer<1), error('1'); end
% %                             freLoading = 1;
% %                         end
% %                         
% % %                         d2inIbin.LU.maxHLayer(luidxPP) = min( d2inIbin.LU.maxL(3,luidxPP), d2inIbin.LU.maxHLayer(luidxPP));                                 
% % %                         d2inIbin.LU.maxHLayer(luidxPP) = min( d2inIbin.LU.maxHLayer(luidxPP));                        
% % %                         d2inIbin.LU.maxHLayer(luidxPP) 
% % %                         curLuLayer(luidxPP)
% % 
% % %                         while any(d2inIbin.LU.maxHLayer(luidxPP) > 1)  % 若当前luidxPP内的托盘 存在非平铺的LU, 更新
% %                             
% %                             % 当前Lu,且堆垛大于1,且相同堆垛最高的
% % %                             a = curLuItemHeight == max(curLuItemHeight(luidxPP)) & luidxPP & curLuLayer > 1;
% %                             
% % %                              m2 = max(curLuLayer(luidxPP)) %最高堆垛的层数
% % %                              idx2 = curLuLayer(:) == m2 %最高堆垛的层数所在index号
% %                              
% % %                              m = max(d2inIbin.LU.maxHLayer(idx2))
% % %                             idx = d2inIbin.LU.maxHLayer(:) == m  %只减最高值应该可有,但低值的再maxHlayer里为高值too
% % %                             a = idx'&luidxPP
% %                             
% % %                             if sum(a) <=0 , 
% % %                                 continue;
% % %                             end
% % %                             d2inIbin.LU.maxHLayer(a) = d2inIbin.LU.maxHLayer(a) - 1;
% %                             
% % %                             d2inIbin.LU.maxHLayer(luidxPP) = d2inIbin.LU.maxHLayer(luidxPP)  - 1;
% % %                             d2inIbin.LU.maxHLayer(d2inIbin.LU.maxHLayer<=1) = 1;  % 修复: 对误减的恢复为1
% %                             
% % %                             d2inIbin.LU.maxHLayer(luidxPP) 
% %                             
% % %                             if any(d2inIbin.LU.maxHLayer(luidxPP) < curLuLayer(luidxPP))  % 若更新后最大层数 < 当前LU对应层数
% % %                                 d2inIbin.LU.maxHLayer(luidxPP)
% % %                                 curLuLayer(luidxPP)
% % %                                 freLoading = 1;
% % %                                 break;
% % %                             else
% % %                                 continue;
% % %                             end
% % 
% % %                         end
% %                         
% %                         if freLoading
% %                             break; 
% %                         end
% %                         
% %                     end
% % %%                    
% % %                     for istrip = nStrip : -1 : 1
% % %                         
% % %                         luidxPP = getLuIdx(istrip,do2inIbin);
% % %                         
% % %                         curLuLayer = getLuLayer(luidxPP,do2inIbin);
% % %                         
% % % %                         d2inIbin.LU.maxHLayer(luidxPP) = min( d2inIbin.LU.maxL(3,luidxPP), d2inIbin.LU.maxHLayer(luidxPP));                                 
% % % %                         d2inIbin.LU.maxHLayer(luidxPP) = min( d2inIbin.LU.maxHLayer(luidxPP));                        
% % %                         d2inIbin.LU.maxHLayer(luidxPP) 
% % %                         
% % %                         if any(d2inIbin.LU.maxHLayer(luidxPP) > 1)  % 若当前luidxPP内的托盘 存在非平铺的LU, 更新
% % %                             d2inIbin.LU.maxHLayer(luidxPP) = d2inIbin.LU.maxHLayer(luidxPP)  - 1;
% % %                             d2inIbin.LU.maxHLayer(d2inIbin.LU.maxHLayer<=1) = 1;  % 修复: 对误减的恢复为1
% % %                             
% % %                             d2inIbin.LU.maxHLayer(luidxPP) 
% % %                             
% % %                             if any(d2inIbin.LU.maxHLayer(luidxPP) < curLuLayer(luidxPP))  % 若更新后最大层数 < 当前LU对应层数
% % %                                 d2inIbin.LU.maxHLayer(luidxPP)
% % %                                 curLuLayer(luidxPP)
% % %                                 fre = 1;
% % %                                 break;
% % %                             else
% % %                                 continue;
% % %                             end
% % %                             
% % %                         else
% % %                             continue; 
% % %                         end
% % %                         
% % %                     end
% %                     
% %                     
% % %% 老版本备份                    
% % %                     istrip= nStrip;
% % %                     
% % %                     luidxPP = getLuIdx(istrip,do2inIbin);
% % %                                         
% % %                     d2inIbin.LU.maxHLayer(luidxPP) = min( d2inIbin.LU.maxL(3,luidxPP), d2inIbin.LU.maxHLayer(luidxPP)) - 1;
% % %                     
% % %                     while all(d2inIbin.LU.maxHLayer(luidxPP)<1)
% % %                         
% % %                         istrip = istrip-1;
% % %                         
% % %                         if istrip==0,  break;  end
% % %                         
% % %                         luidxPP = getLuIdx(istrip,do2inIbin);
% % %                         
% % %                        d2inIbin.LU.maxHLayer(luidxPP) = min( d2inIbin.LU.maxL(3,luidxPP), d2inIbin.LU.maxHLayer(luidxPP)) - 1;
% % %                         
% % %                     end
% % %                     
% % %                     d2inIbin.LU.maxHLayer(d2inIbin.LU.maxHLayer<=1) = 1;  % 修复: 对误减的恢复为1
% %                     if ~freLoading
% %                         break;
% %                     end
% %                     %% 运行主算法及后处理 $5 reRunAlgorithm
% %                     tmpd2 = d2inIbin;
% %                     d2inIbin.LU.maxHLayer(luidxPP) 
% %                     [tmpd2.LU] = setLULWHwithbuff(tmpd2.LU, tmpd2.Veh);          
% %                     do2inIbin = RunAlgorithm(tmpd2,p);        %  do2 = RunAlgorithmPP(d2inIbin,p);  %do2.LU.LU_VehType = ones(size(d2inIbin.LU.ID)) * do2.Veh.order(1); % 针对车型选择,增加变量LU_VehType : 由于Veh内部按体积递减排序,获取order的第一个作为最大值                   
% %                     [do2inIbin.LU,do2inIbin.Item] = setLCwithoutbuff(do2inIbin.LU,do2inIbin.Item);  % do2Array(ibin) = do2; 必须注释，因为是个循环
% %                     do2inIbin.LU.maxHLayer(luidxPP) 
% % 
% %                     % $6 后处理 甩尾平铺一次不够，需要循环再进行甩尾平铺 直到超出一车容量
% %                     if isscalar(unique(do2inIbin.LU.LU_Bin(1,:)))  % 若所有LU在一个车内,对的
% %                     % if max(do2inIbin.LU.LU_Bin(1,:)) == 1  %    do2.LU.LU_VehType = ones(size(do2.LU.ID))*do.Veh.order(1);
% %                         fprintf(1,'       Exsiting 甩尾平铺 in HBinpingpu (do2)...\n');
% %                         
% %                         flagTiledArray(ibin)=2; %2代表甩尾平铺
% %               
% %                         do2Array(ibin) = do2inIbin;                        % todo -> do2 数据不进入d 仅在return2bba中修改    % do2 数据进入d???? return2bba不修改？？？
% %                         
% %                         if ISplotEachPingPuShuaiWei  % plot 甩尾平铺前和甩尾平铺后的bin图
% % %                             plotSolutionT(doinIbin.LU,doinIbin.Veh, 0, 0, 0 , 1 ,3,'甩尾平铺前 Bin');
% %                             plotSolutionT(do2inIbin.LU,do2inIbin.Veh, 0, 0, 0 , 1 ,3,'甩尾平铺后 Bin');  % plot甩尾平铺后的bin
% %                         end
% %                     end                    
% %                 
% %                 end % END OF do2 WHILE 本车甩尾
% %                 
% %                 
% %                 
% % % %                 %% 获取luidxPP 并更新其 d2inIbin.LU.maxHLayer(luidxPP)
% % % %                 while do2inIbin.Bin.isTileNeed == 1 %do2inIbin可能平铺后该bin仍需要平铺,所以有while判断
% % % %                     % $3.4.2 修订d2.LU.maxHLayer (仅对ibin内最后选定的几个strip平铺) TODO $4写的有些复杂,后期简化
% % % %                     % $4.1 从do2inIbin获取 luidxPP ： 某个strip对应的LU逻辑值,
% % % %                     % 循环从本ibin内最后一个strip开始平铺 istrip= nbStrip;
% % % %                     nStrip = numel(do2inIbin.Strip.Weight);  if unique(do2inIbin.Strip.Strip_Bin(2, :)) ~= nStrip,    error('超预期错误');    end
% % % %                     
% % % %                     % 从istrip即最后一个条带开始
% % % %                     istrip= nStrip;
% % % %                     
% % % %                     luidxPP = getLuIdx(istrip,do2inIbin);
% % % %                                         
% % % %                     % $4.2 修订d2.LU.maxHLayer           d2inIbin.LU    maind.LU
% % % %                     % 将本bin内的 某些托盘 (luidxPP) 的允许最大摆放层数 - 1 (先削一层,可能存在同一Strip层数不同case,TODO)
% % % % %                     d2inIbin.LU.maxHLayer(luidxPP)
% % % % %                     d2inIbin.LU.maxL(3,luidxPP)
% % % %                     d2inIbin.LU.maxHLayer(luidxPP) = min( d2inIbin.LU.maxL(3,luidxPP), d2inIbin.LU.maxHLayer(luidxPP)) - 1;
% % % %                     
% % % %                     % $4.2 若当前luidxPP对应Lu的层数均已经为1了, 则需要增加更多的istrip及luidxPP; 再修订d2.LU.maxHLayer
% % % %                     % GET 更新 d2inIbin.LU.maxHLayer(luidxPP) 必须luidxPP的层数>=2层
% % % %                     while all(d2inIbin.LU.maxHLayer(luidxPP)<1)
% % % %                         
% % % %                         istrip = istrip-1;
% % % %                         
% % % %                         if istrip==0,  break;  end
% % % %                         
% % % %                         luidxPP = getLuIdx(istrip,do2inIbin);
% % % %                         
% % % %                        d2inIbin.LU.maxHLayer(luidxPP) = min( d2inIbin.LU.maxL(3,luidxPP), d2inIbin.LU.maxHLayer(luidxPP)) - 1;
% % % %                         
% % % %                     end
% % % %                     % 修复: 对误减的恢复为1
% % % %                     d2inIbin.LU.maxHLayer(d2inIbin.LU.maxHLayer<=1) = 1;
% % % %                     
% % % %                     %% 运行主算法及后处理 $5 reRunAlgorithm
% % % % %                     d2Array(ibin) = do2inIbin;
% % % % %                     d2inIbin.LU.Rotaed(:)
% % % % %                     d2inIbin.LU.Rotaed(:) = 0;
% % % % %                     d2inIbin.LU.Rotaed(:)
% % % % d2inIbin111 = d2inIbin
% % % %                     [d2inIbin111.LU] = setLULWHwithbuff(d2inIbin111.LU, d2inIbin111.Veh);          
% % % % %                     d2inIbin.LU.Rotaed(:)
% % % %                      d2inIbin.LU.LWH
% % % %                     do2inIbin = RunAlgorithm(d2inIbin111,p);        %  do2 = RunAlgorithmPP(d2inIbin,p);  %do2.LU.LU_VehType = ones(size(d2inIbin.LU.ID)) * do2.Veh.order(1); % 针对车型选择,增加变量LU_VehType : 由于Veh内部按体积递减排序,获取order的第一个作为最大值                   
% % % %                     [do2inIbin.LU,do2inIbin.Item] = setLCwithoutbuff(do2inIbin.LU,do2inIbin.Item);  % do2Array(ibin) = do2; 必须注释，因为是个循环
% % % % %                     [d2inIbin.LU,d2inIbin.Item] = setLCwithoutbuff(d2inIbin.LU,d2inIbin.Item);
% % % %                     
% % % % %                     d2inIbin.LU.LWH = do2inIbin.LU.LWH
% % % % %                     d2inIbin.LU.CoordLUBin = do2inIbin.LU.CoordLUBin
% % % % %                     d2inIbin.Item.LWH = do2inIbin.Item.LWH
% % % % %                     d2inIbin.Item.CoordItemBin = do2inIbin.Item.CoordItemBin
% % % %                     
% % % %                     plotSolutionT(do2inIbin.LU,do2inIbin.Veh, 0, 0, 0 , 1 ,3,'甩尾平铺后可能多个bin Bin');
% % % %                     % plot
% % % %                     
% % % %                     
% % % %                     % $6 后处理 甩尾平铺一次不够，需要循环再进行甩尾平铺 直到超出一车容量
% % % %                     if max(do2inIbin.LU.LU_Bin(1,:)) == 1  %    do2.LU.LU_VehType = ones(size(do2.LU.ID))*do.Veh.order(1);
% % % %                         fprintf(1,'       Exsiting 甩尾平铺 in HBinpingpu (do2)...\n');
% % % %                         flagTiledArray(ibin)=2; %2代表甩尾平铺
% % % %                         
% % % % %                         [do2inIbin.LU,do2inIbin.Item] = setLCwithoutbuff(do2inIbin.LU,do2inIbin.Item); 
% % % %                         
% % % %                         do2Array(ibin) = do2inIbin;                        % todo -> do2 数据不进入d 仅在return2bba中修改    % do2 数据进入d???? return2bba不修改？？？
% % % %                         
% % % %                         if ISplotEachPingPuShuaiWei  % plot 甩尾平铺前和甩尾平铺后的bin图
% % % %                             plotSolutionT(doinIbin.LU,doinIbin.Veh, 0, 0, 0 , 1 ,3,'甩尾平铺前 Bin');
% % % %                             plotSolutionT(do2inIbin.LU,do2inIbin.Veh, 0, 0, 0 , 1 ,3,'甩尾平铺后 Bin');  % plot甩尾平铺后的bin
% % % %                         end
% % % %                         
% % % %                     else
% % % %                         break;  %单车甩尾放不下 不再继续甩尾平铺了，到此结束，进入下一辆车甩尾判断
% % % %                     end                    
% % % %                 
% % % %                 end % END OF do2 WHILE 本车甩尾
% %                 
% %             end % END OF ISpingpuShuaiwei
% %             
% %         end% END OF FOR
% %         
% %         %% chk
% %         if any(flagTiledArray)
% %             for ibin=1:nbin
% %                 if flagTiledArray(ibin)==1 %整车平铺
% %                     %                 t3 = struct2table(structfun(@(x) x',d3Array(ibin).LU,'UniformOutput',false));
% %                     %                 to3 = struct2table(structfun(@(x) x',do3Array(ibin).LU,'UniformOutput',false));
% %                     chkLUnewold(d3Array(ibin).LU,do3Array(ibin).LU);
% %                     chktLU(d3Array(ibin).LU);
% %                     chktLU(do3Array(ibin).LU);
% %                     
% %                 end
% %                 if flagTiledArray(ibin)==2 %甩尾平铺
% %                     %                 t2 = struct2table(structfun(@(x) x',d2Array(ibin).LU,'UniformOutput',false));
% %                     %                 to2 = struct2table(structfun(@(x) x',do2Array(ibin).LU,'UniformOutput',false));
% % %                     chkLUnewold(d2Array(ibin).LU,do2Array(ibin).LU);
% %                     chktLU(d2Array(ibin).LU);
% %                     chktLU(do2Array(ibin).LU);
% %                 end
% %                 
% % %                 flagTiledArray;
% % %                 d2Array;
% % %                 do2Array;
% % %                 d3Array;
% % %                 do3Array;
% %             end
% %         end
% % 
% % end