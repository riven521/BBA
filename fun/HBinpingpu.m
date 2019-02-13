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
            do2inIbin = getPartDinThisBin(do,ibin);             % v2  % v1  %    luIdx = do.LU.LU_Bin(1,:) == ibin; %                 stripidx = do.Strip.Strip_Bin(1,:) == ibin; %do.LU.LU_Strip(1,:) == istrip %                 itemidx = do.Item.Item_Bin(1,:) == ibin; %do.LU.LU_Strip(1,:) == istrip %                 do2inIbin.Veh = do.Veh; %                 do2inIbin.Bin = structfun(@(x) x(:,ibin),do.Bin,'UniformOutput',false); %                 do2inIbin.LU = structfun(@(x) x(:,luIdx),do.LU,'UniformOutput',false);                %                 do2inIbin.Strip = structfun(@(x) x(:,stripidx),do.Strip,'UniformOutput',false); %                 do2inIbin.Item = structfun(@(x) x(:,itemidx),do.Item,'UniformOutput',false);

            d3inIbin = d2inIbin;
            do3inIbin = do2inIbin;
            
            %% COMMENT
            %         d2inIbin.Veh = do.Veh;
            %         d2inIbin.Veh = rmfield(d2inIbin.Veh,{'Volume','order'});
            % 2 最后若干/一个strip内的LU
            %         luidx = do.LU.LU_Bin(1,:) == ibin;  %do.LU.LU_Strip(1,:) == istrip
            
            %         d2inIbin.LU = structfun(@(x) x(:,luidx),do.LU,'UniformOutput',false);
            %         d2inIbin.LU.LWH([1,2], d2inIbin.LU.Rotaed ) = flipud(d2inIbin.LU.LWH([1,2], d2inIbin.LU.Rotaed)); %LU.LWH 如旋转,则恢复原形
            %         d2inIbin.LU.PID = d2inIbin.LU.OPID;     d2inIbin.LU.SID = d2inIbin.LU.OSID;  %  d2inIbin.LU.LID = d2inIbin.LU.OLID;
            %         d2inIbin.LU = rmfield(d2inIbin.LU,{'Rotaed','order','LU_Item','DOC','LU_Strip',...
            %             'LU_Bin','CoordLUBin','CoordLUStrip','LU_VehType','OPID','OSID'});
            %
            %         d2inIbin.Par = do.Par;
            %% 2 如果允许全部平铺(可能是非甩尾平铺), 观察本ibin内是否可以全部平铺,如可以,就取消甩尾平铺; 否则,进入甩尾平铺
            if ISpingpuAll==1
                
                % 2.1 确定minmaxLayer：
                % 1 特殊case：如果最后一车 且 一种完全一样的ID的方可1-2-3的递增平铺
                if ibin == numel(bidx) && range(d3inIbin.LU.ID) == 0 
                    minmaxLayer = min(max(d3inIbin.LU.maxHLayer), max(d3inIbin.LU.maxL(3,:))); %maxHLayer：指定层数; maxL(3,:)：计算最大层数;
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
                    do3 = RunAlgorithm(d3inIbin,p);   % do3Array(ibin) = do3; (仅对平铺成功的放入do3Array，其它的放不放无所谓，用不到
                    [do3.LU,do3.Item] = setLCwithoutbuff(do3.LU,do3.Item);               %  margin替换回来
                    
                    % 2.3   当全部平铺没有问题,做后处理
                    if max(do3.LU.LU_Bin(1,:)) == 1
                        fprintf(1,'       Exsiting 整车平铺 in HBinpingpu (do3)...\n');
                        flagTiledArray(ibin)=1;   %1代表整车平铺                       %do3.LU.LU_VehType = ones(size(d3inIbin.LU.ID)) * do3.Veh.order(1); % 针对车型选择,增加变量LU_VehType : 由于Veh内部按体积递减排序,获取order的第一个作为最大值                        
                        
                        do3Array(ibin) = do3;
                                   
                        if ISplotEachPingPuAll == 1  % plot 甩尾平铺前和甩尾平铺后的bin图
                            plotSolutionT(do3inIbin.LU,do3inIbin.Veh, 0, 0);
                            plotSolutionT(do3.LU,do3.Veh, 0, 0);   % plot整车平铺后的bin % do3 修改到d中？？？目前保留到do2Array中，未与d合并
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
                             
                %% 获取luidxPP 并更新其 d2inIbin.LU.maxHLayer(luidxPP)
                %            plotSolutionT(do2inIbin.LU,do2inIbin.Veh);  % 甩尾平铺前的 观察
                while do2inIbin.Bin.isTileNeed(1) == 1 %do2inIbin可能平铺后该bin仍需要平铺,所以有while判断
                    % $3.4.2 修订d2.LU.maxHLayer (仅对ibin内最后选定的几个strip平铺) TODO $4写的有些复杂,后期简化
                    % $4.1 GET luidxPP ： 某个strip对应的LU逻辑值
                    % 循环从本ibin内最后一个strip开始平铺 istrip= nbStrip;
                    nbStrip = numel(do2inIbin.Strip.Weight);  
                        if unique(do2inIbin.Strip.Strip_Bin(2, :)) ~= nbStrip,    error('超预期错误');    end
                    istrip= nbStrip;
                    fi = find(do2inIbin.Strip.Strip_Bin(2,:) >= istrip ); % 同bin的strip序号 and 顺序>=istrip
                    u=unique(do2inIbin.LU.LU_Strip(1,:)); %获取Strip序号的唯一排序值
                    luidxPP = ismember(do2inIbin.LU.LU_Strip(1,:), u(fi)); %%% fi->u(fi) 真正的序号 *********************
                    if ~any(luidxPP),  error('luidxPP全部为空, 不存在u(fi)对应的Lu逻辑判断'); end
                    
                    % $4.2 修订d2.LU.maxHLayer           d2inIbin.LU    maind.LU
                    d2inIbin.LU.maxHLayer(luidxPP) = min( d2inIbin.LU.maxL(3,luidxPP), d2inIbin.LU.maxHLayer(luidxPP)) - 1;
                    
                    % $4.2 若当前luidxPP对应Lu的层数均已经为1了, 则需要增加更多的istrip及luidxPP; 再修订d2.LU.maxHLayer
                    % GET 更新 d2inIbin.LU.maxHLayer(luidxPP) 必须luidxPP的层数>=2层
                    while all(d2inIbin.LU.maxHLayer(luidxPP)<1)
                        istrip = istrip-1;
                        if istrip==0,break;end
                        fi = find( do2inIbin.Strip.Strip_Bin(2,:) >= istrip ); %fi = find( do2inIbin.Strip.Strip_Bin(2,:) == istrip );
                        luidxPP = ismember(do2inIbin.LU.LU_Strip(1,:), u(fi)); %%% fi->u(fi) 真正的序号 *********************
                        if ~any(luidxPP),  error('luidxPP全部为空, 不存在u(fi)对应的Lu逻辑判断'); end
                        if istrip == 0,  error('此bin不存在tileneed,超预期错误');   end
                        d2inIbin.LU.maxHLayer(luidxPP) = min( d2inIbin.LU.maxL(3,luidxPP), d2inIbin.LU.maxHLayer(luidxPP)) - 1;
                    end
                    % 修复: 对误减的恢复为1
                    d2inIbin.LU.maxHLayer(d2inIbin.LU.maxHLayer<=1) = 1;
                    
                    %% 运行主算法及后处理 $5 reRunAlgorithm
                    d2Array(ibin) = d2inIbin;
                    
                    [d2inIbin.LU] = setLULWHwithbuff(d2inIbin.LU, d2inIbin.Veh);                    
                    do2 = RunAlgorithm(d2inIbin,p);        %  do2 = RunAlgorithmPP(d2inIbin,p);  %do2.LU.LU_VehType = ones(size(d2inIbin.LU.ID)) * do2.Veh.order(1); % 针对车型选择,增加变量LU_VehType : 由于Veh内部按体积递减排序,获取order的第一个作为最大值                   
                    [do2.LU,do2.Item] = setLCwithoutbuff(do2.LU,do2.Item);  % do2Array(ibin) = do2; 必须注释，因为是个循环
                    
                    % plot
                    
                    
                    % $6 后处理 甩尾平铺一次不够，需要循环再进行甩尾平铺 直到超出一车容量
                    if max(do2.LU.LU_Bin(1,:)) == 1  %    do2.LU.LU_VehType = ones(size(do2.LU.ID))*do.Veh.order(1);
                        fprintf(1,'       Exsiting 甩尾平铺 in HBinpingpu (do2)...\n');
                        flagTiledArray(ibin)=2; %2代表甩尾平铺
                        
                        do2Array(ibin) = do2;                        % todo -> do2 数据不进入d 仅在return2bba中修改    % do2 数据进入d???? return2bba不修改？？？
                        
                        if ISplotEachPingPuShuaiWei == 1  % plot 甩尾平铺前和甩尾平铺后的bin图
                            plotSolutionT(do2inIbin.LU,do2inIbin.Veh, 0, 0);
                            plotSolutionT(do2.LU,do2.Veh, 0, 0);  % plot甩尾平铺后的bin
                        end
                        
                    else
                        break;  %单车甩尾放不下 不再继续甩尾平铺了，到此结束，进入下一辆车甩尾判断
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
                    chkLUnewold(d2Array(ibin).LU,do2Array(ibin).LU);
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


