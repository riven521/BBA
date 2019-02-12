function [flagTiledArray,do2Array,do3Array] = HBinpingpu(maind,do,p)
global ISpingpuAll ISplotEachPingPu

        % 3.1 初始化5个数据    
        flagTiledArray = zeros(1,length(do.Bin.Weight));  %1代表整车平铺 2代表甩尾平铺

        d2Array(1:length(do.Bin.Weight)) = maind; % fixme 修改maind初始化 容易造成误解
        do2Array(1:length(do.Bin.Weight)) = structfun(@(x) [], do, 'UniformOutput', false);
        d3Array(1:length(do.Bin.Weight)) = maind;
        do3Array(1:length(do.Bin.Weight)) = structfun(@(x) [], do, 'UniformOutput', false); %do;
        
        
        
        %% 3.2 进入bin循环平铺
        bidx = 1:length(do.Bin.Weight); % NOTE: 修改只考虑甩尾平铺到全部Bin纳入考虑, 对非甩尾平铺的进行平铺判断 % bidx = find(do.Bin.isTileNeed);
        % 循环: 每个bin分别尝试平铺(整车平铺和甩尾平铺只能二选一，优先整车平铺)
        for i=1:numel(bidx)
            ibin = bidx(i);
            
            % $1 GET d2/d3 本ibin内待算法计算的输入数据(均来自maind)
            luIdx = do.LU.LU_Bin(1,:) == ibin;
            d2 = getPartdinThisBin(maind,luIdx); %修改成从maind提取IuIdx个输入,而非从运算后的d中提取 %         d2 = getdinThisVeh(do,luIdx)
            d3 = d2;
            
            %% COMMENT
            %         d2.Veh = do.Veh;
            %         d2.Veh = rmfield(d2.Veh,{'Volume','order'});
            % 2 最后若干/一个strip内的LU
            %         luidx = do.LU.LU_Bin(1,:) == ibin;  %do.LU.LU_Strip(1,:) == istrip
            
            %         d2.LU = structfun(@(x) x(:,luidx),do.LU,'UniformOutput',false);
            %         d2.LU.LWH([1,2], d2.LU.Rotaed ) = flipud(d2.LU.LWH([1,2], d2.LU.Rotaed)); %LU.LWH 如旋转,则恢复原形
            %         d2.LU.PID = d2.LU.OPID;     d2.LU.SID = d2.LU.OSID;  %  d2.LU.LID = d2.LU.OLID;
            %         d2.LU = rmfield(d2.LU,{'Rotaed','order','LU_Item','DOC','LU_Strip',...
            %             'LU_Bin','CoordLUBin','CoordLUStrip','LU_VehType','OPID','OSID'});
            %
            %         d2.Par = do.Par;
            %% 3.3 如果允许全部平铺(可能是非甩尾平铺), 观察本ibin内是否可以全部平铺,如可以,就取消甩尾平铺; 否则,进入甩尾平铺
            if ISpingpuAll==1
                if ibin == numel(bidx) && range(d3.LU.ID) == 0 % 最后一车 且 一种完全一样的ID的方可1-2-3的递增平铺
                    minmaxLayer = min(max(d3.LU.maxHLayer), max(d3.LU.maxL(3,:))); %maxHLayer：指定层数; maxL(3,:)：计算最大层数;
                else
                    minmaxLayer = 1;
                end
                iLayer = 1;
                while 1
                    if iLayer > minmaxLayer, break; end
                    d3.LU.maxHLayer(:) = iLayer; %d2内全部LU的层数设定为1 55555 全部平铺的重要条件
                    iLayer=iLayer+1;
                    
                    % $3.3.1 reRunAlgorithm do3是d3运算后的结果
                    d3Array(ibin) = d3;
                    [d3.LU] = setLULWHwithbuff(d3.LU, d3.Veh);
                    do3 = RunAlgorithm(d3,p);   % do3Array(ibin) = do3;
                    
                    % $3.3.2 当全部平铺没有问题,做后处理
                    if max(do3.LU.LU_Bin(1,:)) == 1
                        flagTiledArray(ibin)=1; %1代表整车平铺
                        %do3.LU.LU_VehType = ones(size(d3.LU.ID)) * do3.Veh.order(1); % 针对车型选择,增加变量LU_VehType : 由于Veh内部按体积递减排序,获取order的第一个作为最大值
                        [do3.LU,do3.Item] = setLCwithoutbuff(do3.LU,do3.Item);               %  plot3DBPP(do3,p)
                        do3Array(ibin) = do3;
                        %                 plotSolutionT(do3.LU,do3.Veh);
                        % do3 修改到d中？？？目前保留到do2Array中，未与d合并
                        break;    %  continue;   %不进入下面的甩尾平铺了 加了while不能continue了
                    end
                end
            end
            if flagTiledArray(ibin)==1
                continue;
            end
            %% 3.4 若整车平铺失败 进入甩尾平铺 更为复杂
            
            % 3.4.1 GET do2 本ibin内含输出数据：甩尾平铺使用
            stripidx = do.Strip.Strip_Bin(1,:) == ibin; %do.LU.LU_Strip(1,:) == istrip
            itemidx = do.Item.Item_Bin(1,:) == ibin; %do.LU.LU_Strip(1,:) == istrip
            do2.Veh = do.Veh;
            do2.LU = structfun(@(x) x(:,luIdx),do.LU,'UniformOutput',false);
            do2.Bin = structfun(@(x) x(:,ibin),do.Bin,'UniformOutput',false);
            do2.Strip = structfun(@(x) x(:,stripidx),do.Strip,'UniformOutput',false);
            do2.Item = structfun(@(x) x(:,itemidx),do.Item,'UniformOutput',false);
            
            %% 获取luidxPP 并更新其 d2.LU.maxHLayer(luidxPP)
            % 555 plotSolutionT(do2.LU,do2.Veh);  % 甩尾平铺前的 观察
            while do2.Bin.isTileNeed(1) == 1 %do2内的Bin永远只有1个, 可能平铺后该bin仍需要平铺,所以有while判断
                % $3.4.2 修订d2.LU.maxHLayer (仅对ibin内最后选定的几个strip平铺) TODO $4写的有些复杂,后期简化
                % $4.1 GET luidxPP ： 某个strip对应的LU逻辑值
                % 循环从本ibin内最后一个strip开始平铺 istrip= nbStrip;
                nbStrip = numel(do2.Strip.Weight);
                if unique(do2.Strip.Strip_Bin(2, :)) ~= nbStrip,    error('超预期错误');    end
                istrip= nbStrip;
                fi = find(do2.Strip.Strip_Bin(2,:) >= istrip ); % 同bin的strip序号 and 顺序>=istrip
                u=unique(do2.LU.LU_Strip(1,:)); %获取Strip序号的唯一排序值
                luidxPP = ismember(do2.LU.LU_Strip(1,:), u(fi)); %%% fi->u(fi) 真正的序号 *********************
                if ~any(luidxPP),  error('luidxPP全部为空, 不存在u(fi)对应的Lu逻辑判断'); end
                
                % $4.2 修订d2.LU.maxHLayer           d2.LU    maind.LU
                d2.LU.maxHLayer(luidxPP) = min( d2.LU.maxL(3,luidxPP), d2.LU.maxHLayer(luidxPP)) - 1;
                
                % $4.2 若当前luidxPP对应Lu的层数均已经为1了, 则需要增加更多的istrip及luidxPP; 再修订d2.LU.maxHLayer
                % GET 更新 d2.LU.maxHLayer(luidxPP) 必须luidxPP的层数>=2层
                while all(d2.LU.maxHLayer(luidxPP)<1)
                    istrip = istrip-1;
                    if istrip==0,break;end
                    fi = find( do2.Strip.Strip_Bin(2,:) >= istrip ); %fi = find( do2.Strip.Strip_Bin(2,:) == istrip );
                    luidxPP = ismember(do2.LU.LU_Strip(1,:), u(fi)); %%% fi->u(fi) 真正的序号 *********************
                    if ~any(luidxPP),  error('luidxPP全部为空, 不存在u(fi)对应的Lu逻辑判断'); end
                    if istrip == 0,  error('此bin不存在tileneed,超预期错误');   end
                    d2.LU.maxHLayer(luidxPP) = min( d2.LU.maxL(3,luidxPP), d2.LU.maxHLayer(luidxPP)) - 1;
                end
                % 修复: 对误减的恢复为1
                d2.LU.maxHLayer(d2.LU.maxHLayer<=1) = 1;
                
                %% 运行主算法及后处理 $5 reRunAlgorithm
                d2Array(ibin) = d2;                
                [d2.LU] = setLULWHwithbuff(d2.LU, d2.Veh);
                
                do2 = RunAlgorithm(d2,p);        %  do2 = RunAlgorithmPP(d2,p);  %do2.LU.LU_VehType = ones(size(d2.LU.ID)) * do2.Veh.order(1); % 针对车型选择,增加变量LU_VehType : 由于Veh内部按体积递减排序,获取order的第一个作为最大值                
                
                [do2.LU,do2.Item] = setLCwithoutbuff(do2.LU,do2.Item);  % do2Array(ibin) = do2; 必须注释，因为是个循环
                
                % plot
                if ISplotEachPingPu == 1,     plotSolution(do2,p);       end
                
                % $6 后处理
                if max(do2.LU.LU_Bin(1,:)) == 1  %    do2.LU.LU_VehType = ones(size(do2.LU.ID))*do.Veh.order(1);
                    
                    flagTiledArray(ibin)=2; %2代表甩尾平铺
                    do2Array(ibin) = do2;
                    %                  plotSolutionT(do2.LU,do2.Veh);
                    %                  pause(0.2)
                    
                    % todo -> do2 数据不进入d 仅在return2bba中修改    % do2 数据进入d???? return2bba不修改？？？
                else
                    break;  %单车甩尾放不下 不再继续甩尾平铺了，到此结束，进入下一辆车甩尾判断
                end
                
            end % END OF WHILE
        end% END OF FOR
        
        %% chk
        if any(flagTiledArray)
            for ibin=1:length(flagTiledArray)
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
                
                flagTiledArray;
                d2Array;
                do2Array;
                d3Array;
                do3Array;
            end
        end

end


