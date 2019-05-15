function T = HBinCombine(do,flaggetSmallVeh,do1,flagTiledArray,do2Array,do3Array)
% HBinCombine ==> 对原始 改变车型 平铺 混合间隙 后的结果进行合并

global ISlastVehType ISpingpu parGap

%% 1 ******************获取展示顺序 do数据 -> T
    T = getTableLU(do); chktLU(T);

    % 作图：原始非平铺的图
    %     plotSolutionT(T,struct2table(structfun(@(x) x',d.Veh,'UniformOutput',false)));

%%
% 行1：托盘所在车型号(必须换)    行2：托盘所在车序号(会变,不能换,换就错) 行3：托盘车内安置顺序(必须换) 行4：托盘SID供应商编号(不会变,不用变??)
% 行5：托盘ID型号LID(不会变,不用变?) 行6：托盘堆垛序号ITEM(会变,不能换,换就错,用途?) 行7：托盘零部件编号PID(不会变,不用变?) 增加行8: 展示顺序(必须换)

%% 2 ****************** 针对车型变化 do1数据 获取修订的 output ******************
    if ISlastVehType==1 && flaggetSmallVeh == 1 %如有当允许且车型替换成功
        T1 = getTableLU(do1);     chktLU(T1);

        % 获取 最后一车的逻辑值 
        flaglastLUIdx = T.BINID == max(T.BINID);           % V1 V2等价: % lastVehIdx = max(T.BINID);  %lastVehIdx = max(T{:,'BINID'});     % flaglastLUIdx = T.BINID == lastVehIdx; % lastVehIdx = ibin

        %% 哪些会变化?? fixeme 后期使用时需仔细调整
        % 某个bin内调整,其binID一定不会变化;
        % 其LID/Weight/LWH应该不会变化; SID/PID会变化; 因为OPID OSID OID等原因
        % 某个bin内调整,其BINSEQ,CoordLUBin,LU_VehType一定发生变化 （按bid和binseq排序的） 其ITEMID似乎没用 不返回了把
        % 重点是更新坐标和LU_VehType和BINSEQ，LU_VehType 这几个必定变化(PID/SID需要留意) % LU_VehType   'BINID'   BINSEQ   SID    LID    'ITEMID'    PID  ShowSEQ   'Weight'
        T{flaglastLUIdx,{'CoordLUBin','BINSEQ','LU_VehType'}} = T1{:,{'CoordLUBin','BINSEQ','LU_VehType'}};     chktLU(T);
        
    end

%% 3 ****************** 针对平铺选择 do2（甩尾平铺）/do3Array（整车平铺）数据 获取修订的 output ******************
    if ISpingpu==1 && any(flagTiledArray~=0)  %如有当允许且某个平铺成功（=1或=2）
        
        for ibin=1:length(do2Array) %do*Array 包含所有BIN
            
            if flagTiledArray(ibin)==0,continue; end  % 该ibin未平铺 继续循环
            
            if ISlastVehType==1 && flaggetSmallVeh == 1 && ibin == lastVehIdx, error('最后一车既更换车型,又需要平铺;');   end
            
            if flagTiledArray(ibin)~=0
                
                if flagTiledArray(ibin) == 1    % 该ibin整车平铺成功
                    dd = do3Array(ibin);
                elseif flagTiledArray(ibin) == 2      % 该ibin甩尾平铺成功
                    dd = do2Array(ibin);
                end
                
                % 获取dd->LU的Table格式
                T23 = getTableLU(dd);    chktLU(T23) ;
            end
            
            % 当前table所有属于该bin的逻辑号            
            flagTileLUIdx = T.BINID ==ibin;  % V1" % flagTileLUIdx = T{:,'BINID'}==ibin; (二者等价）
            
            %% CHECK 总表T 和 子表 T23等内部核验
            % 1 是否该ibin内的LU排序是严格递增; 是否T内选中的flagTileLUIdx部分LU属于该ibin且也是严格递增的
            if ~issorted(sort(T23.BINSEQ),'strictascend') || ~issorted(sort(T.BINSEQ(flagTileLUIdx,:)'),'strictascend')
                T.BINSEQ(flagTileLUIdx,:)';
                T23.BINSEQ';
                sort(T23.BINSEQ)';
                sort(T23.BINID)';
                do2Array(ibin).LU.LU_Bin;
                sort(T.BINSEQ(flagTileLUIdx,:))';
                error('1');
            end
            % 2 是否总表T内BINID内是统一binid；是否子表T23内是否同一binid；是否ibin等于T内LU对应的binid号
            if ~isscalar(unique(sort(T.BINID(flagTileLUIdx,:))))  || ~isscalar(unique(sort(T23.BINID)))  || ibin~= unique(sort(T.BINID(flagTileLUIdx,:)))
                sort(T23.BINID)'
                unique(sort(T23.BINID))
                unique(sort(T.BINID(flagTileLUIdx,:)))
                error('11');
            end
            % 3 是否T23内数量等于bin内LU个数
            if sum(flagTileLUIdx) ~= height(T23)
                error('111');
            end
            % 4 是否总表T内初始OPID与返回字表内T23的OPID相同;
            a = T.OPID(flagTileLUIdx,:)';
            b = T23.OPID';
            if  ~isequal(a,b)
                error('1');
            end
            
            %% 哪些会变化?? 讲平铺/gap调整后的bin内的LU的部分属性(坐标,顺序等)返回到总表T中 5555
            % 某个bin内调整,其binID一定不会变化; 其bin的LU_VehType一定不会变化；
            % 其LID/Weight/LWH应该不会变化; SID/PID会变化; 因为OPID OSID OID等原因 idExchange函数
            % 某个bin内调整,其BINSEQ,CoordLUBin一定发生变化 （按bid和binseq排序的）
            % 其LU_Item一定会变化 但似乎没用了 不返回了把
            % 重点是更新坐标 % LU_VehType   'BINID'   BINSEQ   SID    LID    'ITEMID'    PID  ShowSEQ   'Weight'
            
            % 作图：仅平铺的那个BIN的图
            %     plotSolutionT(T23,struct2table(structfun(@(x) x',d.Veh,'UniformOutput',false)));
            
            % fixeme 后期使用时需仔细调整 哪些变量变化的调整
            T{flagTileLUIdx,{'CoordLUBin','BINSEQ','LU_Item'}} = ... %补充增加LU_Item数据切换,虽然用途不大,但不会报chktLU错了.
                T23{:,{'CoordLUBin','BINSEQ','LU_Item'}};

            % 190421 增加
            T.LU_Bin = [T.BINID,T.BINSEQ]; % 190421 增加BINID/BINSEQ的返回LU_Bin,方便作图; 仅看BINSEQ也可
            % if parGap BINSEQ可能变化 -> 重新计算BINSEQ ID不变
            
            % 如果考虑Gap且成功替换Gap,则本bin内的LWH和Rotaed也要替换到主数据T中
            if parGap % && flagGap Gap调整是必须  18-12 此处增加对混装间隙的处理 - 基于Table格式
                T{flagTileLUIdx,{'LWH','Rotaed'}} = T23{:,{'LWH','Rotaed'}};
%                 T{flagTileLUIdx,{'BINSEQ'}} = T23{:,{'BINSEQ'}};
            end
            
            %        sortrows(T.LU_Item)'
            %        if ITEMID相同 其坐标CoordLUBin的长宽必须相同 不同ITEMID的上下重量对比
            %        对table格式LU进行chk 主要是重量
            %        chktLU(T) %仍有无法通过的可能性; 如LU_Item影响不大,建议先注释 TODO
            
        end % END OF BIN
        
    end % END OF PINGPU

%% 4 ****************** 针对获取的T 进行最后返回的output处理 ******************
    %增加 ShowSEQ(按车辆/供应商号/LID区分显示步骤) 和 tblorder -> 给刘强使用
    [T.ShowSEQ, T.tblorder] = getBBASeqTLU(T); 
    %       chktLU(T)  %      上面不通过,猜想是LU_Item未及时调整,在后期甩尾平铺后. TODO

    % NOTE : 在此之前均未改变LU的顺序，改变LU顺序为按照显示顺序ShowSEQ递增
    [~,T.ttt] = sort(T.tblorder);
    T = sortrows(T,'ttt');                               % T = sortrows(T,'BINID')  % T.LID % T.BINID % T.BINSEQ
end







%% V1 V2删除部分似乎可能有用其实无用的注释
% % function T = HBinCombine(do,flaggetSmallVeh,do1,flagTiledArray,do2Array,do3Array)
% % global ISlastVehType ISpingpu parGap
% % %% 1 ******************获取展示顺序 do数据 T=d.LU增加ShowSEQ
% % T = getTableLU(do);
% % chktLU(T)
% %
% % % 作图：原始非平铺的图
% % %     plotSolutionT(T,struct2table(structfun(@(x) x',d.Veh,'UniformOutput',false)));
% %
% % %%
% % % 行1：托盘所在车型号(必须换)    行2：托盘所在车序号(会变,不能换,换就错) 行3：托盘车内安置顺序(必须换) 行4：托盘SID供应商编号(不会变,不用变??)
% % % 行5：托盘ID型号LID(不会变,不用变?) 行6：托盘堆垛序号ITEM(会变,不能换,换就错,用途?) 行7：托盘零部件编号PID(不会变,不用变?) 增加行8: 展示顺序(必须换)
% %
% % %% 2 ****************** 针对车型变化 do1数据 获取修订的 output ******************
% % if ISlastVehType==1 && flaggetSmallVeh == 1 %如有当允许且车型替换成功
% %     T1 = getTableLU(do1);
% %     chktLU(T1)
% %
% %     % 替换T中的最后一车的部分属性 来自T1
% %     lastVehIdx = max(T{:,'BINID'});
% %     flaglastLUIdx = T{:,'BINID'}==lastVehIdx;
% %
% %        %% 哪些会变化??
% %        % 某个bin内调整,其binID一定不会变化;
% %        % 其LID/Weight/LWH应该不会变化; SID/PID会变化; 因为OPID OSID OID等原因
% %        % 某个bin内调整,其BINSEQ,CoordLUBin,LU_VehType一定发生变化 （按bid和binseq排序的） 其ITEMID似乎没用 不返回了把
% %        % 重点是更新坐标和LU_VehType和BINSEQ，LU_VehType 这几个必定变化(PID/SID需要留意) % LU_VehType   'BINID'   BINSEQ   SID    LID    'ITEMID'    PID  ShowSEQ   'Weight'
% %     T{flaglastLUIdx,{'CoordLUBin','BINSEQ','LU_VehType'}} = T1{:,{'CoordLUBin','BINSEQ','LU_VehType'}};
% %     chktLU(T)
% %     %%
% %
% %         % LU_VehType   'BINID'   BINSEQ   SID    LID    'ITEMID'    PID  ShowSEQ   'Weight'
% %                 % %     T{flaglastLUIdx,{'LU_VehType','BINSEQ','ShowSEQ'}} = ...
% %                 % %         T1{:,{'LU_VehType','BINSEQ','ShowSEQ'}};
% %                 % %     T{flaglastLUIdx,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ'}} = ...
% %                 % %         T1{:,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ'}};
% %
% % %      [TLAST_Coord,TLAST_LWH,TLAST_Seq] = getReturnBBA(do1); %% 进行返回处理
% % %     TLAST_Seq1 = finalCheck( [TLAST_Coord,TLAST_LWH,TLAST_Seq],TLUIN_LAST); %参数1：计算； 参数2：原始.
% %
% % %     %由于order改变了,此处仅对最后一个bin的索引进行修改
% % %     lastVehIdx = max(T_Seq{:,'BINID'});
% % %     flaglastLUIdx = T_Seq1{:,'BINID'}==lastVehIdx;
% % %
% % %     % 重点是更新坐标和长宽高
% % %     T_Coord{flaglastLUIdx,:} = TLAST_Coord{:,:};
% % %     T_LWH{flaglastLUIdx,:} = TLAST_LWH{:,:};
% % %     % LU_VehType   'BINID'   BINSEQ   SID    LID    'ITEMID'    PID  ShowSEQ   'Weight'
% % %     T_Seq1{flaglastLUIdx,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ'}} = ...
% % %         TLAST_Seq1{:,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ'}};
% %
% % %     T_Seq = finalCheck([T_Coord,T_LWH,T_Seq],TLUIN); %参数1：计算； 参数2：原始.
% %
% % %     output_CoordLUBin(:,flaglastLUIdx) = TLAST_Coord;
% % %     output_LU_LWH(:,flaglastLUIdx) = TLAST_LWH;
% %
% %     % 部分相关参数需要替换: 行1：托盘所在车型号 行3：托盘车内安置顺序 行8: 展示顺序
% %         %     output_LU_Seq([1,3,4,5,7],flaglastLUIdx) = output_LU_Seq2([1,3,4,5,7],:); %[1,3,4,5,7]表示仅修改这里的几行
% % %     output_LU_Seq([1,3,8],flaglastLUIdx) = TLAST_Seq([1,3,8],:); %[1,3,4,5,7,8]表示仅修改这里的几行
% % %     output_LU_Seq([1,3,4,5,7,8],flaglastLUIdx) = TLAST_Seq([1,3,4,5,7,8],:); %[1,3,4,5,7,8]表示仅修改这里的几行
% %
% %         %     i = [4,5,7]; % 这些应该是不会变的
% %         %     if sum(output_LU_Seq(i,flaglastLUIdx) ~= output_LU_Seq2(i,:) ) >0, error('不会变的变了, 错误'); end
% % end
% %
% % %% 3 ****************** 针对平铺选择 do2（甩尾平铺）/do3Array（整车平铺）数据 获取修订的 output ******************
% % if ISpingpu==1 && any(flagTiledArray~=0)  %如有当允许且某个平铺成功（=1或=2）
% %
% %     for ibin=1:length(do2Array) %do*Array 包含所有BIN
% %         if flagTiledArray(ibin)==0  % 该ibin未平铺 继续循环
% %             continue;
% %         end
% %         if ISlastVehType==1 && flaggetSmallVeh == 1 && ibin == lastVehIdx
% %             error('最后一车既更换车型,又需要平铺;');
% %         end
% %
% %         if flagTiledArray(ibin)==1    % 该ibin整车平铺成功
% %             dd = do3Array(ibin);
% %
% %         % 新增混装间隙优化     %             if parGap==1,   [dd.LU] = getMixedGap(dd.LU, dd.Veh);    end
% %
% %             % 获取LU的Table格式
% %             T23 = getTableLU(dd);   %T23 = getTableLU(do3Array(ibin));
% %             chktLU(T23) ;
% %         end
% %         if flagTiledArray(ibin)==2      % 该ibin甩尾平铺成功
% %             dd = do2Array(ibin);
% %
% %         % 新增混装间隙优化     %             if parGap==1,   [dd.LU] = getMixedGap(dd.LU, dd.Veh);    end
% %
% %             % 获取LU的Table格式
% %             T23 = getTableLU(dd);     % T23 = getTableLU(do2Array(ibin));
% %             chktLU(T23);
% %         end
% %
% %         % 当前table所有属于该bin的逻辑号
% %         flagTileLUIdx = T{:,'BINID'}==ibin;
% %
% %        %% CHECK 1 是否该ibin内的LU排序是严格递增; 是否T内选中的flagTileLUIdx部分LU属于该ibin且也是严格递增的
% %        if ~issorted(sort(T23.BINSEQ),'strictascend') || ~issorted(sort(T.BINSEQ(flagTileLUIdx,:)'),'strictascend')
% %                T.BINSEQ(flagTileLUIdx,:)'
% %                T23.BINSEQ'
% %                sort(T23.BINSEQ)'
% %                sort(T23.BINID)'
% %                do2Array(ibin).LU.LU_Bin
% %                sort(T.BINSEQ(flagTileLUIdx,:))'
% %                error('1');
% %        end
% %        % 2 是否总表T内BINID内是统一binid；是否子表T23内是否同一binid；是否ibin等于T内LU对应的binid号
% %        if ~isscalar(unique(sort(T.BINID(flagTileLUIdx,:))))  || ~isscalar(unique(sort(T23.BINID)))  || ibin~= unique(sort(T.BINID(flagTileLUIdx,:)))
% %             sort(T23.BINID)'
% %             unique(sort(T23.BINID))
% %             unique(sort(T.BINID(flagTileLUIdx,:)))
% %             error('11');
% %        end
% %        % 3 是否T23内数量等于bin内LU个数
% %        if sum(flagTileLUIdx) ~= height(T23)
% %            error('111');
% %        end
% %        % 4 是否总表T内初始OPID与返回字表内T23的OPID相同;
% %        a = T.OPID(flagTileLUIdx,:)';
% %        b = T23.OPID';
% %        if  ~isequal(a,b)
% %            error('1');
% %        end
% %
% %        %% 哪些会变化??
% %        % 某个bin内调整,其binID一定不会变化; 其bin的LU_VehType一定不会变化；
% %        % 其LID/Weight/LWH应该不会变化; SID/PID会变化; 因为OPID OSID OID等原因 idExchange函数
% %        % 某个bin内调整,其BINSEQ,CoordLUBin一定发生变化 （按bid和binseq排序的）
% %        % 其LU_Item一定会变化 但似乎没用了 不返回了把
% %        % 重点是更新坐标 % LU_VehType   'BINID'   BINSEQ   SID    LID    'ITEMID'    PID  ShowSEQ   'Weight'
% % %        sortrows(T.LU_Item)'
% % %        sortrows(T23.LU_Item)'
% % %        T.LU_Item
% % %        T{flagTileLUIdx,{'CoordLUBin','BINSEQ'}} = T23{:,{'CoordLUBin','BINSEQ'}};
% %
% %     % 作图：仅平铺的那个BIN的图
% % %     plotSolutionT(T23,struct2table(structfun(@(x) x',d.Veh,'UniformOutput',false)));
% %
% %         T{flagTileLUIdx,{'CoordLUBin','BINSEQ','LU_Item'}} = ... %补充增加LU_Item数据切换,虽然用途不大,但不会报chktLU错了.
% %             T23{:,{'CoordLUBin','BINSEQ','LU_Item'}};
% %
% %         % 如果考虑Gap且成功替换Gap,则本bin内的LWH和Rotaed也要替换到主数据T中
% %         if parGap % && flagGap Gap调整是必须
% %             T{flagTileLUIdx,{'LWH','Rotaed'}} = ... %补充增加LU_Item数据切换,虽然用途不大,但不会报chktLU错了.
% %                 T23{:,{'LWH','Rotaed'}};
% %         end
% %
% % %        sortrows(T.LU_Item)'
% % %        chktLU(T) %仍有无法通过的可能性; 如LU_Item影响不大,建议先注释 TODO
% %
% %             %% 下面是错的
% %             %        T{flagTileLUIdx,{'LU_VehType','BINSEQ','ShowSEQ'}} = ...
% % %            T23{:,{'LU_VehType','BINSEQ','ShowSEQ'}};
% % %        T{flagTileLUIdx,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ'}} = ...
% % %            T23{:,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ'}};
% % %        T{flagTileLUIdx,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ','CoordLUBin','LWH'}} = ...
% % %            T23{:,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ','CoordLUBin','LWH'}};
% %
% %
% %             %     T{flagTileLUIdx,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ','CoordLUBin','LWH','BINID','ITEMID','Weight'}} = ...
% %             %         T2{:,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ','CoordLUBin','LWH','BINID','ITEMID','Weight'}};
% %
% %             % output_LU_Seq=T{:,{'','BINID',',''ITEMID','Weight'}}'
% %
% %             % % %     [TPP_Coord,TPP_LWH,TPP_Seq]= getReturnBBA(do2Array(ibin)); %% 进行返回处理
% %             % % %     %     [output_CoordLUBin2,output_LU_LWH2,output_LU_Seq2]= getReturnBBA1(do2Array(ibin)); %% 进行返回处理
% %             % % %     TPP_Seq1 = finalCheck([TPP_Coord,TPP_LWH,TPP_Seq],TLUIN_PP1);
% %             % % % %     if  sum(flagTiled(:)~=0)>1
% %             % % % %         sum(flagTiled(:)~=0)
% %             % % % %         1
% %             % % % %     end
% %             % % % %     if flagTiled(ibin)==1
% %             % % % %         if sum(flagTiled(:)==1)>1
% %             % % % %             sum(flagTiled(:)==1)
% %             % % % %             1
% %             % % % %         end
% %             % % % %         TPP_Seq = finalCheck([TPP_Coord,TPP_LWH,TPP_Seq],TLUIN_PP1); %参数1：计算； 参数2：原始.
% %             % % % %     end
% %             % % % %     if flagTiled(ibin)==2
% %             % % % %         if  sum(flagTiled(:)==2)>1
% %             % % % %             sum(flagTiled(:)==2)
% %             % % % %             1
% %             % % % %         end
% %             % % % %         TPP_Seq = finalCheck([TPP_Coord,TPP_LWH,TPP_Seq],TLUIN_PP2); %参数1：计算； 参数2：原始.
% %             % % % %     end
% %             % % %
% %             % % %     % 找出平铺ibin内所有的托盘逻辑值
% %             % % %     flagTileLUIdx = T_Seq1{:,'BINID'} == ibin;
% %             % % %
% %             % % %     % 重点是更新坐标和长宽高
% %             % % %     T_Coord{flagTileLUIdx,:} = TPP_Coord{:,:};
% %             % % %     T_LWH{flagTileLUIdx,:} = TPP_LWH{:,:};
% %             % % %     % LU_VehType   'BINID'   BINSEQ   SID    LID    'ITEMID'    PID  ShowSEQ   'Weight'
% %             % % %     T_Seq1{flagTileLUIdx,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ'}} = ...
% %             % % %         TPP_Seq1{:,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ'}};
% %
% %
% %             %         T_Seq = finalCheck([T_Coord,T_LWH,T_Seq],TLUIN); %参数1：计算； 参数2：原始.
% %
% %             % 部分相关参数需要替换: 行1：托盘所在车型号 行3：托盘车内安置顺序 行8: 展示顺序
% %             %     output_LU_Seq([1,3,8],flagTileLUIdx) = output_LU_Seq2([1,3,8],:);
% %             %     output_LU_Seq([1,3,4,5,7,8],flagTileLUIdx) = output_LU_Seq2([1,3,4,5,7,8],:); %[1,3,4,5,7,8]表示仅修改这里的几行
% %
% %             %     i = [4,5,7]; % 这些应该是不会变的
% %             %     if sum(output_LU_Seq(i,flagTileLUIdx) ~= output_LU_Seq3(i,:) ) >0,
% %             %          output_LU_Seq(i,flagTileLUIdx) ~= output_LU_Seq3(i,:) ;         output_LU_Seq(i,flagTileLUIdx);         output_LU_Seq3(i,:) ;
% %             %         warning('不会变的变了, 错误'); end
% %     end
% % end
% %
% % %% 18-12 此处增加对混装间隙的处理 - 基于Table格式
% %
% % % ****************** 针对车型选择 获取修订的 output ******************
% %
% % T = getShowSeq(T); %增加ShowSEQ tblorder
% %
% % %        chktLU(T)  %      上面不通过,猜想是LU_Item未及时调整,在后期甩尾平铺后. TODO
% %
% %
% % %% table转为结构体后判断是否上轻下重
% % % lu = table2struct(T,'ToScalar',true)
% % % lu = (structfun(@(x) x',lu,'UniformOutput',false));
% %
% %     %% 依据T内最终LU_LUinBin坐标数据判断
% % %     T2 = sortrows(T,{'BINID','BINSEQ'},{'ascend','ascend'})
% %     %     x=T(T.BINID==2&T.ID==1,{'ID','LID','PID','H','Weight','CoordLUBin','BINSEQ','ShowSEQ','ITEMID','ITEMSEQ'})
% % %     x=T2(:,{'ID','LID','PID','H','Weight','X','Y','Z','ITEMID','ITEMSEQ','BINID','BINSEQ','ShowSEQ'})
% %     % if ITEMID相同 其坐标CoordLUBin的长宽必须相同 不同ITEMID的上下重量对比
% %     % 对table格式LU进行chk 主要是重量
% %
% %
% %
% % % NOTE : 在此之前均未改变LU的顺序
% % % 返回BBA数组格式给JAR
% % [~,T.ttt] = sort(T.tblorder);
% % T = sortrows(T,'ttt');
% % % T = sortrows(T,'BINID')
% % % T.LID
% % % T.BINID
% % % T.BINSEQ
% % end
% %
% %
% % %% 注释后
% %
% %
% % %% USELESS
% % % [T_Coord,T_LWH,T_Seq] = getReturnBBA(daBest(bestOne)); %如有多个,返回第一个最优解
% % % T_Seq.tblorder
% % % T_Seq1 = T_Seq;
% %
% % % T.SID = oD.SID;
% % % T.PID = oD.PID;
% % % T.LID = oD.LID;
% % % Tseq = T(:,{'LU_VehType','BINID','BINSEQ','SID','LID','ITEMID','PID','ShowSEQ','Weight','tblorder'});
% %
% % % T_Seq1 = finalCheck([T_Coord,T_LWH,T_Seq],TLUIN); %参数1：计算； 参数2：原始.
% %
% % %  [output_CoordLUBin,output_LU_LWH,output_LU_Seq]= getReturnBBA1(daBest(bestOne)); %% 进行返回处理
