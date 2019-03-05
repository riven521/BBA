%% V5 : 增加可读性 V4 含EID且数据结构基本重写
% 先按SID区分, 再按EID区分, 再给出priority为SID内部顺序
% 重点：63-65 行 tmpM的顺序

function StripPriorityArray = getOrderofLID(SIDord, EIDord, Strip)

%% 目的 -> 获取allPriority, 由多个perPriority组合构成
% 分别对相同SID/EID下的STRIP赋予priority值, 每个下面的priority不重复,都从1开始
% 不同SIDord获取不同的顺序, 从SID/EID的第一个开始

% 初始化
StripPriorityArray = zeros(1,length(SIDord));        

uniSIDord = unique(SIDord);     if ~issorted(uniSIDord,'strictascend'), error('非严格递增'); end

for i=1:length(uniSIDord)
    
    fsid = SIDord==uniSIDord(i);
    
    uniEIDord = unique(EIDord(fsid));  if ~issorted(uniEIDord,'strictascend'), error('非严格递增'); end %本SID下的EID, 也必须是严格递增吧
    
    for j=1:length(uniEIDord)
        
        feid = EIDord==uniEIDord(j);
        
        f = fsid & feid;   if ~any(f), error('此SID和EID下无任何STRIP'); end %后续或许可以改为warning,因为可能不同SID下的EID不同
        
        priority = 1;  % 555 priority: 标量值, 永远递增,最小值=1,最大值为LU个数
        
        perPriority = zeros(1,sum(f));     % perPriority: 数组, 每个SID/EID下的f个正的priority值,需要被赋值        
        
        %% TO DEL 排序用的某些值 改为传递STRIP进来 （如下似乎用不到，可留做参考）
        % UNUSED 获取本SID和EID内的STRIP对应的值tLW;tLL;tL;tID;celltID.
        %         tLW = SLW(:,f);   tLL = SLoadingRateLimit(:,f);  tPured = Spured(:,f);   tSingle = Ssingle(:,f);   %作用未知        
        % FINALLY USED 优先车头摆放的顺序
        %         tnbLU = SnbLU(:,f);     tnbLULID = SnbLULID(:,f); tnbItem = SnbItem(:,f);
        %         tMixed = SisMixed(:,f);     tFull = SHeightfullfull(:,f);   tLR = SLoadingrate(:,f);        % max(tL)  %  min(tL)  % mean(tL)

        %% 初始化
        tID = Strip.LID(:,f);    
        celltID = tID;                   %celltID ： 同一SID/EID下的托盘类型  %Strip.LID实际为LU的ID,拼载判断
        
        if length(tID) == 1    % 如果只有1层  == 等价于  if sum(f)==1
            mID = tID{1};
            isID = all(mID,2);
        else
            [mID, isID] = padcat(tID{:});  
        end
        
         if iscolumn(mID), mID = mID'; end    % mID：每个 f (条带) 对应的ID类型号; isID: 每个 f (条带) 对应的ID类型号是否存在(逻辑判断)
         if iscolumn(isID), isID = isID'; end
        
        % nbIDperStrip = sum(isID,1);
        
        %% 1 获取perPriority: 相同SID/EID下, 所有STRIP的优先顺序; 继而从相邻摆放角度出发, 细化每个内部STRIP的priority
        
        if sum(f)==1 % 1.1 如果该STRIP的只有1个STRIP，且为单纯型STRIP或混合型，均赋值权重为1.
            
            perPriority = 1;     if ~isscalar(isID),    warning('是混合型strip.  也赋值为1, 应该问题不大,是否会出现这个情况要考虑清楚;');   end
            
        else  % 1.2 如果该STRIP的含有多个STRIP, 必须对STRIP进行排序获取order
            
            tmpM = [Strip.nbLU(:,f);  Strip.nbLULID(:,f);  Strip.nbItem(:,f); Strip.isMixed(:,f);  Strip.isHeightFull(:,f);  Strip.loadingrate(:,f); ];
            % 下面的计算有问题，不采用 todo 后期修复. 
%              tmpM = [Strip.nLUID(:,f);  Strip.nLULID(:,f);  Strip.nbItem(:,f); Strip.isMixed(:,f);  Strip.isHeightFull(:,f);  Strip.loadingrate(:,f); ];
          % 非混合优先；托盘类型数递减；托盘LID递减；堆垛数递减（相同托盘数，体积大优先）。 
          [~,order] = sortrows(tmpM',[4,1,2,3,5,6],{'ascend','descend', 'descend','descend','descend','descend'});          if ~isrow(order), order=order'; end
%             [~,order] = sortrows(tmpM',[1,2,3,4,5,6],{'descend', 'descend', 'descend','ascend','descend','descend'});          if ~isrow(order), order=order'; end
              
%             [~,order] = sortrows(tmpM',[4,1,2,5,6],{'ascend','descend', 'descend','descend','descend'});          if ~isrow(order), order=order'; end          

            %% 2 基于1的摆放顺序, 给定最终STRIP的priority到torder
            while any(perPriority==0) %当有SID下面的priority未给全时, 不能跳出循环(未全给原因是找不到相邻的strip了,再基于order给下一个首选strip)
                
                % SIDorder中非0且order中最大的那个作为首选STRIP
                
                % 2.1 找出无相邻( 混合或非混合 )前提下,第一个o及对应的order(o)
                [~,o]=find(perPriority(order) == 0,1,'first');
                
                perPriority(order(o)) = priority;
                
                priority=priority+1;
                
                % 2.2 找出给order(o)位置tLID对应的相邻Strip.
                tnbItem = celltID{:,order(o)};
                
                % NOTE 首次tLID不应该出现混合型STRIP, 但如按高度排序, 是可能出现的; 目前改为按混合排序，非混在前，应该不会出现的了。
                if Strip.isMixed(order(o))
                    warning('首次不希望出现混合型STRIP,若出现不容易找接下来的相邻STRIP');   % 之前为warning
                end
                
                [perPriority,priority] = getAdjPriority(priority,order,perPriority,mID,tnbItem);
                
            end
        end
        
        StripPriorityArray(f) = perPriority;
        
    end % EOF EID 
    
end % EOF SID

% 排错CODE
 if ~all(StripPriorityArray), error('allPriority 存在未赋值列'); end

end

%% 局部函数 %%
%% V2 : 
%函数1: 寻找tLID对应的相邻STRIP,按照一定逻辑
function [SIDorder,priority] = getAdjPriority(priority,order,SIDorder,tID,tLID)
     
      % 如tLID有多个，循环寻找
     for p=1:length(tLID)
        atPID = tLID(p);
        getorderofatLID(atPID);
     end
    
    % 针对单个tPID, 寻找相邻STRIP
    function getorderofatLID(tLID)
        if isscalar(tLID)  %如果只有1个
            %找出是否有对应混合ID
            % 寻找tID中的列: 该列仅有1个LID; 该列对应的SIDorder为0；该列包含（等于）tLID
            [~,nbcol1] = find(tID(:,:)==tLID & sum(~isnan( tID(:,:) ),1) ==1 &  SIDorder==0);
            if isrow(nbcol1),  nbcol1=nbcol1'; end
            % 如有多个nbcol1, 按order的顺序对SIDorder的priority赋值
            if ~isempty(nbcol1)
                [~,b]= ismember(nbcol1,order);
                tmp= [nbcol1,b];
                [a,~] = sortrows(tmp, 2 ,{'ascend'});
                for i1=1:length(nbcol1)
                    SIDorder(a(i1,1)) = priority;
                    priority=priority+1;
                end
            end
            
            % 寻找tID中的列: 该列有2个以上LID; 该列对应的SIDorder为0；该列包含（不等于）tLID
            [~,nbcol2] = find(tID(:,:)==tLID & sum(~isnan( tID(:,:) ),1) ~=1 &  SIDorder==0);
            if ~isempty(nbcol2)
                if ~isscalar(nbcol2)
%                     tID(:,:)
                    % tID(:,:)=tLID
                    warning('nbcol2包含不只一个数: 存在多个相同混合的STRIP; -> 同一个单Strip, 找其混合的其它Strip, 但发现有2个包含该strip的混合strip,无从选择哪个作为相邻strip的矛盾 ');   
%                     error('nbcol2包含不只一个数: 存在多个相同混合的STRIP; -> 同一个单Strip, 找其混合的其它Strip, 但发现有2个包含该strip的混合strip,无从选择哪个作为相邻strip的矛盾 ');   
                end
                SIDorder(nbcol2) = priority;
                priority=priority+1;
            end
            % 寻找nbcol2中除tLID以外的数,并递归找其相邻STRIP,直到无相邻STRIP
            hyLID = tID(:,nbcol2);
            hyLID(hyLID==tLID) =[];
            if ~isempty(hyLID)
                [SIDorder,priority] = getAdjPriority(priority,order,SIDorder,tID,hyLID);
            end
        else
            error('tLID包含多个值');
        end
    end
end


%% V1 函数1: 寻找tLID对应的相邻STRIP,按照一定逻辑
% % function [SIDorder,priority] = getAdjPriority(priority,order,SIDorder,tID,tLID)
% %      % 如tLID有多个，循环寻找
% %      for p=1:length(tLID)
% %         atPID = tLID(p);
% %         getorderofatLID(atPID);
% %      end
% %     
% %     % 针对单个tPID, 寻找相邻STRIP
% %     function getorderofatLID(tLID)
% %         if isscalar(tLID)  %如果只有1个
% %             %找出是否有对应混合ID
% %             % 寻找tID中的列: 该列仅有1个LID; 该列对应的SIDorder为0；该列包含（等于）tLID
% %             [~,nbcol1] = find(tID(:,:)==tLID & sum(~isnan( tID(:,:) ),1) ==1 &  SIDorder==0);
% %             if isrow(nbcol1),  nbcol1=nbcol1'; end
% %             % 如有多个nbcol1, 按order的顺序对SIDorder的priority赋值
% %             if ~isempty(nbcol1)
% %                 [~,b]= ismember(nbcol1,order);
% %                 tmp= [nbcol1,b];
% %                 [a,~] = sortrows(tmp,[2],{'ascend'});
% %                 for i1=1:length(nbcol1)
% %                     SIDorder(a(i1,1)) = priority;
% %                     priority=priority+1;
% %                 end
% %             end
% %             
% %             % 寻找tID中的列: 该列有2个以上LID; 该列对应的SIDorder为0；该列包含（不等于）tLID
% %             [~,nbcol2] = find(tID(:,:)==tLID & sum(~isnan( tID(:,:) ),1) ~=1 &  SIDorder==0);
% %             if ~isempty(nbcol2)
% %                 if ~isscalar(nbcol2),  
% %                     tID(:,:)
% %                     tID(:,:)=tLID
% %                     error('nbcol2包含不只一个数: 存在多个相同混合的STRIP');   
% %                 end
% %                 SIDorder(nbcol2) = priority;
% %                 priority=priority+1;
% %             end
% %             % 寻找nbcol2中除tLID以外的数,并递归找其相邻STRIP,直到无相邻STRIP
% %             hyLID = tID(:,nbcol2);
% %             hyLID(hyLID==tLID) =[];
% %             if ~isempty(hyLID)
% %                 [SIDorder,priority] = getAdjPriority(priority,order,SIDorder,tID,hyLID);
% %             end
% %         else
% %             error('tLID包含多个值');
% %         end
% %     end
% % end


%% %% V3 : 采用cell方式 且分开ITEM和STRIP两类 重点排除order被重复赋值的情况
% % % 先按SID区分, 再给出priority为SID内部顺序
% % % 重点：27-29行tmpM的顺序
% % 
% % % IDorder = getOrderofLID(SIDorder, Strip.isSingleItem, Strip.isAllPured, Strip.nbItem, Strip.isHeightFull,...
% % %                                         Strip.isMixed, Strip.LID, Strip.LW(1:2,:), Strip.loadingrateLimit, Strip.loadingrate);
% %                                     
% % function allPriority = getOrderofLID(SIDord,Ssingle, Spured,SnbItem,SnbLU,SnbLULID, ...
% %     SHeightfullfull, SisMixed, SLID,SLW,SLoadingRateLimit,SLoadingrate)
% % % ID is cell type % ID = [2 3 NaN; 2 1 3; 1 3 NaN; 2 1 NaN]';
% % if ~iscell(SLID),error('输入非cell');end
% % 
% % allPriority = zeros(1,size(SLID,2));
% % uniOrd = unique(SIDord);
% % 
% % % 不同SID获取不同的顺序, 从SID的序号1开始
% % for i=1:length(uniOrd)    
% %     idxSID = SIDord==uniOrd(i);
% %     
% %     % UNUSED 获取本SID内的STRIP对应的值tLW;tLL;tL;tID;celltID.
% %     tLW = SLW(:,idxSID);   tLL = SLoadingRateLimit(:,idxSID);
% %     tID = SLID(:,idxSID);    celltID = tID;
% %     tPured = Spured(:,idxSID);
% %     tSingle = Ssingle(:,idxSID);   %作用未知
% %     
% %     % FINALLY USED 优先车头摆放的顺序
% %     tnbItem = SnbItem(:,idxSID);
% %     tnbLU = SnbLU(:,idxSID);
% %     tnbLULID = SnbLULID(:,idxSID);
% %     
% %     
% %     tMixed = SisMixed(:,idxSID);
% %     tFull = SHeightfullfull(:,idxSID);
% %     tLR = SLoadingrate(:,idxSID);        % max(tL)  %  min(tL)  % mean(tL)
% %     
% %     % 555 priority: 标量值, 永远递增,最小值=1,最大值为LU个数
% %     priority = 1;
% %     SIDpriority = zeros(1,size(tID,2)); % SIDpriority: 数组, 每个SID下的priority值
% %     szRow = cellfun(@(x)size(x,1), tID);
% %     if isscalar(szRow) %如果该STRIP的只有1STRIP，且为单纯型STRIP，赋值为1
% %         SIDpriority = 1;
% %     else
% %         [tID, ~] = padcat(tID{:});   if iscolumn(tID), tID = tID'; end
% %         
% %         % 1 给定不考虑相邻的Strip摆放顺序 ->             
% %         % 改用和Item2Strip类似的以高度优先的排序方式, 其次为LoadingRate
% %         %         tmpM = [tLL;tL;tLW;];
% %         % %         [~,order] = sortrows(tmpM',[1,2,4],{'descend','descend','descend'}); %order is index vector 返回原索引值 如第一个为3，表示原array中第3个目前是第1个
% %         %          tmpM = [tLL;tL;tLW;];[~,order] = sortrows(tmpM',[4,1],{'descend','descend'}); %order is index vector 返回原索引值 如第一个为3，表示原array中第3个目前是第1个
% % 
% %         % 1 单纯 > 满层 > 数量 > LoadingRate
% %         %         tmpM = [tLID; tMixed; tFull; tLL;tL;tLW;];  [~,order] = sortrows(tmpM',[2,3,1,5],{'ascend','descend','descend','descend'}); 
% %         % 2 数量 > 单纯 > 满层 > LoadingRate        
% %         tmpM = [tnbItem; tMixed; tFull; tLR; ];  % tLL; tLW; tID; tPured; tSingle
% %         [~,order] = sortrows(tmpM',[1,2,3,4],{'descend','ascend','descend','descend'}); 
% %         tmpM = [tnbLU; tnbItem; tMixed; tFull; tLR; ];  % tLL; tLW; tID; tPured; tSingle
% %         [~,order] = sortrows(tmpM',[1,2,3,4,5],{'descend', 'descend','ascend','descend','descend'}); 
% %         tmpM = [tnbLU; tnbLULID; tnbItem; tMixed; tFull; tLR; ];  % tLL; tLW; tID; tPured; tSingle
% %         [~,order] = sortrows(tmpM',[1,2,3,4,5,6],{'descend', 'descend', 'descend','ascend','descend','descend'});         
% %         
% % %         tmpM = [tPured;tnbLU; tnbItem; tMixed; tFull; tLR; ];  % tLL; tLW; tID; tPured; tSingle
% % %         [~,order] = sortrows(tmpM',[1,2,3,4,5,6],{'descend', 'descend','descend','ascend','descend','descend'}); 
% %         
% %         % 3 全纯 > 数量 > 单纯 > 满层 > LoadingRate
% %         %         tmpM = [tPured; tLID; tMixed; tFull; tLL;tL;tLW;];  [~,order] = sortrows(tmpM',[1,2,3,4,6],{'descend','descend','ascend','descend','descend'}); 
% %         %         tmpM = [tLID; tMixed; tFull; tLL;tL;tLW;];        [~,order] = sortrows(tmpM',[2,3],{'ascend','descend'});
% %         
% %         if ~isrow(order), order=order'; end
% % 
% %         % 2 基于1的摆放顺序, 给定最终STRIP顺序到torder
% %         while any(SIDpriority==0) %当有SID下面的priority未给全时, 不能跳出循环
% %             % SIDorder中非0且order中最大的那个作为首选STRIP
% %             % 2.1 找出无相邻前提下,第一个o及对应的order(o)
% %             [~,o]=find(SIDpriority(order) == 0,1,'first');
% %             SIDpriority(order(o)) = priority;  
% %             priority=priority+1;
% %             
% %             % 2.2 找出给order(o)位置tLID对应的相邻Strip.
% %             tnbItem = celltID{:,order(o)};
% %             % NOTE 首次tLID不应该出现混合型STRIP, 但如按高度排序, 是可能出现的
% %             [SIDpriority,priority] = getAdjPriority(priority,order,SIDpriority,tID,tnbItem);
% %         end
% %     end
% %     allPriority(idxSID) = SIDpriority;
% % 
% % end
% % 
% % % 防错语句：
% % if any(allPriority==0), error('allPriority 存在未赋值列'); end
% % 
% % end



%% V4 : 含EID且数据结构基本重写
% % % 先按SID区分, 再按EID区分, 再给出priority为SID内部顺序
% % % 重点：27-29行tmpM的顺序                                    
% % % function allPriority = getOrderofLID(SIDord, EIDord, Ssingle, Spured,SnbItem,SnbLU,SnbLULID, ...
% % %     SHeightfullfull, SisMixed, SLID,SLW,SLoadingRateLimit,SLoadingrate)
% % function allPriority = getOrderofLID(SIDord, EIDord, Strip)
% % 
% % %% 目的 -> 获取allPriority, 由多个perPriority组合构成
% % % 分别对相同SID/EID下的STRIP赋予priority值, 每个下面的priority不重复,都从1开始
% % % 不同SIDord获取不同的顺序, 从SID/EID的第一个开始
% % allPriority = zeros(1,length(SIDord));        
% % uniSIDord = unique(SIDord); if ~issorted(uniSIDord,'strictascend'), error('非严格递增'); end
% % for i=1:length(uniSIDord)
% %     fsid = SIDord==uniSIDord(i);
% %     uniEIDord = unique(EIDord(fsid));  if ~issorted(uniEIDord,'strictascend'), error('非严格递增'); end %本SID下的EID, 也必须是严格递增吧
% %     for j=1:length(uniEIDord)
% %         feid = EIDord==uniEIDord(j);
% %         f = fsid & feid;   if ~any(f), error('此SID和EID下无任何STRIP'); end %后续或许可以改为warning,因为可能不同SID下的EID不同
% %         
% %         priority = 1;  % 555 priority: 标量值, 永远递增,最小值=1,最大值为LU个数
% %         perPriority = zeros(1,sum(f));     % perPriority: 数组, 每个SID/EID下的f个正的priority值,需要被赋值
% %         
% %         
% %         %% TO DEL 排序用的某些值 改为传递STRIP进来
% %         % UNUSED 获取本SID和EID内的STRIP对应的值tLW;tLL;tL;tID;celltID.
% % %         tLW = SLW(:,f);   tLL = SLoadingRateLimit(:,f);  tPured = Spured(:,f);   tSingle = Ssingle(:,f);   %作用未知        
% %         % FINALLY USED 优先车头摆放的顺序
% % %         tnbLU = SnbLU(:,f);     tnbLULID = SnbLULID(:,f); tnbItem = SnbItem(:,f);
% % %         tMixed = SisMixed(:,f);     tFull = SHeightfullfull(:,f);   tLR = SLoadingrate(:,f);        % max(tL)  %  min(tL)  % mean(tL)
% % 
% %         %% 0 初始化     tID = SLID(:,f);    celltID = tID;
% %         tID = Strip.LID(:,f);    celltID = tID;  %Strip.LID实际为LU的ID,拼载判断
% %         if length(tID) == 1
% %             mID = tID{1};
% %             isID = all(mID,2);
% %         else
% %             [mID, isID] = padcat(tID{:});  
% %         end
% %          if iscolumn(mID), mID = mID'; end;  if iscolumn(isID), isID = isID'; end
% %         
% %         % nbIDperStrip = sum(isID,1);
% %         
% %         %% 1 获取perPriority: 相同SID/EID下, 所有STRIP的优先顺序; 继而从相邻摆放角度出发, 细化每个内部STRIP的priority
% %         % 1.1 如果该STRIP的只有1个STRIP，且为单纯型STRIP或混合型，均赋值权重为1.
% %         if sum(f)==1
% %             perPriority = 1;     if ~isscalar(isID), warning('是混合型strip也赋值为1, 应该问题不大,是否会出现这个情况要考虑清楚;'); end
% %             % 1.2 如果该STRIP的含有多个STRIP, 必须对STRIP进行排序获取order
% %         else
% %             %             tmpM = [tnbLU; tnbLULID; tnbItem; tMixed; tFull; tLR; ];  % tLL; tLW; tID; tPured; tSingle
% % 
% %             tmpM = [Strip.nbLU(:,f);  Strip.nbLULID(:,f);  Strip.nbItem(:,f); Strip.isMixed(:,f);  Strip.isHeightFull(:,f);  Strip.loadingrate(:,f); ];
% %             [~,order] = sortrows(tmpM',[4,1,2,3,5,6],{'ascend','descend', 'descend', 'descend','descend','descend'});          if ~isrow(order), order=order'; end
% % %             [~,order] = sortrows(tmpM',[1,2,3,4,5,6],{'descend', 'descend', 'descend','ascend','descend','descend'});          if ~isrow(order), order=order'; end
% % %             [~,order] = sortrows(tmpM',[4],{'descend'});          if ~isrow(order), order=order'; end
% %             
% %             %% 2 基于1的摆放顺序, 给定最终STRIP的priority到torder
% %             while any(perPriority==0) %当有SID下面的priority未给全时, 不能跳出循环(未全给原因是找不到相邻的strip了,再基于order给下一个首选strip)
% %                 % SIDorder中非0且order中最大的那个作为首选STRIP
% %                 % 2.1 找出无相邻前提下,第一个o及对应的order(o)
% %                 [~,o]=find(perPriority(order) == 0,1,'first');
% %                 perPriority(order(o)) = priority;
% %                 priority=priority+1;
% %                 
% %                 % 2.2 找出给order(o)位置tLID对应的相邻Strip.
% %                 tnbItem = celltID{:,order(o)};
% %                 % NOTE 首次tLID不应该出现混合型STRIP, 但如按高度排序, 是可能出现的
% %                 if Strip.isMixed(order(o))
% %                     warning('首次不希望出现混合型STRIP,若出现不容易找接下来的相邻STRIP'); 
% %                 end
% %                 [perPriority,priority] = getAdjPriority(priority,order,perPriority,mID,tnbItem);
% %             end
% %         end
% %         allPriority(f) = perPriority;
% %     end
% % end
% % 
% % % 排错CODE
% %  if ~all(allPriority), error('allPriority 存在未赋值列'); end
% % 
% % end
