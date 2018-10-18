%% V3 : 采用cell方式 且分开ITEM和STRIP两类 重点排除order被重复赋值的情况
% 先按SID区分, 再给出priority为SID内部顺序
% 重点：27-29行tmpM的顺序

% IDorder = getOrderofLID(SIDorder, Strip.isSingleItem, Strip.isAllPured, Strip.nbItem, Strip.isHeightFull,...
%                                         Strip.isMixed, Strip.LID, Strip.LW(1:2,:), Strip.loadingrateLimit, Strip.loadingrate);
                                    
function allPriority = getOrderofLID(SIDord,Ssingle, Spured,SnbItem,SnbLU,SnbLULID, ...
    SHeightfullfull, SisMixed, SLID,SLW,SLoadingRateLimit,SLoadingrate)
% ID is cell type % ID = [2 3 NaN; 2 1 3; 1 3 NaN; 2 1 NaN]';
if ~iscell(SLID),error('输入非cell');end

allPriority = zeros(1,size(SLID,2));
uniOrd = unique(SIDord);

% 不同SID获取不同的顺序, 从SID的序号1开始
for i=1:length(uniOrd)    
    idxSID = SIDord==uniOrd(i);
    
    % UNUSED 获取本SID内的STRIP对应的值tLW;tLL;tL;tID;celltID.
    tLW = SLW(:,idxSID);   tLL = SLoadingRateLimit(:,idxSID);
    tID = SLID(:,idxSID);    celltID = tID;
    tPured = Spured(:,idxSID);
    tSingle = Ssingle(:,idxSID);   %作用未知
    
    % FINALLY USED 优先车头摆放的顺序
    tnbItem = SnbItem(:,idxSID);
    tnbLU = SnbLU(:,idxSID);
    tnbLULID = SnbLULID(:,idxSID);
    
    
    tMixed = SisMixed(:,idxSID);
    tFull = SHeightfullfull(:,idxSID);
    tLR = SLoadingrate(:,idxSID);        % max(tL)  %  min(tL)  % mean(tL)
    
    % 555 priority: 标量值, 永远递增,最小值=1,最大值为LU个数
    priority = 1;
    SIDpriority = zeros(1,size(tID,2)); % SIDpriority: 数组, 每个SID下的priority值
    szRow = cellfun(@(x)size(x,1), tID);
    if isscalar(szRow) %如果该STRIP的只有1STRIP，且为单纯型STRIP，赋值为1
        SIDpriority = 1;
    else
        [tID, ~] = padcat(tID{:});   if iscolumn(tID), tID = tID'; end
        
        % 1 给定不考虑相邻的Strip摆放顺序 ->             
        % 改用和Item2Strip类似的以高度优先的排序方式, 其次为LoadingRate
        %         tmpM = [tLL;tL;tLW;];
        % %         [~,order] = sortrows(tmpM',[1,2,4],{'descend','descend','descend'}); %order is index vector 返回原索引值 如第一个为3，表示原array中第3个目前是第1个
        %          tmpM = [tLL;tL;tLW;];[~,order] = sortrows(tmpM',[4,1],{'descend','descend'}); %order is index vector 返回原索引值 如第一个为3，表示原array中第3个目前是第1个

        % 1 单纯 > 满层 > 数量 > LoadingRate
        %         tmpM = [tLID; tMixed; tFull; tLL;tL;tLW;];  [~,order] = sortrows(tmpM',[2,3,1,5],{'ascend','descend','descend','descend'}); 
        % 2 数量 > 单纯 > 满层 > LoadingRate        
        tmpM = [tnbItem; tMixed; tFull; tLR; ];  % tLL; tLW; tID; tPured; tSingle
        [~,order] = sortrows(tmpM',[1,2,3,4],{'descend','ascend','descend','descend'}); 
        tmpM = [tnbLU; tnbItem; tMixed; tFull; tLR; ];  % tLL; tLW; tID; tPured; tSingle
        [~,order] = sortrows(tmpM',[1,2,3,4,5],{'descend', 'descend','ascend','descend','descend'}); 
        tmpM = [tnbLU; tnbLULID; tnbItem; tMixed; tFull; tLR; ];  % tLL; tLW; tID; tPured; tSingle
        [~,order] = sortrows(tmpM',[1,2,3,4,5,6],{'descend', 'descend', 'descend','ascend','descend','descend'});         
        
%         tmpM = [tPured;tnbLU; tnbItem; tMixed; tFull; tLR; ];  % tLL; tLW; tID; tPured; tSingle
%         [~,order] = sortrows(tmpM',[1,2,3,4,5,6],{'descend', 'descend','descend','ascend','descend','descend'}); 
        
        % 3 全纯 > 数量 > 单纯 > 满层 > LoadingRate
        %         tmpM = [tPured; tLID; tMixed; tFull; tLL;tL;tLW;];  [~,order] = sortrows(tmpM',[1,2,3,4,6],{'descend','descend','ascend','descend','descend'}); 
        %         tmpM = [tLID; tMixed; tFull; tLL;tL;tLW;];        [~,order] = sortrows(tmpM',[2,3],{'ascend','descend'});
        
        if ~isrow(order), order=order'; end

        % 2 基于1的摆放顺序, 给定最终STRIP顺序到torder
        while any(SIDpriority==0) %当有SID下面的priority未给全时, 不能跳出循环
            % SIDorder中非0且order中最大的那个作为首选STRIP
            % 2.1 找出无相邻前提下,第一个o及对应的order(o)
            [~,o]=find(SIDpriority(order) == 0,1,'first');
            SIDpriority(order(o)) = priority;  
            priority=priority+1;
            
            % 2.2 找出给order(o)位置tLID对应的相邻Strip.
            tnbItem = celltID{:,order(o)};
            % NOTE 首次tLID不应该出现混合型STRIP, 但如按高度排序, 是可能出现的
            [SIDpriority,priority] = getAdjPriority(priority,order,SIDpriority,tID,tnbItem);
        end
    end
    allPriority(idxSID) = SIDpriority;

end

% 防错语句：
if any(allPriority==0), error('allPriority 存在未赋值列'); end

end

%% 局部函数 %%
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
                [a,~] = sortrows(tmp,[2],{'ascend'});
                for i1=1:length(nbcol1)
                    SIDorder(a(i1,1)) = priority;
                    priority=priority+1;
                end
            end
            
            % 寻找tID中的列: 该列有2个以上LID; 该列对应的SIDorder为0；该列包含（不等于）tLID
            [~,nbcol2] = find(tID(:,:)==tLID & sum(~isnan( tID(:,:) ),1) ~=1 &  SIDorder==0);
            if ~isempty(nbcol2)
                if ~isscalar(nbcol2),  error('nbcol2包含不只一个数: 存在多个相同混合的STRIP');   end
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