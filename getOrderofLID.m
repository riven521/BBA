%% V3 : 采用cell方式 且分开ITEM和STRIP两类 重点排除order被重复赋值的情况
% 先按SID区分, 再给出priority为SID内部顺序
% 重点：27-29行tmpM的顺序
function allPriority = getOrderofLID(SIDord,ID,LW,LoadingRateLimit,Loadingrate)
% ID is cell type % ID = [2 3 NaN; 2 1 3; 1 3 NaN; 2 1 NaN]';

if ~iscell(ID),error('输入非cell');end

allPriority = zeros(1,size(ID,2));
uniOrd = unique(SIDord);
for i=1:length(uniOrd)
    idxSID = SIDord==uniOrd(i);
    
    tLW = LW(:,idxSID);   tLL = LoadingRateLimit(:,idxSID);     tL = Loadingrate(:,idxSID);
    
    tID = ID(:,idxSID);   celltID = tID;
    
    priority = 1;
    SIDorder = zeros(1,size(tID,2));
    szRow = cellfun(@(x)size(x,1), tID);
    if isscalar(szRow)
        SIDorder = 1;
    else
        [tID, ~] = padcat(tID{:});   if iscolumn(tID), tID = tID'; end
        
        % 1 给定不考虑相邻的Strip摆放顺序
        tmpM = [tLL;tL;tLW;];
            %         [~,order] = sortrows(tmpM',[1,2,4],{'descend','descend','descend'}); %order is index vector 返回原索引值 如第一个为3，表示原array中第3个目前是第1个
        % 改用和Item2Strip类似的以高度优先的排序方式, 其次未LoadingRate
        [~,order] = sortrows(tmpM',[4,1],{'descend','descend'}); %order is index vector 返回原索引值 如第一个为3，表示原array中第3个目前是第1个
        if ~isrow(order), order=order'; end
        
        % 2 基于1的摆放顺序, 给定最终STRIP顺序到torder
        while any(SIDorder==0)
            % SIDorder中非0且order中最大的那个作为首选STRIP
            % 2.1 找出无相邻前提下,第一个o及对应的order(o)
            [~,o]=find(SIDorder(order) == 0,1,'first');
            SIDorder(order(o)) = priority;  priority=priority+1;
            
            % 2.2 找出给order(o)位置tLID对应的相邻Strip.
            tLID = celltID{:,order(o)};
            % 首次tLID不应该出现混合型STRIP
            [SIDorder,priority] = getAdjPriority(priority,order,SIDorder,tID,tLID);
        end
    end
    allPriority(idxSID) = SIDorder;
end
end

% 寻找tLID对应的相邻STRIP,按照一定逻辑
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
                if ~isscalar(nbcol2),  error('nbcol2包含不只一个数');   end
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
            error('eeee');
        end
    end
end


%% V2 : 采用cell方式  获取STRIP内部包含ID后的排列顺序 BUG: 未考虑单strip多LUID,后续非该LUID内的strip行为
% % function order = getOrderofID(ID)
% % % ID is cell type % ID = [2 3 NaN; 2 1 3; 1 3 NaN; 2 1 NaN]';
% %
% % if ~iscell(ID),error('输入非cell');end
% % sz = cellfun(@(x)size(x,1), ID);
% %
% % % 如cell为等长,则直接转换;否则调用padcat
% % if numel(unique(sz)) == 1
% %     ID = cell2mat(ID);
% % else
% %     % 将cell内的矩阵转换为相等行数 padcat
% %     [ID, ~] = padcat(ID{:});
% % end
% %
% % % 初始化
% % ID = sort(ID,1); %先对ID按列排序->最小的放前面
% % priority=0;
% % order = zeros(1,size(ID,2));
% %
% % for t=1:max(ID(:)) %从PID/SID的第1个直到最后1个,逐次寻找并给定顺序
% %     idx3 = find(sum(~isnan( ID(:,:) ),1) >=3); %找出第t行有t的列idx 且非NaN值有多于1个列其次赋权
% %     if ~isempty(idx3), warning('三个以上同在一个strip,需要关注');end
% %
% %     idx1 = find(ID(1,:)==t & sum(~isnan( ID(:,:) ),1) ==1); %找出第t行有t的列idx 且非NaN值仅有1个的列优先赋权
% %     idx2 = find(ID(1,:)==t & sum(~isnan( ID(:,:) ),1) ~=1); %找出第t行有t的列idx 且非NaN值有多于1个列其次赋权
% %
% %     if ~isempty(idx1)
% %         priority=priority+1;
% %         order(idx1) = priority;
% %     end
% %     if ~isempty(idx2)
% %         priority=priority+1;
% %         order(idx2) = priority;
% %     end
% % end
% %
% % end


%% V1 : 直接采用矩阵方式
% % % function z = getOrderofID(ID)
% % % t = ID;
% % %  ss = sum(t,1);  %每个STRIP内包含的SID个数
% % %
% % %  z = zeros(1,size(t,2));
% % %  k=1;
% % %  for i=1:numel(ss)
% % %     if ss(i) ==1
% % %         if i>1 && ss(i-1) ==1 && find(t(:, i)==1) ~= find(t(:, i-1)==1) %判断当前与前一个strip是否属于同样SID
% % %              k=k+1;
% % %         end
% % %         z(i) = k;
% % % %     elseif  i>1 && ss(i) >1 %只要遇到STRIP包含2个及以上的STRIP时，更新顺序
% % % %         k=k+1;
% % % %         z(i)=k;
% % % %         k=k+1;
% % %     elseif ss(i) >1 %只要遇到STRIP包含2个及以上的STRIP时，更新顺序
% % %         k=k+1;
% % %         z(i)=k;
% % %         k=k+1;
% % %  end
% % % end


%  z = zeros(1,size(t,2));
%  k=1;
%  for i=1:numel(ss)
%     if ss(i) ==1
%         if i>1 && ss(i-1) ==1 && find(t(:, i)==1) ~= find(t(:, i-1)==1) %判断当前与前一个strip是否属于同样SID
%             k=k+1;
%         end
%         z(i) = k;
%     elseif ss(i) >1 %只要遇到STRIP包含2个及以上的STRIP时，更新顺序
%         k=k+1;
%         z(i)=k;
%         k=k+1;
%     end
%  end
% end