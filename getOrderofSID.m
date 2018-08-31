%% V3 : 采用cell方式 且分开ITEM和STRIP两类 重点排除order被重复赋值的情况
% SID: 按照SID由小到大的顺序排列
function order = getOrderofSID(ID)
% ID is cell type % ID = [2 3 NaN; 2 1 3; 1 3 NaN; 2 1 NaN]';

if ~iscell(ID),error('输入非cell');end
szRow = cellfun(@(x)size(x,1), ID);

    if (max(szRow)==min(szRow) && ~isscalar(szRow)) %所有STRIP均为单一型
        mID = cell2mat(ID);
        order = mID;
        return;
    else %有STRIP为混合型-但不存在三个以上的混合
        if max(szRow)>=3
            error('同一STRIP混合型存在三个及以上SID');
        end
        warning('同一STRIP混合型存在两个SID');
        
        [mID, ~] = padcat(ID{:}); % 否则调用padcat:将cell内的矩阵转换为相等行数
        order = zeros(1,size(mID,2));
        priority=1;
        
        SID = 1;  %其实SID,从1开始
        
        while 1
            if all(order)
                 mID
                 order
                break;
            end %如果order全部大于0,表示已经全部排序
            
            % idx1(单strip为单一型)可能有多个 idx2(单strip为混合型) 只能有一个
            [~,nbcol1] = find(mID(:,:)==SID & sum(~isnan( mID(:,:) ),1) ==1); %找出mID中包含SID的列,且该列仅有1个
            [~,nbcol2] = find( mID(:,:)==SID & sum(~isnan( mID(:,:) ),1) ~=1); %找出mID中包含SID的列,且该列不只1个
            
            % 单一型strip的idx1非空 且 其order未全部被安排
            if ~isempty(nbcol1) &&  any(order(nbcol1)==0)
                if any(order(nbcol1)~=0)
                    error('单strip的order存在重复赋值可能！');
                end
                order(nbcol1) = priority;
                priority=priority+1;
                
                % 更新下一个SID: 非混合型按递增排序 
                nextSID = SID + 1; %逐个SID判断,SID中间不能非连续
            end
            
            % 混合型strip的idx2非空 且 其order未全部被安排
            if ~isempty(nbcol2) &&  any(order(nbcol2)==0)
                if any(diff(mID(:,nbcol2))~=1)  %如果存在混合型SID，且混合的非连续值, 只能是错误
                    error('混合型strip的存在混合非连续值！');
                end
                
                
                %将已赋值的order排除在外,即order以首次为准, 不可重复安排
                if any(order(nbcol2)~=0)
                    flagRemove=any(order(nbcol2)~=0);
                    nbcol2(flagRemove)=[];
                end
                
                
                if numel(nbcol2)>1
                    error('混合型strip且未赋值的存在多个！');
                end
                order(nbcol2) = priority;
                priority=priority+1;
                
                % 更新下一个SID: 找出idx2中排除本SID后的另一个SID 不能为多个
                a = mID(:, nbcol2);
                a(a(:)==SID)=[];
                if ~isscalar(a),   error('同一STRIP混合型存在三个及以上SID！');    end
                nextSID = a;
            end
            
            % 防错措施 - 不应该有
            if isempty(nbcol1) && isempty(nbcol2)
                error('不存在包含该SID的strip！');
            end
            
            SID = nextSID;            
            
        end
    end
end

% 初始化
% % mID = sort(mID,1); %先对ID按列排序->最小的放前面
% % lastTID = -1;
% % nbID = max(mID(:));
% % typeID = zeros(1,nbID);
% % tmpidx=zeros(1,numel(szRow));
% % while 1
% %     if all(typeID), break; end %如果aID全部=1,表示已经全部排序
% %     idx3 = find(sum(~isnan( mID(:,:) ),1) >=3); %找出第t行有t的列idx 且非NaN值有多于1个列其次赋权
% %     if ~isempty(idx3), warning('三个以上同在一个strip,需要关注');end
% %     
% % 
% %     1
% %     if ~isempty(idx1)
% %         priority=priority+1;
% %         order(idx1) = priority;
% %         typeID(tID) = 1; %tID已经使用
% %     end
% %     if ~isempty(idx2)
% %         priority=priority+1;
% %         order(idx2) = priority;
% %         typeID(tID) = 1; %tID已经使用
% %     end    
% %     if isempty(idx1) && isempty(idx2)
% %         tmptID = find(typeID(:) == 0); %第一个满足条件的序号返回
% %         if ~isempty(tmptID)
% %             for i=1:numel(tmptID)
% %               tID = find(typeID(:) == 0,i); 
% %                if tID ~= lastTID,
% %                    lastTID = tID;
% %                    continue;
% %                end
% %             end
% %         end
% %     end
% %     if ~isempty(idx2)
% %         tID =  mID(2,idx2); %使用本strip对应的多余的1个
% %     elseif ~isempty(idx1)
% %         tID = find(typeID(:) == 0,1); %第一个满足条件的序号返回
% %     end
% %     
% % end
% % 
% % end


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