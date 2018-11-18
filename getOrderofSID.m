%% V4 : SID is cell方式 且分开ITEM和STRIP两类 重点排除order被重复赋值的情况
function SIDallPriority = getOrderofSID(SID)

    if ~iscell(SID),error('输入非cell');end
    
    [mSID, isSID] = padcat(SID{:});     
    if iscolumn(mSID), mSID = mSID'; end 
    if iscolumn(isSID), isSID = isSID'; end % 调用padcat: 将cell内的矩阵转换为相等行数
    nbSIDperStrip = sum(isSID,1);     %nbSIDperStrip: 向量,每个strip内SID的个数     nbSIDperStrip = cellfun(@(x)size(x,1), SID);

    % CASE 1: 所有STRIP均为单一型 或仅有1个Strip的单一型(isSID是一维scalar=1)或混合型(isSID是一维列向量)
    if isvector(isSID) || isscalar(isSID)
        % mSID = cell2mat(SID);
        SIDallPriority = mSID;  return;

    % CASE 2: 存在包含多个SID的STRIP
    else %有STRIP为混合型-但不存在三个以上的混合
        if max(nbSIDperStrip)>=3,  error('同一STRIP混合存在三个及以上SID');    end
        if max(nbSIDperStrip)==1,  error('STRIP均非混合型,应该已经被排除了');    end
        warning('同一STRIP混合型存在两个SID');
        
        uniSID = unique(mSID(~isnan(mSID)));          if ~issorted(uniSID,'strictascend'), error('本BIN内STRIP的SID号未按指定SID(严格)递增排序.'); end
        priority=1;  
        SIDallPriority = zeros(1,length(nbSIDperStrip));
        
        for idxSID=1:length(uniSID)
         
            fsid = any(ismember(mSID,uniSID(idxSID)),1);     % 包含本sid的trip逻辑判断
            fprio = SIDallPriority==0;                                       % 未赋权重值逻辑判断
            fmix = nbSIDperStrip>1;                                         % 混合型逻辑判断
            
            f1 = fsid & fprio & ~fmix;       % 非混合
            f2 = fsid & fprio & fmix;          % 混合
            
            % 单一型strip的idx非空 且 其priority未全部被安排 的存在
            if any(f1)
                SIDallPriority(f1) = priority;
                priority=priority+1;
            end
            
            % 混合型strip的idx非空 且 其priority未全部被安排 的存在
            if any(f2)
                if sum(f2)>1, error('同一SID,其出现在多个混合型strip! '); end
                if ~issorted(mSID(:,f2),'strictascend'), error('同一SID,其混合strip内SID是非递增连续值(跨SID排序)'); end % TODO , 排除mSID(:,f2)包含NaN的情况
                SIDallPriority(f2) = priority;
                priority=priority+1;
                
                % 排错CODE
                nextSID = mSID(:,f2);  nextSID(nextSID(:)==uniSID(idxSID))=[];   if ~isscalar(nextSID) || isempty(nextSID),   error('同一STRIP混合型存在三个及以上SID或为空！');    end
                if uniSID(idxSID+1) ~= nextSID, error('如果存在多个SID,下一个SID必须是和本SID混合的供应商,可能是非连续的问题'); end
            end
            
            % 排错CODE
            if ~any(f1) && ~any(f2), warning('不存在包含该SID的strip！,可能是该SID数量过少,或已在混合case中安排'); end
        end
        
        % 排错CODE
        if ~all(SIDallPriority),  error('还有未权重的SID值'); end
        
    end
end

            

%         while ~done
%             if all(SIDallPriority),  done=false;   continue;  end   %如果order全部大于0,表示已经全部排序
   
            % V2: idx1(单strip为单一型)可能有多个 idx2(单strip为混合型) 只能有一个
%             [~,nbcol1] = find(mID(:,:)==uniSID(SID) & sum(~isnan( mID(:,:) ),1) ==1 & allPriority==0 ); %找出mID中包含SID的列,且该列仅有1个
%             [~,nbcol2] = find( mID(:,:)==uniSID(SID) & sum(~isnan( mID(:,:) ),1) ~=1 & allPriority==0); %找出mID中包含SID的列,且该列不只1个
%             
            % V1: idx1(strip为单一型)可能有多个 idx2(strip为混合型) 只能有一个
%             [~,nbcol1] = find(mSID(:,:)==idxSID & sum(~isnan( mSID(:,:) ),1) ==1 & SIDallPriority==0 ); %找出mID中包含SID的列,且该列仅有1个
%             [~,nbcol2] = find( mSID(:,:)==idxSID & sum(~isnan( mSID(:,:) ),1) ~=1 & SIDallPriority==0); %找出mID中包含SID的列,且该列不只1个
             
%             % 单一型strip的idx1非空 且 其order未全部被安排
%             if ~isempty(nbcol1)
%                 SIDallPriority(nbcol1) = priority;
%                 priority=priority+1;                
%                 % 更新下一个SID: 非混合型按递增排序 
%                 nextSID = idxSID + 1; %逐个SID判断,SID中间不能非连续
%             end
            
%             % 混合型strip的idx2非空 且 其order未全部被安排
%             if ~isempty(nbcol2)
%                 if any(diff(mSID(:,nbcol2))~=1)  %如果存在混合型SID，且混合的非连续值, 只能是错误
%                     error('混合型strip的存在混合非连续值！');
%                 end                               
%                 if numel(nbcol2)>1
%                     error('混合型strip且未赋值的存在多个！');
%                 end
%                 SIDallPriority(nbcol2) = priority;
%                 priority=priority+1;
            
%                 % 更新下一个SID: 找出idx2中排除本SID后的另一个SID 不能为多个
%                 hyLID = mSID(:, nbcol2);
%                 hyLID(hyLID(:)==idxSID)=[];
%                 if ~isscalar(hyLID),   error('同一STRIP混合型存在三个及以上SID！');    end
%                 nextSID = hyLID;

            
%             % 防错措施 - 不应该有
%             if isempty(nbcol1) && isempty(nbcol2)
%                 warning('不存在包含该SID的strip！,可能是该SID数量过少,已在混合case中安排');
%                 if ~exist('nextSID'),
%                     error('不存在nextSID');
%                     nextSID=idxSID; end;
%                 nextSID=nextSID+1;
%             end
            
%             idxSID = nextSID;
            



%% 
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


%% %% V3 : 采用cell方式 且分开ITEM和STRIP两类 重点排除order被重复赋值的情况
% % % SID: 按照SID由小到大的顺序排列
% % function allPriority = getOrderofSID(ID)
% % % ID is cell type % ID = [2 3 NaN; 2 1 3; 1 3 NaN; 2 1 NaN]';
% % 
% % if ~iscell(ID),error('输入非cell');end
% % szRow = cellfun(@(x)size(x,1), ID);
% % 
% %     if (max(szRow)==min(szRow) && ~isscalar(szRow)) %所有STRIP均为单一型
% %         mID = cell2mat(ID);
% %         allPriority = mID;
% %         return;
% %     else %有STRIP为混合型-但不存在三个以上的混合
% %         if max(szRow)>=3
% %             error('同一STRIP混合型存在三个及以上SID');
% %         end
% %         warning('同一STRIP混合型存在两个SID');
% %         
% %         
% %         [mID, ~] = padcat(ID{:}); % 否则调用padcat:将cell内的矩阵转换为相等行数        
% %         uniSID = unique(mID(~isnan(mID)));
% %         allPriority = zeros(1,size(mID,2));
% %         priority=1;        
% %         SID = 1;  %其实SID,从1开始        
% %         while 1
% %             if all(allPriority)
% %                  %mID
% %                  %allPriority
% %                 break;
% %             end %如果order全部大于0,表示已经全部排序
% %             
% %             % V2: idx1(单strip为单一型)可能有多个 idx2(单strip为混合型) 只能有一个
% % %             [~,nbcol1] = find(mID(:,:)==uniSID(SID) & sum(~isnan( mID(:,:) ),1) ==1 & allPriority==0 ); %找出mID中包含SID的列,且该列仅有1个
% % %             [~,nbcol2] = find( mID(:,:)==uniSID(SID) & sum(~isnan( mID(:,:) ),1) ~=1 & allPriority==0); %找出mID中包含SID的列,且该列不只1个
% % %             
% %             % V1: idx1(单strip为单一型)可能有多个 idx2(单strip为混合型) 只能有一个
% %             [~,nbcol1] = find(mID(:,:)==SID & sum(~isnan( mID(:,:) ),1) ==1 & allPriority==0 ); %找出mID中包含SID的列,且该列仅有1个
% %             [~,nbcol2] = find( mID(:,:)==SID & sum(~isnan( mID(:,:) ),1) ~=1 & allPriority==0); %找出mID中包含SID的列,且该列不只1个
% %             
% %             % 单一型strip的idx1非空 且 其order未全部被安排
% %             if ~isempty(nbcol1)
% %                 allPriority(nbcol1) = priority;
% %                 priority=priority+1;                
% %                 % 更新下一个SID: 非混合型按递增排序 
% %                 nextSID = SID + 1; %逐个SID判断,SID中间不能非连续
% %             end
% %             
% %             % 混合型strip的idx2非空 且 其order未全部被安排
% %             if ~isempty(nbcol2)
% %                 if any(diff(mID(:,nbcol2))~=1)  %如果存在混合型SID，且混合的非连续值, 只能是错误
% %                     error('混合型strip的存在混合非连续值！');
% %                 end                               
% %                 if numel(nbcol2)>1
% %                     error('混合型strip且未赋值的存在多个！');
% %                 end
% %                 allPriority(nbcol2) = priority;
% %                 priority=priority+1;
% %             
% %                 % 更新下一个SID: 找出idx2中排除本SID后的另一个SID 不能为多个
% %                 hyLID = mID(:, nbcol2);
% %                 hyLID(hyLID(:)==SID)=[];
% %                 if ~isscalar(hyLID),   error('同一STRIP混合型存在三个及以上SID！');    end
% %                 nextSID = hyLID;
% %             end
% %             
% %             % 防错措施 - 不应该有
% %             if isempty(nbcol1) && isempty(nbcol2)
% %                 warning('不存在包含该SID的strip！,可能是该SID数量过少,已在混合case中安排');
% %                 if ~exist('nextSID'),
% %                     error('不存在nextSID');
% %                     nextSID=SID; end;
% %                 nextSID=nextSID+1;
% %             end
% %             
% %             SID = nextSID;
% %             
% %         end
% %     end
% % end

%% 
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