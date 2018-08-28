%% V3 : 采用cell方式 且
function order = getOrderofID(ID)
% ID is cell type % ID = [2 3 NaN; 2 1 3; 1 3 NaN; 2 1 NaN]';

if ~iscell(ID),error('输入非cell');end
sz = cellfun(@(x)size(x,1), ID);

% 如cell为等长,则直接cell2mta转换; 且直接输出order
if  (max(sz)==min(sz)) 
    mID = cell2mat(ID)
    order = sz;
    return;
else
    % 否则调用padcat:将cell内的矩阵转换为相等行数 padcat
    [mID, ~] = padcat(mID{:});
end

% 初始化
mID = sort(mID,1); %先对ID按列排序->最小的放前面
priority=0;
order = zeros(1,size(mID,2));

tID = 1;
lastTID = -1;
nbID = max(mID(:));
typeID = zeros(1,nbID);
tmpidx=zeros(1,numel(sz));
while 1
    if all(typeID), break; end %如果aID全部=1,表示已经全部排序
    idx3 = find(sum(~isnan( mID(:,:) ),1) >=3); %找出第t行有t的列idx 且非NaN值有多于1个列其次赋权
    if ~isempty(idx3), warning('三个以上同在一个strip,需要关注');end
    
    % 找到该列所有的
    [~,idx1] = find(mID(:,:)==tID & sum(~isnan( mID(:,:) ),1) ==1); %找出第t行有t的列idx 且非NaN值仅有1个的列优先赋权
    tmpidx(idx1) = 1;
    [~,idx2] = find( mID(:,:)==tID & sum(~isnan( mID(:,:) ),1) ~=1 & tmpidx ==0); %找出第t行有t的列idx 且非NaN值有多于1个列其次赋权
    tmpidx(idx2) = 1;
    
    if ~isempty(idx1)
        priority=priority+1;
        order(idx1) = priority;
        typeID(tID) = 1; %tID已经使用
    end
    if ~isempty(idx2)
        priority=priority+1;
        order(idx2) = priority;
        typeID(tID) = 1; %tID已经使用
    end    
    if isempty(idx1) && isempty(idx2)
        tmptID = find(typeID(:) == 0); %第一个满足条件的序号返回
        if ~isempty(tmptID)
            for i=1:numel(tmptID)
              tID = find(typeID(:) == 0,i); 
               if tID ~= lastTID,
                   lastTID = tID;
                   continue;
               end
            end
        end
    end
    if ~isempty(idx2)
        tID =  mID(2,idx2); %使用本strip对应的多余的1个
    elseif ~isempty(idx1)
        tID = find(typeID(:) == 0,1); %第一个满足条件的序号返回
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