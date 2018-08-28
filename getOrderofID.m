%% V3 : ����cell��ʽ ��
function order = getOrderofID(ID)
% ID is cell type % ID = [2 3 NaN; 2 1 3; 1 3 NaN; 2 1 NaN]';

if ~iscell(ID),error('�����cell');end
sz = cellfun(@(x)size(x,1), ID);

% ��cellΪ�ȳ�,��ֱ��cell2mtaת��; ��ֱ�����order
if  (max(sz)==min(sz)) 
    mID = cell2mat(ID)
    order = sz;
    return;
else
    % �������padcat:��cell�ڵľ���ת��Ϊ������� padcat
    [mID, ~] = padcat(mID{:});
end

% ��ʼ��
mID = sort(mID,1); %�ȶ�ID��������->��С�ķ�ǰ��
priority=0;
order = zeros(1,size(mID,2));

tID = 1;
lastTID = -1;
nbID = max(mID(:));
typeID = zeros(1,nbID);
tmpidx=zeros(1,numel(sz));
while 1
    if all(typeID), break; end %���aIDȫ��=1,��ʾ�Ѿ�ȫ������
    idx3 = find(sum(~isnan( mID(:,:) ),1) >=3); %�ҳ���t����t����idx �ҷ�NaNֵ�ж���1������θ�Ȩ
    if ~isempty(idx3), warning('��������ͬ��һ��strip,��Ҫ��ע');end
    
    % �ҵ��������е�
    [~,idx1] = find(mID(:,:)==tID & sum(~isnan( mID(:,:) ),1) ==1); %�ҳ���t����t����idx �ҷ�NaNֵ����1���������ȸ�Ȩ
    tmpidx(idx1) = 1;
    [~,idx2] = find( mID(:,:)==tID & sum(~isnan( mID(:,:) ),1) ~=1 & tmpidx ==0); %�ҳ���t����t����idx �ҷ�NaNֵ�ж���1������θ�Ȩ
    tmpidx(idx2) = 1;
    
    if ~isempty(idx1)
        priority=priority+1;
        order(idx1) = priority;
        typeID(tID) = 1; %tID�Ѿ�ʹ��
    end
    if ~isempty(idx2)
        priority=priority+1;
        order(idx2) = priority;
        typeID(tID) = 1; %tID�Ѿ�ʹ��
    end    
    if isempty(idx1) && isempty(idx2)
        tmptID = find(typeID(:) == 0); %��һ��������������ŷ���
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
        tID =  mID(2,idx2); %ʹ�ñ�strip��Ӧ�Ķ����1��
    elseif ~isempty(idx1)
        tID = find(typeID(:) == 0,1); %��һ��������������ŷ���
    end
    
end

end


%% V2 : ����cell��ʽ  ��ȡSTRIP�ڲ�����ID�������˳�� BUG: δ���ǵ�strip��LUID,�����Ǹ�LUID�ڵ�strip��Ϊ
% % function order = getOrderofID(ID)
% % % ID is cell type % ID = [2 3 NaN; 2 1 3; 1 3 NaN; 2 1 NaN]';
% % 
% % if ~iscell(ID),error('�����cell');end
% % sz = cellfun(@(x)size(x,1), ID);
% % 
% % % ��cellΪ�ȳ�,��ֱ��ת��;�������padcat
% % if numel(unique(sz)) == 1
% %     ID = cell2mat(ID);
% % else
% %     % ��cell�ڵľ���ת��Ϊ������� padcat
% %     [ID, ~] = padcat(ID{:});
% % end
% % 
% % % ��ʼ��
% % ID = sort(ID,1); %�ȶ�ID��������->��С�ķ�ǰ��
% % priority=0;
% % order = zeros(1,size(ID,2));
% % 
% % for t=1:max(ID(:)) %��PID/SID�ĵ�1��ֱ�����1��,���Ѱ�Ҳ�����˳��
% %     idx3 = find(sum(~isnan( ID(:,:) ),1) >=3); %�ҳ���t����t����idx �ҷ�NaNֵ�ж���1������θ�Ȩ
% %     if ~isempty(idx3), warning('��������ͬ��һ��strip,��Ҫ��ע');end
% %     
% %     idx1 = find(ID(1,:)==t & sum(~isnan( ID(:,:) ),1) ==1); %�ҳ���t����t����idx �ҷ�NaNֵ����1���������ȸ�Ȩ
% %     idx2 = find(ID(1,:)==t & sum(~isnan( ID(:,:) ),1) ~=1); %�ҳ���t����t����idx �ҷ�NaNֵ�ж���1������θ�Ȩ
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


%% V1 : ֱ�Ӳ��þ���ʽ
% % % function z = getOrderofID(ID)
% % % t = ID;
% % %  ss = sum(t,1);  %ÿ��STRIP�ڰ�����SID����
% % % 
% % %  z = zeros(1,size(t,2)); 
% % %  k=1;
% % %  for i=1:numel(ss)
% % %     if ss(i) ==1
% % %         if i>1 && ss(i-1) ==1 && find(t(:, i)==1) ~= find(t(:, i-1)==1) %�жϵ�ǰ��ǰһ��strip�Ƿ�����ͬ��SID
% % %              k=k+1;
% % %         end     
% % %         z(i) = k;
% % % %     elseif  i>1 && ss(i) >1 %ֻҪ����STRIP����2�������ϵ�STRIPʱ������˳��                
% % % %         k=k+1;
% % % %         z(i)=k;
% % % %         k=k+1;
% % %     elseif ss(i) >1 %ֻҪ����STRIP����2�������ϵ�STRIPʱ������˳��      
% % %         k=k+1;
% % %         z(i)=k;
% % %         k=k+1;
% % %  end
% % % end


%  z = zeros(1,size(t,2)); 
%  k=1;
%  for i=1:numel(ss)
%     if ss(i) ==1
%         if i>1 && ss(i-1) ==1 && find(t(:, i)==1) ~= find(t(:, i-1)==1) %�жϵ�ǰ��ǰһ��strip�Ƿ�����ͬ��SID
%             k=k+1;
%         end
%         z(i) = k;
%     elseif ss(i) >1 %ֻҪ����STRIP����2�������ϵ�STRIPʱ������˳��        
%         k=k+1;
%         z(i)=k;
%         k=k+1;
%     end
%  end
% end