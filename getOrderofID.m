 % �жϽṹ��s�Ƿ�ȫ����ͬ�еľ��󹹳�
function z = getOrderofID(ID)
t = ID;
 ss = sum(t,1);  %ÿ��STRIP�ڰ�����SID����
 
%   for i=1:size(t,1)
%      if sum(ss( find(t(i,:))  ) > 1) > 1
%          error('��ͬһ��SID��2��������STRIP����');
%      end
%  end
% if  any(sum(t)>2)
%     error('��ͬһ��strip����3�������ϸ�SID');
% end

 z = zeros(1,size(t,2)); 
 k=1;
 for i=1:numel(ss)
    if ss(i) ==1
        if i>1 && ss(i-1) ==1 && find(t(:, i)==1) ~= find(t(:, i-1)==1) %�жϵ�ǰ��ǰһ��strip�Ƿ�����ͬ��SID
             k=k+1;
        end     
        z(i) = k;
%     elseif  i>1 && ss(i) >1 %ֻҪ����STRIP����2�������ϵ�STRIPʱ������˳��                
%         k=k+1;
%         z(i)=k;
%         k=k+1;
    elseif ss(i) >1 %ֻҪ����STRIP����2�������ϵ�STRIPʱ������˳��      
        k=k+1;
        z(i)=k;
        k=k+1;
 end
end


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