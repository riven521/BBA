 % 判断结构体s是否全部相同列的矩阵构成
function z = getOrderofID(ID)
t = ID;
 ss = sum(t,1);  %每个STRIP内包含的SID个数
 
%   for i=1:size(t,1)
%      if sum(ss( find(t(i,:))  ) > 1) > 1
%          error('有同一个SID被2个及以上STRIP包括');
%      end
%  end
% if  any(sum(t)>2)
%     error('有同一个strip包括3个及以上各SID');
% end

 z = zeros(1,size(t,2)); 
 k=1;
 for i=1:numel(ss)
    if ss(i) ==1
        if i>1 && ss(i-1) ==1 && find(t(:, i)==1) ~= find(t(:, i-1)==1) %判断当前与前一个strip是否属于同样SID
             k=k+1;
        end     
        z(i) = k;
%     elseif  i>1 && ss(i) >1 %只要遇到STRIP包含2个及以上的STRIP时，更新顺序                
%         k=k+1;
%         z(i)=k;
%         k=k+1;
    elseif ss(i) >1 %只要遇到STRIP包含2个及以上的STRIP时，更新顺序      
        k=k+1;
        z(i)=k;
        k=k+1;
 end
end


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