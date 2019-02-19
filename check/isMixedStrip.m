%% 函数3: 判断STRIP是否混合型isMixedStrip
function [isMixed,isMixedSID,isMixedEID] = isMixedStrip(Strip)

    n = length(Strip.Weight);
    [isMixed,isMixedSID,isMixedEID] = deal(zeros(1,n));
        
    % 循环判断Strip是否为不同LU.ID的混合型
    for i=1:n
        
         if numel(Strip.LID{i}) > 1
             isMixed(i) = 1;             
         else
             isMixed(i) = 0;
         end

         if numel(Strip.SID{i}) > 1
             isMixedSID(i) = 1;
         else
             isMixedSID(i) = 0;
         end
         
         if numel(Strip.EID{i}) > 1
             isMixedEID(i) = 1;
         else
             isMixedEID(i) = 0;
         end
         
    end
end
