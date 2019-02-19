%% ����3: �ж�STRIP�Ƿ�����isMixedStrip
function [isMixed,isMixedSID,isMixedEID] = isMixedStrip(Strip)

    n = length(Strip.Weight);
    [isMixed,isMixedSID,isMixedEID] = deal(zeros(1,n));
        
    % ѭ���ж�Strip�Ƿ�Ϊ��ͬLU.ID�Ļ����
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
