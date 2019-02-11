function [TF] = isrepeated(id, fatherid)
%isrepeated �ж���ͬһfatherid��, id�Ƿ��ظ�
%   return a logical scalar

if nargin < 2,  fatherid = ones(1:length(id));   end
validateattributes(id,{'double'},{'integer','vector'});
validateattributes(fatherid,{'double'},{'integer','vector'});

TF = false;
% �ع�id���� ����fatherid��ͬ��Ϊһ������
ufid = unique(fatherid);
for i=1:numel(ufid)
    id1 = id(fatherid == ufid(i));
    id2 = id(fatherid ~= ufid(i));
    interid = intersect(id1,id2);
    if ~isempty(interid)
        TF = true;
         return;
    end
end

end










% function [TF] = issimplified(pshape)
% % ISSIMPLIFIED Determine if a polyshape is simplified
% %
% % TF = ISSIMPLIFIED(pshape) turns a logical array whose elements are true 
% % if the corresponding element of pshape is a well-defined polygon.
% %
% % See also polyshape, simplify, rmslivers
% 
% % Copyright 2016-2017 The MathWorks, Inc.
% 
% n = polyshape.checkArray(pshape);
% 
% TF = false(n);
% for i=1:numel(pshape)
%     if pshape(i).isEmptyShape()
%         TF(i) = true;
%     elseif pshape(i).SimplifyState >= 0
%         TF(i) = logical(pshape(i).SimplifyState);
%     else
%         [~, canBeSimplified] = checkAndSimplify(pshape(i), false);
%         TF(i) = ~canBeSimplified;
%     end
% end
