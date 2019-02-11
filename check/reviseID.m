function [id] = reviseID(id, fatherid)
%reviseID 在同一fatherid下,确定重复时，改变当下id与其他重复的数字,
%   return a non-repeated id

ufid = unique(fatherid);
for i=1:numel(ufid)
    id1 = id(fatherid == ufid(i));
    id2 = id(fatherid ~= ufid(i));
    interid = intersect(id1,id2);
    if ~isempty(interid) %存在重复，需要更正
         id(id1 == interid) = id(id1 == interid) + max(id(id1 ~= interid));  % 不允许和同一SID下其它ID重复
         if isrepeated(id,fatherid)
             id = reviseID(id,fatherid);
         end
    end
end

end