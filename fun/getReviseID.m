function [id] = getReviseID(id, fatherid)
%getReviseID 在同一fatherid下,确定重复时，改变当下id与其他重复的数字,
%   return a non-repeated id

ufid = unique(fatherid);
for i=1:numel(ufid)
    id1 = id(fatherid == ufid(i));
    id2 = id(fatherid ~= ufid(i));
    interid = intersect(id1,id2);
    if ~isempty(interid) %存在重复，需要更正
        if ~isempty(id(id1 ~= interid))  %如果非空，增加max值
            id(id1 == interid) = id(id1 == interid) + max(id(id1 ~= interid));  % 不允许和同一SID下其它ID重复
        else  %如果空:表明均是交集，增加1肯定没问题
            id(id1 == interid) = id(id1 == interid) + 1;
        end
         if isrepeated(id,fatherid)
             id = getReviseID(id,fatherid);
         end
    end
end

end