function [id] = getReviseID(id, fatherid)
%getReviseID ��ͬһfatherid��,ȷ���ظ�ʱ���ı䵱��id�������ظ�������,
%   return a non-repeated id

ufid = unique(fatherid);
for i=1:numel(ufid)
    id1 = id(fatherid == ufid(i));
    id2 = id(fatherid ~= ufid(i));
    interid = intersect(id1,id2);
    if ~isempty(interid) %�����ظ�����Ҫ����
        if ~isempty(id(id1 ~= interid))  %����ǿգ�����maxֵ
            id(id1 == interid) = id(id1 == interid) + max(id(id1 ~= interid));  % �������ͬһSID������ID�ظ�
        else  %�����:�������ǽ���������1�϶�û����
            id(id1 == interid) = id(id1 == interid) + 1;
        end
         if isrepeated(id,fatherid)
             id = getReviseID(id,fatherid);
         end
    end
end

end