function oS = getReorderStruct(o,S)
oS = struct;
if isstruct(S) && isrow(o)
    % ��̬���鲢����
    fields = fieldnames(S);
    for idx = 1:length(fields) 
        oS.(fields{idx})(:,o) = S.(fields{idx});   %    aField = S.(fields{idx});
    end
else
    error('���벻�ǽṹ��');
end
end
