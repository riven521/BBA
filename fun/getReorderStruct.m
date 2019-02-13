function oS = getReorderStruct(o,S)
oS = struct;
if isstruct(S) && isrow(o)
    % 动态数组并排序
    fields = fieldnames(S);
    for idx = 1:length(fields) 
        oS.(fields{idx})(:,o) = S.(fields{idx});   %    aField = S.(fields{idx});
    end
else
    error('输入不是结构体');
end
end
