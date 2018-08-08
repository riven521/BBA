 % 判断结构体s是否全部相同列的矩阵构成
function oS = reorderStruct(o,S)
os = struct;
if isstruct(S) && isrow(o)
    fields = fieldnames(S);
for idx = 1:length(fields)
%    aField = S.(fields{idx});
   oS.(fields{idx})(:,o) = S.(fields{idx});
end

%     sCol = struct2cell( structfun(@(x) size(x,2), s, 'UniformOutput', false));
%     nCol = numel(unique(cell2mat( sCol )));
%     if nCol == 1
%         flag = true;
%     else
%         flag = false; end
else
    error('输入不是结构体');
end
end
