% 判断结构体s是否全部相同列的矩阵构成
function flag = isSameCol(s)
if isstruct(s)
    sCol = struct2cell( structfun(@(x) size(x,2), s, 'UniformOutput', false));
    nCol = numel(unique(cell2mat( sCol )));
    if nCol == 1
        flag = true;
    else
        flag = false; end
else
    error('输入不是结构体');
end
end

