% �жϽṹ��s�Ƿ�ȫ����ͬ�еľ��󹹳�
function flag = isSameCol(s)
if isstruct(s)
    sCol = struct2cell( structfun(@(x) size(x,2), s, 'UniformOutput', false));
    nCol = numel(unique(cell2mat( sCol )));
    if nCol == 1
        flag = true;
    else
        flag = false; end
else
    error('���벻�ǽṹ��');
end
end

