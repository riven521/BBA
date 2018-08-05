%% 行列转换检验
function [LU,Veh] = initCheck(LU,Veh)
% 判断小于0
fLU = structfun(@(x) any( x(:) < 0), LU, 'UniformOutput', false);
fVeh = structfun(@(x) any( x(:) < 0), LU, 'UniformOutput', false);
if any(cell2mat(struct2cell(fLU))) || any(cell2mat(struct2cell(fVeh)))
    error('输入数据存在小于0的数'); end
        
nVeh = length(Veh.ID);
nLU = length(LU.ID);
% 判断Veh向量是否行向量; 矩阵列是否一致; 不一致进行旋转
fields  = fieldnames(Veh);
for idx = 1:length(fields)
    aField = Veh.(fields{idx});
    if ~ismatrix(aField) || ~isnumeric(aField), error('输入数据存在错误'); end
    if isvector(aField) && ~strcmp('LWH',fields{idx}) && ~strcmp('buff',fields{idx})
        if ~isrow(aField)
            Veh.(fields{idx}) = aField';    end
        if length(aField) ~= nVeh
            error('Veh数据输入存在错误');   end
    else % ismatrix
        if size(aField,1) == nVeh && size(aField,2) ~= nVeh
            Veh.(fields{idx}) = aField';
            if size(aField',2)  ~= nVeh
                error('Veh数据输入存在错误');   end
        end
    end
end
% 判断LU向量是否行向量; 矩阵列是否一致
fields  = fieldnames(LU);
for idx = 1:length(fields)
    aField = LU.(fields{idx});
    if ~ismatrix(aField) || ~isnumeric(aField), error('输入数据存在错误'); end
    if isvector(aField)
        if ~isrow(aField)
            LU.(fields{idx}) = aField';     end
        if length(aField) ~= nLU
            error('LU数据输入存在错误');   end
    else % ismatrix
        if size(aField,1) == nLU && size(aField,2) ~= nLU
            LU.(fields{idx}) = aField';    end
        if size(LU.(fields{idx}),2)  ~= nLU
            error('LU数据输入存在错误');   end
    end
end
end