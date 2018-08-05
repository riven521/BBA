%% ����ת������
function [LU,Veh] = initCheck(LU,Veh)
% �ж�С��0
fLU = structfun(@(x) any( x(:) < 0), LU, 'UniformOutput', false);
fVeh = structfun(@(x) any( x(:) < 0), LU, 'UniformOutput', false);
if any(cell2mat(struct2cell(fLU))) || any(cell2mat(struct2cell(fVeh)))
    error('�������ݴ���С��0����'); end
        
nVeh = length(Veh.ID);
nLU = length(LU.ID);
% �ж�Veh�����Ƿ�������; �������Ƿ�һ��; ��һ�½�����ת
fields  = fieldnames(Veh);
for idx = 1:length(fields)
    aField = Veh.(fields{idx});
    if ~ismatrix(aField) || ~isnumeric(aField), error('�������ݴ��ڴ���'); end
    if isvector(aField) && ~strcmp('LWH',fields{idx}) && ~strcmp('buff',fields{idx})
        if ~isrow(aField)
            Veh.(fields{idx}) = aField';    end
        if length(aField) ~= nVeh
            error('Veh����������ڴ���');   end
    else % ismatrix
        if size(aField,1) == nVeh && size(aField,2) ~= nVeh
            Veh.(fields{idx}) = aField';
            if size(aField',2)  ~= nVeh
                error('Veh����������ڴ���');   end
        end
    end
end
% �ж�LU�����Ƿ�������; �������Ƿ�һ��
fields  = fieldnames(LU);
for idx = 1:length(fields)
    aField = LU.(fields{idx});
    if ~ismatrix(aField) || ~isnumeric(aField), error('�������ݴ��ڴ���'); end
    if isvector(aField)
        if ~isrow(aField)
            LU.(fields{idx}) = aField';     end
        if length(aField) ~= nLU
            error('LU����������ڴ���');   end
    else % ismatrix
        if size(aField,1) == nLU && size(aField,2) ~= nLU
            LU.(fields{idx}) = aField';    end
        if size(LU.(fields{idx}),2)  ~= nLU
            error('LU����������ڴ���');   end
    end
end
end