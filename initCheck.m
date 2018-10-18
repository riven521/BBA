%% ����ת������
function [LU,Veh] = initCheck(LU,Veh)
% �ж�С��0
fLU = structfun(@(x) any( x(:) < 0), LU, 'UniformOutput', false);
fVeh = structfun(@(x) any( x(:) < 0), LU, 'UniformOutput', false);
if any(cell2mat(struct2cell(fLU))) || any(cell2mat(struct2cell(fVeh)))
    error('�������ݴ���С��0����'); end

nVeh = length(Veh.ID);
nLU = length(LU.ID);
%% *************** Veh�����ж� ***************
% �ж�Veh�����Ƿ�������; �������Ƿ�һ��; ��һ�½�����ת
fields  = fieldnames(Veh);
% ��Veh������ÿ���ֶηֱ��ж�, ÿ���ֶε����������ǳ�����nVeh
for idx = 1:length(fields)
    aField = Veh.(fields{idx});
    % �ж�1: ��������ֵ��
    if ~ismatrix(aField) || ~isnumeric(aField), error('�������ݴ��ڴ���(����ֵ�;���)'); end
    % �ж�2: �����粻��nVeh,����ת��(�����������,���ڿ����Ǿ������͵�LWH��buff�ֿ��ж�)
    if isvector(aField) && ~strcmp('LWH',fields{idx}) && ~strcmp('buff',fields{idx}) %��������������Ҳ��Ǻ�������ֵ(�˴����жϳ��������LWH)
        if ~isrow(aField) %����������ݷ�������,�޸�
            Veh.(fields{idx}) = aField';         end
                                    %         if length(aField) ~= nVeh %���ж�vector�ſ�����length�ж�
                                    %             error('Veh����������ڴ���1');            end
    else % ismatrix %�˴��жϳ���������Ƿ�д�� LWH buff��
        if size(aField,1) == nVeh && size(aField,2) ~= nVeh %������������==������ �� ����~=������
            Veh.(fields{idx}) = aField'; %����ת��
        end
        if size(aField,1) == nVeh && size(aField,2) == nVeh
            
            if nVeh==3
                Veh.(fields{idx}) = aField';
                warning('���ڴ���3�����������,δ�����ж�,������ж����������Ƿ���ȷ');
                fprintf('�����ĵĳ����Ϊ: %d \n ',Veh.(fields{idx}));
                tmpVeh = Veh.(fields{idx});
                if ~(tmpVeh(1)*1.5 <= tmpVeh(2)) %�����һ�������Ŀ�ȵ�1.5�����ڳ����ĳ���,����Ϊ����߷ŷ���
                    Veh.(fields{idx}) = aField;
                end             
                
            else
                warning('����3���������Ƿ���,�������ڷǳ����ʱ����');  
            end
            
        end        
    end
    % �ж�3: ת�������������nVeh
    if size(Veh.(fields{idx}),2)  ~= nVeh %5555 ��Ҫ�ж�
        error('Veh���ݴ������������ڳ��������,�������2');
    end
end 
%% *************** LU�����ж� ***************
% �ж�LU�����Ƿ�������; �������Ƿ�һ��; ��һ�½�����ת
fields  = fieldnames(LU);
for idx = 1:length(fields)
    aField = LU.(fields{idx});
    % �ж�1: ��������ֵ��
    if ~ismatrix(aField) || ~isnumeric(aField), error('�������ݴ��ڴ���(����ֵ�;���)'); end
    % �ж�2: �����粻��nLU,����ת��    (�����������,���ڿ����Ǿ������͵�LWH��buff��margin�ֿ��ж�)
    if isvector(aField) && ~strcmp('LWH',fields{idx}) && ~strcmp('maxL',fields{idx}) && ~strcmp('margin',fields{idx})%&& ~strcmp('buff',fields{idx} ) %��������������Ҳ��Ǻ�������ֵ
        if ~isrow(aField)
            LU.(fields{idx}) = aField';     end
                                            %         if length(aField) ~= nLU %���ж�vector�ſ�����length�ж�
                                            %             error('LU����������ڲ���ԭ�µĴ���1');   end
    else % ismatrix %�˴��ж����̳�����Ƿ�д��
        if size(aField,1) == nLU && size(aField,2) ~= nLU
            LU.(fields{idx}) = aField';   end
        if size(aField,1) == nLU && size(aField,2) == nLU 

            if nLU==3
                LU.(fields{idx}) = aField';
                warning('���ڴ���3�����̵����,���������');    
                fprintf('���̵ĵĳ����Ϊ: %d \n ',LU.(fields{idx}));
                uniID = unique(LU.ID);
                % ����ID - ��Ӧ����LW/isRota/yID/maxL  etc... ���������ת��
                for iLU = 1:length(uniID)
                    tmpM = LU.(fields{idx});
                    tmpVeh =tmpM(1:2,LU.ID==uniID(iLU))';    %��ȡ��ͳ�,��Ҫ��
                    if isempty(tmpVeh),   error('��Ӧ�ó��ֵĴ���');        end
                    if numel(unique(tmpVeh,'rows')) > 2 
                        LU.(fields{idx}) = aField;
                        warning('ԭ���̳���ߵȽ�����ת�� \n');    
                        break;
                    end
                end                
            else
                 warning('����3���������Ƿ���,�������ڷǳ����ʱ����');  
            end
            
        end
    end
    % �ж�3: ת�������������nLU    
    if size(LU.(fields{idx}),2)  ~= nLU
        size(LU.(fields{idx}),2)
        nLU
        error('LU���ݴ����������������������,���ڲ���ԭ�µĴ���2');   end
end 

end