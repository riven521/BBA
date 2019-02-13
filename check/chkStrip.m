%% chkStrip Strip�ṹ��˶�
%% Form
%    chkStrip(Strip)

function chkStrip(Strip)
     if ~isstruct(Strip), error('NOT a STRUCT!'); end
     if isstruct(Strip)
     T = struct2table(structfun(@(x) x', Strip,'UniformOutput',false));  end
    
     % ��ά���жϣ�Strip������ Strip��LID���ܳ�2��
     % if any(~iscell(T.LID(fmix))),   error('���Strip��cell��LID������Ϊ2'); end
     fmix = T.isMixed == 1;
     if any(T.nbLID(fmix) >=3) || any(T.nbLID(fmix) <= 1)
         error('���Strip��LID������Ϊ2'); end     
     if any(T.nbLID(~fmix) ~= 1)
         error('�ǻ��Strip��LID������Ϊ1'); end
     
     % �жϣ�Width Not Full �п���ʱ Mixed
     fwid = T.isWidthFull == 0;
     if any(fwid&fmix), 
         fwid
         fmix
         warning('��Ȳ�����Strip����Ϊ���Strip'); end

end

% TO �޸�ʹ��
function deepCheck(LU,Veh)
    % �߶�Լ�� ����
    if min(Veh.LWH(3, :)) < max(LU.LWH(3, : )),     error('����: �������̸߶� ���� �����Ϳ��ø߶����� \n');     end

    % ����Լ�� ���� (���̳������)
    flag1 = LU.isRota == 1;  %����rota��LU���
    flag0 = LU.isRota == 0;  %������rota��LU���
    % ���Բ�������ת��LU
    if min(Veh.LWH(1, :)) < max(LU.LWH(1,flag0)), error('����: �������̿�� ���� �����Ϳ��ÿ������ \n');    end
    if min(Veh.LWH(2, :)) < max(LU.LWH(2,flag0)), error('����: �������̳��� ���� �����Ϳ��ó������� \n');    end
    % ����������ת��LU
    if ~all(~flag1) && min(Veh.LWH(1:2)) < max(min(LU.LWH(1:2, flag1)))  %flag1����ȫ0 �������������
        error('����: ����������̱ߵ����ֵ ���� �����ͳ���� \n');
    end
    if ~all(~flag1) && max(Veh.LWH(1:2)) < max(max(LU.LWH(1:2, flag1 ))) %flag1����ȫ0 �������������
        error('����: ����������ߵ����ֵ ���� �����ͳ���� \n');
    end
end