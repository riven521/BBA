%% GCHECKINPUT ��Ҫ����:��������daת��+��������da�˶�
%% Form
%    [d] = GcheckInput(d,ParaArray)
%% Descripttion
%    �������ݼ���
%
%% Inputs
%    d                      (1,1)    �ṹ�壺�û���������
%    ParaArray          (1,1)    �ṹ�壺�㷨���Բ���
%
%% Outputs
%    d                     (1,1)
%

function [d] = GcheckInput(d,ParaArray)
%  ����:�����(row);  ����:��������(coloum);

initCheck();
inputExchange(); %333
inputLUIDExchange(); %555
deepCheck(); %333
getLUIDArray(); %% ���㣺LU�����������
printInput();

    %% Ƕ�׺���
    %% initCheck �������ݳ�������    
    function initCheck()
        % �����ж�����: d.Veh
        if numel(d.Veh.LWH(:))~=numel(d.Veh.BUFF(:)) || numel(d.Veh.LWH(:)) ~= 3 ||...
                ~isscalar(d.Veh.Weight)
            error('������da.BinArray����,����Ԥ�� \n');
        end
        if any(d.Veh.LWH(:)<=0) || any(d.Veh.Weight(:)<=0) || any(d.Veh.BUFF(:)<0)
            error('������da.BinArray�����쳣����,��������� \n');
        end
        % �����ж�����: d.LU
        if numel(d.LU.ID(:))~=numel(d.LU.Weight(:)) || numel(d.LU.BUFF(:)) ~= 2 || ...
                numel(d.LU.ID(:))*3~=numel(d.LU.LWH(:)) || ~ismatrix(d.LU.LWH)
                        %             numel(d.LU.ID(:))~=numel(d.LU.Weight(:)) 
                        %             numel(d.LU.BUFF(:)) ~= 2 
                        %             numel(d.LU.ID(:))*3~=numel(d.LU.LWH(:)) 
                        %             ~ismatrix(d.LU.LWH)           
            error('������da.LUArray����,����Ԥ�� \n');
        end
        if any(d.LU.LWH(:)<=0) || any(d.LU.Weight(:)<=0)  || any(d.LU.BUFF(:)<0) || ...
                any(d.LU.isRota(:)<0) || any(d.LU.isRota(:)>1)
            error('������da.BinArray�����쳣����,��������� \n');
        end
    end

    %% Input����ת��
    function inputExchange()
        % 1 ���������ж�ת��+�����ж�ת��
        % d.LU.ID,d.LU.Weight: ������row������        
        % d.Veh.LWH: ������column������
        % d.Veh.BUFF: ������column������
        if ~isrow(d.LU.ID),              d.LU.ID=d.LU.ID';              end
        if ~isrow(d.LU.Weight),              d.LU.Weight=d.LU.Weight';              end
        if ~iscolumn(d.Veh.LWH),    d.Veh.LWH=d.Veh.LWH';    end
        if ~iscolumn(d.Veh.BUFF),   d.Veh.BUFF=d.Veh.BUFF';   end        
        % d.LU.BUFF: ������column������,�̶�ת��Ϊmatrix����,����Ϊ��������,����Ϊ3;
        if ~iscolumn(d.LU.BUFF),    d.LU.BUFF=d.LU.BUFF';    end
        % d.LU.LWH: ������matrix����,����Ϊ��������,����Ϊ3; δ����3��3������
        if size(d.LU.LWH,1)~=3,     d.LU.LWH=d.LU.LWH';     end
        
        % ������û�и���ֵ
        if ~isfield(d.LU,'isRota')
            d.LU.isRota = ones(size(d.LU.ID));
        end
            
        % ��whichRotation��0��1ʱ���ж�, ��Ϊ2ʱ���򲻱�
        if ParaArray.whichRotation == 1 % ȫ��������ת
            d.LU.isRota = ones(size(d.LU.isRota));
        elseif ParaArray.whichRotation == 0    % ȫ����ֹ��ת
            d.LU.isRota = zeros(size(d.LU.isRota));
        elseif ParaArray.whichRotation == 2    % ����������ת
            1;
        end

        % 2 Input���Ӽ�϶BUFF���feasible��LU��BIN�ĳ����ת��
        d.Veh.LWH = d.Veh.LWH - d.Veh.BUFF;
        d.LU.BUFF = [d.LU.BUFF; 0]; %�û�������BuffΪÿ���������ӵĳߴ�(�ܼ�϶)
        d.LU.BUFF = repmat(d.LU.BUFF,1,numel(d.LU.ID));
        d.LU.LWH = d.LU.LWH + d.LU.BUFF;
        %TODO: �˴����Ӽ�϶ΪȨ��֮�ʣ�δ����rotation��ı仯�����ڿ������㷨�����Ӽ�϶
        
        % 3 Ĭ�Ͻ�LUȫ������Horizontal������ת��ǰ�᣺��LU������ת��
        % NOTE: �˴������1: Horizontal�����LWH���Ƿ�Rotaed���
        % NOTE : ֱ���滻��ԭʼORIGINAL �� LWH
        [d.LU.Rotaed]= placeItemHori(d.LU.LWH,d.LU.isRota,1); %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��
        d.LU.LWH = getRotaedLWH(d.LU.LWH, d.LU.Rotaed, d.LU.BUFF); 
        
        
        % 3 Input����LUArray�Ƿ�Rotation���ж�flag TOBE DELE
%         d.LU.RotaFlag = zeros(1,numel(d.LU.ID)); TOBE DELE
    end

    %% LUID����ת��: d.LU.ID ת��Ϊ��1��ʼ������� ������������ID��Ϣ
    function inputLUIDExchange()
        uniLUID = unique(d.LU.ID);
        nLUid = size(uniLUID,2);
        tmpLUID=d.LU.ID; %�м����
        for i=1:nLUid
            tmpLUID(d.LU.ID(:)==uniLUID(i)) = i;
        end
        d.LU.ID=tmpLUID;
    end

    %% �����ж�
    function deepCheck()        
        % ���������ȡ
        dBin = unique(d.Veh.LWH','rows')'; %��ȡ����ͬ��������
         dLU = d.LU.LWH;
        dLUisRota = d.LU.isRota;
        dLUid = d.LU.ID;
        nLUid = size(unique(dLUid),2);
        % �߶�Լ�� ����
        if any(dBin(3) < dLU(3,:)),     error('����: �������̸߶� ���� �����Ϳ��ø߶����� \n');     end
        
        % ����Լ�� ���� (���̳������)
        flag1 = dLUisRota == 1;  %����rota��LU���
        flag0 = dLUisRota == 0;  %������rota��LU���
        % ���Բ�������ת��LU
        if any(dBin(1) < dLU(1,flag0)),  error('����: �������̿�� ���� �����Ϳ��ÿ������ \n');    end
        if any(dBin(2) < dLU(2,flag0)),  error('����: �������̳��� ���� �����Ϳ��ó������� \n');    end
        % ����������ת��LU
        if ~all(~flag1) && min(dBin(1:2)) < max(min(dLU(1:2, flag1)))  %flag1����ȫ0 �������������
            error('����: ����������̱ߵ����ֵ ���� �����ͳ���� \n');
        end
        if ~all(~flag1) && max(dBin(1:2)) < max(max(dLU(1:2, flag1 ))) %flag1����ȫ0 �������������
            error('����: ����������ߵ����ֵ ���� �����ͳ���� \n');
        end

                            % ����Լ�� ���� (���̳������) TOBE DELE
                    %         if ParaArray.whichRotation == 1
                    %             if min(dBin(1:2)) < max(min(dLU(1:2,:))) || max(dBin(1:2)) < max(max(dLU(1:2,:)))
                    %                 error('����: �������̳���� ���� �����ͳ���� \n');
                    %             end
                    %         else %��׼rotation
                    %             if any(dBin(1) < dLU(1,:)),  error('����: �������̿�� ���� �����Ϳ��ÿ������ \n');    end
                    %             if any(dBin(2) < dLU(2,:)),  error('����: �������̳��� ���� �����Ϳ��ó������� \n');    end
                    %         end

        % ����ID�����ʣ����ͣ���Ӧ����Լ��  ����555
        for iLU = 1:nLUid
            tmp = dLU(1:2,dLUid==iLU)';   %��ȡ��ͳ�,��Ҫ��
            if numel(unique(tmp,'rows')) > 2,  error('����: ��������ID������ͬ ���䳤��ͬ���� \n');    end
        end
    end
    
    %%  ��ȡLUID�����������(ͬ����ID��������������Ƿ����ת)
    function getLUIDArray()
%         printstruct(d);
        d.LUID.ID = unique(d.LU.ID);
        nLUID = numel(d.LUID.ID);        

        for iID = 1:nLUID
        d.LUID.Weight(iID) = sum(d.LU.Weight .* (d.LU.ID == d.LUID.ID(iID)) );
        d.LUID.Volume(iID) = sum(prod(d.LU.LWH) .* (d.LU.ID == d.LUID.ID(iID)) );
        d.LUID.isRota(iID) =  unique(d.LU.isRota(d.LU.ID == d.LUID.ID(iID))); if ~isscalar(d.LUID.isRota(iID)), error('��������'); end
        end

    end

   %% ��ӡInput����
    function printInput()
        fprintf('������ֻ��һ������ ��=%1.0f ��=%1.0f ��=%1.0f \n', unique(d.Veh.LWH','rows')');
        fprintf('������ֻ��һ������ ���϶=%1.0f ����϶=%1.0f �߼�϶=%1.0f \n', unique(d.Veh.BUFF','rows')');
        fprintf('�������� %d ����Ʒ,����߷ֱ�Ϊ \n',numel(d.LU.ID(:)));
        fprintf('%1.1f %1.1f %1.1f \n',d.LU.LWH);
    end

end





