%% GCHECKINPUT ��Ҫ����:��������daת��+��������da�˶�
%% Form
%    [da] = GcheckInput(da,ParaArray)
%% Descripttion
%    �������ݼ���
%
%% Inputs
%    da                      (1,1)    �ṹ�壺�û���������
%    ParaArray          (1,1)    �ṹ�壺�㷨���Բ���
%
%% Outputs
%    da                     (1,1)
%

function [da] = GcheckInput(da,ParaArray)
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
        % �����ж�����: da.BinArray
        if numel(da.BinArray.LWH(:))~=numel(da.BinArray.BUFF(:)) || numel(da.BinArray.LWH(:)) ~= 3 ||...
                ~isscalar(da.BinArray.Weight)
            error('������da.BinArray����,����Ԥ�� \n');
        end
        if any(da.BinArray.LWH(:)<=0) || any(da.BinArray.Weight(:)<=0) || any(da.BinArray.BUFF(:)<0)
            error('������da.BinArray�����쳣����,��������� \n');
        end
        % �����ж�����: da.LUArray
        if numel(da.LUArray.ID(:))~=numel(da.LUArray.Weight(:)) || numel(da.LUArray.BUFF(:)) ~= 2 || ...
                numel(da.LUArray.ID(:))*3~=numel(da.LUArray.LWH(:)) || ~ismatrix(da.LUArray.LWH)
                        %             numel(da.LUArray.ID(:))~=numel(da.LUArray.Weight(:)) 
                        %             numel(da.LUArray.BUFF(:)) ~= 2 
                        %             numel(da.LUArray.ID(:))*3~=numel(da.LUArray.LWH(:)) 
                        %             ~ismatrix(da.LUArray.LWH)           
            error('������da.LUArray����,����Ԥ�� \n');
        end
        if any(da.LUArray.LWH(:)<=0) || any(da.LUArray.Weight(:)<=0)  || any(da.LUArray.BUFF(:)<0) || ...
                any(da.LUArray.isRota(:)<0) || any(da.LUArray.isRota(:)>1)
            error('������da.BinArray�����쳣����,��������� \n');
        end
    end

    %% Input����ת��
    function inputExchange()
        % 1 ���������ж�ת��+�����ж�ת��
        % da.LUArray.ID,da.LUArray.Weight: ������row������        
        % da.BinArray.LWH: ������column������
        % da.BinArray.BUFF: ������column������
        if ~isrow(da.LUArray.ID),              da.LUArray.ID=da.LUArray.ID';              end
        if ~isrow(da.LUArray.Weight),              da.LUArray.Weight=da.LUArray.Weight';              end
        if ~iscolumn(da.BinArray.LWH),    da.BinArray.LWH=da.BinArray.LWH';    end
        if ~iscolumn(da.BinArray.BUFF),   da.BinArray.BUFF=da.BinArray.BUFF';   end        
        % da.LUArray.BUFF: ������column������,�̶�ת��Ϊmatrix����,����Ϊ��������,����Ϊ3;
        if ~iscolumn(da.LUArray.BUFF),    da.LUArray.BUFF=da.LUArray.BUFF';    end
        % da.LUArray.LWH: ������matrix����,����Ϊ��������,����Ϊ3; δ����3��3������
        if size(da.LUArray.LWH,1)~=3,     da.LUArray.LWH=da.LUArray.LWH';     end
        
        % ������û�и���ֵ
        if ~isfield(da.LUArray,'isRota')
            da.LUArray.isRota = ones(size(da.LUArray.ID));
        end
            
        % ��whichRotation��0��1ʱ���ж�, ��Ϊ2ʱ���򲻱�
        if ParaArray.whichRotation == 1 % ȫ��������ת
            da.LUArray.isRota = ones(size(da.LUArray.isRota));
        elseif ParaArray.whichRotation == 0    % ȫ����ֹ��ת
            da.LUArray.isRota = zeros(size(da.LUArray.isRota));
        elseif ParaArray.whichRotation == 2    % ����������ת
            1;
        end

        % 2 Input���Ӽ�϶BUFF���feasible��LU��BIN�ĳ����ת��
        da.BinArray.LWH = da.BinArray.LWH - da.BinArray.BUFF;
        da.LUArray.BUFF = [da.LUArray.BUFF; 0]; %�û�������BuffΪÿ���������ӵĳߴ�(�ܼ�϶)
        da.LUArray.BUFF = repmat(da.LUArray.BUFF,1,numel(da.LUArray.ID));
        da.LUArray.LWH = da.LUArray.LWH + da.LUArray.BUFF;
        %TODO: �˴����Ӽ�϶ΪȨ��֮�ʣ�δ����rotation��ı仯�����ڿ������㷨�����Ӽ�϶
        
        % 3 Ĭ�Ͻ�LUȫ������Horizontal������ת��ǰ�᣺��LU������ת��
        % NOTE: �˴������1: Horizontal�����LWH���Ƿ�Rotaed���
        % NOTE : ֱ���滻��ԭʼORIGINAL �� LWH
        [da.LUArray.Rotaed]= placeItemHori(da.LUArray.LWH,da.LUArray.isRota,1); %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��
        da.LUArray.LWH = getRotaedLWH(da.LUArray.LWH, da.LUArray.Rotaed, da.LUArray.BUFF); 
        
        
        % 3 Input����LUArray�Ƿ�Rotation���ж�flag TOBE DELE
%         da.LUArray.RotaFlag = zeros(1,numel(da.LUArray.ID)); TOBE DELE
    end

    %% LUID����ת��: da.LUArray.ID ת��Ϊ��1��ʼ������� ������������ID��Ϣ
    function inputLUIDExchange()
        uniLUID = unique(da.LUArray.ID);
        nLUid = size(uniLUID,2);
        tmpLUID=da.LUArray.ID; %�м����
        for i=1:nLUid
            tmpLUID(da.LUArray.ID(:)==uniLUID(i)) = i;
        end
        da.LUArray.ID=tmpLUID;
    end

    %% �����ж�
    function deepCheck()        
        % ���������ȡ
        dBin = unique(da.BinArray.LWH','rows')'; %��ȡ����ͬ��������
         dLU = da.LUArray.LWH;
        dLUisRota = da.LUArray.isRota;
        dLUid = da.LUArray.ID;
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
%         printstruct(da);
        da.LUIDArray.ID = unique(da.LUArray.ID);
        nLUID = numel(da.LUIDArray.ID);        

        for iID = 1:nLUID
        da.LUIDArray.Weight(iID) = sum(da.LUArray.Weight .* (da.LUArray.ID == da.LUIDArray.ID(iID)) );
        da.LUIDArray.Volume(iID) = sum(prod(da.LUArray.LWH) .* (da.LUArray.ID == da.LUIDArray.ID(iID)) );
        da.LUIDArray.isRota(iID) =  unique(da.LUArray.isRota(da.LUArray.ID == da.LUIDArray.ID(iID))); if ~isscalar(da.LUIDArray.isRota(iID)), error('��������'); end
        end

    end

   %% ��ӡInput����
    function printInput()
        fprintf('������ֻ��һ������ ��=%1.0f ��=%1.0f ��=%1.0f \n', unique(da.BinArray.LWH','rows')');
        fprintf('������ֻ��һ������ ���϶=%1.0f ����϶=%1.0f �߼�϶=%1.0f \n', unique(da.BinArray.BUFF','rows')');
        fprintf('�������� %d ����Ʒ,����߷ֱ�Ϊ \n',numel(da.LUArray.ID(:)));
        fprintf('%1.1f %1.1f %1.1f \n',da.LUArray.LWH);
    end

end





