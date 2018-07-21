function [da] = GcheckInput(da,ParaArray)
% ��Ҫ����:��������daת��+��������da�˶�
% Input --- 
% Output ---
%  ����:�����(row); 
%  ����:��������(coloum);
% initCheck();
% inputExchange();
% inputLUIDExchange();
% deepCheck();
% printInput();

initCheck();
inputExchange(); %333
inputLUIDExchange(); %555
deepCheck(); %333
printInput();

    %% �����ж�
    function initCheck()
        % �����ж�����: da.BinArray
        if numel(da.BinArray.LWH(:))~=numel(da.BinArray.BUFF(:)) || numel(da.BinArray.LWH(:)) ~= 3 || ~isscalar(da.BinArray.Weight)
            error('������da.BinArray����,����Ԥ�� \n');
        end
        % �����ж�����: da.LUArray
        if numel(da.LUArray.ID(:))~=numel(da.LUArray.Weight(:)) || numel(da.LUArray.BUFF(:)) ~= 2 || numel(da.LUArray.ID(:))*3~=numel(da.LUArray.LWH(:)) || ~ismatrix(da.LUArray.LWH)
            error('������da.LUArray����,����Ԥ�� \n');
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
        
        % 2 Input���Ӽ�϶BUFF���feasible��LU��BIN�ĳ����ת��
        da.BinArray.LWH = da.BinArray.LWH - da.BinArray.BUFF;
        da.LUArray.BUFF = [da.LUArray.BUFF; 0]; %�û�������BuffΪÿ���������ӵĳߴ�(�ܼ�϶)
        da.LUArray.BUFF = repmat(da.LUArray.BUFF,1,numel(da.LUArray.ID));
        da.LUArray.LWH = da.LUArray.LWH + da.LUArray.BUFF;
        %TODO: �˴����Ӽ�϶ΪȨ��֮�ʣ�δ����rotation��ı仯�����ڿ������㷨�����Ӽ�϶
        
        % 3 Input����LUArray�Ƿ�Rotation���ж�flag
        da.LUArray.RotaFlag = zeros(1,numel(da.LUArray.ID));
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
        dLUid = da.LUArray.ID;
        nLUid = size(unique(dLUid),2);
        % �߶�Լ�� ����
        if any(dBin(3) < dLU(3,:)),     error('����: �������̸߶� ���� �����Ϳ��ø߶����� \n');     end
        % ����Լ�� ����
        if ParaArray.whichRotation == 1
            if min(dBin(1:2)) < max(min(dLU(1:2,:))) || max(dBin(1:2)) < max(max(dLU(1:2,:)))
                error('����: �������̳���� ���� �����ͳ���� \n');
            end
        else %��׼rotation
            if any(dBin(1) < dLU(1,:)),  error('����: �������̿�� ���� �����Ϳ��ÿ������ \n');    end
            if any(dBin(2) < dLU(2,:)),  error('����: �������̳��� ���� �����Ϳ��ó������� \n');    end
        end
        % ����ID�����ʣ����ͣ���Ӧ����Լ��  ����555
        for iLU = 1:nLUid
            tmp = dLU(1:2,dLUid==iLU)';   %��ȡ��ͳ�,��Ҫ��
            if numel(unique(tmp,'rows')) > 2,  error('����: ��������ID������ͬ ���䳤��ͬ���� \n');    end
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





