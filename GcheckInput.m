function [da] = GcheckInput(da,ParaArray)
% �ȶ��������ݽ���ת��,������ݽ���check
% ����ת��
    if size(da.LUArray.ID,1)~=1 && size(da.LUArray.ID,1)>size(da.LUArray.ID,2) %������
        da.LUArray.ID=da.LUArray.ID';
    end
    if size(da.BinArray.LWH,2)~=1 && size(da.BinArray.LWH,2)>size(da.BinArray.LWH,1) %������
        da.BinArray.LWH=da.BinArray.LWH';
    end
    if size(da.LUArray.LWH,2)~=3  %ȷ���û���������ά������;
        da.LUArray.LWH=da.LUArray.LWH';        
    end
    da.LUArray.LWH = da.LUArray.LWH';  %ת����ά������Ϊ��ά������
% �����϶��ת��
    if size(da.BinArray.BUFF,2)~=1  %BinArray.BUFF����Ϊ������
        da.BinArray.BUFF=da.BinArray.BUFF';
    end
    if size(da.LUArray.BUFF,2)~=1   %LUArray.BUFF����Ϊ������
        da.LUArray.BUFF=da.LUArray.BUFF';
    end

    
% ����ת��: da.LUArray.ID ת��Ϊ��1��ʼ������� ������������ID��Ϣ
    nLUid = size(unique(da.LUArray.ID),2);
    uniLUID = unique(da.LUArray.ID);
    tmpLUID=da.LUArray.ID; %�м����
    for i=1:nLUid
        tmpLUID(da.LUArray.ID(:)==uniLUID(i)) = i;
    end
    da.LUArray.ID=tmpLUID;
% ����ת��2: ���Ӽ�϶���ת��
    da.LUArray.BUFF = [da.LUArray.BUFF;0];
    da.BinArray.LWH = da.BinArray.LWH - da.BinArray.BUFF;
    da.LUArray.LWH = da.LUArray.LWH + da.LUArray.BUFF;
    

% ����Check
    nDim = size(da.LUArray.LWH,1); 
%     for i=1:nDim
%         fprintf('%1.0f', da.LUArray.LWH(i,:));
%         fprintf('\n');
%     end
    if nDim ~=3, error('����������ά����,����Ԥ�� \n'); end
    dLU = da.LUArray.LWH(1:nDim,:);
    dLUid = da.LUArray.ID;
    nLUid = size(unique(dLUid),2);
    tmpbinDataMatrix = da.BinArray.LWH(1:nDim,:);
    dBin = unique(tmpbinDataMatrix','rows')';
    nItem = size(dLU,2);    
    if size(dBin,2)==1
        fprintf('������ֻ��һ������ ��=%1.0f ��=%1.0f ��=%1.0f \n', dBin);
        fprintf('�������� %d ����Ʒ,����߷ֱ�Ϊ \n',nItem);
        fprintf('%1.0f %1.0f %1.0f \n',dLU);
    else
        error('�������ж������,�������� \n');
    end
    if numel(da.BinArray.BUFF)~=3, error('���������ͼ�϶��������3��,�������� \n'); end
    if numel(da.LUArray.BUFF)~=3, error('���������ͼ�϶��������3��,�������� \n'); end

    %% �߶�Լ��
    if any(dBin(3) < dLU(3,:))
        error('����: �������̸߶� ���� �����Ϳ��ø߶����� \n');
    end 
    %% ����Լ��
    if ParaArray.whichRotation == 1
        if min(dBin(1:2)) < max(min(dLU(1:2,:))) || max(dBin(1:2)) < max(max(dLU(1:2,:)))
            error('����: �������̳���� ���� �����ͳ���� \n');
        end
    else %��׼rotation 
        if any(dBin(1) < dLU(1,:))
            error('����: �������̿�� ���� �����Ϳ��ÿ������ \n');
        end
        if any(dBin(2) < dLU(2,:))
            error('����: �������̳��� ���� �����Ϳ��ó������� \n');
        end
    end
    %% ����ID�����ʣ����ͣ��볤��Լ��  
    for iLU = 1:nLUid
%         dLULW = dLU(1:2,:); %��ȡ��ͳ� ��Ҫ��
        tmp = dLU(1:2,dLUid==iLU)'; %��ȡ
        if numel(unique(tmp,'rows')) > 2
            error('����: ��������ID��ͬ ���䳤��ͬ���� \n');
        end
    end    
end





