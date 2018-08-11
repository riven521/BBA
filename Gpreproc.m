%% GPREPROC ��Ҫ����:��������daת��+��������da�˶�
%% Form
%    [d] = Gpreproc(d)
%% Descripttion
%    ��������Ԥ����
%
%% Inputs
%    d                      (1,1)    �ṹ�壺�û���������
%
%% Outputs
%    d                     (1,1)
%

function [LU,Veh] = Gpreproc(LU,Veh,pwhichSortItemOrder)
    % 1 ID ת��Ϊ��1��ʼ������� ������������ID����Ϣ
    if isfield(LU, 'ID'),  LU.ID = idExchange(LU.ID); end
    if isfield(LU, 'PID'),  LU.PID = idExchange(LU.PID); end
    if isfield(LU, 'PID'),  LU.SID = idExchange(LU.SID); end
    if isfield(LU, 'PID'),  LU.UID = idExchange(LU.UID); end
    
    % 2 Input���Ӽ�϶BUFF���feasible��LU��BIN�ĳ����ת��
% %     Veh.LWH = Veh.LWH - Veh.buff;
% %     LU.LWH(1,:) =  LU.LWH(1,:) +  LU.buff(1,: ) + LU.buff(2,: );
% %     LU.LWH(2,:) =  LU.LWH(2,:) +  LU.buff(3,: ) + LU.buff(4,: );
    
    % V1 : buff: ���̼�ļ�϶
%     LU.buff = [LU.buff; 0]; %�û�������BuffΪÿ���������ӵĳߴ�(�ܼ�϶)
%     LU.buff = repmat(LU.buff,1,numel(LU.ID));
    LU.LWH = LU.LWH + LU.buff;
    %TODO: �˴����Ӽ�϶ΪȨ��֮�ʣ�δ����rotation��ı仯�����ڿ������㷨�����Ӽ�϶

    % 3 Ĭ�Ͻ�LUȫ������Horizontal������ת��ǰ�᣺��LU������ת��
    % NOTE: �˴������1: Horizontal�����LWH���Ƿ�Rotaed���
    % NOTE : ֱ���滻��ԭʼORIGINAL �� LWH
%        LU.LWH

    if pwhichSortItemOrder ==1
        [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,1); %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��        
    elseif pwhichSortItemOrder ==2
        [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,0); %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��
    elseif pwhichSortItemOrder ==3
        [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,Veh.LWH(3,1)); %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��  3����϶��С����      
    end
     LU.LWH = getRotaedLWH(LU.LWH, LU.Rotaed, LU.buff);
%         LU.LWH
        
    % 4 Veh�������->С   Ĭ��˳��
    [~,order] = sortrows(Veh.volume', [1],{'descend'});    
    Veh = structfun(@(x) x(:,order),Veh,'UniformOutput',false);
    
%     LUID = getLUIDArray(LU); %% ���㣺LU����������� ��ʱ����

     printInput();
%% Ƕ�׺���:
function printInput()
    fprintf('������ֻ��һ������ ��=%1.0f ��=%1.0f ��=%1.0f \n', unique(Veh.LWH','rows')');
    fprintf('������ֻ��һ������ ���϶=%1.0f ����϶=%1.0f �߼�϶=%1.0f \n', unique(Veh.buff','rows')');
    fprintf('�������� %d ����Ʒ,����߷ֱ�Ϊ \n',numel(LU.ID(:)));
    fprintf('%1.1f %1.1f %1.1f \n',LU.LWH);
end

end


function exID = idExchange(ID)
        uniID = unique(ID);
        exID=ID; %�м����
        for i=1:length(uniID)
            exID(ID(:)==uniID(i)) = i;
        end
end


%%  ��ȡLUID�����������(ͬ����ID��������������Ƿ����ת)
function LUID = getLUIDArray(LU)
LUID.ID = unique(LU.ID);
for iID = 1:length(LUID.ID)
    LUID.Weight(iID) = sum(LU.Weight .* (LU.ID == LUID.ID(iID)) );
    LUID.Volume(iID) = sum(prod(LU.LWH) .* (LU.ID == LUID.ID(iID)) );
    LUID.isRota(iID) =  unique(LU.isRota(LU.ID == LUID.ID(iID)));
    LUID.maxL(iID) =  unique(LU.maxL(LU.ID == LUID.ID(iID)));
    LUID.yID(iID) =  unique(LU.yID(LU.ID == LUID.ID(iID)));
    if ~isscalar(LUID.isRota(iID))||~isscalar(LUID.maxL(iID))||~isscalar(LUID.yID(iID)), error('��������'); end
end
end
        

    %% 
% %     function inputExchange()
% % % %         % 1 ���������ж�ת��+�����ж�ת��
% % % %         % d.LU.ID,d.LU.Weight: ������row������        
% % % %         % d.Veh.LWH: ������column������
% % % %         % d.Veh.BUFF: ������column������
% % % %         if ~isrow(d.LU.ID),              d.LU.ID=d.LU.ID';              end
% % % %         if ~isrow(d.LU.Weight),              d.LU.Weight=d.LU.Weight';              end
% % % %         if ~iscolumn(d.Veh.LWH),    d.Veh.LWH=d.Veh.LWH';    end
% % % %         if ~iscolumn(d.Veh.BUFF),   d.Veh.BUFF=d.Veh.BUFF';   end        
% % % %         % d.LU.BUFF: ������column������,�̶�ת��Ϊmatrix����,����Ϊ��������,����Ϊ3;
% % % %         if ~iscolumn(d.LU.BUFF),    d.LU.BUFF=d.LU.BUFF';    end
% % % %         % d.LU.LWH: ������matrix����,����Ϊ��������,����Ϊ3; δ����3��3������
% % % %         if size(d.LU.LWH,1)~=3,     d.LU.LWH=d.LU.LWH';     end
% %         
% %         % ������û�и���ֵ
% % %         if ~isfield(d.LU,'isRota')
% % %             d.LU.isRota = ones(size(d.LU.ID));
% % %         end
% %             
% %         % ��whichRotation��0��1ʱ���ж�, ��Ϊ2ʱ���򲻱�
% % %         if ParaArray.whichRotation == 1 % ȫ��������ת
% % %             d.LU.isRota = ones(size(d.LU.isRota));
% % %         elseif ParaArray.whichRotation == 0    % ȫ����ֹ��ת
% % %             d.LU.isRota = zeros(size(d.LU.isRota));
% % %         elseif ParaArray.whichRotation == 2    % ����������ת
% % %             1;
% % %         end
% % 
% %         % 2 Input���Ӽ�϶BUFF���feasible��LU��BIN�ĳ����ת��
% %         d.Veh.LWH = d.Veh.LWH - d.Veh.BUFF;
% %         d.LU.BUFF = [d.LU.BUFF; 0]; %�û�������BuffΪÿ���������ӵĳߴ�(�ܼ�϶)
% %         d.LU.BUFF = repmat(d.LU.BUFF,1,numel(d.LU.ID));
% %         d.LU.LWH = d.LU.LWH + d.LU.BUFF;
% %         %TODO: �˴����Ӽ�϶ΪȨ��֮�ʣ�δ����rotation��ı仯�����ڿ������㷨�����Ӽ�϶
% %         
% %         % 3 Ĭ�Ͻ�LUȫ������Horizontal������ת��ǰ�᣺��LU������ת��
% %         % NOTE: �˴������1: Horizontal�����LWH���Ƿ�Rotaed���
% %         % NOTE : ֱ���滻��ԭʼORIGINAL �� LWH
% %         [d.LU.Rotaed]= placeItemHori(d.LU.LWH,d.LU.isRota,1); %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��
% %         d.LU.LWH = getRotaedLWH(d.LU.LWH, d.LU.Rotaed, d.LU.BUFF); 
% %         
% % end

    