function [LU,Veh] = cpuLULWH(LU,Veh)
% cpuLULWH ��Ҫ����:
% ���� LU ��LWH
        
    %% 3 ����LU����margin��LU��LWH������LU��Rotaed���
    % GET LU's LWH with margin and Rotaed or not
    % 2 Input���Ӽ�϶BUFF���feasible��LU��BIN�ĳ����ת��
    LU.LWH(1,:) =  LU.LWH(1,:) +  LU.margin(1,: ) + LU.margin(2,: ); %��ȣ����ң�
    LU.LWH(2,:) =  LU.LWH(2,:) +  LU.margin(3,: ) + LU.margin(4,: ); %���ȣ����£�
    
    [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,Veh.LWH(1,1));  %�ڶ���������  3��VEH�������Ұڷŵķ�϶��С����
    LU.LWH = getRotaedLWH(LU.LWH, LU.Rotaed, LU.margin);
    
    % 3 Ĭ�Ͻ�LUȫ������Horizontal������ת��ǰ�᣺��LU������ת��
    % NOTE: �˴������1: Horizontal�����LWH���Ƿ�Rotaed���
    % NOTE : ֱ���滻��ԭʼORIGINAL �� LWH
    % LU.LWH

%     if pwhichSortItemOrder ==1
%         [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,1); %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��        
%     elseif pwhichSortItemOrder ==2
%         [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,0); %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��
%     elseif pwhichSortItemOrder ==3 %Ĭ�ϴ�ѡ��
%         [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,Veh.LWH(1,1)); %�ڶ���������  3��VEH�������Ұڷŵķ�϶��С����
%     end
end
    