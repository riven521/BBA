function [LU] = setLULWHwithbuff(LU,Veh)
% setLULWHwithbuff ==> �޶�LU��LWH���� ��1������margin��2��ȷ��������תRotaed��
        
    % 1 UPDATE : ����LU��Rotaed���( LU��Ż����� ȴ���ھ���߽�Ŀ�� �ĸ�Сѡ�ĸ�
    if any(LU.Rotaed~=0)
        error('run�˺���ǰ������ת״̬��Ӧ�У��������º�����ȡ'); % ����У�����Ӧ��Ҳ������Ϊ�����ò����ò���
    end
    [LU.Rotaed]= getLuRotaedBBA(LU.LWH,LU.isRota,LU.margin,Veh.LWH(1,1));  %�ڶ���������  3��VEH�������Ұڷŵķ�϶��С����
    
    % 2 UPDATE : ��margin���Һ�Rotaed��LU��LWH   
    [LU.LWH,LU.OLWH] = LWHbuffer(LU.LWH, LU.margin, LU.Rotaed);   % V1:    LU.LWH = getRotaedLWH(LU.LWH, LU.Rotaed, LU.margin);

    % LU.OLWH = LWHunbuffer(LU.LWH, LU.margin, LU.Rotaed);


    % v1 ���ݣ�����������ʽȷ��LU�Ƿ���ת�� Ĭ�Ͻ�LUȫ������Horizontal������ת��ǰ�᣺��LU������ת��
    % NOTE: �˴������1: Horizontal�����LWH���Ƿ�Rotaed���
    % NOTE : ֱ���滻��ԭʼORIGINAL �� LWH

%     if pwhichSortItemOrder ==1
%         [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,1); %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��        
%     elseif pwhichSortItemOrder ==2
%         [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,0); %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��
%     elseif pwhichSortItemOrder ==3 %Ĭ�ϴ�ѡ��
%         [LU.Rotaed]= placeItemHori(LU.LWH,LU.isRota,Veh.LWH(1,1)); %�ڶ���������  3��VEH�������Ұڷŵķ�϶��С����
%     end
end



