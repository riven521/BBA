function [d] = getRandDa(n)
% ��ͨ����:���������������
      rng(1);
%      rng('default');
                                                ID = randi([11,18],1,n); %ID ������
    LWH = zeros(n,3);
    isRota = zeros(1,n);
    [uniID] = unique(ID);
    nUniId =  length(uniID);
    isRotaId = randi([0,1], 1, nUniId);
    for i=1: nUniId
    idx = find(ID(:)==uniID(i));
                                                LWH(idx,1) = randi([100,1000]);
                                                LWH(idx,2) = randi([100,1000]);
                                                isRota(1,idx) = isRotaId(i);
        for j=1:length(idx)
                                                LWH(idx(j),3) = randi([500,1000]);%randi([1000,1500]);
        end            
    end
    
    d.LU.ID = ID;
    d.LU.isRota = isRota; %isRota: �ж�ĳ����LU�Ƿ������ת,Ĭ�Ͼ�����ת,0��ʾ��������ת
    
    d.LU.LWH = LWH';
                                                d.Veh.LWH = [2400   5000   2400];%[2400   12750   2400];
    d.Veh.BUFF = [50;50;0];
    d.LU.BUFF = [100;100];%[2;2];
                    %     d.Veh.LWH = d.Veh.LWHREAL - d.Veh.BUFF;
                    %     d.LU.LWH = d.LU.LWHREAL + d.LU.BUFF;
    d.Veh.Weight = 1000;
    d.LU.Weight = d.LU.LWH(1,:);
end


% function [d] = getRandDa(n)
% % ��ͨ����:���������������
%      rng(222121);
% %      rng('default');
%                                                 ID = randi([11,13],1,n); %ID������
%     LWH = zeros(n,3);
%     [uniID] = unique(ID);
%     for i=1:length(uniID)
%     idx = find(ID(:)==uniID(i));
%                                                 LWH(idx,1) = randi([5,15]);
%                                                 LWH(idx,2) = randi([5,15]);
%         for j=1:length(idx)
%                                                 LWH(idx(j),3) = randi([1,8]);
%         end            
%     end
%     d.LU.ID = ID;
%     d.LU.LWH = LWH';
%                                                 d.Veh.LWH = [53; 55; 10];
%     d.Veh.BUFF = [0;0;0];
%     d.LU.BUFF = [0.1;0.1];%[2;2];
% %     d.Veh.LWH = d.Veh.LWHREAL - d.Veh.BUFF;
% %     d.LU.LWH = d.LU.LWHREAL + d.LU.BUFF;
%     d.Veh.Weight = 1000;
%     d.LU.Weight = d.LU.LWH(1,:);
% end