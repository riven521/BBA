function [d] = getRandDa(n,m)
%     rng(10000);
%     rng('default');
    s = rng;     save('srng','s');
%     load('srng');    rng(s);
    
    Par.maxHeavey = 18;
%     Par.LUbuff = [0;0]; %�Ժ�����
    
%% Veh
    Veh.ID = randperm(m);   %Veh ������
    Veh.LWH = zeros(3,m);
    
%     maxSeg = 1; %�������� 3
%     Veh.yID = zeros(maxSeg,m);
%     Veh.xID = zeros(1,m);
%     Veh.buff = zeros(3,m);

    for i=1: m
        Veh.Weight(i) = randi([1000,1500]);    %���� randi([1000,1500]);
        Veh.LWH(1,i) = randi([2400,2400]); %��� randi([2200,2400]);
        Veh.LWH(2,i) = randi([12750,12750]); %��/Height��  randi([5000,6000]);
        Veh.LWH(3,i) = randi([2400,2400]); %�߶�  randi([2200,2400]); 
        
%         Veh.LWH(1,2) = 1800; %��� randi([2200,2400]);
%         Veh.LWH(2,2) = 7050; %��/Height��  randi([5000,6000]);
%         Veh.LWH(3,2) = 1800;%�߶�  randi([2200,2400]); 

%         nYID = randi([1,maxSeg]);
%         nXID = 1;
%         for j=1:nYID
%             if j==1
%                 Veh.yID(j,i) = round(Veh.LWH(2,i)/nYID);
%             elseif j==nYID
%                 Veh.yID(j,i) = Veh.LWH(2,i);
%             else
%                 Veh.yID(j,i) = Veh.yID(j-1,i)+ round(Veh.LWH(2,i)/nYID);
%             end
%         end
%         for j=1:nXID
%             Veh.xID(j,i) = Veh.LWH(2,i);
%         end
%     end

    % remove extra rows
%     emptyYID =  all(Veh.yID ==0,2) == true;
%     Veh.yID(emptyYID,:) = [];    
    
    
%% LU   
%     LU.buff = zeros(3,n); %�Ժ�����
%     LU.maxL = zeros(1,n);
%     LU.yID = zeros(1,n);
%     LU.xID = zeros(1,n);
%     LU.UID = zeros(1,n);
%     LU.isH = zeros(1,n);

    LU.ID = randi([11,13],1,n); %ID ������ 
    LU.LID = LU.ID;

    LU.LWH = zeros(3,n);
    LU.isRota = zeros(1,n);

    LU.margin = zeros(4,n); 

    LU.PID = zeros(1,n);
    LU.SID = zeros(1,n);
    LU.Weight = zeros(1,n);


    [uniID] = unique(LU.ID);    %   isRotaId = randi([0,1], 1, nUniId);
    for i=1: length(uniID)
        % ��LUID���
        idx = find(LU.ID(:)==uniID(i));       
        LU.LWH(1,idx) = randi([1486-200,1486+200]); %��� randi([1000,1300]);
        LU.LWH(2,idx) = randi([765-200,765+200]); %���� randi([700,850]); 
        LU.isRota(idx) = randi([1,1]);            %�Ƿ���ת
%         LU.HightL(idx) = randi([1,1]);       %���ɶѶ����
%         LU.maxL(idx) = randi([1,4]);
%         LU.yID(idx) = randi([0,nYID]); % TODO ���ڸ��복��һ��
%         LU.xID(idx) = randi([0,nXID]);
%         LU.UID(idx(j)) = randi([300,300]);
        LU.margin(:,idx) = randi([50,50]);  %��������
        for j=1:length(idx)
            LU.LWH(3,idx(j)) = randi([852-50,852+50]); %�߶�250,1150
            LU.PID(idx(j)) = randi([100,100]); %100,103
            LU.SID(idx(j)) = randi([200, 200]); %200, 203
            LU.EID(idx(j)) = randi([200, 200]); %200, 203
            
            
            LU.Weight(idx(j)) = randi([10,50]);
            
%             if LU.Weight(idx(j)) >= Par.maxHeavey
%                 LU.isH(idx(j)) = 1;
%             end
        end
    end

%%
    d.LU = LU;
    d.Par = Par;
    d.Veh = Veh;
    
end

%     d.LU.ID = ID;
%     d.LU.isRota = isRota; %isRota: �ж�ĳ����LU�Ƿ������ת,Ĭ�Ͼ�����ת,0��ʾ��������ת    
%     d.LU.LWH = LWH';
%         d.LU.buff = [2;2]; 
%         d.LU.Weight = d.LU.LWH(1,:);

    %     Veh.LWH = [2400   5000   2400]';%[2400   12750   2400];
    %     d.Veh.Weight = 1000;

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