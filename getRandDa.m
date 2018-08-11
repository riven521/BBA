function [d] = getRandDa(n,m)
%     rng(1);
    rng('default');
    
    Par.maxHeavey = 18;
%     Par.LUbuff = [0;0]; %�Ժ�����
    
% Veh
    Veh.ID = randperm(m);   %Veh ������
    Veh.buff = zeros(3,m);
    maxSeg = 1; %�������� 3

    Veh.LWH = zeros(3,m);
    Veh.yID = zeros(maxSeg,m);
    Veh.xID = zeros(1,m);

    for i=1: m
        Veh.Weight(i) = randi([1000,1500]);    %���� randi([1000,1500]);
        Veh.LWH(1,i) = randi([10,10]); %��� randi([2200,2400]);
        Veh.LWH(2,i) = randi([20,20]); %��/Height��  randi([5000,6000]);
        Veh.LWH(3,i) = randi([10,10]); %�߶�  randi([2200,2400]); 
        Veh.volume(i) = Veh.LWH(1,i) * Veh.LWH(2,i) * Veh.LWH(3,i); % ���
        
        nYID = randi([1,maxSeg]);
        nXID = 1;
        for j=1:nYID
            if j==1
                Veh.yID(j,i) = round(Veh.LWH(2,i)/nYID);
            elseif j==nYID
                Veh.yID(j,i) = Veh.LWH(2,i);
            else
                Veh.yID(j,i) = Veh.yID(j-1,i)+ round(Veh.LWH(2,i)/nYID);
            end
        end
        for j=1:nXID
            Veh.xID(j,i) = Veh.LWH(2,i);
        end
    end

    % remove extra rows
    emptyYID =  all(Veh.yID ==0,2) == true;
    Veh.yID(emptyYID,:) = [];    
    
    
% LU   
    LU.ID = randi([11,11],1,n); %ID ������
    LU.buff = zeros(3,n); %�Ժ�����    

    LU.LWH = zeros(3,n);
    LU.isRota = zeros(1,n);
%     LU.maxL = zeros(1,n);
%     LU.yID = zeros(1,n);
%     LU.xID = zeros(1,n);
    LU.margin = zeros(4,n);

    LU.PID = zeros(1,n);
%     LU.UID = zeros(1,n);
    LU.SID = zeros(1,n);
    LU.Weight = zeros(1,n);
    LU.isH = zeros(1,n);

    [uniID] = unique(LU.ID);    %   isRotaId = randi([0,1], 1, nUniId);
    for i=1: length(uniID)
        % ��LUID���
        idx = find(LU.ID(:)==uniID(i));
        LU.LWH(1,idx) = randi([4,5]); %��� randi([1000,1300]);
        LU.LWH(2,idx) = randi([2,6]); %���� randi([700,850]); 
        LU.isRota(idx) = randi([0,1]);            %�Ƿ���ת
%         LU.maxL(idx) = randi([1,4]);
%         LU.yID(idx) = randi([0,nYID]); % TODO ���ڸ��복��һ��
%         LU.xID(idx) = randi([0,nXID]);
        LU.margin(:,idx) = randi([1,4]);
        for j=1:length(idx)
            LU.LWH(3,idx(j)) = randi([1,3]); %�߶�250,1150
            LU.PID(idx(j)) = randi([100,101]); %100,103
            LU.SID(idx(j)) = randi([200, 200]); %200, 203
%             LU.UID(idx(j)) = randi([300,300]);
            LU.Weight(idx(j)) = randi([10,20]);
            if LU.Weight(idx(j)) >= Par.maxHeavey
                LU.isH(idx(j)) = 1;
            end
        end
    end



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