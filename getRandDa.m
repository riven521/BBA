function [d] = getRandDa(n,m)
%     rng(10000);
%     rng('default');
    s = rng;     save('srng','s');
%     load('srng');    rng(s);
    
    Par.maxHeavey = 18;
%     Par.LUbuff = [0;0]; %以后无用
    
%% Veh
    Veh.ID = randperm(m);   %Veh 类型数
    Veh.LWH = zeros(3,m);
    
%     maxSeg = 1; %车辆段数 3
%     Veh.yID = zeros(maxSeg,m);
%     Veh.xID = zeros(1,m);
%     Veh.buff = zeros(3,m);

    for i=1: m
        Veh.Weight(i) = randi([1000,1500]);    %重量 randi([1000,1500]);
        Veh.LWH(1,i) = randi([2400,2400]); %宽度 randi([2200,2400]);
        Veh.LWH(2,i) = randi([12750,12750]); %长/Height度  randi([5000,6000]);
        Veh.LWH(3,i) = randi([2400,2400]); %高度  randi([2200,2400]); 
        
%         Veh.LWH(1,2) = 1800; %宽度 randi([2200,2400]);
%         Veh.LWH(2,2) = 7050; %长/Height度  randi([5000,6000]);
%         Veh.LWH(3,2) = 1800;%高度  randi([2200,2400]); 

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
%     LU.buff = zeros(3,n); %以后无用
%     LU.maxL = zeros(1,n);
%     LU.yID = zeros(1,n);
%     LU.xID = zeros(1,n);
%     LU.UID = zeros(1,n);
%     LU.isH = zeros(1,n);

    LU.ID = randi([11,13],1,n); %ID 类型数 
    LU.LID = LU.ID;

    LU.LWH = zeros(3,n);
    LU.isRota = zeros(1,n);

    LU.margin = zeros(4,n); 

    LU.PID = zeros(1,n);
    LU.SID = zeros(1,n);
    LU.Weight = zeros(1,n);


    [uniID] = unique(LU.ID);    %   isRotaId = randi([0,1], 1, nUniId);
    for i=1: length(uniID)
        % 与LUID相关
        idx = find(LU.ID(:)==uniID(i));       
        LU.LWH(1,idx) = randi([1486-200,1486+200]); %宽度 randi([1000,1300]);
        LU.LWH(2,idx) = randi([765-200,765+200]); %长度 randi([700,850]); 
        LU.isRota(idx) = randi([1,1]);            %是否旋转
%         LU.HightL(idx) = randi([1,1]);       %最大可堆垛层数
%         LU.maxL(idx) = randi([1,4]);
%         LU.yID(idx) = randi([0,nYID]); % TODO 后期改与车型一致
%         LU.xID(idx) = randi([0,nXID]);
%         LU.UID(idx(j)) = randi([300,300]);
        LU.margin(:,idx) = randi([50,50]);  %左上右下
        for j=1:length(idx)
            LU.LWH(3,idx(j)) = randi([852-50,852+50]); %高度250,1150
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
%     d.LU.isRota = isRota; %isRota: 判断某类型LU是否可以旋转,默认均可旋转,0表示不可以旋转    
%     d.LU.LWH = LWH';
%         d.LU.buff = [2;2]; 
%         d.LU.Weight = d.LU.LWH(1,:);

    %     Veh.LWH = [2400   5000   2400]';%[2400   12750   2400];
    %     d.Veh.Weight = 1000;

% function [d] = getRandDa(n)
% % 普通函数:随机产生输入数据
%      rng(222121);
% %      rng('default');
%                                                 ID = randi([11,13],1,n); %ID类型数
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