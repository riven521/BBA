%% GET BIN ������� isTileNeed

%% ����
function   [Bin,LU] = cpuBin(Bin,Strip,Item,LU,Veh)
%% ��ʼ��
    sz = size(Bin.Weight);
    Bin.Binarea = ones(sz)*-1; 
    Bin.BinareaLimit = ones(sz)*-1; 
    Bin.Itemarea =  ones(sz)*-1; 
    Bin.loadingrate =  ones(sz)*-1; 
    Bin.loadingrateLimit =  ones(sz)*-1; 
    
    Bin.isTileNeed =  zeros(sz); 

%% 0: ����LU_Bin and BIN��PID,LID,SID
% �ɻ�ϵ�LU.DOC����LU_BIN, ����BIN�ڰ�����PID,LID,SID������ 1808����

% % % nbLU = size(LU.LWH,2);
% % % LU.LU_Bin = [zeros(1,nbLU);zeros(1,nbLU)];
% % % for iLU=1:nbLU
% % %     theStrip = LU.LU_Strip(1,iLU); %iLU���ڵڼ���Item
% % %     LU.LU_Bin(1,iLU)= Strip.Strip_Bin(1,theStrip);
% % % end

LU.DOC=[LU.DOC; LU.LU_Bin];
nBin = size(Bin.LW,2);
for iBin=1:nBin
    tmp = LU.DOC([1,2,3], LU.DOC(10,:) == iBin);
    Bin.PID(:,iBin) = num2cell(unique(tmp(1,:))',1);
    Bin.LID(:,iBin) = num2cell(unique(tmp(2,:))',1);
    Bin.SID(:,iBin) = num2cell(unique(tmp(3,:))',1);
end
    
%% 1: ����binװ����
% ItemloadingrateLimit - ÿ��bin��Item�������/ÿ��binȥ��ʣ���ߺ�������
% Itemloadingrate - ÿ��bin��Item�������/ÿ��bin���������
Bin = computeLoadingRate2DBin(Bin,Item,Veh); 

%% 2: Bin.isTileOneNeed �ж�Bin�Ƿ�ȫ����Ҫƽ�̵�1��
% 1 �ܳ���С�ڳ�����1/4
% % f = Bin.LW(2,:) >= 0.75*Veh.LWH(2,1); %���г���ʣ�೤�� >= 3/4 ����
% % Bin.isTileNeed(f) = 1;

%% 3: Bin.isTileNeed �ж�Bin�Ƿ���Ҫƽ��
% V2: ���ǿ��/�߶�Լ��, �ҿ����Ѿ������Strip, �����к�������strip
nbBin=length(Bin.Weight);
for ibin=1:nbBin
    % 2.1 �Ҵ�ibin�еĿ�Ȳ�����strip    (��Ȳ���: �����϶ > ����Item�Ŀ��)
    fS = Strip.Strip_Bin(1,:) == ibin & Strip.isWidthFull==0 ;
    if any(fS)
        % itemidx:fS�е�item
        fiS = find(fS);
        itemidx = ismember(Item.Item_Strip(1,:), fiS); %itemidx:fiS��Strip��Ӧ��Item�߼�ֵ luidx = ismember(LU.LU_Strip(1,:), fiS);
%         Item.HLayer(itemidx)
%         Item.HLayer(itemidx)>1
        %�������Strip��Ӧ��Item��>1���,��ƽ��.  ��ibin: isTileNeed
        if any(Item.HLayer(itemidx)>1)
            Bin.isTileNeed(ibin) = 1;
        end        
    end
    
     % 2.2 �Ҵ�ibin�еĸ߶Ȳ�����strip     (�߶Ȳ���: �����϶ > ����Item�Ŀ��)
    fS = Strip.Strip_Bin(1,:) == ibin & Strip.isHeightFull==0 ;
    if any(fS)
        % itemidx:fS�е�item
        fiS = find(fS);
        itemidx = ismember(Item.Item_Strip(1,:), fiS); % luidx = ismember(LU.LU_Strip(1,:), fiS);
         %�������Strip��Ӧ��Item��>1���,��ƽ��.  ��ibin: isTileNeed
        if any(Item.HLayer(itemidx)>1)
            Bin.isTileNeed(ibin) = 1;
        end
    end    
end

% V1: �����ǿ��/�߶�Լ��, �������Ѿ�����
% iBin = Strip.Strip_Bin(1,~Strip.isWidthFull);
% Bin.isTileNeed(iBin) = 1;
% iBin = Strip.Strip_Bin(1,~Strip.isHeightFull);
% Bin.isTileNeed(iBin) = 1;
% Bin.isTileNeed
end

%% �ֲ����� %%

%% ����1: computeLoadingRate2DBin
function Bin = computeLoadingRate2DBin(Bin,Item,Veh)
    % ��ʼ��
    nBin = size(Bin.LW,2);
    % ����ÿ��Bin��װ����
    BinWidth = Veh.LWH(1,1);
    BinHeight = Veh.LWH(2,1);
    BinArea = BinWidth .* BinHeight;
    %ÿ��Bin�Ŀ������ = �����߶�*�������
    Bin.Binarea = repmat(BinArea,1,nBin);
    %ÿ��Bin �����޿������ = ���(binʹ�ÿ��=�������-binʣ����) *�߶�(binʹ�ø߶�=�����߶�-binʣ��߶�)
    Bin.BinareaLimit = (BinWidth - Bin.LW(1,:)) .* (BinHeight - Bin.LW(2,:));

    a = Item.LWH;
    b = Item.Item_Bin;
    for iBin =1:nBin
        %ÿ��Bin��װ�����
        Bin.Itemarea(iBin)= sum(a(1, (b(1,:)==iBin)) .* a(2, (b(1,:)==iBin)));
    end
    %ÿ��bin��װ�ر���
    Bin.loadingrate =  Bin.Itemarea ./ Bin.Binarea;
    %ÿ��bin������װ�ر���
    Bin.loadingrateLimit =  Bin.Itemarea ./ Bin.BinareaLimit;
end
