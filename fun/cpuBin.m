%% GET BIN ������� isTileNeed

%% ����
function   [Bin] = cpuBin(Bin,Strip,Item,LU,Veh)
%% ��ʼ��
    sz = size(Bin.Weight);
    nBin = size(Bin.LW,2);
%     Bin.Binarea = ones(sz)*-1; 
%     Bin.BinareaLimit = ones(sz)*-1; 
%     Bin.Itemarea =  ones(sz)*-1; 
%     Bin.loadingrate =  ones(sz)*-1; 
%     Bin.loadingrateLimit =  ones(sz)*-1; 
    Bin.isTileNeed =  zeros(sz); 

%% ����1: V2 ����BIN��PID,LID,SID
    t = struct2table(structfun(@(x) x',LU,'UniformOutput',false));
    
    for iBin=1:nBin
        f = t.LU_Bin(:,1) == iBin;   
        Bin.LID(:,iBin) = {unique(t.ID(f))};           % NOTE: Bin���LID��LU��ID
        %         Item.LID(:,iItem) = {unique(t.LID(f))};
        Bin.SID(:,iBin) = {unique(t.SID(f))};
        Bin.EID(:,iBin) = {unique(t.EID(f))};
        Bin.PID(:,iBin) = {unique(t.PID(f))};
    end
    

    
%% 1: ����binװ���� Ŀǰ��;���� ��ʱע��
% loadingrateLimit - ÿ��bin��Item�������/ÿ��binȥ��ʣ���ߺ�������
% loadingrate - ÿ��bin��Item�������/ÿ��bin���������
% [Bin.loadingrate,Bin.loadingrateLimit] = computeLoadingRate2DBin(Bin,Item,Veh);

    % ���ú��� ����binװ����
    % tmpBin = computeLoadingRate2DBin1(Bin,Item,Veh);  if any(Bin.loadingrate~= tmpBin.loadingrate),     error('�汾��������ͬ');  end %���ڿ�ɾ

%% 2: Bin.isTileOneNeed �ж�Bin�Ƿ�ȫ����Ҫƽ�̵�1��
% 1 �ܳ���С�ڳ�����1/4
% % f = Bin.LW(2,:) >= 0.75*Veh.LWH(2,1); %���г���ʣ�೤�� >= 3/4 ����

%% 3: Bin.isTileNeed �ж�Bin�Ƿ���Ҫ˦βƽ�� ��HBinpingpuʹ��: ����������Ȼ�߶Ȳ�������Ҫ˦βƽ��
% ����Item.HLayer Strip.isHeightFull Strip.isWidthFull �ж�
Bin.isTileNeed = computeisTileNeedofBin(nBin,Strip,Item);
end




%% �ֲ����� %%
function TF = computeisTileNeedofBin(n,Strip,Item)

% V3: ���ǿ��/�߶�Լ��(��˦β������,����˦βƽ��), �ҿ����Ѿ������Strip, �����к�������strip
[TF] = deal(zeros(1,n));

for ibin=1:n
    
    % 2.1 �Ҵ�ibin�еĿ�Ȳ��� OR �߶Ȳ�����strip    (��Ȳ���: �����϶ > ����Item�Ŀ��)  (�߶Ȳ���: �����϶ > ����Item�Ŀ��)
    fS = (Strip.Strip_Bin(1,:) == ibin & Strip.isWidthFull==0) | (Strip.Strip_Bin(1,:) == ibin & Strip.isHeightFull==0); %...
                                                        %         | (Strip.Strip_Bin(1,:) == ibin & Strip.isHeightBalance==0);
    if any(fS)
        % itemidx:fS�е�item
        fiS = find(fS);
        itemidx = ismember(Item.Item_Strip(1,:), fiS); %itemidx:fiS��Strip��Ӧ��Item�߼�ֵ luidx = ismember(LU.LU_Strip(1,:), fiS);
        %�������Strip��Ӧ��Item��>1���,��ƽ��.  ��ibin: isTileNeed
        if any(Item.HLayer(itemidx)>1)
            TF(ibin) = 1;
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


%% ����1: v1 computeLoadingRate2DBin
function Bin = computeLoadingRate2DBin1(Bin,Item,Veh)
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

%% ����1: v2 computeLoadingRate2DBin
function [loadingrate,loadingrateLimit] = computeLoadingRate2DBin(Bin,Item,Veh)
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
    loadingrate =  Bin.Itemarea ./ Bin.Binarea;
    %ÿ��bin������װ�ر���
    loadingrateLimit =  Bin.Itemarea ./ Bin.BinareaLimit;
end

%% ����Ϊע��

%% 0: % �ɻ�ϵ�LU.DOC����LU_BIN, ����BIN�ڰ�����PID,LID,SID������ 1808����

% % % nbLU = size(LU.LWH,2);
% % % LU.LU_Bin = [zeros(1,nbLU);zeros(1,nbLU)];
% % % for iLU=1:nbLU
% % %     theStrip = LU.LU_Strip(1,iLU); %iLU���ڵڼ���Item
% % %     LU.LU_Bin(1,iLU)= Strip.Strip_Bin(1,theStrip);
% % % end

%% V1 ����LU_Bin and BIN��PID,LID,SID
% % LU.DOC=[LU.DOC; LU.LU_Bin];
% % nBin = size(Bin.LW,2);
% % for iBin=1:nBin
% %     tmp = LU.DOC([1,2,3], LU.DOC(10,:) == iBin);
% %     Bin.PID(:,iBin) = num2cell(unique(tmp(1,:))',1);
% %     Bin.LID(:,iBin) = num2cell(unique(tmp(2,:))',1);
% %     Bin.SID(:,iBin) = num2cell(unique(tmp(3,:))',1);
% % end