%% GET BIN �������
% 1 Strip.isMixed % 1����ϲ㣻 0��������

%% ����
function   [Bin] = cpuBin(Bin,Item,Veh)
%% ��ʼ��
    sz = size(Bin.Weight);
    Bin.Binarea = ones(sz)*-1; 
    Bin.Itemarea =  ones(sz)*-1; 
    Bin.Itemloadingrate =  ones(sz)*-1; 
    Bin.ItemloadingrateLimit =  ones(sz)*-1; 
    
%% 0: ����binװ����
% ItemloadingrateLimit - ÿ��bin��Item�������/ÿ��binȥ��ʣ���ߺ�������
% Itemloadingrate - ÿ��bin��Item�������/ÿ��bin���������
Bin = computeLoadingRate2DBin(Bin,Item,Veh); 

%% 1:

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
