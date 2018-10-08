%% GET BIN 相关属性 isTileNeed

%% 函数
function   [Bin,LU] = cpuBin(Bin,Strip,Item,LU,Veh)
%% 初始化
    sz = size(Bin.Weight);
    Bin.Binarea = ones(sz)*-1; 
    Bin.BinareaLimit = ones(sz)*-1; 
    Bin.Itemarea =  ones(sz)*-1; 
    Bin.loadingrate =  ones(sz)*-1; 
    Bin.loadingrateLimit =  ones(sz)*-1; 
    
    Bin.isTileNeed =  zeros(sz); 

%% 0: 计算LU_Bin and BIN的PID,LID,SID
% 由混合的LU.DOC新增LU_BIN, 计算BIN内包含的PID,LID,SID等数据 1808新增

% % % nbLU = size(LU.LWH,2);
% % % LU.LU_Bin = [zeros(1,nbLU);zeros(1,nbLU)];
% % % for iLU=1:nbLU
% % %     theStrip = LU.LU_Strip(1,iLU); %iLU属于第几个Item
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
    
%% 1: 计算bin装载率
% ItemloadingrateLimit - 每个bin内Item的体积和/每个bin去除剩余宽高后的总体积
% Itemloadingrate - 每个bin内Item的体积和/每个bin可用总体积
Bin = computeLoadingRate2DBin(Bin,Item,Veh); 

%% 2: Bin.isTileOneNeed 判断Bin是否全部需要平铺到1层
% 1 总长度小于车长的1/4
% % f = Bin.LW(2,:) >= 0.75*Veh.LWH(2,1); %所有车的剩余长度 >= 3/4 车长
% % Bin.isTileNeed(f) = 1;

%% 3: Bin.isTileNeed 判断Bin是否需要平铺
% V2: 考虑宽度/高度约束, 且考虑已经单层的Strip, 可运行后向增加strip
nbBin=length(Bin.Weight);
for ibin=1:nbBin
    % 2.1 找处ibin中的宽度不满的strip    (宽度不满: 横向间隙 > 单个Item的宽度)
    fS = Strip.Strip_Bin(1,:) == ibin & Strip.isWidthFull==0 ;
    if any(fS)
        % itemidx:fS中的item
        fiS = find(fS);
        itemidx = ismember(Item.Item_Strip(1,:), fiS); %itemidx:fiS个Strip对应的Item逻辑值 luidx = ismember(LU.LU_Strip(1,:), fiS);
%         Item.HLayer(itemidx)
%         Item.HLayer(itemidx)>1
        %如果所有Strip对应的Item有>1层的,则平铺.  即ibin: isTileNeed
        if any(Item.HLayer(itemidx)>1)
            Bin.isTileNeed(ibin) = 1;
        end        
    end
    
     % 2.2 找处ibin中的高度不满的strip     (高度不满: 竖向间隙 > 单个Item的宽度)
    fS = Strip.Strip_Bin(1,:) == ibin & Strip.isHeightFull==0 ;
    if any(fS)
        % itemidx:fS中的item
        fiS = find(fS);
        itemidx = ismember(Item.Item_Strip(1,:), fiS); % luidx = ismember(LU.LU_Strip(1,:), fiS);
         %如果所有Strip对应的Item有>1层的,则平铺.  即ibin: isTileNeed
        if any(Item.HLayer(itemidx)>1)
            Bin.isTileNeed(ibin) = 1;
        end
    end    
end

% V1: 仅考虑宽度/高度约束, 不考虑已经单层
% iBin = Strip.Strip_Bin(1,~Strip.isWidthFull);
% Bin.isTileNeed(iBin) = 1;
% iBin = Strip.Strip_Bin(1,~Strip.isHeightFull);
% Bin.isTileNeed(iBin) = 1;
% Bin.isTileNeed
end

%% 局部函数 %%

%% 函数1: computeLoadingRate2DBin
function Bin = computeLoadingRate2DBin(Bin,Item,Veh)
    % 初始化
    nBin = size(Bin.LW,2);
    % 计算每个Bin的装载率
    BinWidth = Veh.LWH(1,1);
    BinHeight = Veh.LWH(2,1);
    BinArea = BinWidth .* BinHeight;
    %每个Bin的可用体积 = 车辆高度*车辆宽度
    Bin.Binarea = repmat(BinArea,1,nBin);
    %每个Bin 的有限可用体积 = 宽度(bin使用宽度=车辆宽度-bin剩余宽度) *高度(bin使用高度=车辆高度-bin剩余高度)
    Bin.BinareaLimit = (BinWidth - Bin.LW(1,:)) .* (BinHeight - Bin.LW(2,:));

    a = Item.LWH;
    b = Item.Item_Bin;
    for iBin =1:nBin
        %每个Bin的装载体积
        Bin.Itemarea(iBin)= sum(a(1, (b(1,:)==iBin)) .* a(2, (b(1,:)==iBin)));
    end
    %每个bin的装载比率
    Bin.loadingrate =  Bin.Itemarea ./ Bin.Binarea;
    %每个bin的有限装载比率
    Bin.loadingrateLimit =  Bin.Itemarea ./ Bin.BinareaLimit;
end
