function [Res1_LUBeBinMatrix,Res2_CoordLUBin,Res3_LWHRota,Res4_DrawSeq,da] = ...
    BBA_Main(LUID,LULWH,BINLWH,PARA_whichRotationHori,LUBUFF,BINBUFF)
% 输入六个参数 LUID,LULWH,BINLWH,PARA,LUBUFF,BINBUFF
%  LUID - 行向量(n列) 托盘类型 相同数字表明同一类型,允许堆垛
%  LULWH - 矩阵(3行*n列) 托盘长宽高 画图使用该值
%  BINLWH - 列向量(3行) 车辆长宽高 目前仅考虑单车型 TODO 后期增加到多车型按顺序使用
%  PARA_whichRotationHori - 标量(1*1) 允许rotation,但考虑rotation的差异
%  LUBUFF - 列向量(2行) 托盘间长宽的总间隙(2倍的单边长宽间隙) 可用托盘长宽高=每个托盘的实际长宽高+增加的buff
%  BINBUFF - 列向量(3行) 车辆的长宽高的总间隙(2倍的单边长宽间隙,1倍的高度间隙) 可用车型长宽高=车型的实际长宽高-BINBUFF
% 输出四个参数 Res1_LUBeBinMatrix,Res2_CoordLUBin,Res3_LWHRota,Res4_DrawSeq
%  Res1_LUBeBinMatrix - 矩阵(2*n列)
%  Res2_CoordLUBin - double矩阵(3*n列)
%  Res3_LWHRota - 矩阵(*n列)
%  Res4_DrawSeq - double矩阵(2*n列)

%% Call 本函数:
% clear;close all; format long g; format bank; %NOTE 不被MATLAB CODE 支持
% rng('default');rng(1); % NOTE 是否随机的标志
% LUArray = struct('ID',[],'LWH',[],...
%     'weight',[],'Lbuffer',[],'Wbuffer',[],'Type',[],'Material',[]);
% ItemArray = struct('ID',[],'LWH',[],...
%     'weight',[],'Lbuffer',[],'Wbuffer',[]);
% StripArray = struct('ID',[],'LWH',[],... %ONLY LW
%     'weight',[],'Lbuffer',[],'Wbuffer',[]);
% BinArray = struct('ID',[],'LWH',[],...
%    'Capacity',[],'Lbuffer',[],'Wbuffer',[],'Hbuffer',[]);
% da = struct('LUArray',LUArray,'ItemArray',ItemArray,'StripArray',StripArray,'BinArray',BinArray);
%% 1参数初始化
% whichStripH 1 best 2 first 3 next; whichBinH 1 best; TODO 增加其它分批方式
% whichSortItemOrder 1 长高递减 2 最短边递减; 
% whichRotation 1:允许rotation 0:禁止
% rotation组合 1 1 2 1 0 (1 1 2 1 1 )(1 1 2 1 2) % 非rotation组合 1 1 1 0 0 （2/3 1 1 0 0）
% whichRotationHori 0:在安置顺序时按FBS_{RG}方式; 1：New/NoNew按Horizon方式 2：New/NoNew按Vertical方式
% ParaArray = struct('whichStripH',1,'whichBinH',1,'whichSortItemOrder',2,...
%     'whichRotation',1,'whichRotationHori',0,'timeLimit',100,'ub0',10);
ParaArray = struct('whichStripH',1,'whichBinH',1,'whichSortItemOrder',2,...
    'whichRotation',1,'whichRotationHori',0,'timeLimit',100,'ub0',10);
%% 2结构体da赋值
if nargin ~=0 %如果参数不为空
    da.LUArray.ID = LUID; 
    da.LUArray.LWH = LULWH;   %LWHREAL 真实尺寸    
    da.BinArray.LWH = BINLWH; %LWHREAL 真实尺寸    
    % 增加间隙 -
    da.LUArray.BUFF = LUBUFF; %BUFF 托盘LU的间隙
    da.BinArray.BUFF = BINBUFF; %BUFF 车辆BIN的间隙
    % 增加重量 -
    da.BinArray.Weight = 1000;
    da.LUArray.Weight = ones(size(LULWH,1),1);
    % 多参数测试
    ParaArray.whichStripH = PARA_whichRotationHori(1);
    ParaArray.whichBinH = PARA_whichRotationHori(2);
    ParaArray.whichSortItemOrder = PARA_whichRotationHori(3);
    ParaArray.whichRotation = PARA_whichRotationHori(4);
    ParaArray.whichRotationHori = PARA_whichRotationHori(5);
    ParaArray.whichRotationAll = PARA_whichRotationHori(6);
else
%     da.LUArray.ID = [1 1 2 2];
%     da.BinArray.LWH = [5;20;4]';
%     da.LUArray.LWH = [2 2 3 3; 5 5 6 6; 4 4 4 4];
%     da.LUArray.BUFF = [0,0];
%     da.BinArray.BUFF = [0,0,0];
%     da.BinArray.Weight = 1000;
%     da.LUArray.Weight = da.LUArray.LWH(1,:);
    
% %         da.LUArray.ID = [5236934585 4 5236934585 5236934585];
% %     da.BinArray.LWH = [8;10;4];  %[6; 10; 4];
% %     da.LUArray.LWH = [3 2 3 3; 6 4 6 6; 2 2 1 2];
% %     
% %     da.LUArray.BUFF = [0,0]';
% %     da.BinArray.BUFF = [0,0,0]'; 
            
         % 增加重量 -
%          da.BinArray.Weight = 1000;
%     da.LUArray.Weight = da.LUArray.LWH(1,:);
    %% 产生随机算例
     n=15;da=getRandDa(n);save('rndDa.mat','da');
     load('rndDa.mat')

    %%
    %     load insLU3.mat;   % load ins.mat;
    %     da.LUArray.ID = LUid;%     da.LUArray.LWH = w;%     da.BinArray.LWH = W;
    %     clear LUid w W;%     error('No input! ');
end
clc;
%% 检验Input输入数据
printstruct(da);
da = GcheckInput(da,ParaArray);
%% 启发式: LU到Item的算法
printstruct(da);
[da] = HLUtoItem(da,ParaArray);
%% 启发式：Item到Strip的算法
printstruct(da);
[da] = HItemToStrip(da,ParaArray);
%% 计算装载率 
da = computeLoadingRate(da);
    function da = computeLoadingRate(da)
        % printstruct(da);
        % 初始化
        nStrip = size(da.StripArray.LW,2);
        da.StripArray.Stripvolume = zeros(1,nStrip);
        da.StripArray.Itemvolume = zeros(1,nStrip);
        da.StripArray.Itemloadingrate = zeros(1,nStrip);
        
        % 计算每个strip的装载率
        da.StripArray.Stripvolume = da.StripArray.LW(2,:)*da.BinArray.LWH(1,1);      %每个strip的可用体积
        a = da.ItemArray.LWH;
        b = da.ItemArray.itemBeStripMatrix;
        for iStrip =1:nStrip
            %每个strip的装载体积
            da.StripArray.Itemvolume(iStrip)= sum(a(1, (b(1,:)==iStrip)) .* a(2, (b(1,:)==iStrip)));
        end
        %每个strip的装载比率
        da.StripArray.Itemloadingrate =  da.StripArray.Itemvolume ./ da.StripArray.Stripvolume;
    end
%% 启发式：Strip到Bin的算法
printstruct(da);
[da] = HStripToBin(da,ParaArray); %todo CHECK CHECK CHECK
%% Item到bin的信息获取:
printstruct(da);
[da] = HItemToBin(da); 
printstruct(da);
%% 修正输出结果（增加LWHRota:旋转后的LWH,考虑减去buffer因素)
da.LUArray.LWHRota = da.LUArray.LWH;
nLU = numel(da.LUArray.LURotaFlag);
for iLU=1:nLU
    if da.LUArray.LURotaFlag(iLU)
        da.LUArray.LWHRota(1,iLU)=da.LUArray.LWH(2,iLU);
        da.LUArray.LWHRota(2,iLU)=da.LUArray.LWH(1,iLU);
    end
end
printstruct(da);
%% 返回输出结果(原始顺序) 输出4个参数
% 参数1 - 行1:Bin序号；行2：该bin内顺序
Res1_LUBeBinMatrix=da.LUArray.LUBeBinMatrix;

% 参数2 - LU在Bin内的坐标
% 增加间隙-增加CoordLUBinWithBuff变量
da.LUArray.CoordLUBinWithBuff = da.LUArray.CoordLUBin + da.LUArray.BUFF./2;
Res2_CoordLUBin=da.LUArray.CoordLUBinWithBuff; %Res2_CoordLUBin：DOUBLE类型: Lu的xyz值 TTTTTTTTTT

% 参数3 - LU的长宽高(旋转后)
% 增加间隙-修订LWHRota为减小Buffer后的实际数据变量
da.LUArray.LWHRota = da.LUArray.LWHRota - da.LUArray.BUFF;
Res3_LWHRota=da.LUArray.LWHRota;  %Res3_LWHRota：DOUBLE LU的长宽高（旋转后）

% 参数4 - 行1：LU在Item的位置；行2：LU的ID类型
LUBeItemArray=da.LUArray.LUBeItemArray(1,:);
LUID=da.LUArray.ID;
Res4_DrawSeq = [LUBeItemArray; LUID];

% Res5_BinLWH = da.BinArray.LWH; %减去Buffer后实际可用的长宽高
%% 返回输出结果(安放顺序)
% % [~,x]=sort(Res1_LUBeBinMatrix(2,:));
% % Res1_LUBeBinMatrix=Res1_LUBeBinMatrix(:,x)
% % LUBeItemArray=LUBeItemArray(1,x)
% % Res2_CoordLUBin=Res2_CoordLUBin(:,x)
% % Res3_LWHRota=Res3_LWHRota(:,x)
% % fprintf('本算例计算全部完成 \n');
%% 画图
if nargin == 0
%     printstruct(da);
%     plot3DBPP(da,ParaArray);
% 以下修订纯为画图使用
% plot2DBPP(da,ParaArray);hold on;
da.ItemArray.LWH = da.ItemArray.LWH - da.LUArray.BUFF(:,1:size(da.ItemArray.LWH,2));
% plot2DBPP(da,ParaArray);hold on;
da.ItemArray.CoordItemBin = da.ItemArray.CoordItemBin + da.LUArray.BUFF(:,1:size(da.ItemArray.LWH,2))/2;
plot2DBPP(da,ParaArray);
else
%     % 以下修订纯为画图使用
% plot2DBPP(da,ParaArray);hold on;
da.ItemArray.LWH = da.ItemArray.LWH - da.LUArray.BUFF(:,1:size(da.ItemArray.LWH,2));
% plot2DBPP(da,ParaArray);hold on;
da.ItemArray.CoordItemBin = da.ItemArray.CoordItemBin + da.LUArray.BUFF(:,1:size(da.ItemArray.LWH,2))/2;
plot2DBPP(da,ParaArray);
end

% printstruct(da);
% mcc -W 'java:BBA_Main,Class1,1.0' -T link:lib BBA_Main.m -d '.\new'
%% 计算下届
% lb = computerLB(da); fprintf('LB = %d \n', lb); %以某个bin类型为准
end





%% ********************** 下面是ts算法的代码 暂时不用 ****************

% [ub,x,b] = HnextFit(ItemArray,BinArray);
% [ub,x,b] = HnextFit_origin(ItemArray,BinArray);
% disp(b');
% disp(x);
% fprintf('UB = %d \n', ub);
% [ub,x,b] = TSpack(d,n,w,W,lb,timeLimit, ub0,x,b,whichH);

function [toReturn,x,b] = TSpack(d,n,w,W,lb,timeLimit, ub0,x,b,whichH)
[nb,x,b] = heur(d,n,w,W,x,b,whichH);
ub0 = nb;

if nb == lb
    toReturn = nb;
    fprintf('best lb = %d ', nb);
end

%/* initial (trivial) solution */
cnb = n;
cb = 1:n;
cx = zeros(d,n);
cw = w;

%/* external loop */
D = 1; tt = 0.0;
toReturn = nb;

end

function [nb,x,b] = heur(d,n,w,W,x,b,whichH)
nb = -1;
which = (d-2)*100 + whichH; % which 为0或100
if which == 0
    [nb,x,b]  = HnextFit(n,w,W,x,b,n+1); %/* first heuristic for 2d bin packing */
    %          disp(x);
    %          disp(b);
elseif which == 100
    [nb,x,b]  = HHnextFit(n,w,W,x,b); %/* first heuristic for 3d bin packing */
end
end

    function [ub,px,pb]  = HnextFit(ItemArray,BinArray)
        % Initialize
        d = size(ItemArray.LWH,1)-1;
        n = size(ItemArray.LWH,2);
        nn = n + 1;
        w = ItemArray.LWH(1:d,:);
        W = BinArray.LWH(1:d,:);
        x = zeros(d,n); b = zeros(n,1); bNb = zeros(n,1);
        
        %/* sort the items */
        % sortD = size(w,1);%获取需要排序的维度
        [~,ord] = sort(w(d,:),'descend');%对w进行排序,只需要它的顺序ord;按第d行排序（高度)
        pw = w(:,ord);
        px = x;
        pb = (b+999); % 0 + 999
        pbNb = bNb;
        
        %/* next fit packing */
        % binLeftArray(1,ub) ： wleft
        % binLeftArray(2,ub) :  hleft
        nBin = n;
        binLeftArray = repmat(W,1,nBin);  %初始
        ub = 1;
        for i=1:n
            %     if (binLeftArray(1,ub) == W(1)) & (binLeftArray(2,ub) == W(2)) %如果是空bin
            if pbNb(ub) == 0   %如果是空bin
                if (pw(1,i) <= binLeftArray(1,ub)) && (pw(2,i) <= binLeftArray(2,ub)) %如果宽高都不超标
                    px(1,i) = 0; px(2,i) = 0;
                    binLeftArray(1,ub) = binLeftArray(1,ub) - pw(1,i);
                    binLeftArray(2,ub) = binLeftArray(2,ub) - pw(2,i);
                    pbNb(ub) = pbNb(ub) + 1;
                else
                    error('EEE');
                end
            else               %如果不是空bin
                if pw(1,i) <= binLeftArray(1,ub)  %如果i的宽满足当前bin的剩余宽度，剩余高度应该不变
                    px(1,i) = W(1) - binLeftArray(1,ub);
                    px(2,i) = W(2) - binLeftArray(2,ub) - pw(2,i);     %高度为????
                    binLeftArray(1,ub) = binLeftArray(1,ub) - pw(1,i);
                    binLeftArray(2,ub) = binLeftArray(2,ub);
                    pbNb(ub) = pbNb(ub) + 1;
                else
                    if pw(2,i)  <= binLeftArray(2,ub)  %如果i的高满足当前bin的剩余高度
                        px(1,i) = 0;
                        px(2,i) = binLeftArray(2,ub);
                        binLeftArray(1,ub) = W(1) - pw(1,i);
                        binLeftArray(2,ub) = binLeftArray(2,ub) - pw(2,i);
                        pbNb(ub) = pbNb(ub) + 1;
                    else  %如果i的高不能满足当前bin的剩余高度
                        ub = ub + 1;
                        px(1,i) = 0;   px(2,i) = 0;
                        pbNb(ub) = pbNb(ub) + 1;
                        binLeftArray(1,ub) = binLeftArray(1,ub) - pw(1,i);
                        binLeftArray(2,ub) = binLeftArray(2,ub) - pw(2,i);
                    end
                end
            end
            pb(i) = ub-1;
        end
        
        
        %原始的
        x = px(:,ord);
        b = pb(ord);
        bNb = pbNb(ord);
        ub = ub +1;
    end

    function [ub,px,pb]  = HnextFit_origin(ItemArray,BinArray)
        % Initialize
        d = size(ItemArray.LWH,1);
        if d==3
            d=d-1;
        end
        n = size(ItemArray.LWH,2);
        nn = n + 1;
        w = ItemArray.LWH(1:d,:);
        W = BinArray.LWH(1:d,:);
        x = zeros(d,n); b = zeros(n,1); bNb = zeros(n,1);
        
        %/* sort the items */
        sortD = size(w,1);%获取需要排序的维度
        [~,ord] = sort(w(sortD,:),'descend');%对w进行排序,只需要它的顺序ord
        pw = w(:,ord);
        ord;
        
        px = zeros(size(x,1),size(x,2));
        pb = 999*ones(size(b,1),size(b,2));
        %/* next fit packing */
        hleft = W(2) - pw(2,1);
        wleft = W(1);
        ub = 0; hcurr = 0;
        for i=1:n  %从第一个item开始安置
            if pw(1,i) <= wleft  %如果item的w 比wleft小，安置item到本bin本层：更新x1值，更新wleft；hleft不变
                px(1,i) = W(1) - wleft;
                wleft = wleft - pw(1,i);
            else    %否则往上一层安排。
                if pw(2,i) <= hleft  %如果item的h 比hleft小：表明bin高度充足，安置item到上一曾：更新坐标hleft，更新hcurr，wleft？ 更新坐标x值，更新wleft
                    hcurr = W(2) - hleft; %安排在同一层，所以hcurr不变，也等于pw(2,1)(但在其他bin就不对了，所以用hcurr)
                    hleft = hleft - pw(2,i);
                else  %如果放不下，开新bin，更新hcurr；更新hleft；更新数量ub+1（如果达到nn值，跳出）；更新坐标x值0，更新wleft
                    hcurr = 0;    %安排在新的bin，所以hcurr为0;
                    hleft = W(2) - pw(2,i);
                    if (ub+1 == nn)
                        break;
                    end
                    ub = ub + 1;
                end
                % 无论放在上层或开新bin，更新x1为0；更新wleft为W(1)-此item的宽w
                px(1,i) = 0;
                wleft = W(1) - pw(1,i);
            end
            % 此处统一更新x1值，即高度值，为hcurr；统一更新b值=ub；
            px(2,i) = hcurr;
            pb(i) = ub;
        end
        
        %原始的
        x = px(:,ord);
        b = pb(ord);
        ub = ub +1;
    end


    function [ nb ] = HHnextFit(n,w,W,x,b)
        
        nb = 1;
        
    end

    function [ lb ] = computerLB(da)
        sum1 = sum(prod(da.ItemArray.LWH,1));
        % todo 增加判断是否所有的BinArray中所有的bin是相同的 如果是 则继续执行
        sum2 = prod(da.BinArray.LWH(:,1));
        lb = ceil(sum1/sum2);
        if lb <=0, error('EEE');end
    end

%% ************************* 下面是注释代码  ************************

%% OLD lower计算代码
% % function [ lb ] = lower(d,n,w,W,whichL)
% % if whichL == 0
% %     lb = 1;
% % elseif whichL == 1
% %      sum1 = sum(prod(w,1));
% %      sum2 = prod(W);
% %      lb = ceil(sum1/sum2);
% % end
% %      if lb <=0, error('EEE');end
% % end

%% 结构体的三种strip算法

% % %% function [StripSolutionSort] = HnextFitDH(da)
% % function [StripSolutionSort] = HnextFitDH(da)
% % % 输入: da
% % % 输出: StripSolutionSort
% % %% 提取单类型bin,二维item数据
% % % nDim nItem nBin
% % % itemDataMatrix uniBinDataMatrix
% % nDim = size(da.ItemArray.LWH,1);  if nDim ==3, nDim = nDim-1;end
% % itemDataMatrix = da.ItemArray.LWH(1:nDim,:);
% % tmpbinDataMatrix = da.BinArray.LWH(1:nDim,:);
% % uniBinDataMatrix = unique(tmpbinDataMatrix','rows')';
% % nItem = size(itemDataMatrix,2);  nBin = nItem;
% % if size(uniBinDataMatrix,2)==1
% %     fprintf('本算例只有一个箱型 宽=%1.0f 长=%1.0f  \n', uniBinDataMatrix);
% %     fprintf('本算例有 %d 个物品,其宽长分别为 \n',nItem);
% %     fprintf('%1.0f %1.0f \n',itemDataMatrix);
% % else
% %     error('本算例有多个箱型,超出期望 \n');
% % end
% % %% 输出StripSolutionSort初始化
% % % StripSolutionSort: stripWidth stripDataMatrix stripBeItemArray
% % % StripSolutionSort: itemOrd itemDataMatrixSort itemCoordMatrixSort itemBeLevelMatrixSort
% % StripSolutionSort.stripWidth = uniBinDataMatrix (1,1); %只需要strip的宽度,dim1为宽度
% % StripSolutionSort.stripDataMatrix = zeros(2,nBin);  %dim2-长(高)度(以最高的计算)
% % StripSolutionSort.stripDataMatrix(1,:) = StripSolutionSort.stripWidth; %dim1-宽度剩余 ;
% % StripSolutionSort.stripBeItemArray = zeros(1,nBin); %某个strip包含多少个item,具体编号不计算
% % 
% %  %/* sort the items */
% % [~,itemOrd] = sort(itemDataMatrix(nDim,:),'descend'); %对w进行排序,只需要它的顺序ord;按第nDim行排序（长/高度)
% % StripSolutionSort.itemOrd = itemOrd;
% % StripSolutionSort.itemDataMatrixSort = itemDataMatrix(:,itemOrd);
% % StripSolutionSort.itemCoordMatrixSort = zeros(nDim,nItem);
% % StripSolutionSort.itemBeStripMatrixSort = zeros(2,nItem); %dim1:属于第几个level dim2:属于该level第几个排放
% % %% NF循环
% % iLevel = 1; iItem = 1;
% % %/* next fit packing */
% % while 1
% %     if iItem > nItem, break; end
% %     % 不同条件下的选择：如果当前item的宽<=当前strip的当前level的宽
% %     flag = StripSolutionSort.itemDataMatrixSort(1,iItem) <= StripSolutionSort.stripDataMatrix(1,iLevel);
% %     if ~isempty(flag)
% %         thisLevel = iLevel;
% %         [StripSolutionSort] = insertItemToStrip(thisLevel,iItem,StripSolutionSort);
% %         iItem = iItem + 1;            
% %     else
% %         iLevel = iLevel + 1;% 如果宽度不满足，则level升级
% %     end
% % end
% % printstruct(StripSolutionSort);
% % end
% % 
% % %% function [StripSolutionSort] = HfirstFitDH(da)
% % function [StripSolutionSort] = HfirstFitDH(da)
% % % 输入: da
% % % 输出: StripSolutionSort
% % %% 提取单类型bin,二维item数据
% % % nDim nItem nBin
% % % itemDataMatrix uniBinDataMatrix
% % nDim = size(da.ItemArray.LWH,1);  if nDim ==3, nDim = nDim-1;end
% % itemDataMatrix = da.ItemArray.LWH(1:nDim,:);
% % tmpbinDataMatrix = da.BinArray.LWH(1:nDim,:);
% % uniBinDataMatrix = unique(tmpbinDataMatrix','rows')';
% % nItem = size(itemDataMatrix,2);  nBin = nItem;
% % if size(uniBinDataMatrix,2)==1
% %     fprintf('本算例只有一个箱型 宽=%1.0f 长=%1.0f  \n', uniBinDataMatrix);
% %     fprintf('本算例有 %d 个物品,其宽长分别为 \n',nItem);
% %     fprintf('%1.0f %1.0f \n',itemDataMatrix);
% % else
% %     error('本算例有多个箱型,超出期望 \n');
% % end
% % %% 输出StripSolutionSort初始化
% % % StripSolutionSort: stripWidth stripDataMatrix stripBeItemArray
% % % StripSolutionSort: itemOrd itemDataMatrixSort itemCoordMatrixSort itemBeLevelMatrixSort
% % StripSolutionSort.stripWidth = uniBinDataMatrix (1,1); %只需要strip的宽度,dim1为宽度
% % StripSolutionSort.stripDataMatrix = zeros(2,nBin);  %dim2-长(高)度(以最高的计算)
% % StripSolutionSort.stripDataMatrix(1,:) = StripSolutionSort.stripWidth; %dim1-宽度剩余 ;
% % StripSolutionSort.stripBeItemArray = zeros(1,nBin); %某个strip包含多少个item,具体编号不计算
% % 
% %  %/* sort the items */
% % [~,itemOrd] = sort(itemDataMatrix(nDim,:),'descend'); %对w进行排序,只需要它的顺序ord;按第nDim行排序（长/高度)
% % StripSolutionSort.itemOrd = itemOrd;
% % StripSolutionSort.itemDataMatrixSort = itemDataMatrix(:,itemOrd);
% % StripSolutionSort.itemCoordMatrixSort = zeros(nDim,nItem);
% % StripSolutionSort.itemBeStripMatrixSort = zeros(2,nItem); %dim1:属于第几个level dim2:属于该level第几个排放
% % %% FF循环
% % iLevel = 1; iItem = 1;
% % %/* next fit packing */
% % while 1
% %     if iItem > nItem, break; end
% %     % 不同条件下的选择：如果find宽度足够的多个level,并安置在第一个遇到的 唯一区别是thisLevel的获取
% %     flag = find(StripSolutionSort.stripDataMatrix(1,1:iLevel) >= StripSolutionSort.itemDataMatrixSort(1,iItem));
% %     if ~isempty(flag)
% %         thisLevel = flag(1);
% %         [StripSolutionSort] = insertItemToStrip(thisLevel,iItem,StripSolutionSort);
% %         iItem = iItem + 1;
% %     else
% %         iLevel = iLevel + 1;% 如果宽度不满足，则level升级
% %     end
% % end
% % printstruct(StripSolutionSort);
% % end
% % 
% % %% function [StripSolutionSort] = HbestFitDH(da)
% % function [StripSolutionSort] = HbestFitDH(da)
% % % 输入: da
% % % 输出: StripSolutionSort
% % %% 提取单类型bin,二维item数据
% % % nDim nItem nBin
% % % itemDataMatrix uniBinDataMatrix
% % nDim = size(da.ItemArray.LWH,1);  if nDim ==3, nDim = nDim-1;end
% % itemDataMatrix = da.ItemArray.LWH(1:nDim,:);
% % tmpbinDataMatrix = da.BinArray.LWH(1:nDim,:);
% % uniBinDataMatrix = unique(tmpbinDataMatrix','rows')';
% % nItem = size(itemDataMatrix,2);  nBin = nItem;
% % if size(uniBinDataMatrix,2)==1
% %     fprintf('本算例只有一个箱型 宽=%1.0f 长=%1.0f  \n', uniBinDataMatrix);
% %     fprintf('本算例有 %d 个物品,其宽长分别为 \n',nItem);
% %     fprintf('%1.0f %1.0f \n',itemDataMatrix);
% % else
% %     error('本算例有多个箱型,超出期望 \n');
% % end
% % %% 输出StripSolutionSort初始化
% % % StripSolutionSort: stripWidth stripDataMatrix stripBeItemArray
% % % StripSolutionSort: itemOrd itemDataMatrixSort itemCoordMatrixSort itemBeLevelMatrixSort
% % StripSolutionSort.stripWidth = uniBinDataMatrix (1,1); %只需要strip的宽度,dim1为宽度
% % StripSolutionSort.stripDataMatrix = zeros(2,nBin);  %dim2-长(高)度(以最高的计算)
% % StripSolutionSort.stripDataMatrix(1,:) = StripSolutionSort.stripWidth; %dim1-宽度剩余 ;
% % StripSolutionSort.stripBeItemArray = zeros(1,nBin); %某个strip包含多少个item,具体编号不计算
% % 
% %  %/* sort the items */
% % [~,itemOrd] = sort(itemDataMatrix(nDim,:),'descend'); %对w进行排序,只需要它的顺序ord;按第nDim行排序（长/高度)
% % StripSolutionSort.itemOrd = itemOrd;
% % StripSolutionSort.itemDataMatrixSort = itemDataMatrix(:,itemOrd);
% % StripSolutionSort.itemCoordMatrixSort = zeros(nDim,nItem);
% % StripSolutionSort.itemBeStripMatrixSort = zeros(2,nItem); %dim1:属于第几个level dim2:属于该level第几个排放
% % %% FF循环
% % iLevel = 1; iItem = 1;
% % %/* next fit packing */
% % while 1
% %     if iItem > nItem, break; end
% %     % 不同条件下的选择：如果 find宽度足够的多个level,并安置在最小剩余宽度的 
% %     flag = find(StripSolutionSort.stripDataMatrix(1,1:iLevel) >= StripSolutionSort.itemDataMatrixSort(1,iItem));
% %     if ~isempty(flag)
% %         % 唯一与FF区别从这到thisLevel的计算（选中满足条件且最小的
% %         tepMin = StripSolutionSort.stripDataMatrix(1,1:iLevel);
% %         tepMin = min(tepMin(flag));
% %         thisLevel = find(StripSolutionSort.stripDataMatrix(1,1:iLevel)==tepMin);
% %         if length(thisLevel)>1
% %             thisLevel = thisLevel(1);
% %         end 
% %         
% %         [StripSolutionSort] = insertItemToStrip(thisLevel,iItem,StripSolutionSort);
% %         iItem = iItem + 1;
% %     else
% %         iLevel = iLevel + 1;% 如果宽度不满足，则level升级
% %     end
% % end
% % printstruct(StripSolutionSort);
% % end
% % 
% % 
% % 
% % 

%% 非结构体的HfirstFitDH2算法
% % function [stripLeftMatrix,pbelongMatrix,pitemMatrix,pcoordMatrix,ord]  = HfirstFitDH2(ItemArray,BinArray)
% % % 输入参数初始化
% % nDim = size(ItemArray.LWH,1);
% % if nDim ==3, nDim = nDim-1;end
% % nItem = size(ItemArray.LWH,2);
% % nBin = nItem;
% % % nn = n + 1;
% % itemMatrix = ItemArray.LWH(1:nDim,:);
% % binMatrix = BinArray.LWH(1:nDim,:);
% % % 输出参数初始化
% % coordMatrix = zeros(nDim,nItem);
% % stripWidth = binMatrix(1,1); %只需要strip的宽度,dim1为宽度
% % stripLeftMatrix = [stripWidth*(ones(1,nBin));zeros(1,nBin);zeros(1,nBin)];%初始化strip: dim1-宽度剩余 ; dim2-长(高)度(以最高的计算); dim3-该strip包含的item个数
% % belongMatrix = zeros(2,nItem); %dim1:属于第几个level dim2:属于该level第几个排放
% % 
% % %/* sort the items */
% % [~,ord] = sort(itemMatrix(nDim,:),'descend');%对w进行排序,只需要它的顺序ord;按第d行排序（高度)
% % pitemMatrix = itemMatrix(:,ord);
% % pcoordMatrix = coordMatrix;
% % pbelongMatrix = belongMatrix;
% % iLevel = 1; iItem = 1;  
% % 
% % %%
% % %/* first fit packing */
% % while 1
% %     if iItem > nItem, break; end
% %     % find宽度足够的多个level,并安置在第一个遇到的 唯一区别从这到thisLevel的计算 + 后面的thisLevel的替换
% %     findLevelArray = find(stripLeftMatrix(1,1:iLevel) >= pitemMatrix(1,iItem));
% %     if findLevelArray
% %         thisLevel = findLevelArray(1);
% %         [pcoordMatrix,stripLeftMatrix,pbelongMatrix] = insertItemToStrip2(thisLevel,iItem,pitemMatrix,stripWidth,pcoordMatrix,stripLeftMatrix,pbelongMatrix);
% %         iItem = iItem + 1;        
% %     % 如果宽度不满足，则level升级
% %     else
% %         iLevel = iLevel + 1;
% %     end
% % end
% %      pcoordMatrix
% %      stripLeftMatrix
% %      pbelongMatrix
% % end

%% 非结构体的HbesttFitDH2算法
% % function [stripLeftMatrix,pbelongMatrix,pitemMatrix,pcoordMatrix,ord]  = HbestFitDH2(ItemArray,BinArray)
% % % 输入参数初始化
% % nDim = size(ItemArray.LWH,1);
% % if nDim ==3, nDim = nDim-1;end
% % nItem = size(ItemArray.LWH,2);
% % nBin = nItem;
% % % nn = n + 1;
% % itemMatrix = ItemArray.LWH(1:nDim,:);
% % binMatrix = BinArray.LWH(1:nDim,:);
% % % 输出参数初始化
% % coordMatrix = zeros(nDim,nItem);
% % stripWidth = binMatrix(1,1); %只需要strip的宽度,dim1为宽度
% % stripLeftMatrix = [stripWidth*(ones(1,nBin));zeros(1,nBin);zeros(1,nBin)];%初始化strip: dim1-宽度剩余 ; dim2-长(高)度(以最高的计算); dim3-该strip包含的item个数
% % belongMatrix = zeros(2,nItem); %dim1:属于第几个level dim2:属于该level第几个排放
% % 
% % %/* sort the items */
% % [~,ord] = sort(itemMatrix(nDim,:),'descend');%对w进行排序,只需要它的顺序ord;按第d行排序（高度)
% % %         ord = 1:nItem;    % 此语句目的不对items排序
% % pitemMatrix = itemMatrix(:,ord);
% % pcoordMatrix = coordMatrix;
% % pbelongMatrix = belongMatrix;
% % iLevel = 1; iItem = 1;  
% % 
% % %%
% % %/* best fit packing */
% % while 1
% %     if iItem > nItem, break; end
% %     % find宽度足够的多个level,并安置在最小剩余宽度的 
% %     findLevelArray = find(stripLeftMatrix(1,1:iLevel) >= pitemMatrix(1,iItem));
% %     if findLevelArray
% % %         唯一与FF区别从这到thisLevel的计算（选中满足条件且最小的
% %         tepMin = stripLeftMatrix(1,1:iLevel);
% %         tepMin = min(tepMin(findLevelArray));
% %         thisLevel = find(stripLeftMatrix(1,1:iLevel)==tepMin);
% %         if length(thisLevel)>1
% %             thisLevel = thisLevel(1);
% %         end
% %         [pcoordMatrix,stripLeftMatrix,pbelongMatrix] = insertItemToStrip2(thisLevel,iItem,pitemMatrix,stripWidth,pcoordMatrix,stripLeftMatrix,pbelongMatrix);
% %         iItem = iItem + 1;
% %         
% %     % 如果宽度不满足，则level升级
% %     else
% %         iLevel = iLevel + 1;
% %     end
% % end
% %      pcoordMatrix
% %      stripLeftMatrix
% %      pbelongMatrix
% % 
% % end

%% 非结构体的HbestFitBinDH算法
% % function [pbelongItemBinMatrix,pbelongStripBinMatrix,pcoordItemBinMatrix,binLeftMatrix ] = HbestFitBinDH(stripLeftMatrix,pbelongMatrix,pitemMatrix,pcoordMatrix,ItemArray,BinArray)
% % % 输入参数初始化
% % nDim = size(ItemArray.LWH,1);
% % if nDim == 3, nDim = nDim-1;end
% % nItem = size(ItemArray.LWH,2);
% % nBin = nItem;
% % 
% % nStrip = sum(stripLeftMatrix(3,:)>0); %具体使用的Strip的数量
% % % nn = n + 1;
% % binMatrix = BinArray.LWH(1:nDim,:);
% % 
% % stripWidth = binMatrix(1,1); %只需要strip的宽度,dim1为宽度
% % % 输出参数初始化
% % pbelongItemBinMatrix = zeros(2,nItem); % dim1:序号item在某个bin dim2:进入顺序
% % pbelongStripBinMatrix = zeros(2,nBin); % dim1:序号strip在某个bin dim2:进入顺序
% % pcoordItemBinMatrix = zeros(nDim,nItem); %坐标值
% % binLeftMatrix = [binMatrix(1,1)*(ones(1,nBin));binMatrix(2,1)*ones(1,nBin);zeros(1,nBin);zeros(1,nBin)];
% % %初始化bin: dim1-bin宽度剩余 ; dim2-bin长(高)度(555剩余）; dim3-该bin包含的item个数; dim4-该bin包含的strip个数;
% % 
% % %/* sort the strips by 长(高) 默认已经是按这个顺序 无需再行排序 */
% % 
% % % Best fit 
% % iStrip=1;iBin=1;
% % while 1
% %     if iStrip > nStrip, break; end
% %     findBinArray = find(binLeftMatrix(2,1:iBin) >= stripLeftMatrix(2,iStrip));
% %     if findBinArray
% %         tepMin = binLeftMatrix(2,1:iBin);
% %         tepMin = min(tepMin(findBinArray)); % 555 check
% %         thisBin = find(binLeftMatrix(2,1:iBin)==tepMin);
% %         if length(thisBin)>1
% %             thisBin = thisBin(1);
% %         end
% %         %更新strip归属信息
% %         pbelongStripBinMatrix(1,iStrip) = thisBin;
% %         binLeftMatrix(4,thisBin) = binLeftMatrix(4,thisBin) + 1; %本bin下第几次安置strip
% %         pbelongStripBinMatrix(2,iStrip) = binLeftMatrix(4,thisBin);
% %         
% %         %获取本iStrip内的item序号, 并更新Item归属信息
% %         idxItemStrip = find(pbelongMatrix(1,:)==iStrip);
% %         pbelongItemBinMatrix(1,idxItemStrip) = thisBin;    %第几个bin
% %         
% %         %更新bin内信息
% %         binLeftMatrix(1,thisBin) = min(binLeftMatrix(1,thisBin),stripLeftMatrix(1,iStrip)); %所有剩余宽度的最小值
% %         binLeftMatrix(2,thisBin) = binLeftMatrix(2,thisBin) - stripLeftMatrix(2,iStrip); %更新剩余高度
% %         binLeftMatrix(3,thisBin) = binLeftMatrix(3,thisBin) + length(idxItemStrip); %本bin下合计几个item
% %         
% %         
% %         %更新xy坐标信息 x不变 y通过bin高度-bin剩余高度-本次strip高度
% %         pcoordItemBinMatrix(1,idxItemStrip) = pcoordMatrix(1,idxItemStrip);
% %         pcoordItemBinMatrix(2,idxItemStrip) = binMatrix(2,1) - (binLeftMatrix(2,thisBin) + stripLeftMatrix(2,iStrip));
% %      
% %         iStrip = iStrip + 1;
% %     else
% %         iBin = iBin + 1;
% %     end    
% % end
% % 
% % % 增加更新pbelongItemBinMatrix中访问顺序的步骤
% % for iItem=1:nItem
% %     tmp = find(pbelongItemBinMatrix(1,:)==iItem);
% %     if isempty(tmp),  break;   end
% %     for i=1:length(tmp)
% %         pbelongItemBinMatrix(2,tmp(i)) = pbelongItemBinMatrix(2,tmp(i)) + i;
% %     end
% % end
% % 
% % pbelongStripBinMatrix
% % pbelongItemBinMatrix
% % binLeftMatrix
% % pcoordItemBinMatrix
% % end
