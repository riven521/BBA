function [output_CoordLUBin,output_LU_LWH,output_LU_Seq] = ...
    BBA_Main(LUID,LULWH,VEHID,VEHLWH,varargin)  %前4个必须
% Inputs：LUID: 种类少 ; LULID: 种类多
%   LUID	                (1,n)   托盘类型 相同数字表明同一类型,允许堆垛
%   cpuLUVeh后其unique后必须是从1开始严格递增
%   LULWH                (3,n)   托盘宽长高
%   VEHID                 (1,m)  车型编号
%   VEHLWH              (3,m)   车型宽长高（考虑多车型）
%   ------------------------------------------------------
%   LUSID                   (1,n)   托盘供应商编号
%   cpuLUVeh后其unique后必须是从1开始严格递增
%   LUPID                   (1,n)   托盘零部件编号
%   cpuLUVeh后其unique后必须是从1开始严格递增
%   LUISROTA            (1,n)  托盘是否允许旋转
%   LUMARGIN          (4,n)   托盘间margin(1-4左右上下)  可用托盘长宽高=每个托盘的实际长宽高+增加的margin
%   LUWEIGHT           (1,n)  托盘重量
%   VEHWEIGHT        (1,m)  车型最大承载重量
%   LULID                   (1,n)   托盘类型编号 暂无规定？
%   LUINDEX            (1,n） 托盘索引号-刘强专用
%   sort后必须是从1开始严格递增
%   LUEID                  (1,n） 托盘EP LOCATION - Milkrun版本参数
%   cpuLUVeh后其unique后必须是从1开始严格递增
% Outputs
%   output_CoordLUBin      (3,n)    每个LU的X,Y,Z
%   output_LU_LWH            (3,n)    每个LU的宽长高（旋转后的：实际值）
%   output_LU_Seq             (8,n)    行1: LU在某个BIN内；行2: LU在该BIN内的安放顺序 。。。

%% Initialize Global Variable
% clear;close all; format long g; format bank; %NOTE 不被MATLAB CODE 支持
% rng('default');rng(1); % NOTE 是否随机的标志
close all;

global ISshuaiwei ISdiagItem ISstripbalance ISreStripToBin ISisNonMixed ISisMixTile  
global ISisNonMixedLU ISisMixTileLU
global parBalance

% 全局变量1： 平铺开关
%   ISpingpu : 是否甩尾平铺 ISpingpuall：是否整车平铺  二选一 若平铺：必须甩尾平铺，可不整车平铺；但不可整车平铺，不甩尾平铺。
global ISpingpu ISpingpuAll
% 全局变量2： 混装间隙开关
%   parGap ： 是否允许间隙调整（仅在甩尾平铺或整车平铺成功的车内进行）；parMulipleGap：是否在间隙调整过程中允许多次调整（必须有）
global  parGap parMulipleGap
% 全局变量3： 改变车型开关
% ISlastVehType：将最后不满的一车改为更小的车型
global ISlastVehType


% 全局变量5： 版本控制
% verMilkRun：1 milkrun版本（包含多一个EID输入） 0 非milkrun；
global verMilkRun
% 全局变量6： 作图开关
global ISplotBBA  ISplotEachPingPu ISplotStrip ISplotPause ISplotShowType ISplotShowGapAdjust % plotStrip
% 全局变量7： 已删除
% global ISisGpreprocLU1 ISsItemAdjust ISreStripToBinMixed ISplotSolu

% 开关2
parGap = 1  % 是否允许主函数的间隙调整
parMulipleGap = 1 % 是否允许间隙递归多次调整

parBalance = 8/30;
% ISisNonMixedLU    1 LU可以形成满垛 0 必定有非满垛生成 LU排序依据
% ISisMixTileLU          1 当isNonMixed=0时, 将非满垛对应的LU赋值为1（结合LU排序生成ITEM知识）
% ISisGpreprocLU1    1 与ISisNonMixedLU取值有关
% ISstripbalance        1 调用高度均衡开关



% 开关 + Gpreproc 的V2版本 修复业务3问题(即ITEM非满垛且一层的均衡问题)
ISstripbalance = 1     % 555：有了图形好看, 堆垛均衡使用 同一Strip非混合且高度不均衡且LU层数差异值>1时操作 （方法：对应LU的最大层数递减; 如无法）
% ISisGpreprocLU1 = 1 % 必须1; 配合ISstripbalance使用 1表示对同一水平strip内ITEM可能有2个以上的判断为ISisNonMixedLU=1->堆垛均衡使用 0 表示正常判断
% ISstripbalance=0 即不均衡时,可以ISisGpreprocLU1=0; 表示LU使劲高度堆,更多的ISisNonMixedLU=0.

%
% DEL ISisMixedStrip = 1 可以删除了 % 1表示依据LU.ID判断是否混合 0依据LU.LID判断 NOTE: 所有STRIP
% ITEM均为ID替换LID
% ISsItemAdjust = 1              % 暂时不用 用途忘记了
% ISreStripToBinMixed = 1   %车头优先非AllPure类型, 再考虑优先LU数量排序参数 默认为1 应该可以删除的参数

ISplotBBA = 1 % 是否显示LU/Strip/Bin的结果（均已排序）
ISplotShowGapAdjust = 0 % 是否显示Gap调整过程

ISplotShowType = 3 % 1 LID 2 PID 3 ID 作图的颜色标记选项 4 SID 5EID
        % ISplotSolu = 0

                            ISplotStrip = 0             % 每次Run algorithm 生成Strip就显示结果 看细节 后期替换为同一
                            ISplotEachPingPu = 0 % 每次Main 平铺时 生成Strip就显示结果 看细节 后期替换为同一

ISplotPause = 0.0;  %-0.05 % plot间隔时间

ISdiagItem = 0  % 默认为 0 吧 为1 总有些过于低的被认为Item高度满层, check原因吧

% 下面还不完整, 可能要调 目前全部为1
ISisNonMixedLU = 1 % 555: 优先非混合LU形成ITEM, 图好看许多 必须有 默认为 1
ISisMixTileLU = 1       % 555: 优先混合LU的单纯ITEM部分来形成ITEM, 图好看许多 必须有 默认为 1

ISisNonMixed = 1    % 555: 优先非混合Item形成STRIP, 图好看许多 必须有 默认为 1
ISisMixTile  = 1         % 555: 优先混合Item的单纯Strip部分来形成STRIP, 图好看许多 必须有 默认为 1 但可能出现混合现象

ISreStripToBin = 1   % 车头优先LU数量排序参数 默认为1 必须

ISshuaiwei = 1         % 555 : 宽度和高度不满, 甩尾   ******  该参数需要和下面的pingpu结合使用 不甩尾 平铺无法进行*******


ISlastVehType = 0   % 555: 最后一车的调整, 与其它无关, 暂不考虑

% Milkrun版本 特殊全局变量
verMilkRun = 0  % 555: 默认不是MilkRun版本(MilkRun是9个参数的版本)

    % MR的参数9: EP LOCATION增加, 多的变量的自动增加
     if nargin > 1 && length(varargin) < 9 %'LUEID',varargin{9});
          verMilkRun  = 0;   varargin{9} = ones(1,length(LUID));
     elseif nargin > 1
           verMilkRun = 1; % 9个输入参数为milkrun版本
     end
     
if verMilkRun == 1
ISpingpu = 0         % 555 : 宽度和高度不满, 且层数>1, 平铺. 可能有问题 (在于平铺后与ISisNonMixed矛盾)
ISpingpuAll = 0       %555: 所有均平铺, 只要该车辆放得下; 若放不下, 考虑上面甩尾平铺问题
else
ISpingpu = 1          % 555 : 宽度和高度不满, 且层数>1, 平铺. 可能有问题 (在于平铺后与ISisNonMixed矛盾)
ISpingpuAll = 1       %555: 所有均平铺, 只要该车辆放得下; 若放不下, 考虑上面甩尾平铺问题
end

%% Initialize Parameter Variable
nAlg = 1;
pA(nAlg) = ParameterInitialize('whichStripH', 3,...
                             'whichBinH',3, ...
                             'whichSortItemOrder',3, ... 
                             'whichRotation',2, ...
                             'whichRotationHori', 1);
                         
%% Initialize Data Structure  - d
if nargin < 1 % Randome Generate
    n=32; m=1;                                          % 16需要注意 250 srng1
    d = DataRandomInitialize(n,m);        % 0 默认值; >0 随机产生托盘n个算例 仅在直接允许BBA时采用
    save( strcat( '.\useless\', strcat('GoodIns',num2str(n))), 'd');          %     load .\new\GoodIns200.mat;
else
    d = DataInitialize( ...
            'LUID', LUID,...
            'LULWH',LULWH, ...
            'VEHID',VEHID,...
            'VEHLWH',VEHLWH,...
            'LUSID',varargin{1},...
            'LUPID',varargin{2},...
            'LUISROTA',varargin{3},...
            'LUMARGIN',varargin{4},...
            'LUWEIGHT',varargin{5},...
            'VEHWEIGHT',varargin{6},...
            'LULID',varargin{7},...
            'LUINDEX',varargin{8},...
            'LUEID',varargin{9}); 
end

%% 1：确保LU/Veh结构体的列数相同 且 非负的矩阵
if ~isSameCol(d.LU) || ~isSameCol(d.Veh)
    warning('LU或Veh的列数不同');
end
structfun(@(x) validateattributes(x,{'double'},{'nonnegative','2d'}), d.LU, 'UniformOutput', false);
structfun(@(x) validateattributes(x,{'double'},{'nonnegative','2d'}), d.Veh, 'UniformOutput', false);

%% 2:  在此基础上进行数据行列调整，数据修改（LWH, ID号等）
fprintf(1,'check and preprocessing Lu and Veh processing ... ');

d = chkInput(d);
if ~isSameCol(d.LU) || ~isSameCol(d.Veh)
    errorr('LU或Veh的列数不同');
end
%% 3：对新增后将进行运算的数据进行check Gpreproc=cpuLUVeh d是需要保留的
[d.LU,d.Veh] = cpuLUVeh(d.LU,d.Veh);
d = chkInput(d);

% ************* 保留输入且处理过后的数据 ***********
maind = d; % 主要的输入数据保留 （NOTE : ********** maind的LWH必须为不含margin的 ********）

if verMilkRun == 0 && ~isscalar(unique(d.LU.SID))
    error('目前为非MilkRun版本,不允许有多家供应商编号');
end
if verMilkRun == 1 && (isscalar(unique(d.LU.SID)) && isscalar(unique(d.LU.EID)))
    error('目前为MilkRun版本,不允许只有单个供应商编号/单个EP LOCATION编号');
end            % if verMilkRun == 1 && (~isscalar(unique(d.LU.SID)) && ~isscalar(unique(d.LU.EID))),  error('目前为MilkRun测试版本,只允许多个供应商编号 或 多个EP LOCATION编号,不能同时存在');  end
%% 没有属性的临时增加 - CAN DEL
%     n = numel(d.LU.Weight);

%     if ~isfield(d.LU, 'Index') %% 给个初始进入顺序 暂时没有用 等同刘强所需index参数
%         error('No index of LU');  %d.LU.Index = 1:n;
%     end
    
%     % 平铺使用属性
%     if ~isfield(d.LU, 'maxL')
%         d.LU.maxL(1,:) =  floor(d.Veh.LWH(1,1)./d.LU.LWH(1,:));
%         d.LU.maxL(2,:) =  floor(d.Veh.LWH(2,1)./d.LU.LWH(2,:));
%         d.LU.maxL(3,:) =  floor(d.Veh.LWH(3,1)./d.LU.LWH(3,:));   %具体每个托盘LU的高度的最大层数
%     end
    
%     if ~isfield(d.LU, 'maxHLayer'),     d.LU.maxHLayer = d.LU.maxL(3,:); end% maximum given height layer
                     
%% Simulate - All ALGORITHM

fprintf(1,'\nRunning the simulation...\n');

for iAlg = 1:nAlg
    %% 1 运行主算法 输入d 输出 do
    
    % 1.1 修订LU的LWH数据 包含margin 考虑Rotaed
    d = maind;
    [d.LU] = setLULWHwithbuff(d.LU, d.Veh);
    
    % 1.2 主函数RunAlgorithm
    do = RunAlgorithm(d,pA(iAlg));   %获取可行解结构体 ( 在此之前获取固定LU方向和固定LU大小的LWH）   % 可删 - 已在RunAlgorithm中声明 do.LU.LU_VehType = ones(size(do.LU.ID)) * do.Veh.order(1); % 针对车型选择,增加变量LU_VehType : 由于Veh内部按体积递减排序,获取order的第一个作为最大值

    % 1.3 修订LU和Item的LWH/Coord数据 删除margin
    [do.LU,do.Item] = setLCwithoutbuff(do.LU,do.Item);   
    
    % 1.4 CHECK 输入和输出的LU
    chkLUnewold(d.LU,do.LU);  %新老数据对比 todo 完善新老数据对比 预防算法内部错误
    chktLU(do.LU); %核验是否数据是上轻下重；是否Z高度与ItemSEQ一致；是否ITEMID与XY坐标一致; todo :完善结果核对函数
    
    % 1.5 % dA(iAlg)=do; % 仅在多算法版本使用    
              %     plotSolutionT(do.LU,do.Veh);
              % plotSolution(do,pA(iAlg)); %尽量不用
    
    %% 2 运行车型调整算法,不改变d 获取d1和do1, flaggetSmallVeh
    if ISlastVehType
        [flaggetSmallVeh,do1] = HBinChange(maind,do,pA(iAlg));  
    end

    %% 3 运行堆垛平铺（含正常平铺和甩尾平铺）算法,不改变d 来自maind 
    % 获取 d2Array(甩尾平铺） d3Array（整车平铺） 的输入数据，每个是一个bin
    % 获取 do2Array do3Array  的输出数据，每个是一个bin
    % 获取 flagTiledArray 逻辑判定，1 本bin甩尾成功 0 不成功
    if ISpingpu==1
        [flagTiledArray,do2Array,do3Array] = HBinpingpu(maind,do,pA(iAlg));            
    end
    
    %% 4 运行间隙调整算法 输出do3
    if parGap == 1
        
        if ISpingpu==1 && any(flagTiledArray==1) %（仅允许甩尾平铺和整车平铺成功的车进行）todo 开放给所有bin
            idxArray = find(flagTiledArray~=0);
            for i=1:length(idxArray)
                if flagTiledArray(idxArray(i)) == 1
                    dd = do3Array(idxArray(i)); 
                elseif flagTiledArray(idxArray(i)) == 2
                    dd = do2Array(idxArray(i)); 
                end
                
                [dd.LU] = HBinGapAdjust(dd.LU, dd.Veh);
                
                if flagTiledArray(idxArray(i)) == 1
                    do3Array(idxArray(i)) = dd;
                elseif flagTiledArray(idxArray(i)) == 2
                    do2Array(idxArray(i)) = dd;
                end                
            end            
        end
    end
    
end
   
%% Simulate - CHOOSE BEST ONE
% 555 算法首先判断并排除bin内相同类型托盘不相邻的解 TODO 数据的CHECK
%  flagA(iAlg) =  isAdjacent(dA(iAlg));           % 算法判断是否相同类型托盘相邻摆放 +
% % dA = dA(1,logical(flagA));
% % pA = pA(1,logical(flagA));

% TODO 从多次算法结果中选出从必定bin内相邻的最优结果 - NOTE: 采用单参数时无需考虑
% if isempty(dA), error('本算例内所有解都存在托盘不相邻的情况 \n'); end
% [daBest,paBest] = getbestsol(dA,pA);  %可不考虑d1 -> 解的优劣与最后一个bin关系不大
%%%% Return length(parMax) 个 solutions to BBA
% if isempty(daBest), error('本算例内未找出最优解返回BBA \n'); end
% bestOne = 1;
                                %%%% dA = do = daBest(1) = daBest(bestOne) %isequal(do,d1)
%% POST PROCESSING
% 基于do do1 do2 do3 flagTiledArray flaggetSmallVeh 等结果数据进行处理
% 基于table格式

if ISlastVehType && ISpingpu
    T = HBinCombine(do,flaggetSmallVeh,do1,flagTiledArray,do2Array,do3Array);  
elseif ISlastVehType && ~ISpingpu
    T = HBinCombine(do,flaggetSmallVeh,do1);
elseif ~ISlastVehType && ISpingpu
    T = HBinCombine(do,[],[],flagTiledArray,do2Array,do3Array);
elseif ~ISlastVehType && ~ISpingpu
    T = HBinCombine(do);
end

%% return必须的参数
% 行1：托盘所在车型号(必须换)    行2：托盘所在车序号(会变,不能换,换就错) 行3：托盘车内安置顺序(必须换) 行4：托盘SID供应商编号(不会变,不用变??)
% 行5：托盘ID型号LID(不会变,不用变?) 行6：托盘堆垛序号ITEM(会变,不能换,换就错,用途?) 行7：托盘零部件编号PID(不会变,不用变?) 增加行8: 展示顺序(必须换)

output_CoordLUBin=T.CoordLUBin';
output_LU_LWH=T.LWH';
% V2 get output_LU_Seq
if verMilkRun == 1 
output_LU_Seq=T{:,{'LU_VehType','BINID','BINSEQ','OSID','LID','ITEMID','OPID','ShowSEQ','Weight','Index','OEID'}}'; %增加返回行10: LuIndex来自刘强只要是数字就可以
else % 带回去OPID OSID OEID 等 其实用途也不大 %用途大的是
output_LU_Seq=T{:,{'LU_VehType','BINID','BINSEQ','OSID','LID','ITEMID','OPID','ShowSEQ','Weight','Index'}}'; %增加返回行10: LuIndex来自刘强只要是数字就可以
end
% V1
% output_LU_Seq=T{:,{'LU_VehType','BINID','BINSEQ','OSID','LID','ITEMID','OPID','ShowSEQ','Weight'}}'; % ITEMID意义不大    output_LU_Seq([2,3,5,8],:)

% plot
if ISplotBBA
    %     plotSolutionBBA(output_CoordLUBin,output_LU_LWH,output_LU_Seq,do); 
    V = struct2table(structfun(@(x) x',d.Veh,'UniformOutput',false));
    plotSolutionT(T,V); %按table格式做图
end

    % if  ISplotSolu 
    % %      plotSolution(daBest(bestOne),paBest(bestOne)); %尽量不用 包含plotStrip 不包含单车型作图
    % %      if max(do2.LU.LU_Bin(1,:)) == 1
    % %      plotSolution(do2,pA(iAlg));
    % %      end
    % %        if flaggetSmallVeh,   plotSolution(do1,paBest(bestOne));   end %尽量不用 包含plotStrip 仅包含单车型作图
    % end

    
%% 常用语句
% T.Properties.VariableNames
% % 简单查看某个bin内的部分子表
% tmpT = sortrows(T,{'ITEMID'})
% x=T(T.BINID==2&T.ID==1,{'ID','LID','PID','H','Weight','CoordLUBin','BINSEQ','ShowSEQ','ITEMID','ITEMSEQ'})
% sortrows(x,{'CoordLUBin'})
% 剔除展示顺序 % output_LU_Seq = output_LU_Seq(1:7,:);
% whos

clearvars -except output*
fprintf(1,'Simulation done.\n');

% mcc -W 'java:BBA_Main,Class1,1.0' -T link:lib BBA_Main.m -d '.\new'
% mcc -W 'java:BBA_MR_Main,Class1,1.0' -T link:lib BBA_MR_Main.m -d '.\new'

% printstruct(do,'sortfields',1,'PRINTCONTENTS',0);    printstruct(do.Veh);
% do = rmfield(do, {'Veh', 'LU'});
% pcode 'H*.m'.
% lu = table2struct(T,'ToScalar',true)
% lu = (structfun(@(x) x',lu,'UniformOutput',false));
end %END MAIN




%% ********************** 下面是基本无用的代码 暂时不用 ****************



%% dd: 需要调整Gap的Bin
% function [flagGap,T23LU] = getMixedGap1(dd)
% % Initilize output
% flagGap=0;
% T23LU = struct2table(structfun(@(x) x',dd.LU,'UniformOutput',false));
% 
% % find and sort 'fidx' (fidx: 有Gap的Strip的Index)
% fidx=find(dd.Strip.isGapBalance==0);  %fidx: 存在混装间隙的strip的序号
% [~,ff]=sort(dd.Strip.Strip_Bin(2,fidx)); %
% fidx=fidx(ff); %从最小（最靠近车顶）的进入顺序进行, 进入bin内是12在前,11在后
% 
% for i=1:length(fidx)
%     flagGap = 0;
%     idxs = fidx(i);  % strip12 找到对应的LU LU_Strip 第12个 它的LWH和它的CoordLUBin CoordLUStrip
%     
%     %% 1: pgLUinStrip: strip对应的pg 和 pgLUinStrip: 此strip内lu对应的pg -> pgBlanksinStrip:消减后的pg
%     %pgLUinStrip: strip内的LU的多边形,可能是不规则的或是多个的
%     pgLUinStrip = pgStripLU(idxs,T23LU);    % pgLUinStrip.NumRegions
%     
%     %pgStrip: 此strip的多边形,一定是个矩阵
%     stripWidth = dd.Veh.LWH(1, unique(T23LU.LU_VehType));     %strip宽度
%     stripHeight = dd.Strip.LW(2,idxs);                                       %strip高度
%     [x,y]=boundary(pgLUinStrip);
%     pgStrip = polyshape(pgRectangle(min(x),min(y),stripWidth,stripHeight));  if size(pgStrip.Vertices,1) ~= 4, error('strip不是矩形'); end
%     
%     %pgBlanksinStrip: strip内的剩余空白区域,每个区域必须是矩阵
%     pgBlanksinStrip = regions(subtract(pgStrip,pgLUinStrip)); % 获取本strip内的剩余可用区域    
%     %pgBlanksinStrip排序, （目前: 面积小的在前面） todo: y小的在前面
%     pgBlanksinStrip = sortregions(pgBlanksinStrip,'area','ascend');
%     %                 plot(pgBlanksinStrip)
%     %                 hold on
%     %                 plot(pgLUinStrip)
%     
%     %% 2: idxLUArray：找到idxs后一个strip内的LU,依据LU的大小判定能否放下
%     % Find nextStrip
%     nextStrip = dd.Strip.Strip_Bin(2,idxs) + 1;                                 %此strip在bin内的进入顺序 + 1：即下一个strip的序号
%     nextStripIdx = find(dd.Strip.Strip_Bin(2,:) == nextStrip);
%     
%     % Get and sort idxLUArray by area of LU in nextStrip  %更新LU的顺序，从大面积到到小面积
%     idxLUArray = find(T23LU.LU_Strip(:,1)==nextStripIdx);                   %在混装间隙下一个strip内包含的LU序号 
%     
%     idxLUareaArray = T23LU.LWH(idxLUArray,1).*T23LU.LWH(idxLUArray,2);    % areas
%     [~, ord] = sort(idxLUareaArray,'descend');
%     idxLUArray=idxLUArray(ord);
%     
%     % TODO 找出能放的下的region，并用某种规则（横竖）放入.
%     % 逐个LU尝试放入pgBlandsinStrip
%     for j=1:length(idxLUArray)
%         % 从第j个LU（LU应该排序，从面积大的）开始, 尝试放入到上面的strip去
%         l = idxLUArray(j);
%         
%         for i = 1:length(pgBlanksinStrip)
%             % 2.1 get pgBlank
%             pgBlank= pgBlanksinStrip(i);
%             
%             % 2.2.1 get pgLU
%             % Find pgBlank's 顶点坐标 x and y) 
%             [originx,originy] = boundary(pgBlank);
%             originLUx=min(originx);
%             originLUy=min(originy);
%             
%             % 2.2.2 get pgLU 不旋转 基于pgBlank's Region
%             w=T23LU.LWH(l,1);
%             h=T23LU.LWH(l,2);
%             pgLU = polyshape(pgRectangle(originLUx,originLUy,w,h));
%             
%             % 2.2.3 get pgLURota 旋转 基于pgBlank's Region
%             if T23LU.isRota(l)
%                 wRota = T23LU.LWH(l,2);
%                 hRota = T23LU.LWH(l,1);
%                 pgLURota = polyshape(pgRectangle(originLUx,originLUy,wRota,hRota));
%             end
%             
%             % 2.3 if pgLU/pgLURota can be assigned into pgBlank
%             [xLU,yLU] = boundary(pgLU);
%             [xLURota,yLURota] = boundary(pgLURota);
%             xLU
%             yLU
%             %% 判定是否可以放下的重要条件
%             % 555: if all of pgLU's boundary points belong to pgBlank
%             flagLUx =  all(isinterior(pgBlank,xLU([1,4]),yLU([1,4])));
%             flagLU =  all(isinterior(pgBlank,xLU,yLU));
%             flagLURotax =  all(isinterior(pgBlank,xLURota([1,4]),yLURota([1,4])));
%             flagLURota =  all(isinterior(pgBlank,xLURota,yLURota));
%             
%             % 3 如果可以替换,
%             if flagLU || flagLUx
% %                 [x,y]=boundary(pgBlank);
% %                 originLUx = min(xLU)
% %                 originLUy = min(yLU)
%                 %TODO : update 第l个LU的strip序号/进入顺序.CoordLUStrip坐标等
%                 T23LU.CoordLUBin(l,1) = originLUx;
%                 T23LU.CoordLUBin(l,2) = originLUy;
%                 T23LU.LWH(l,1)=w;
%                 T23LU.LWH(l,2)=h;
%                 flagGap=1;
%                 break
%             end
%             
%             if flagLURota || flagLURotax
% %                 [x,y]=boundary(pgBlank);
% %                 originLUx = min(xLU);
% %                 originLUy = min(yLU);
%                 %TODO : update 第l个LU的strip序号/进入顺序.CoordLUStrip坐标等
%                 T23LU.CoordLUBin(l,1) = originLUx;
%                 T23LU.CoordLUBin(l,2) = originLUy;
%                 T23LU.LWH(l,1)=h;
%                 T23LU.LWH(l,2)=w;
%                 T23LU.Rotaed(l)=~T23LU.Rotaed(l);
%                 flagGap=1;
%                 break
%             end
%             
%             %     area(pgLU)                         area(pgBlank)
%             
%             %                         figure
%             %                         plot(pgLU)
%             %                         hold on
%             %                         plot(pgBlank)
%             
%             %                     area(pg12)
%             %                     area(pg1)
%             
%         end
%         if flagGap,   break;    end
%     end
%     if flagGap,   break;    end
%     
% end
% end






function lastd = getdinLastVeh(tmpd)
    % tmpd中的Bin是排序后的, 从最小的开始试
    tmpusedVehIdx = max(tmpd.LU.LU_Bin(1,:)); %tmpusedVehIdx: 最后一个Bin的index值
    flagusedLUIdx = tmpd.LU.LU_Bin(1,:)==tmpusedVehIdx; % flagused: 找出最后一个Bin对应的LUindex值
    if isSameCol(tmpd.LU)
        % 获取仅最后一个Bin的输入数据
        lastd.LU = structfun(@(x) x(:,flagusedLUIdx),tmpd.LU,'UniformOutput',false);  %仅取最后一辆车内的LU
        lastd.LU.LWH([1,2], lastd.LU.Rotaed ) = flipud(lastd.LU.LWH([1,2], lastd.LU.Rotaed)); %LU.LWH 如旋转,则恢复原形
        lastd.LU = rmfield(lastd.LU,{'Rotaed','order','LU_Item','DOC','LU_Strip','LU_Bin','CoordLUBin','maxL','CoordLUStrip'}); 
        lastd.Par = tmpd.Par;
    else
        error('不能使用structfun');
    end
end

%% **** 算法指标选择最优解 ****
function [daMax,parMax] = getbestsol(DaS,Par)

    % 如果仅有一个可行解, 直接返回;
    if size(DaS,2)==1    %仅当dA有多次时采用,目前参数锁定,
        daMax = DaS(1); parMax = Par(1);
        return
    end
    
%获取评价指标和对应参数
for r=1:length(DaS)
    resLoadingRateBin(r) = mean(DaS(r).Bin.loadingrate); %bin的装载率均值最大 Itemloadingrate ItemloadingrateLimit
    resLoadingRateStripLimit(r) = mean(DaS(r).Strip.loadingrateLimit); %strip的limit装载率最大 Itemloadingrate ItemloadingrateLimit
    resLoadingRateBinLimit(r) = mean(DaS(r).Bin.loadingrateLimit); %bin的limit装载率最大 Itemloadingrate ItemloadingrateLimit
    resLoadingRateStrip(r) = mean(DaS(r).Strip.loadingrate); %strip的装载率最大 Itemloadingrate ItemloadingrateLimit    
%     Par(r);
end

%% 算法选择最优的解给用户
% maxresBin=max(resLoadingRateBinLimit(1,idxStrip)); %找出idxStrip中的最大bin
% if ~all(ismember(idxBin,idxStrip)),   error('not all member of bin in strip'); end %错误有可能出现 
%% 1 maxresBin代表常规车辆的平均装载率,物品总量一定,bin越多,该值越小,解越差,此最大值是必须
%% 2 maxresStrip代表常规Strip的平均装载率,物品总量一定,strip宽度一定,高度越高,该值越小,解越差,此最大值不一定时必须（因为看起来不好看）
%% 但该值在不同bin高度时可能有影响，且该值好时，人为看起来可能并不好
%% 3 maxresStripLimit代表特殊Strip的平均装载率,物品总量一定,strip内部宽度越大,间隙越大,值越小,此最大值几乎是必须
%% 该值好时，人为看起来可能好（strip内部间隙小）；但不一定时最优（还有可能相同托盘不在一起）
%% 4 maxresBinLimit代表特殊Bin的平均装载率,物品总量一定?? 对特殊情况有用,待观察
idxBin=find(resLoadingRateBin==max(resLoadingRateBin)); %取
idxStripLimit=find(resLoadingRateStripLimit==max(resLoadingRateStripLimit));
idxBinLimit=find(resLoadingRateBinLimit==max(resLoadingRateBinLimit));
idxStrip=find(resLoadingRateStrip==max(resLoadingRateStrip));
%% 5 找出idxStrip和idxBin两者的交集
% % if isempty(intersect(idxBin,idxStrip))
idx =idxBin;
if isempty(idx), error('idxBin为空 '); end %错误几乎不可能出现
idx0 =intersect(idx,idxStripLimit);
if ~isempty(idx0),
    idx = idx0; 
else
    warning('idx0 is empty');
end
% if isempty(idx), error('idxBin and idxStripLimit 的交集为空 '); end %错误几乎不可能出现
idx1 = intersect(idx,idxBinLimit);
if ~isempty(idx1),  
%      idx = idx1; 
else
    warning('idx1 is empty');
end
idx2 = intersect(idx,idxStrip);
if ~isempty(idx2),  
%     idx = idx2; 
else
    warning('idx2 is empty');
end

%% 将idx剩余的返回到主函数
if ~isempty(idx)
    for tmpidx=1:length(idx)
        daMax(tmpidx) = DaS(idx(tmpidx));
        parMax(tmpidx) = Par(idx(tmpidx));
    end
end

end % END OF ALL





        
function plotSolution(do,par)
%% 画图
% V3 margin 提前到RunAlgorithm运行后就执行:
% plot2DBPP(do,par);
plot3DBPP(do,par);

        % V1 buff version
        % do.Item.LWH = do.Item.LWH - do.LU.buff(:,1:size(do.Item.LWH,2));
        % do.Item.LWH(1,:) = do.Item.LWH(1,:) - ( do.LU.margin(1, 1:size(do.Item.LWH,2) ) + do.LU.margin(2,: )); 
        % do.Item.LWH(2,:) = do.Item.LWH(2,:) - (do.LU.margin(3,: ) + do.LU.margin(4,: )); 
        % do.Item.CoordItemBin = do.Item.CoordItemBin + do.LU.buff(:,1:size(do.Item.LWH,2))/2;

        % V2 margin version
        % 作图前更新LU ITEM的Coord和LW; 更新ITEM同时更新LU
        % [do.LU,do.Item] = setLCwithoutbuff(do.LU,do.Item);
        

end


    %% ************ 判断是否相同类型托盘相邻摆放
%     function flag = isAdjacent(d)
%         flag = 1;
%         printstruct(d);
%         % 每个bin中找出各类型ID所在Strip是否相邻
%         nBin = size(d.Bin.LW,2);
%         for iBin = 1:nBin
%             t = [d.Item.LID; d.Item.Item_Strip; d.Item.Item_Bin ];
%             tiBin = t( : , t(4,:) == iBin );
%             nIdType = unique(tiBin(1,:)); %nIdType: 本iBin内包含的LU的ID类型
%             for iId = 1:nIdType
%                 tiId = tiBin( : , tiBin(1,:) == iId );
%                 nIdStrip = unique(tiId(2,:)); %nIdStrip: 本iBin及本iID下包含的Strip的序号
%                 % 判断排序后的放入本ID类型的Strip序号是否相邻
%                 if ~all(diff(sort(nIdStrip))==1)
%                     flag = 0;
%                 end
%             end                    
%         end
%     end
    
% % %             % ns - 本bin内strip个数及顺序
% % %             ns = d.Strip.stripBeBinMatrix(2,d.Strip.stripBeBinMatrix(1,:) == iBin);
% % %             % ni - 本bin内item内LU类型及顺序
% % %             d.Item.Item_Bin(1,:) == iBin
% % %             ni = d.Item.LID(d.Item.Item_Bin(1,:) == iBin);
% % %             [a,b] = find(d.Item.LID(d.Item.Item_Bin(1,:) == iBin));
% % %             ni_uni = unique(ni);
% % %             for ini = 1:length(ni_uni)
% % % %                 d.Item.
% % % %                 d.Item.Item_Strip(:,
% % %             end
% % %             nStrip = length(ns);
% % %             % i,j is adjacent strips(levels)
% % %             for iStrip = 1:nStrip
% % %                 for jStrip = (iStrip+1):(nStrip-1)
% % %                 [is] = find(ns==iStrip); %第3个strip放第1层
% % %                 [js] = find(ns==jStrip); %第1个strip放第2层
% % %                 LUIDInis = d.Item.LID(1,(d.Item.Item_Strip(1,:)==is))
% % %                 LUIDInjs = d.Item.LID(1,(d.Item.Item_Strip(1,:)==js))
% % %                 
% % %                 end
% % %             end
% % %         end
% % %         



%% ********************** 下面是ts算法的代码 暂时不用 ****************

%%

    %% 注释
    
%         while (all(do2.Strip.isHeightFull(fi) == 1) && all(do2.Strip.isWidthFull(fi) == 1))
% %                  || all(d2.LU.maxHLayer(luidxPP)==1)% 目前仅能对最后一个strip调整, 或增加最后一个Strip内的Lu的maxHLayer全部为1
%              
%             istrip = istrip-1;
%             fi = find( do2.Strip.Strip_Bin(2,:) >= istrip ); 
%             luidxPP = ismember(do2.LU.LU_Strip(1,:), u(fi)); %%% fi->u(fi) 真正的序号 ********************* 
%             if ~any(luidxPP),  error('luidxPP全部为空, 不存在u(fi)对应的Lu逻辑判断'); end
%             if istrip == 1,  error('此bin不存在tileneed,超预期错误');   end
%             d2.LU.maxHLayer(luidxPP)
%         end
        
        
        % 部分LU(luidxPP),修改其maxHLayer
%         u=unique(do2.LU.LU_Strip(1,:)); %获取Strip序号的唯一排序值
%         luidxPP = ismember(do2.LU.LU_Strip(1,:), u(fi)); %%% fi->u(fi) 真正的序号 *********************
%         if ~any(luidxPP),  error('luidxPP全部为空, 不存在u(fi)对应的Lu逻辑判断'); end
%         do2.LU.LU_Item(1,luidxPP)


%                     d2.LU.tmpHLayer = zeros(size(d2.LU.maxHLayer))
%                     d2.LU.tmpHLayer(luidxPP) = do2.Item.HLayer(do2.LU.LU_Item(1,luidxPP))
%                     dlu=[ d2.LU.maxL(3,luidxPP);
%                      d2.LU.maxHLayer(luidxPP);
%                      d2.LU.tmpHLayer(luidxPP)];
%                      min(dlu)
%         d2.LU.maxHLayer(luidxPP) = min(dlu) - 1;

% % % [ub,x,b] = HnextFit(Item,Veh);
% % % [ub,x,b] = HnextFit_origin(Item,Veh);
% % % disp(b');
% % % disp(x);
% % % fprintf('UB = %d \n', ub);
% % % [ub,x,b] = TSpack(d,n,w,W,lb,timeLimit, ub0,x,b,whichH);
% % 
% % function [toReturn,x,b] = TSpack(d,n,w,W,lb,timeLimit, ub0,x,b,whichH)
% % [nb,x,b] = heur(d,n,w,W,x,b,whichH);
% % ub0 = nb;
% % 
% % if nb == lb
% %     toReturn = nb;
% %     fprintf('best lb = %d ', nb);
% % end
% % 
% % %/* initial (trivial) solution */
% % cnb = n;
% % cb = 1:n;
% % cx = zeros(d,n);
% % cw = w;
% % 
% % %/* external loop */
% % D = 1; tt = 0.0;
% % toReturn = nb;
% % 
% % end
% % 
% % function [nb,x,b] = heur(d,n,w,W,x,b,whichH)
% % nb = -1;
% % which = (d-2)*100 + whichH; % which 为0或100
% % if which == 0
% %     [nb,x,b]  = HnextFit(n,w,W,x,b,n+1); %/* first heuristic for 2d bin packing */
% %     %          disp(x);
% %     %          disp(b);
% % elseif which == 100
% %     [nb,x,b]  = HHnextFit(n,w,W,x,b); %/* first heuristic for 3d bin packing */
% % end
% % end
% % 
% %     function [ub,px,pb]  = HnextFit(Item,Veh)
% %         % Initialize
% %         d = size(Item.LWH,1)-1;
% %         n = size(Item.LWH,2);
% %         nn = n + 1;
% %         w = Item.LWH(1:d,:);
% %         W = Veh.LWH(1:d,:);
% %         x = zeros(d,n); b = zeros(n,1); bNb = zeros(n,1);
% %         
% %         %/* sort the items */
% %         % sortD = size(w,1);%获取需要排序的维度
% %         [~,ord] = sort(w(d,:),'descend');%对w进行排序,只需要它的顺序ord;按第d行排序（高度)
% %         pw = w(:,ord);
% %         px = x;
% %         pb = (b+999); % 0 + 999
% %         pbNb = bNb;
% %         
% %         %/* next fit packing */
% %         % binLeftArray(1,ub) ： wleft
% %         % binLeftArray(2,ub) :  hleft
% %         nBin = n;
% %         binLeftArray = repmat(W,1,nBin);  %初始
% %         ub = 1;
% %         for i=1:n
% %             %     if (binLeftArray(1,ub) == W(1)) & (binLeftArray(2,ub) == W(2)) %如果是空bin
% %             if pbNb(ub) == 0   %如果是空bin
% %                 if (pw(1,i) <= binLeftArray(1,ub)) && (pw(2,i) <= binLeftArray(2,ub)) %如果宽高都不超标
% %                     px(1,i) = 0; px(2,i) = 0;
% %                     binLeftArray(1,ub) = binLeftArray(1,ub) - pw(1,i);
% %                     binLeftArray(2,ub) = binLeftArray(2,ub) - pw(2,i);
% %                     pbNb(ub) = pbNb(ub) + 1;
% %                 else
% %                     error('EEE');
% %                 end
% %             else               %如果不是空bin
% %                 if pw(1,i) <= binLeftArray(1,ub)  %如果i的宽满足当前bin的剩余宽度，剩余高度应该不变
% %                     px(1,i) = W(1) - binLeftArray(1,ub);
% %                     px(2,i) = W(2) - binLeftArray(2,ub) - pw(2,i);     %高度为????
% %                     binLeftArray(1,ub) = binLeftArray(1,ub) - pw(1,i);
% %                     binLeftArray(2,ub) = binLeftArray(2,ub);
% %                     pbNb(ub) = pbNb(ub) + 1;
% %                 else
% %                     if pw(2,i)  <= binLeftArray(2,ub)  %如果i的高满足当前bin的剩余高度
% %                         px(1,i) = 0;
% %                         px(2,i) = binLeftArray(2,ub);
% %                         binLeftArray(1,ub) = W(1) - pw(1,i);
% %                         binLeftArray(2,ub) = binLeftArray(2,ub) - pw(2,i);
% %                         pbNb(ub) = pbNb(ub) + 1;
% %                     else  %如果i的高不能满足当前bin的剩余高度
% %                         ub = ub + 1;
% %                         px(1,i) = 0;   px(2,i) = 0;
% %                         pbNb(ub) = pbNb(ub) + 1;
% %                         binLeftArray(1,ub) = binLeftArray(1,ub) - pw(1,i);
% %                         binLeftArray(2,ub) = binLeftArray(2,ub) - pw(2,i);
% %                     end
% %                 end
% %             end
% %             pb(i) = ub-1;
% %         end
% %         
% %         
% %         %原始的
% %         x = px(:,ord);
% %         b = pb(ord);
% %         bNb = pbNb(ord);
% %         ub = ub +1;
% %     end
% % 
% %     function [ub,px,pb]  = HnextFit_origin(Item,Veh)
% %         % Initialize
% %         d = size(Item.LWH,1);
% %         if d==3
% %             d=d-1;
% %         end
% %         n = size(Item.LWH,2);
% %         nn = n + 1;
% %         w = Item.LWH(1:d,:);
% %         W = Veh.LWH(1:d,:);
% %         x = zeros(d,n); b = zeros(n,1); bNb = zeros(n,1);
% %         
% %         %/* sort the items */
% %         sortD = size(w,1);%获取需要排序的维度
% %         [~,ord] = sort(w(sortD,:),'descend');%对w进行排序,只需要它的顺序ord
% %         pw = w(:,ord);
% %         ord;
% %         
% %         px = zeros(size(x,1),size(x,2));
% %         pb = 999*ones(size(b,1),size(b,2));
% %         %/* next fit packing */
% %         hleft = W(2) - pw(2,1);
% %         wleft = W(1);
% %         ub = 0; hcurr = 0;
% %         for i=1:n  %从第一个item开始安置
% %             if pw(1,i) <= wleft  %如果item的w 比wleft小，安置item到本bin本层：更新x1值，更新wleft；hleft不变
% %                 px(1,i) = W(1) - wleft;
% %                 wleft = wleft - pw(1,i);
% %             else    %否则往上一层安排。
% %                 if pw(2,i) <= hleft  %如果item的h 比hleft小：表明bin高度充足，安置item到上一曾：更新坐标hleft，更新hcurr，wleft？ 更新坐标x值，更新wleft
% %                     hcurr = W(2) - hleft; %安排在同一层，所以hcurr不变，也等于pw(2,1)(但在其他bin就不对了，所以用hcurr)
% %                     hleft = hleft - pw(2,i);
% %                 else  %如果放不下，开新bin，更新hcurr；更新hleft；更新数量ub+1（如果达到nn值，跳出）；更新坐标x值0，更新wleft
% %                     hcurr = 0;    %安排在新的bin，所以hcurr为0;
% %                     hleft = W(2) - pw(2,i);
% %                     if (ub+1 == nn)
% %                         break;
% %                     end
% %                     ub = ub + 1;
% %                 end
% %                 % 无论放在上层或开新bin，更新x1为0；更新wleft为W(1)-此item的宽w
% %                 px(1,i) = 0;
% %                 wleft = W(1) - pw(1,i);
% %             end
% %             % 此处统一更新x1值，即高度值，为hcurr；统一更新b值=ub；
% %             px(2,i) = hcurr;
% %             pb(i) = ub;
% %         end
% %         
% %         %原始的
% %         x = px(:,ord);
% %         b = pb(ord);
% %         ub = ub +1;
% %     end
% % 
% % 
% %     function [ nb ] = HHnextFit(n,w,W,x,b)
% %         
% %         nb = 1;
% %         
% %     end



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

% % %% function [StripSolutionSort] = HnextFitDH(d)
% % function [StripSolutionSort] = HnextFitDH(d)
% % % 输入: d
% % % 输出: StripSolutionSort
% % %% 提取单类型bin,二维item数据
% % % nDim nItem nBin
% % % itemDataMatrix uniBinDataMatrix
% % nDim = size(d.Item.LWH,1);  if nDim ==3, nDim = nDim-1;end
% % itemDataMatrix = d.Item.LWH(1:nDim,:);
% % tmpbinDataMatrix = d.Veh.LWH(1:nDim,:);
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
% % %% function [StripSolutionSort] = HfirstFitDH(d)
% % function [StripSolutionSort] = HfirstFitDH(d)
% % % 输入: d
% % % 输出: StripSolutionSort
% % %% 提取单类型bin,二维item数据
% % % nDim nItem nBin
% % % itemDataMatrix uniBinDataMatrix
% % nDim = size(d.Item.LWH,1);  if nDim ==3, nDim = nDim-1;end
% % itemDataMatrix = d.Item.LWH(1:nDim,:);
% % tmpbinDataMatrix = d.Veh.LWH(1:nDim,:);
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
% % %% function [StripSolutionSort] = HbestFitDH(d)
% % function [StripSolutionSort] = HbestFitDH(d)
% % % 输入: d
% % % 输出: StripSolutionSort
% % %% 提取单类型bin,二维item数据
% % % nDim nItem nBin
% % % itemDataMatrix uniBinDataMatrix
% % nDim = size(d.Item.LWH,1);  if nDim ==3, nDim = nDim-1;end
% % itemDataMatrix = d.Item.LWH(1:nDim,:);
% % tmpbinDataMatrix = d.Veh.LWH(1:nDim,:);
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
% % function [stripLeftMatrix,pbelongMatrix,pitemMatrix,pcoordMatrix,ord]  = HfirstFitDH2(Item,Veh)
% % % 输入参数初始化
% % nDim = size(Item.LWH,1);
% % if nDim ==3, nDim = nDim-1;end
% % nItem = size(Item.LWH,2);
% % nBin = nItem;
% % % nn = n + 1;
% % itemMatrix = Item.LWH(1:nDim,:);
% % binMatrix = Veh.LWH(1:nDim,:);
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
% % function [stripLeftMatrix,pbelongMatrix,pitemMatrix,pcoordMatrix,ord]  = HbestFitDH2(Item,Veh)
% % % 输入参数初始化
% % nDim = size(Item.LWH,1);
% % if nDim ==3, nDim = nDim-1;end
% % nItem = size(Item.LWH,2);
% % nBin = nItem;
% % % nn = n + 1;
% % itemMatrix = Item.LWH(1:nDim,:);
% % binMatrix = Veh.LWH(1:nDim,:);
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
% % function [pbelongItemBinMatrix,pbelongStripBinMatrix,pcoordItemBinMatrix,binLeftMatrix ] = HbestFitBinDH(stripLeftMatrix,pbelongMatrix,pitemMatrix,pcoordMatrix,Item,Veh)
% % % 输入参数初始化
% % nDim = size(Item.LWH,1);
% % if nDim == 3, nDim = nDim-1;end
% % nItem = size(Item.LWH,2);
% % nBin = nItem;
% % 
% % nStrip = sum(stripLeftMatrix(3,:)>0); %具体使用的Strip的数量
% % % nn = n + 1;
% % binMatrix = Veh.LWH(1:nDim,:);
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
% %         tepMin = min(tepMin(findBinArray)); % 555 chk
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
%% Call 本函数:
% clear;close all; format long g; format bank; %NOTE 不被MATLAB CODE 支持
% rng('default');rng(1); % NOTE 是否随机的标志
% LU = struct('ID',[],'LWH',[],...
%     'weight',[],'Lbuffer',[],'Wbuffer',[],'Type',[],'Material',[]);
% Item = struct('ID',[],'LWH',[],...
%     'weight',[],'Lbuffer',[],'Wbuffer',[]);
% Strip = struct('ID',[],'LWH',[],... %ONLY LW
%     'weight',[],'Lbuffer',[],'Wbuffer',[]);
% Veh = struct('ID',[],'LWH',[],...
%    'Capacity',[],'Lbuffer',[],'Wbuffer',[],'Hbuffer',[]);
% d = struct('LU',LU,'Item',Item,'Strip',Strip,'Veh',Veh);

% 1参数初始化
% whichStripH 1 best 2 first 3 next; whichBinH 1 best; TODO 增加其它分批方式
% whichSortItemOrder 1 长高递减 2 最短边递减; 
% whichRotation 1:允许rotation 0:禁止
% rotation组合 1 1 2 1 0 (1 1 2 1 1 )(1 1 2 1 2) % 非rotation组合 1 1 1 0 0 （2/3 1 1 0 0）
% whichRotationHori 0:在安置顺序时按FBS_{RG}方式; 1：New/NoNew按Horizon方式 2：New/NoNew按Vertical方式
% ParaArray = struct('whichStripH',1,'whichBinH',1,'whichSortItemOrder',2,...
%     'whichRotation',1,'whichRotationHori',0,'timeLimit',100,'ub0',10);
% % ParaArray = struct('whichStripH',1,'whichBinH',1,'whichSortItemOrder',2,...
% %     'whichRotation',1,'whichRotationHori',0,'whichRotationAll',1,'whichRotationBin',1,'timeLimit',100,'ub0',10);

% % nAlg = 1;
% % for i = 3:3 %1-3 best first next均可 设为3: 不允许前面小间隙放其它东西 因为一旦允许, 会大概率违背相邻约束
% %     for j=3:3 %0-3 排序: 0: Vert；1: Hori; 2:error  3:按缝隙最小排序   Gpreproc 此处替代HItemToStrip函数中的物品摆放
% %         for k=2:2 %0-2 默认0 不可旋转 1全部可旋转 2: 按人为设置是否允许Rotation 
% %             for TLUl=1:1 % 已无用 :  % 0-2 0已取消 保留1-2 RotaHori 1hori 2 vert 555 横放不了会纵放，不允许；纵放后不会横放（放不下）；
% %                 for m=3:3 %1-3 best first next均可 选用的best fit 是否改位NEXT FIT 1002日改为m=3
% %                 % pA nAlg 
% %                 pA(nAlg) = ParameterInitialize( ...
% %                              'whichStripH', i,...
% %                              'whichBinH',m, ...
% %                              'whichSortItemOrder',j, ... 
% %                              'whichRotation',k, ...
% %                              'whichRotationHori', TLUl);
% %                  nAlg=nAlg+1;
% %                 end
% %             end
% %         end
% %     end
% % end
% % nAlg = nAlg - 1;
