V12 180710 
1 增加间隙，输入参数5和6
增加第五个输入 LUBUFF = [1,1]; %长宽间隙
增加第六个输入 BINBUFF = [1,1,1]; 
2 增加GcheckInput中对BUFF的数据判断
% 输入间隙的转换
% 输入转换2: 增加间隙后的转换
3 

V0730
重要更新可能bug；允许个性化LU旋转设置
1 增加函数：getLUIDArray 获取LUID类型相关数据(同类型ID的体积，重量，是否可旋转)
2 增加函数：getITEMIDArray 获取ItemID类型相关数据(同类型ItemID的体积，面积，重量，item是否可旋转)
3 修改函数：Main()->ParameterInitialize 初始化不同参数后的算法计算
4 删除：ParaArray.whichRotationAll ParaArray.whichRotationBin LU和BIN颠倒

大幅度修订：
5 增加/修改有关旋转的配置：
	5.1 LU和ITEM增加isRota：LU/ITEM是否允许ROTATE的标记（statistical）默认1均可旋转,0表示不可以旋转
	5.2 LU和ITEM增加Rotaed：LU/ITEM是否做了ROTATE的标记（dynamic）
	5.3 LU和ITEM的LWH：LU/ITEM做了ROTATE后的LWH（dynamic）（NOTE:不保留原始LWH）
	5.4 whichRotation = 1/0 全部允许/禁止ROTATE的标记用途不大；依据LU/ITEM的ISROTA判断
	5.5 增加函数: placeItemHori (按照Hori/Vert方式摆放,返回是否需要对LU/Item旋转)+ getRotaedLWH（依据返回值对LU/ITEM+BUFF进行更新）

V0731-1
1 增加代码：HStripToBin-> 增加first/next 安排strip到bin
V0731-2
1 增加代码：BBA_MAIN-> 注释格式化
2 增加函数:BBA_MAIN-> d = DataInitialize( varargin )

V0803-1
1 修改代码: BBA_MAIN-> (输入改用varargin)
2 增加变量: StripWeight重量增加到HItemToStrip中
3 增加变量: BinWeight重量增加到HStripToBin中
4 增加bin重量：不超过车辆最大重量（best选bin过程是否考虑重量,目前仅考虑维度）

V0803-2
函数无变化,仅命名做了调整
1 find and replace 所有文件 替换**Array -> **
2 find and replace 所有文件 替换**BE-- -> **_--
3 find and replace 所有文件 替换da -> d

V0803-3
函数无变化,仅命名做了调整
1 find and replace paArray -> pA  daArray -> dA flagArray -> flagA
1 find and replace 4个output
1 find and replace
V0803-4
1 HLUToItem: 更新函数输入输出,方面了解哪些参数需要被使用
V0803-5
1 HStripToBin: 更新函数输入输出,方面了解哪些参数需要被使用

V0805-1
1 d初始化
2 所有子函数不传递d,传递d.LU等具体子struc
3 完成除LU->ITEM STRIP->BIN

1 检验输入更改 
2 Main函数 BBA_Main(LUID,LULWH,BINID,BINLWH); 前4参数必须,其它参数可有可无; 如无则自己生成随机数据; 如有则再swith中获取
3 polyshape使用(画LU/Item的图) 待改进
4 555 完成DataInitiliza初始化 增加较多输入参数

V0805-2
1 增加随机生成的车型按volume从大到小排序
2 数据初始化即开始: [d.LU, d.Veh] = initCheck(d.LU,d.Veh); 
3 更改HItemToStrip输入输出

V0805-3
1 增加多车型首选第一个车型（默认应该是体积最大车型）- 修改Veh.LWH相关参数
2  

TODO:
FPRINT
2 test main error 画图
3 modifyStripWithOneItem(d);修改

2 修改局部函数HItemToStrip
为多车型（随机生成，车型按从大到小排序，首先安排第一种大车型）

1 完成所有的局部函数
2 


TODO
1：用户输入和对应算法参数，以及运算结果，返回BBA结果全部保存到MAT文件，方便后期核验和重复算法。
2：增加matlab画图功能，可视化每一步
3：改输入参数为非nested structure,并不在内部取structure（局部函数）
4：初始化每个struct结构，采用initilize方式
5：采用ployshape函数画图，和表示strip代替系数

6：STRIP给高度（最大的堆垛高度3/最小4）

sign
Y = sign(x) 返回与 x 大小相同的数组 Y，其中 Y 的每个元素是：
    1，前提是 x 的对应元素大于 0。
    0，前提是 x 的对应元素等于 0。
    -1，前提是 x 的对应元素小于 0。
false

v = v(:).'; % Make sure v is a row vector.
lh ~= zeros(size(lh))
ones(length(r),1);
r = find(x)
idxs(r,:)
find(unvisited,1);
visited = nextpt; unvisited(startour) = 0; % update unvisited points
size/find/reshape/error/double/prod/round/drawnow
idx/true


1 同一类型的托盘，前提：一辆车可以放下该类型托盘；可否分开放入多辆车？
2 同一供应商的托盘，前提：一辆车可以放下该类型托盘；可否分开放入多辆车？
3 最小的托盘的长宽/最大的托盘长宽（是否超过车型的宽度）
4 每家供应商的托盘类型最多数量（）
5 重物的定义

5 增加数据核对函数（重要部分）




1 whichRotation == 1 等 whichRotationHori = 0 等调整 DONE
2 Rotation旋转等变化 DONE 使用到placeItemHori的地方都有getRotaedLWH -> 即在获取必须旋转的标记后，对LWH进行调整，包括对BUFF的调整；
	在HItemToStrip增加如下：调整LU的Rotaed更新->后期无需更新

    % LUArray旋转相关,及时更新    
    nbItem=length(d.ItemArray.Rotaed);
    % 循环每个item
    for idxItem=1:nbItem
        flagThisItem = (d.LUArray.LUBeItemArray(1,:)==idxItem );
        % 对应位置LU.Rotaed更新
        if d.ItemArray.Rotaed(idxItem)
            d.LUArray.Rotaed(flagThisItem) = ~d.LUArray.Rotaed(flagThisItem);
            % 对应位置LU.LWH更新
            d.LUArray.LWH(1, flagThisItem) = d.ItemArray.LWH(1, idxItem);
            d.LUArray.LWH(2, flagThisItem) = d.ItemArray.LWH(2, idxItem);
        end
    end

3 Main函数对LU的返回 OK
4 PLOT问题 OK



        if ParaArray.whichRotation == 1 %TODO 
            da.LUArray.isRota = ones(size(da.LUArray.isRota));
