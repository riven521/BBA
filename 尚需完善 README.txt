TODO :

2 材料增加考虑（通过托盘ID号辨别 DONE）

1 重量增加考虑
5 重货判定及位置

3 托盘Buffer增加考虑（长宽）- 可作为托盘间间隙的距离*2
4 车辆buffer增加考虑（长宽高） - 可作为边支LU间隙的调节器


%长宽高间隙

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

V0801-1
1 修改代码: BBA_MAIN-> (输入改用varargin)
2 增加变量: StripWeight重量增加到HItemToStrip中
3 增加变量: BinWeight重量增加到HStripToBin中
4 增加bin重量：不超过车辆最大重量（best选bin过程是否考虑重量,目前仅考虑维度）


TODO
1：用户输入和对应算法参数，以及运算结果，返回BBA结果全部保存到MAT文件，方便后期核验和重复算法。
2：增加matlab画图功能，可视化每一步
3：改输入参数为非nested structure,并不在内部取structure


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
