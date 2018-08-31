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
2 [thisLevel,iLevel] = getThisLevel(iItem,iLevel,sItem, Strip, p); 增加输出iLevel,否则NEXTFIT 变为 FIRSTFIT

V0805-4
1 小的 更改

V0806-1
1 Strip.Strip_Item: 增加该变量，行1：该strip内包括Item数量；并增加行2，放入该Strip内包括Item种类

V0808-1
1 更新LU的排序，tmpLUMatrix = [LU.SID; LU.ID; LU.PID; LU.LWH]; sortrows(tmpLUMatrix',[1, 2, 3, 6],{'ascend','ascend','ascend','descend'}); 
2 更新ITEM进入Strip后的排序; ITEM按高度排序

V0808-2
1 更新HLUtoItem,原先算是BUG，采用NEXT FIT方式安排LU到ITEM

V0808-3
1 ITEM.ID -> ITEM.LID
2 增加isSameCol函数，判断struct的矩阵内列数全部相同
3 增加reorderStruct函数，对结构体S进行重新排序
4 Item旋转后，LUArray也要旋转相关,及时更新  已由后续处理转为嵌入到rotateItem
5 所有用到structfun地方全部做了struct的判断，并增加reorderStruct的顺序处理

V0809-1
1 ITEM 新增PID （多行n列） 0为未出现，其它为出现次数
2 STRIP 新增LID/PID/SID/UID等 （多行n列）
3 BIN 新增LID/PID/SID/UID等 （多行n列） - 但不好计算。

V0809-2
1 Gpreproc增加pwhichSortItemOrder 判断初始物品如何摆放，全部横放或全部竖放
2 修改striporder：优先SID，其次？ 增加fff函数

V0809-3
1 放弃HItemToStrip中对ITEM的顺序
2 

V0811-1
1 放弃HItemToStrip的rotateItem及对应update：右侧剩余间隙很小，但Hori改Rota可以放也不采用，因为会增加strip高度
同时，必须修改getThisLevel中的可旋转的flag，取消旋转条件允许的旋转。
即：在同一strip内，严禁任何ITEM再进行旋转；所有旋转在数据初始化完成，并给予ITEM排序给方案，中间不允许再旋转。1旋转就会增加STRIP高度。
2 在Gpreproc增加pwhichSortItemOrder=3->placeItemHori=3 判断初始物品如何摆放，全部横放或全部竖放(1/2),按横竖摆放右侧剩余间隙大小决定横竖放（3）
3 修改HItemToStrip对ITEM的排序，SID 第一，其次ITEM的长度和宽度，再其次LUID，确保相同ID摆放一起（STRIP2BIN确保类别越单纯的越在前端）

V0811-2
1 大幅度修改主函数输出（参数1：Coord 参数2：LWH 参数3：LU显示顺序，方便合并输出）
2 大幅度修改HItemToBin(LU,Item,Strip)；
3 确保零部件相邻：在LU2ITEM排序保障

V0811-3 ：实现margin约束
1 保障LU BUFF 为0
2 增加LU margin 4个值全为2
3 增加LU margin -> LU.LWH
4 旋转 getRotaedLWH -> 参数3由buff变为margin 确保旋转时margin随之变化
5 作图(包括main输出)时更新LU ITEM的Coord; LU ITEM的LW;
6 [d.LU,d.Item] = updateItemMargin(d.LU,d.Item);
7 return回bb的主程序李增加 updateItemMargin
8 555 重大bug修订：unique问题, updateItemMargin.m

V0811-4: 返回刘强的版本
1 修改输入输出
2 修改plot画图
3 修改best旋转标准

-------------------------------------------------------------
存在BUG的版本：
V0820-1:
1 在HItemToBin中 增加LU.LU_Strip
2 LU2Item 获取Item.PID2等
3 Item2Strip 获取Strip.PID2等 获取LU.LU_Strip第一行
V0822-1:
1 除LU的ID外，其它ID修改全部为V2版基于cell的SID/PID;
V0829-1:
2 待修改bug-确保strip内相连的LU,必须后面是单纯的LU
同一SID下区分LID单纯，单纯后的必须相连，且混合LID可能非连续！！！！！！！！！！！！！！！
-------------------------------------------------------------
bug修订的版本：
V0831-1:
1 修改bug-确保strip的排序, 总体按strip高度和loadingrate排序,内部按strip相邻约束排序
2 除Lu外，其它的PID/SID/LID全部采用cell方式,注意和patcat函数结合使用
3 HItemToStrip的order简单获取； HStripToBin的order复杂获取：
    3.1 getOrderofSID获取strip对应的sid的order为第一优先；
    3.2 getOrderofLID获取strip再sid内的order（改为priority）




修订Vehicle问题（最后一个不满Vehicle改用小）



[d.LU,d.Item,d.ItemID] = HLUtoItem(d.LU,d.Veh); 
[d.Item,d.Strip] = HItemToStrip(d.LU,d.Item,d.Veh,p);
 d = computeLoadingRateStrip(d);
 [d.Strip,d.Bin]= HStripToBin(d.Strip,d.Veh,d.LU,p);
  [d.LU,d.Item] = HItemToBin(d.LU,d.Item,d.Strip);
  d = computeLoadingRate2DBin(d);



V0811-1 ：实现指定位置摆放


TODO:
0 plot strip时更新LW和Coord，确保含margin的显示正确
1 从getbestsol由指标判断->最好变为由摆放顺序和算法直接获取最优，无需指标判断。
2 增加同BIN/STRIP包含多个指标SID/PID等如何计算的问题？
3 isAdjacent为基础增加更多解的check判断函数
4 getStriporder：zs zl 等STRIP排序确保相同SID在一个里面




1 增加strip新（LUID类别，LU剩余宽度是否可放另一个小LU，即是否满level）
2 每次算法对strip摆放都要尝试hori和vert两种方式，分别判断不同类型LU如何摆放为好，暂不考虑部分vert，部分hori在同一level情形，（可能与其它非空strip拼载可用考虑） -> 改为先确固定hori和vert摆放，后期不调整的方式（依赖顺序）

1 FPRINT
2 test main error （modifyStripWithOneItem(d);修改）
3 车型约束5-边界约束4-摆放相邻约束1-摆放位置约束4

TODO:
1 分别每种LUID分别做Hori=1/2，而非全文统一 NEXTFIT 4 STRIPH
2 供应商ID增加, 确保同供应商放在一起（按供应商排序，由指定顺序，从小到大排序）（LU按供应商排序，再按高度排序）
2 供应商在bin内排序要一起，而非strip; strip先按供应商排序,后按高度排序
3 逻辑确保多供应商时没有摆放位置要求
4 strip在bin内局部优化，也可以bin之间调整，确保满足要求即可。
5 供应商之间衔接时，注意可以一起堆垛并连接的strip高度变异最小，可以固定strip选LUID
6 暂且规定不同供应商即使LUID相同，也不能堆垛；否则无法获取堆垛后ITEM的SID等。后期可局部优化。

1 对于可变高度托盘的堆垛，也需要考虑BEST FIRST NEXT

TODO
1：用户输入和对应算法参数，以及运算结果，返回BBA结果全部保存到MAT文件，方便后期核验和重复算法。
2：增加matlab画图功能，可视化每一步
3：改输入参数为非nested structure,并不在内部取structure（局部函数） DONE
4：初始化每个struct结构，采用initilize方式
5：采用ployshape函数画图，和表示strip代替系数

6：STRIP给高度（最大的堆垛高度3/最小4）


1：研发投入，如名门项目落地，从项目分成中抵扣；（5/5）10w起；
0：0研发投入，如名门项目落地，在既定项目分成比例提高20%；（7/3）；

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
