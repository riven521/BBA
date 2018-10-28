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
    3.3 依据strip高度和loadingrate先排序，后按strip相邻约束排序
4 新增3个m文件

V0831-2:
1 修订getOrderofSID: 做了比较大的简化调整，参考getOrderofLID函数, 测试了几次没有问题。但不保险
--------------------------------------------------------------------------------------
V0916-1:
1 修订initCheck中的3个托盘的问题（转置）
2 修订initCheck中的1或2个托盘的问题（转置）
3 修订initCheck中的3个车辆的问题（转置）
    if ~(tmpVeh(1)*1.5 <= tmpVeh(2))
V0916-2:    
1 修订Vehicle问题（最后一个不满Vehicle改用小）
难点：从大vehicle改小vehicle, LU->ITEM要重新算(高度有差); ITEM->STRIP要重新算（宽度有差）; STRIP->BIN要重新算（长度有差）。无需重新CHECKINPUT.
2 修改getReturnBBA函数 -> 由嵌套函数放入局部函数
    a. 增加LU_VehType=daMax.LU.LU_VehType;
    b. 增加第1行,但是第9个放入的: tmpShow[9]
3 完成Vehicle修订问题
    a. main函数中94行增加% ****************** 针对车型选择 获取 smalld ： 最后一车放入更小车型后的输出
    b. POST PROCESSING中 增加  % ****************** 针对车型选择 获取修订的 output ：修订3个output输出
    
--------------------------------------------------------------------------------------
V1001-1:
0 修订updatemargin,在作图之前完成
1 增加plot3dbpp功能 + 结合plotcube函数
2 增加Item和Strip的多个属性计算
3 增加甩尾功能(待完善)

V1001-2:
1 增加cpu函数(计算Item,Strip,Bin等属性计算)

V1001-3:
1 修复getITEMorder:bug: 出现nbcol2错误 -> 原因在于对Item优先高度排序,致使混合Strip过多且重复.
% 目前顺序 : 1: SID ; 2: isNonMixed; 3: Longth/Height; 4:Width; 5: LID; 6: Height
[~,order] = sortrows(tmpItem',[1,6, 4, 3, 2, 5 ],{'ascend','descend','descend','descend','descend','descend'});  

V1002-1:
1 增加plot3dstrip:在strip阶段进行作图;以Item为基础

V1002-2:
1 修正Bin内顺序. 目前Strip顺序固定后,进入不同Bin的结果会变差; 准备修改为每次进入新Bin,排除已安排的strip,余下strip重新进bin

V1002-3:
1 增加量大车头方案2个, 函数名HreStripToBin和HreStripToEachBin
2 LU.DOC的计算全部移动到cpu***函数中
3 Main函数:m修改为3

V1002-4:
1 甩尾单独放到函数HStripSW中

V1003-1:
1 新增:cpuStrip中单独增加Strip.isHeightFull的判断,依据|maxH-minH| >= 最高值的1/3也为full 
2 新增:增加Lu的HeightL属性,其堆垛高度不能超过给定值

V1003-2:
1 LU.maxL: 在当前车型下最大长宽高层次计算 Gpreproc函数

V1004-1:
1 清理d数据中无用的变量和属性, 规范化命名
2 cpuStrip内增加LU.LU_Strip, LU.CoordLUStrip的计算;
3 增加对Strip的Item或LU作图 plot3DStrip
4 对d的数据梳理了一遍

V1004-2:
1 采用返回给BBA的数据画图(3维) DONE
2 cpuBin内增加Bin是否需要平铺的判断条件 isPing

V1005-1(里程碑):
1 增加平铺功能 DONE (BUG MAYBE) main函数增加类似车型的平铺函数
2 cpuStrip修订Strip.isHeightFull判断条件,(过低的Strip也视为非Full)
3 cpuBin修订:修改Bin.isTileNeed函数 -> 增加即使not full但Item的层数仅为1次case 并不需要Tile
4 return bba 数据包含平铺

V1007-1
1 增加global全局变量 - 不同约束策略的考虑
2 cpuItem:Item高度满层定义为二选一

V1008-1
1: 增加plot

V1009-1
1: Main函数增加全部平铺要求(在车辆允许情况下),否则按甩尾平铺处理
2: 预处理Gpreproc种增加LU.nbLID计算
3: 增加COMP_STRUCT函数
4: cpuItem中增加Item.ISisMixTile:对Item.isNonMixed=0的单独对待并后期处理.
5: 修改cpuStripnbItem, 增加Strip.nbLU的计算; 在S2B排序时, 将LU数量多的放最前面.
6: cpustrip增加Strip.isHeightBalance高度均衡计算, 并修改Strip.isHeightFull计算模式
6: cpustrip增加Strip.lowestHeight计算
7: 修改甩尾排序方式(HStripSW)
8: 增加11个问题算例及目前效果

V1012-1
1: Copy-of 增加若干李浩算例
2: 高度增加对角线和所有堆垛最高值

V1012-2
1: 整合量大车头函数
2: 高度增加对角线和所有堆垛最高值


V1014-1
1: 增加ISisNonMixedLU = 1 ISisMixTileLU = 0 对LU堆垛进行排序
2: Gpreproc中增加LU.isNonMixed/isMixedTile的计算
3: 修改getLUorder, 增加LU.isNonMixed/isMixedTile排序

V1014-2
1: BBA-MAIN 返回BBA的输出值 增加中断后的LID新排序  output_LU_Seq增加第八行, 输出顺序.
2: 修改函数1: 量大车头方案2HreStripToBin(Bin,Strip,Item,LU,Veh,p) 确保正确
3: 修改cpuStripnbItem(Strip,Item,LU): 保证strip计算正确性. 但目前还无法确保对所有f都正确.

V1018-1
1: 修改函数1: 量大车头方案2HreStripToBin(Bin,Strip,Item,LU,Veh,p) 修改确保正确的bug 对
2: 增加nbLULID在LU和Strip结构体
3: 修改函数HStripToBin, 增加Strip.nbLULID对Strip排序
4: 修改cpuStripnbItem, 增加对Strip.nbLULID计算

V1019-1
1: 增加LU.isShuaiWei, 在HStripSW中计算, 目的在plot时可判断顺序
2: 增加甩尾与否: plotSolutionBBA m8 展示甩尾
3: 更新main函数平铺和最后一车对getReturnBBA后得到的结果的更新

V1019-2
1: 修改Item的高度满层定义, 

V1027-1
1: getReturnBBA(daMax)增加数据的Table格式（结构体->表格）
2： LU.ID = idExchange(LU.ID); 使得内部
V1027-2
1: getReturnBBA返回Table非矩阵
2: 增加finalCheck补充检查，并依靠order从输入数据获取SID/PID/LID等不变值；LWH也可以获取，但海从输出结果获取了
3: 车型选择增加updateItemMargin，修改margin
4: 测试多次运行未发现崩溃现象（但有bug出现，还是要多测试才行）

里程碑:
V1028-1
1 BBA-Main修改成从maind提取IuIdx个输入,而非从运算后的d中提取
2 BBA-Main修改99 没有属性的临时增加
3 BBA-Main增加 checkLU
4 BBA-Main增加 224 d2Array d3Array  do2Array do3Array flagTiled
5 BBA-Main增加 getdinThisVeh 从原始输入数据maind旋转d2 d3等数据
6 BBA-Main增加 getTableLU 将d.LU转换table
6 BBA-Main增加 T1 T23等向T的某个bin内的变化修改，不排序等
6 BBA-Main增加 599 T = getShowSeq(T) 获取排序的Seq

TODO:
MAXL等确实属性的调整；
OPID OSID等属性的注意
平铺和车型在最后一车的更改，只能以某个为准，目前是平铺 LASTVEHICLE 测验
返回值准确性的确定???
PLOT的修改???
返回数据的简单CHECK
甩尾排序还有问题, 可能.
展示顺序还有问题, 即相邻Strip,前后相同LID,可能
对方崩溃问题

平铺all strip排序有问题 DONE 但未更新
循环多次运行 崩溃否？ 否 DONE
************甩尾顺序综合属性选择***********
各种strip属性记录
各种其它难点
PPT书写
FIGURE返回测试

NOTE:
1: Item.Item_Strip(1,:) == iStrip 避免==右侧为i索引顺序. 因为Item_Strip的第一行可能并非从1开始,也非连续

TODO:
1: LU是否满托的计算 DONE
2: 混合Strip的非混合(且满层?)的优先orderItem
3: HItemToStrip增加改变ITEM order的函数, Item没问题, 但LU作图有问题.
CPUSTRIP 五次注释 移到I2S中; CPUBIN和Item2Bin也可能有问题, 明天排查

QA:
Q1: I2S: order, 对同ID的LU,以高度递减排序(同STRIP混合可以), 异Strip混合有问题. 容易形成最后不满维尾垛结合其它满垛的问题.
A : 

TODO:
1 1/4车头的平铺问题(修订Item.isHeightFull修订) DONE
4 修订数量多的在车头(基于LU数量而非Item数量) BASIC DONE
1 增加平铺功能 DONE
2 甩尾后展示顺序功能
3 最高限制和最多层次,哪个条件先达到按哪个做（或给选择，按高度或按层数） DONE
4 实现指定位置摆放
5 修订莫名BUG ncol2错误; DONE
6 内存溢出错误
7 甩尾功能完善 + 测试

TODO:
0 plot strip时更新LW和Coord，确保含margin的显示正确 DONE 
1 从getbestsol由指标判断->最好变为由摆放顺序和算法直接获取最优，无需指标判断。 DONE
2 增加同BIN/STRIP包含多个指标SID/PID等如何计算的问题？ DONE cell patcat
3 isAdjacent为基础增加更多解的check判断函数  
4 getStriporder：zs zl 等STRIP排序确保相同SID在一个里面 DONE 

TODO:
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

-----------------------------------------------------------------------
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
