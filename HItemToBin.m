function [da] = HItemToBin(da)
% 重要函数:Item放入Bin中
%% 获取item与bin的关系 itemBeBinMatrix  + 坐标 CoordItemBin
%% 获取LU与bin的关系 LUBeBinMatrix  + 坐标 CoordLUBin
%% 获取LU与item的关系 LURotaFlag
%% 初始化
nDim = size(da.ItemArray.LWH,1);  
nItem = size(da.ItemArray.LWH,2); %具体使用的Item的数量
nLU= size(da.LUArray.LWH,2); %具体使用的LU的数量
uniBinDataMatrix = unique((da.BinArray.LWH(1:nDim,:))','rows')';
%% 结构体提取
    LUBeItemArray = da.LUArray.LUBeItemArray;
stripBeBinMatrix = da.StripArray.stripBeBinMatrix;
itemBeStripMatrix = da.ItemArray.itemBeStripMatrix;
itemRotaFlag = da.ItemArray.itemRotaFlag;
CoordItemStrip = da.ItemArray.CoordItemStrip;
% LWStripSort = da.StripArray.LW(:,da.StripArray.striporder);
LWStrip = da.StripArray.LW;
LWHLU = da.LUArray.LWH;
    
itemBeBinMatrix=zeros(2,nItem);
CoordItemBin=zeros(2,nItem);
    LUBeBinMatrix=zeros(2,nLU);
    CoordLUBin=zeros(2,nLU);
    LURotaFlag=zeros(1,nLU);
    
for iBin=1:max(stripBeBinMatrix(1,:))
    tmpItemSeq=1;  %每个bin内顺序从1开始
    tmpLUSeq=1;
    nbStrip=numel(find(stripBeBinMatrix(1,:)==iBin));       
%         tmpLWStrip = LWStrip(:,stripBeBinMatrix(1,:)==iBin); %此bin内剩余宽度和长(高)度
          tmpStrip = []; %计算CoordItemBin长(高)度使用
    for iStrip=1:nbStrip %从安放的第一个Strip开始逐步放置
        [~,thisStrip] = find(stripBeBinMatrix(1,:)==iBin & stripBeBinMatrix(2,:)==iStrip );%此bin内第iStrip个strip
            
        if isempty(thisStrip) || numel(thisStrip)>1, error('意外错误');  end
        nbItem=numel(find(itemBeStripMatrix(1,:)==thisStrip));        
        for iItem=1:nbItem %同一strip内 循环取值iItem
            [~,thisItem] = find(itemBeStripMatrix(1,:)==thisStrip & itemBeStripMatrix(2,:)==iItem );%此Strip内第iItem个item
            if isempty(thisItem) || numel(thisItem)>1, error('意外错误');  end
            % 更新itemBeBinMatrix 555
            itemBeBinMatrix(1,thisItem)=iBin; % 本thisItem所在iBin的序号
            itemBeBinMatrix(2,thisItem)=tmpItemSeq; % 本thisItem在本iBin内的顺序 从1开始
            tmpItemSeq=tmpItemSeq+1;     
            % 更新CoordItemBin(按strip更新而非Item) 555 HAVE checked?
            CoordItemBin(1,thisItem) = CoordItemStrip(1,thisItem); % 55 本thisItem的x坐标不变
                % iStrip = 1： 首个strip高度一定为0；iStrip=2：高度一定不为0
                CoordItemBin(2,thisItem) = sum(LWStrip(2,tmpStrip)); % 555 错误原因2：thisStrip顺序不确定 % 555 LWStripSort错误原因:iStrip每个bin重新取样;LWStripSort是所有的Strip %             CoordItemBin(2,thisItem) = sum(tmpLWStrip(2,1:iStrip-1)); 
            CoordItemBin(3,thisItem) = 0; %Item均从0开始
            
            % 增加对LU的更新
            nbLU=numel(find(LUBeItemArray(1,:)==thisItem));
%                 tmpLWLU = LWHLU(:,LUBeItemArray(1,:)==thisItem); %计算
                tmpLU = []; %计算CoordLUBin高度使用
            for iLU=1:nbLU
                [~,thisLU] = find(LUBeItemArray(1,:)==thisItem & LUBeItemArray(2,:)==iLU);
                if isempty(thisLU) || numel(thisLU)>1, error('意外错误');  end
                % 更新LURotaFlag 555
                LURotaFlag(1,thisLU)=itemRotaFlag(1,thisItem);
                % 更新LUBeBinMatrix 555
                LUBeBinMatrix(1,thisLU)=iBin;
                LUBeBinMatrix(2,thisLU)=tmpLUSeq;
                tmpLUSeq=tmpLUSeq+1;               
                % 更新CoordLUBin 555 HAVE checked
                CoordLUBin(1,thisLU) = CoordItemBin(1,thisItem);
                CoordLUBin(2,thisLU) = CoordItemBin(2,thisItem);                
                    CoordLUBin(3,thisLU) = sum(LWHLU(3,tmpLU));  %555 错误原因同上% CoordLUBin(3,thisLU) = sum(tmpLWLU(3,1:iLU-1));
                    tmpLU = [tmpLU thisLU];
            end            
        end        
        tmpStrip = [tmpStrip thisStrip];        
    end
end

%% 结构体赋值
da.ItemArray.itemBeBinMatrix=itemBeBinMatrix;
da.ItemArray.CoordItemBin=CoordItemBin;
    da.LUArray.LUBeBinMatrix= LUBeBinMatrix;
    da.LUArray.CoordLUBin = CoordLUBin;
    da.LUArray.LURotaFlag = LURotaFlag;
% printstruct(da);
% 输出主要结果:获得从1开始每个bin包含的item数据
% CoordItemBin itemBeBinMatrix
for iBin = 1:max(da.ItemArray.itemBeBinMatrix(1,:))
    [~,idx] = find(da.ItemArray.itemBeBinMatrix(1,:)==iBin); %本iBin下的item索引号
    idxSeq = da.ItemArray.itemBeBinMatrix(2,idx); %本iBin内item放入顺序Seq
    fprintf('bin 的宽+长+高为: ' );
    fprintf(' %d  ',uniBinDataMatrix);
    fprintf('\n');
    fprintf('bin %d 的剩余宽+剩余长为:  ',iBin);fprintf('\n');
    fprintf('( %d ) ',da.BinSArray.LW(:,iBin));fprintf('\n');
    fprintf('\n');
    
    fprintf('bin %d 包含 original item 索引号{顺序}(长宽)[旋转标志]{坐标}为  \n  ',iBin);
    fprintf('%d ',idx);fprintf('\n');
    fprintf('{%d} ',idxSeq);fprintf('\n');
    fprintf(' (%d %d %d) ', da.ItemArray.LWH(1:nDim,idx));fprintf('\n');
    fprintf(' [%d]     ', da.ItemArray.itemRotaFlag(:,idx));fprintf('\n');
    fprintf(' {%d %d %d} ', da.ItemArray.CoordItemBin(:,idx));fprintf('\n');
    fprintf('\n');
    
    [~,idxLU] = find(da.LUArray.LUBeBinMatrix(1,:)==iBin); %本iBin下的item索引号
    fprintf('bin %d 包含 original LU 索引号{顺序}[item序号](长宽)[旋转标志]{坐标}为  \n  ',iBin);
    idxLUSeq = da.LUArray.LUBeBinMatrix(2,idxLU); %本iBin内item放入顺序Seq
    idxLUItem = da.LUArray.LUBeItemArray(1,idxLU);
    fprintf('%d ',idxLU);fprintf('\n');    
    fprintf('{%d} ',idxLUSeq);fprintf('\n');
    fprintf('[%d] ',idxLUItem);fprintf('\n');
    fprintf(' (%d %d %d) ', da.LUArray.LWH(1:nDim,idxLU));fprintf('\n');
    fprintf(' [%d]     ', da.LUArray.LURotaFlag(:,idxLU));fprintf('\n');
    fprintf(' {%d %d %d} ', da.LUArray.CoordLUBin(:,idxLU));fprintf('\n');
    fprintf('\n');    

    % 按安放顺序展示
    [~,x]=sort(da.LUArray.LUBeBinMatrix(2,idxLU));    
    idxLUSeq = idxLUSeq(x); %本iBin内item放入顺序Seq
    idxLUItem = idxLUItem(x);
    fprintf('%d ',idxLU);fprintf('\n');
    fprintf('{%d} ',idxLUSeq);fprintf('\n');
    fprintf('[%d] ',idxLUItem);fprintf('\n');
    fprintf(' (%d %d %d) ', da.LUArray.LWH(1:nDim,idxLU(x)));fprintf('\n');
    fprintf(' [%d]     ', da.LUArray.LURotaFlag(:,idxLU(x)));fprintf('\n');
    fprintf(' {%d %d %d} ', da.LUArray.CoordLUBin(:,idxLU(x)));fprintf('\n');
    fprintf('\n');
end
end
