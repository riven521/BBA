function [d] = HItemToBin(d)
% 重要函数:Item放入Bin中
%% 获取item与bin的关系 Item_Bin  + 坐标 CoordItemBin
%% 获取LU与bin的关系 LU_Bin  + 坐标 CoordLUBin
%% 获取LU与item的关系 LURotaFlag

%% 初始化
nDim = size(d.Item.LWH,1);  
nItem = size(d.Item.LWH,2); %具体使用的Item的数量
nLU= size(d.LU.LWH,2); %具体使用的LU的数量

%% 结构体提取
LU_Item = d.LU.LU_Item;
stripBeBinMatrix = d.Strip.Strip_Bin;
Item_Strip = d.Item.Item_Strip;

itemRotaed = d.Item.Rotaed;
% % itemRotaFlag = d.Item.itemRotaFlag;
CoordItemStrip = d.Item.CoordItemStrip;

LWStrip = d.Strip.LW; % LWStripSort = d.Strip.LW(:,d.Strip.striporder);
LWHLU = d.LU.LWH;
    
Item_Bin=zeros(2,nItem);
CoordItemBin=zeros(2,nItem);

LU_Bin=zeros(2,nLU);
CoordLUBin=zeros(2,nLU);
LURotaed=zeros(1,nLU);
    
% 循环每个bin
for iBin=1:max(stripBeBinMatrix(1,:))
    tmpItemSeq=1;  %每个bin内顺序从1开始
    tmpLUSeq=1;
    tmpStrip = []; %计算CoordItemBin长(高)度使用 %       tmpLWStrip = LWStrip(:,stripBeBinMatrix(1,:)==iBin); %此bin内剩余宽度和长(高)度
    nbStrip=numel(find(stripBeBinMatrix(1,:)==iBin));    %nbStrip：该iBin内strip数量    
    % 循环bin内每个strip
    for iStrip=1:nbStrip %从安放的第一个Strip开始逐步放置
        [~,thisStrip] = find(stripBeBinMatrix(1,:)==iBin & stripBeBinMatrix(2,:)==iStrip ); %此bin内第iStrip个strip      
        if ~isscalar(thisStrip), error('意外错误');  end
        nbItem=numel(find(Item_Strip(1,:)==thisStrip));
        % 循环strip每个item
        for iItem=1:nbItem %同一strip内 循环取值iItem
            [~,thisItem] = find(Item_Strip(1,:)==thisStrip & Item_Strip(2,:)==iItem );%此Strip内第iItem个item
            if ~isscalar(thisItem), error('意外错误');  end
            
            % 更新itemBeBinMatrix 555
            Item_Bin(1,thisItem)=iBin;  % 本thisItem所在iBin的序号
            Item_Bin(2,thisItem)=tmpItemSeq;  % 本thisItem在本iBin内的顺序 从1开始
            tmpItemSeq=tmpItemSeq+1;     
            % 更新CoordItemBin(按strip更新而非Item) 555 HAVE checked?
            CoordItemBin(1,thisItem) = CoordItemStrip(1,thisItem); % 55 本thisItem的x坐标不变
                % iStrip = 1： 首个strip高度一定为0；iStrip=2：高度一定不为0
            CoordItemBin(2,thisItem) = sum(LWStrip(2,tmpStrip)); % 555 错误原因2：thisStrip顺序不确定 % 555 LWStripSort错误原因:iStrip每个bin重新取样;LWStripSort是所有的Strip %             CoordItemBin(2,thisItem) = sum(tmpLWStrip(2,1:iStrip-1)); 
            CoordItemBin(3,thisItem) = 0; %Item均从0开始
            
            % 增加对LU的更新
            tmpLU = []; %计算CoordLUBin高度使用 % tmpLWLU = LWHLU(:,LU_Item(1,:)==thisItem); %计算
            nbLU=numel(find(LU_Item(1,:)==thisItem));            
            % 循环item每个LU
            for iLU=1:nbLU
                [~,thisLU] = find(LU_Item(1,:)==thisItem & LU_Item(2,:)==iLU);
                if ~isscalar(thisLU), error('意外错误');  end
                
                % 更新LURotaed 555
                LURotaed(1,thisLU)=itemRotaed(1,thisItem);
                % 更新LUBeBinMatrix 555
                LU_Bin(1,thisLU)=iBin;
                LU_Bin(2,thisLU)=tmpLUSeq;
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
d.Item.Item_Bin=Item_Bin;
d.Item.CoordItemBin=CoordItemBin;

d.LU.LU_Bin= LU_Bin;
d.LU.CoordLUBin = CoordLUBin;

%NOTE: 此处把LU在ITEM阶段的旋转赋值给LU了   
%GcheekInput放入的对比的是与原始值，此次要替换为再一次相对原始值.
%% % 下面不需要了, 因为LUArray.Rotaed已经再HItemToStrip更新过了
% d.LU.Rotaed = LURotaed;


% printstruct(d);
% 输出主要结果:获得从1开始每个bin包含的item数据
% CoordItemBin Item_Bin
for iBin = 1:max(d.Item.Item_Bin(1,:))
    [~,idx] = find(d.Item.Item_Bin(1,:)==iBin); %本iBin下的item索引号
    idxSeq = d.Item.Item_Bin(2,idx); %本iBin内item放入顺序Seq
    fprintf('bin 的宽+长+高为: ' );
    fprintf(' %d  ',d.Veh.LWH);
    fprintf('\n');
    fprintf('bin %d 的剩余宽+剩余长为:  ',iBin);fprintf('\n');
    fprintf('( %d ) ',d.Bin.LW(:,iBin));fprintf('\n');
    fprintf('\n');
    
    fprintf('bin %d 包含 original item 索引号{顺序}(长宽)[旋转标志]{坐标}为  \n  ',iBin);
    fprintf('%d ',idx);fprintf('\n');
    fprintf('{%d} ',idxSeq);fprintf('\n');
    fprintf(' (%d %d %d) ', d.Item.LWH(1:nDim,idx));fprintf('\n');
    fprintf(' [%d]     ', d.Item.Rotaed(:,idx));fprintf('\n');
    fprintf(' {%d %d %d} ', d.Item.CoordItemBin(:,idx));fprintf('\n');
    fprintf('\n');
    
    [~,idxLU] = find(d.LU.LU_Bin(1,:)==iBin); %本iBin下的item索引号
    fprintf('bin %d 包含 original LU 索引号{顺序}[item序号](长宽)[旋转标志]{坐标}为  \n  ',iBin);
    idxLUSeq = d.LU.LU_Bin(2,idxLU); %本iBin内item放入顺序Seq
    idxLUItem = d.LU.LU_Item(1,idxLU);
    fprintf('%d ',idxLU);fprintf('\n');
    fprintf('{%d} ',idxLUSeq);fprintf('\n');
    fprintf('[%d] ',idxLUItem);fprintf('\n');
    fprintf(' (%d %d %d) ', d.LU.LWH(1:nDim,idxLU));fprintf('\n');
    fprintf(' [%d]     ', d.LU.Rotaed(:,idxLU));fprintf('\n');
    fprintf(' {%d %d %d} ', d.LU.CoordLUBin(:,idxLU));fprintf('\n');
    fprintf('\n');

    % 按安放顺序展示
    % %     [~,x]=sort(d.LU.LU_Bin(2,idxLU));
    % %     idxLUSeq = idxLUSeq(x); %本iBin内item放入顺序Seq
    % %     idxLUItem = idxLUItem(x);
    % %     fprintf('%d ',idxLU);fprintf('\n');
    % %     fprintf('{%d} ',idxLUSeq);fprintf('\n');
    % %     fprintf('[%d] ',idxLUItem);fprintf('\n');
    % %     fprintf(' (%d %d %d) ', d.LU.LWH(1:nDim,idxLU(x)));fprintf('\n');
    % %     fprintf(' [%d]     ', d.LU.LURotaFlag(:,idxLU(x)));fprintf('\n');
    % %     fprintf(' {%d %d %d} ', d.LU.CoordLUBin(:,idxLU(x)));fprintf('\n');
    % %     fprintf('\n');
end
end
