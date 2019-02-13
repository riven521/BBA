function [LU,Item] = HItemToBin(LU,Item,Strip)
% 重要函数:Item放入Bin中
%% 获取item与bin的关系 Item_Bin  + 坐标 CoordItemBin
%% 获取LU与bin的关系 LU_Bin  + 坐标 CoordLUBin

%% 初始化
%nDim = 3;
nItem = size(Item.LWH,2); %具体使用的Item的数量
nLU= size(LU.LWH,2); %具体使用的LU的数量
    
% 新增变量初始化
Item.Item_Bin=zeros(2,nItem);
Item.CoordItemBin=zeros(2,nItem);
LU.LU_Bin=zeros(2,nLU);
LU.CoordLUBin=zeros(2,nLU); % LU.LU_Strip=zeros(2,nLU); %% NOTE: 在plot3DBPP时 出错, 发现此处赋值无用, so 删除该语句
[LU.leftA,LU.bottomA,LU.rightA,LU.topA] = deal(zeros(1,nLU));

iiStrip = 0;
% 循环每个bin 5555555555 非常重要的函数 55555555555555
for iBin=1:max(Strip.Strip_Bin(1,:))
    tmpItemSeq=1;  %每个bin内顺序从1开始
    tmpLUSeq=1;
    tmpStrip = []; %计算CoordItemBin长(高)度使用 %       tmpLWStrip = Strip.LW(:,Strip.Strip_Bin(1,:)==iBin); %此bin内剩余宽度和长(高)度
    nbStrip=numel(find(Strip.Strip_Bin(1,:)==iBin));    %nbStrip：该iBin内strip数量    
    % 循环bin内每个strip
    for iStrip=1:nbStrip %从安放的第一个Strip开始逐步放置
        iiStrip = iiStrip+1;
        tmpLUSeqinStrip=1;
        [~,thisStrip] = find(Strip.Strip_Bin(1,:)==iBin & Strip.Strip_Bin(2,:)==iStrip ); %此bin内第iStrip个strip        

        %%
        % 平铺要求 1:单个strip至少有多个 2: 不能是混合的 3: Item的层数要均匀,相差不能>=2
%         if Strip.isSingleItem(thisStrip) ~= 1 && Strip.isMixed(thisStrip) ~= 1
% %             flagItem = Item.Item_Strip(1,:) == thisStrip;
% %             idxItem = find(flagItem)
% %             Item.Layer(flagItem)
% %             difLayer = diff( Item.Layer(flagItem))
% %             di = find(abs(difLayer)>=2);
% %              if ~isempty(di) && ~isscalar(di), error('EEE'); end
% %              flagItem(di)
% %              flagItem(di+1)
% %             if abs (diff( Item.Layer(flagItem)) ) > 1 % 如果Item的层数 >= 2 则需要平铺
% %                 LU.LU_Item(1, : )
% %                 LU.LU_Item(2, : )
% %                 Item.Layer(flagItem)
% %                 for i=1:sum(flagItem)
% %                     
% %                 end
%             end
%             [~,y] = find(Item.Item_Strip(1,:) == thisStrip); %y是Item索引号,当前strip内
%             for i=1:length(y) % 每个Item的循环
%                 
%                 sum(LU.LU_Item(1,:) == y(i))
%                 Item.Item_Strip(2,flagItem)
%                 Item.LWH(:,flagItem)
% 
%             end
%             Strip.isFull(thisStrip)
%             Strip.isMixed(thisStrip)
%             Strip.isSingleItem(thisStrip)
%             Strip.loadingrateLimit(thisStrip)
%         end
        
        % 平铺要求 1:单个strip只有1个, 但层数过多, 可以周边平铺
        
        
        % 平铺要求 1:单个strip只有多个, 后面的可以往前平铺
        
 %%       
        if ~isscalar(thisStrip), 
            error('意外错误'); 
        end
        nbItem=numel(find(Item.Item_Strip(1,:)==thisStrip));
        % 循环strip每个item,  从第一个顺序开始逐步获取ITEM在BIN内的坐标, 依据ITEM_STRIP的顺序
        for iItem=1:nbItem %同一strip内 循环取值iItem
            [~,thisItem] = find(Item.Item_Strip(1,:)==thisStrip & Item.Item_Strip(2,:)==iItem );    %此Strip内第iItem个item
            if ~isscalar(thisItem), error('意外错误');  end
            
            % 更新itemBeBinMatrix 555
            Item.Item_Bin(1,thisItem)=iBin;  % 本thisItem所在iBin的序号
            Item.Item_Bin(2,thisItem)=tmpItemSeq;  % 本thisItem在本iBin内的顺序 从1开始
            tmpItemSeq=tmpItemSeq+1;     
            % 更新CoordItemBin(按strip更新而非Item) 555 HAVE checked?
            Item.CoordItemBin(1,thisItem) = Item.CoordItemStrip(1,thisItem); % 55 本thisItem的x坐标不变
                % iStrip = 1： 首个strip高度一定为0；iStrip=2：高度一定不为0
            Item.CoordItemBin(2,thisItem) = sum(Strip.LW(2,tmpStrip)); % 555 错误原因2：thisStrip顺序不确定 % 555 LWStripSort错误原因:iStrip每个bin重新取样;LWStripSort是所有的Strip %             Item.CoordItemBin(2,thisItem) = sum(tmpLWStrip(2,1:iStrip-1)); 
            Item.CoordItemBin(3,thisItem) = 0; %Item均从0开始
            
            % 增加对LU的更新
            tmpLU = [];             %计算CoordLUBin高度使用 % tmpLWLU = LU.LWH(:,LU.LU_Item(1,:)==thisItem); 
            nbLU=numel(find(LU.LU_Item(1,:)==thisItem));                  
            % 循环item每个LU, 从第一个顺序开始逐步获取LU在BIN内的坐标
            for iLU=1:nbLU                
                [~,thisLU] = find(LU.LU_Item(1,:)==thisItem & LU.LU_Item(2,:)==iLU);
                if ~isscalar(thisLU), error('意外错误');  end
                
                                            % 更新LU_Strip % NOTE: 在plot3DBPP时 出错, 发现此处赋值无用, so 删除该语句
                                            %                 LU.LU_Strip(1,thisLU)=iiStrip;
                                            %                 LU.LU_Strip(2,thisLU)=tmpLUSeqinStrip;

                tmpLUSeqinStrip=tmpLUSeqinStrip+1;
                                            % 更新LURotaed 555  % LURotaed(1,thisLU)=Item.Rotaed(1,thisItem);
                % 更新LU_Bin 555
                LU.LU_Bin(1,thisLU)=iBin;
                LU.LU_Bin(2,thisLU)=tmpLUSeq;
                tmpLUSeq=tmpLUSeq+1;
                
                % 更新CoordLUBin 555 HAVE checked
                LU.CoordLUBin(1,thisLU) = Item.CoordItemBin(1,thisItem);
                LU.CoordLUBin(2,thisLU) = Item.CoordItemBin(2,thisItem);                
                LU.CoordLUBin(3,thisLU) = sum(LU.LWH(3,tmpLU));  %555 错误原因同上% LU.CoordLUBin(3,thisLU) = sum(tmpLWLU(3,1:iLU-1));
                tmpLU = [tmpLU thisLU];
                
                % 新增：LU的四个坐标
                LU.leftA(thisLU) = LU.CoordLUBin(1,thisLU);
                LU.bottomA(thisLU) = LU.CoordLUBin(2,thisLU);
                LU.rightA(thisLU) = LU.leftA(thisLU) + LU.LWH(1,thisLU);
                LU.topA(thisLU) = LU.bottomA(thisLU) + LU.LWH(2,thisLU);        
                
            end
        end
        tmpStrip = [tmpStrip thisStrip];
    end
end

            %NOTE: 此处把LU在ITEM阶段的旋转赋值给LU了   
            %GcheekInput放入的对比的是与原始值，此次要替换为再一次相对原始值.
            % LU.Rotaed = LURotaed; % 下面不需要了, 因为LUArray.Rotaed已经再HItemToStrip更新过了


% printstruct(d);
%     printscript();
    
% 输出主要结果:获得从1开始每个bin包含的item数据
% CoordItemBin Item_Bin
% %     function printscript()
% %         nBin = max(Item.Item_Bin(1,:));
% %         for iBin = 1: nBin
% %             [~,idx] = find(Item.Item_Bin(1,:)==iBin); %本iBin下的item索引号
% %             idxSeq = Item.Item_Bin(2,idx); %本iBin内item放入顺序Seq
% %             fprintf('bin 的宽+长+高为: ' );
% %             fprintf(' %d  ',Veh.LWH);
% %             fprintf('\n');
% %             fprintf('bin %d 的剩余宽+剩余长为:  ',iBin);fprintf('\n');
% %             fprintf('( %d ) ',Bin.LW(:,iBin));fprintf('\n');
% %             fprintf('\n');
% %             
% %             fprintf('bin %d 包含 original item 索引号{顺序}(长宽)[旋转标志]{坐标}为  \n  ',iBin);
% %             fprintf('%d ',idx);fprintf('\n');
% %             fprintf('{%d} ',idxSeq);fprintf('\n');
% %             fprintf(' (%d %d %d) ', Item.LWH(1:nDim,idx));fprintf('\n');
% %             fprintf(' [%d]     ', Item.Rotaed(:,idx));fprintf('\n');
% %             fprintf(' {%d %d %d} ', Item.CoordItemBin(:,idx));fprintf('\n');
% %             fprintf('\n');
% %             
% %             [~,idxLU] = find(LU.LU_Bin(1,:)==iBin); %本iBin下的item索引号
% %             fprintf('bin %d 包含 original LU 索引号{顺序}[item序号](长宽)[旋转标志]{坐标}为  \n  ',iBin);
% %             idxLUSeq = LU.LU_Bin(2,idxLU); %本iBin内item放入顺序Seq
% %             idxLUItem = LU.LU_Item(1,idxLU);
% %             fprintf('%d ',idxLU);fprintf('\n');
% %             fprintf('{%d} ',idxLUSeq);fprintf('\n');
% %             fprintf('[%d] ',idxLUItem);fprintf('\n');
% %             fprintf(' (%d %d %d) ', LU.LWH(1:nDim,idxLU));fprintf('\n');
% %             fprintf(' [%d]     ', LU.Rotaed(:,idxLU));fprintf('\n');
% %             fprintf(' {%d %d %d} ', LU.CoordLUBin(:,idxLU));fprintf('\n');
% %             fprintf('\n');
% %             
% %             % 按安放顺序展示
% %             % %     [~,x]=sort(LU.LU_Bin(2,idxLU));
% %             % %     idxLUSeq = idxLUSeq(x); %本iBin内item放入顺序Seq
% %             % %     idxLUItem = idxLUItem(x);
% %             % %     fprintf('%d ',idxLU);fprintf('\n');
% %             % %     fprintf('{%d} ',idxLUSeq);fprintf('\n');
% %             % %     fprintf('[%d] ',idxLUItem);fprintf('\n');
% %             % %     fprintf(' (%d %d %d) ', LU.LWH(1:nDim,idxLU(x)));fprintf('\n');
% %             % %     fprintf(' [%d]     ', LU.LURotaFlag(:,idxLU(x)));fprintf('\n');
% %             % %     fprintf(' {%d %d %d} ', LU.CoordLUBin(:,idxLU(x)));fprintf('\n');
% %             % %     fprintf('\n');
% %         end
% %     end
end
